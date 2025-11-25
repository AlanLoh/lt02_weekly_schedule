#! /usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Alan Loh"
__copyright__ = "Copyright 2025, exoschedule"
__credits__ = ["Alan Loh"]
__maintainer__ = "Alan"
__email__ = "alan.loh@obspm.fr"
__status__ = "Production"
__all__ = [
    "bookings_from_vcr",
    "booking_difference",
    "write_booking_file",
    "get_exoplanet_dict",
    "get_star_dict",
    "obsid_to_timerange",
    "get_source_directories",
    "get_source_exposure_time",
    "update_exposure_time",
    "update_minimal_elevation",
    "clear_exposure_time",
    "get_imaging_sources"
]


from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
import astropy.units as u
from typing import List, Tuple
import numpy as np
import operator
import os
import glob
import json
import subprocess
import copy
import logging
log = logging.getLogger("exoschedule")

from nenupy.schedule import Schedule, ReservedBlock
from nenupy.astro.target import FixedTarget

from exoschedule import (
    MIN_DURATION,
    MID_DAY_START_HOUR,
    MID_DAY_STOP_HOUR,
    KP_CODE,
    REMOTE_SERVER,
    SPECIAL_CASES,
    IMAGING_SOURCES
)


# ============================================================= #
# --------------------- bookings_from_vcr --------------------- #
def bookings_from_vcr(
        vcr_current_booking: str,
        start_time: Time,
        stop_time: Time,
        merge_slots: bool = True,
        only_kp_slots: bool = False,
        min_dt: TimeDelta = TimeDelta(60, format="sec")
    ) -> List[Tuple[Time, Time]]:
    """Computes a list of free booking slots from the current VCR booking schedule.

    The rules are:
    * No booking less than MIN_DURATION (by default 1 hour)
    * No booking in the UTC interval MID_DAY_START_HOUR -- MID_DAY_STOP_HOUR unless Saturday or Sunday

    Parameters
    ----------
    vcr_current_booking : str
        VCR booking file (that can be downloaded via 'https://gui-nenufar.obs-nancay.fr/' > 'Booking' > 'Current Booking.csv')
    start_time : Time
        Start time at which a new booking will be considered
    stop_time : Time
        Stop time at which a new booking will be considered
    merge_slots : bool
        If the vcr_current_booking already contains booking slots already assigned to the KP_CODE, merge them with the new ones.
        Otherwise consider them as already booked. By default True
    only_kp_slots : bool
        Restrict the booking to the already booked KP slots. If set to True, merge_slots argument is ignored. By default False
    min_dt : TimeDelta, optional
        Minimal time resolution of the booking time definition, by default TimeDelta(60, format="sec")

    Returns
    -------
    List[Tuple[Time, Time]]
        List of bookings in the shape [(start_0, stop_0), (start_1, stop_1), ...].
    """
    
    # Get the reserved slots currently scheduled within the VCR
    booked_slots = ReservedBlock.from_VCR(vcr_current_booking)
    if only_kp_slots:
        log.info(f"Limiting the booking to already booked {KP_CODE} slots.")
        # Only keeps slots belonging to the current KP
        booked_slots = booked_slots.get(name=KP_CODE, operation=operator.eq)
        # Produce the list of start and stop times and ignore the rest
        return [(blk.time_min, blk.time_max) for blk in booked_slots if blk.is_within(start_time, stop_time)]
    elif merge_slots:
        # Remove the slots for which the KP is already scheduled
        booked_slots = booked_slots.get(name=KP_CODE, operation=operator.ne)

    # Generate an instance of Schedule at the min_dt resolution and 
    # in between start and stop times.
    schedule = Schedule(time_min=start_time, time_max=stop_time, dt=min_dt)

    # Insert the reserved slots from the VCR and extract the free slots indices
    schedule.insert(booked_slots)
    free_slots = schedule.freeSlots

    free_slot_indices = np.arange(schedule.size)[schedule.freeSlots]

    # Avoid mid-day
    start_hours = np.array([dt_slot.hour for dt_slot in schedule.starts.datetime])
    week_end_mask = np.array([time_datetime.weekday() < 5 for time_datetime in schedule.starts.datetime], dtype=bool) # Saturday is 5, Sunday is 6
    mid_day_mask = (start_hours >= MID_DAY_START_HOUR) * (start_hours < MID_DAY_STOP_HOUR)
    free_slot_indices[(mid_day_mask * week_end_mask)[free_slots]] = False

    # Group consecutive indices
    free_groups_indices = np.split(free_slot_indices, np.where(np.diff(free_slot_indices) != 1)[0] + 1)

    # Return the list of free booking slots
    free_bookings = []
    for free_group in free_groups_indices:
        if len(free_group) == 1: continue
        group_start = schedule.starts[free_group[0]]
        group_stop = schedule.stops[free_group[-1]]

        group_duration = group_stop - group_start

        if (group_duration - MIN_DURATION) < (- 1 * u.s):
            log.info(f"Booking < min_duration ({MIN_DURATION.sec}s): {group_start.isot} - {group_stop.isot} --> skipped.")
            continue

        free_bookings.append( (group_start, group_stop) )

    return free_bookings

# ============================================================= #
# -------------------- booking_difference --------------------- #
def booking_difference(booking_1: List[Tuple[Time, Time]], booking_2: List[Tuple[Time, Time]]) -> List[Tuple[Time, Time]]:

    booking_1 = copy.copy(booking_1)

    for s2, e2 in booking_2:

        tmp = []
        idx_to_remove = []
        for i, (s1, e1) in enumerate(booking_1):
            # If the second booking is unrelated to the first, skip
            if (e2 < s1) or (s2 > e1):
                continue

            # If the start is before, align them, same for end
            if s2 < s1:
                s2 = s1
            if e2 > e1:
                e2 = e1
            
            # Remove the previous booking and create as many more
            idx_to_remove.append(i)

            boundaries = np.array([s1, s2, e2, e1])
            boundaries, counts = np.unique(boundaries, return_counts=True)
            for start, stop in zip(*[iter(sorted(boundaries[counts == 1]))]*2):
                tmp.append((start, stop))
        
        for idx in sorted(idx_to_remove)[::-1]:
            del booking_1[idx]

        for booking in tmp:
            booking_1.append(booking)

    return sorted(booking_1)


# ============================================================= #
# -------------------- write_booking_file --------------------- #
def write_booking_file(booking_starts_stops: List[Tuple[Time, Time]], output_path: str = "", file_prefix: str = "") -> None:
    """Write a NenuFAR booking file that can be directly uploaded to the Virtual Control Room.
    The name of the booking file is automatic and w

    Parameters
    ----------
    booking_starts_stops : List[Tuple[Time, Time]]
        List of bookings (tuples of booking_start / booking_stop).
    output_path : str, optional
        Path where the booking file will be saved, by default "".
    file_prefix : str, optional
        Prefix of the name so that the resulting file is "<file_prefix><kp_code>_booking.csv", by default "".
    """
    file_name = os.path.join(output_path, f"{file_prefix}{KP_CODE.lower()}_booking.csv")
    with open(file_name, "w") as wfile:
        for booking in booking_starts_stops:
            # Only write time at the second precision
            booking[0].precision = 0
            booking[1].precision = 0
            wfile.write(f"{booking[0].iso},{booking[1].iso},{KP_CODE},nenupy Weekly Schedule\n")
    log.info(f"{file_name} written.")


# ============================================================= #
# --------------------- _get_source_dict ---------------------- #
def _get_source_dict(file_name: str) -> dict:
    """Generic function to be called by get_exoplanet_dict() and get_star_dict().

    Parameters
    ----------
    file_name : str
        Name of the json file to read.

    Returns
    -------
    dict
        The json file loaded into a dictionnary
    """

    this_directory = os.path.dirname(__file__)
    source_list_json = os.path.join(this_directory, file_name)

    with open(source_list_json, "r") as rfile:
        src_dict = json.load(rfile)
    
    return src_dict


# ============================================================= #
# -------------------- get_exoplanet_dict --------------------- #
def get_exoplanet_dict() -> dict:
    """Return the dictionnary of all exoplanets listed in exoplanets.json

    Returns
    -------
    dict
        Dict of exoplanets
    """
    return _get_source_dict("exoplanets.json")


# ============================================================= #
# ----------------------- get_star_dict ----------------------- #
def get_star_dict() -> dict:
    """Return the dictionnary of all stars listed in stars.json

    Returns
    -------
    dict
        Dict of stars
    """
    return _get_source_dict("stars.json")


# ============================================================= #
# --------------------- obsid_to_timerange -------------------- #
def obsid_to_timerange(observation_repository: str) -> Tuple[Time, Time]:
    """Convert a NenuFAR observation repository name to the start / stop times of this observation.
    It basically consists of parsing the observation OBSID YYYYMMDD_hhmmss_YYYYMMDD_hhmmss which is the start of the repo name.

    Parameters
    ----------
    observation_repository : str
        Name of the observation repository

    Returns
    -------
    Tuple[Time, Time]
        Tuple of (observation_start, observation_stop)
    """
    obs_name = os.path.basename(observation_repository)
    obs_id = obs_name[:31] # YYYYMMDD_hhmmss_YYYYMMDD_hhmmss

    start_day, start_hour, stop_day, stop_hour = obs_id.split("_")
    start_day_iso = f"{start_day[:4]}-{start_day[4:6]}-{start_day[6:8]}"
    start_hour_iso = ":".join([start_hour[i * 2: i * 2 + 2] for i in range(3)])
    stop_day_iso = f"{stop_day[:4]}-{stop_day[4:6]}-{stop_day[6:8]}"
    stop_hour_iso = ":".join([stop_hour[i * 2: i * 2 + 2] for i in range(3)])

    return Time(f"{start_day_iso}T{start_hour_iso}"), Time(f"{stop_day_iso}T{stop_hour_iso}")


# ============================================================= #
# ------------------ find_source_directories ------------------ #
def _find_via_ssh(data_path: str) -> List[str]:
    log.info(f"Searching for observations in '{data_path}' via SSH to {REMOTE_SERVER}...")
    commands = [
        f"ssh {REMOTE_SERVER} 'find {data_path} -mindepth 3 -maxdepth 3 -type d'"
    ]
    proc = subprocess.Popen("; ".join(commands), shell=True, stdout=subprocess.PIPE, executable="/bin/bash")
    stdout, stderr = proc.communicate()
    exitcode = proc.returncode
    observation_directories =  stdout.decode("utf-8").split("\n")[:-1]
    return observation_directories

def _find_via_glob(data_path: str) -> List[str]:
    log.info(f"Searching for observations in '{data_path}' locally...")
    directories = []
    year_paths = glob.glob(os.path.join(data_path, "*"))
    for year_path in sorted(year_paths):
        month_paths = glob.glob(os.path.join(year_path, "*"))
        for month_path in sorted(month_paths):
            observations = glob.glob(os.path.join(month_path, "*"))
            for observation in sorted(observations):
                directories.append(observation)
    return directories

def find_directories(*data_paths: str) -> List[str]:
    """Limit datapath to /databf/<instr>/<kp>. Then the function will look for obs dir stored after <year>/<month>/*"""
    observations_by_path = []
    for path in data_paths:
        if os.path.isdir(path):
            # The data are local
            observations_by_path.append(
                _find_via_glob(data_path=path)
            )
        else:
            # Remote connexion using REMOTE_SERVER
            observations_by_path.append(
                _find_via_ssh(data_path=path)
            )
    return sum(observations_by_path, [])


# ============================================================= #
# ------------------ get_source_directories ------------------- #
def get_source_directories(source_name: str, observation_directories: List[str]) -> List[str]:
    """Filter the observation direcotries by source name.
    This uses the SPECIAL_CASES variable where several alternative source name may be defined.

    Parameters
    ----------
    source_name : str
        The name of the source as it should appear in the JSON file
    observation_directories : List[str]
        List of observation directories such as returned by find_directories()

    Returns
    -------
    List[str]
        List of source directories
    """
    source_name = source_name.upper()

    # Filter out directories containing the source name
    directories = []
    for obs_dir in observation_directories:
        src_known_names = SPECIAL_CASES.get(source_name, {}).get("include", [source_name])
        not_the_src = SPECIAL_CASES.get(source_name, {}).get("exclude", [])
        if any(name in obs_dir for name in src_known_names) and not any(other in obs_dir for other in not_the_src):
            directories.append(obs_dir)
    log.debug(f"{len(directories)} observations found for source '{source_name}'.")

    return directories

# ============================================================= #
# ----------------- get_source_exposure_time ------------------ #
def get_source_exposure_time(source_name: str, observation_directories: List[str]) -> TimeDelta:
    """Get the total amount of exposure time spent on source_name.
    This will look into all repo within data_path and search for any directory containing source_name.

    Parameters
    ----------
    source_name : str
        The name of the source as it should appear in the JSON file
    observation_directories : List[str]
        List of observation directories such as returned by find_directories()

    Returns
    -------
    TimeDelta
        The exposure time on source_name
    """
    directories = get_source_directories(
        source_name=source_name,
        observation_directories=observation_directories
    )

    # Computing the cumulative exposure time
    exposure_time = TimeDelta(0, format="sec")
    for observation in directories:
        start, stop = obsid_to_timerange(observation)
        exposure_time += stop - start

    return exposure_time


# ============================================================= #
# ---------------- get_source_last_observation ---------------- #
def get_source_last_observation(source_name: str, observation_directories: List[str]) -> Time:
    directories = get_source_directories(
        source_name=source_name,
        observation_directories=observation_directories
    )

    if len(directories) == 0:
        # The source has never been observed
        log.warning(f"{source_name} has never been observed.")
        default_time = Time("2000-01-01T00:00:00")
        default_time.precision = 0
        return default_time

    directories = sorted(directories)
    latest_directory = directories[-1]
    observation_time, _ = obsid_to_timerange(latest_directory)

    observation_time.precision = 0

    return observation_time


# ============================================================= #
# -------------------- _modify_json_files --------------------- #
def _modify_json_files(source_dict: dict) -> None:
    """Find out which JSON file corresponds with the source_dict and save it.

    Parameters
    ----------
    source_dict : dict
        Dictionnary of sources
    """
    first_source = list(source_dict.keys())[0]
    # Find out which of the json file this dict corresponds with
    this_directory = os.path.dirname(__file__)
    if first_source in list(get_star_dict().keys()):            
        json_file = os.path.join(this_directory, "stars.json")
    else:
        json_file = os.path.join(this_directory, "exoplanets.json")
    
    log.info(f"Updating '{json_file}'...")
    
    with open(json_file, "w") as wfile:
        json.dump(source_dict, wfile)


# ============================================================= #
# ------------------- update_exposure_time -------------------- #
def update_exposure_time(source_dict: dict, *data_paths: str, save: bool = False) -> dict:
    """Update the exposure time of every source in source_dict with observations found in datapath.
    For every source get_source_exposure_time() is called: the data repository names are parsed to figure out the exposure time from the obsid.

    Parameters
    ----------
    source_dict : dict
        The source dictionnary (such as returned by get_exoplanet_dict() for instance)
    data_paths : str
        The path(s) where all the data repositories are
    save : bool, optional
        If True the JSON files are updated with the new exposure times, by default False

    Returns
    -------
    dict
        The source dictionnary with the exposure time updated
    
    Example
    -------
    .. code-block:: python

        >>> from exoschedule.functions import get_exoplanet_dict, update_exposure_time
        >>> ep = get_exoplanet_dict()
        >>> ep = update_exposure_time(ep, "/databf/nenufar/LT02", "/databf/nenufar/ES02", save=True)
    """

    obs_directories = find_directories(*data_paths)

    for src in source_dict:
        dt = get_source_exposure_time(source_name=src, observation_directories=obs_directories)
        last_obs = get_source_last_observation(source_name=src, observation_directories=obs_directories)
        previous_exposure_time = source_dict[src]["exposure_time_hours"]
        current_exposure_time = dt.sec / 3600
        source_dict[src]["exposure_time_hours"] = current_exposure_time
        source_dict[src]["last_observation"] = last_obs.isot
        if current_exposure_time != previous_exposure_time:
            log.info(f"{src}: DT={(current_exposure_time - previous_exposure_time)}h added.")

    if save:
        _modify_json_files(source_dict)
    
    return source_dict


# ============================================================= #
# ----------------- update_minimal_elevation ------------------ #
def update_minimal_elevation(source_dict: dict, save: bool = False, minimal_elevation: u.Quantity = None) -> dict:
    """Update the minimal observing elevation in the source dictionnary.
    For every source, its visibility at Nançay is computed, namely its maximal and minimal altitude (a_max, and a_min).
    if a_max - 10 deg < 20 deg: minimal_elevation = 20 deg
    elif a_max - 10 deg > 60 deg: minimal_elevation = 60 deg
    else minimal_elevation = a_max - 10 deg

    Parameters
    ----------
    source_dict : dict
        The source dictionnary (such as returned by get_exoplanet_dict() for instance)
    save : bool, optional
        If True the JSON files are updated with the new exposure times, by default False
    minimal_elevation : astropy.units.Quantity, optional
        If a value is given, every source will be updated with this value, by default None

    Returns
    -------
    dict
        The source dictionnary with the minimal_elevation updated
    """

    for src in source_dict:

        if minimal_elevation is None:
            skycoord = SkyCoord(
                source_dict[src]["ra"],
                source_dict[src]["dec"],
                unit=(u.hourangle, u.deg)
            )
            target = FixedTarget(skycoord)
            max_elev = target.elevation_max.deg
            minimal_elev = max_elev - 10
            if minimal_elev < 20:
                minimal_elev = 20
            elif minimal_elev > 60:
                minimal_elev = 60
        else:
            minimal_elev = minimal_elevation

        source_dict[src]["minimal_elevation"] = minimal_elev

    if save:
        _modify_json_files(source_dict)
    
    return source_dict


# ============================================================= #
# -------------------- clear_exposure_time -------------------- #
def clear_exposure_time(source_dict: dict, save: bool = False) -> dict:
    """Set every source's exposure_time listed in source_dict to 0.

    Parameters
    ----------
    source_dict : dict
        Disctionnary of sources.
    save : bool, optional
        If True the JSON files are updated with the new exposure times, by default False

    Returns
    -------
    dict
        The source dictionnary with the exposure time updated
    """
    for src in source_dict:
        source_dict[src]["exposure_time_hours"] = 0

    if save:
        _modify_json_files(source_dict)

    return source_dict

# ============================================================= #
# -------------------- get_imaging_sources -------------------- #
def get_imaging_sources(all_sources: dict, imaging_sources: List[str] = IMAGING_SOURCES) -> dict:

    # Make sure the sources listed in IMAGING_SOURCES are known and present in all sources
    for source in imaging_sources:
        if source not in all_sources:
            raise ValueError(f"Source '{source}' is not in the known source list.")

    log.info(f"Imaging sources considered are: {imaging_sources}.")

    # Filter out the dictionnary
    return {source: all_sources[source] for source in imaging_sources}


# # ============================================================= #
# # -------------------- write_booking_file --------------------- #
# def cumulative_observation_time(nenufar_path: str = "/databf/nenufar/LT02/") -> TimeDelta:
#     n_observations = 0
#     cumulative_time = TimeDelta(0, format="sec")
#     year_paths = glob.glob(os.path.join(nenufar_path, "*"))
#     for year_path in year_paths:
#         print(f"Working on year {year_path}...")
#         month_paths = glob.glob(os.path.join(year_path, "*"))
#         for month_path in month_paths:
#             print(f"\tmonth {month_path}...")
#             observations = glob.glob(os.path.join(month_path, "*"))
#             n_observations += len(observations)
#             for observation in observations:
#                 obs_name = os.path.basename(observation)
#                 obs_id = obs_name[:31] # YYYYMMDD_hhmmss_YYYYMMDD_hhmmss
#                 start_day, start_hour, stop_day, stop_hour = obs_id.split("_")
#                 start_day_iso = f"{start_day[:4]}-{start_day[4:6]}-{start_day[6:8]}"
#                 start_hour_iso = ":".join([start_hour[i * 2: i * 2 + 2] for i in range(3)])
#                 stop_day_iso = f"{stop_day[:4]}-{stop_day[4:6]}-{stop_day[6:8]}"
#                 stop_hour_iso = ":".join([stop_hour[i * 2: i * 2 + 2] for i in range(3)])
#                 cumulative_time += Time(f"{stop_day_iso}T{stop_hour_iso}") -  Time(f"{start_day_iso}T{start_hour_iso}")
#     print(f"Number of observartions in {nenufar_path}: {n_observations}.")
#     print(f"Observation time: {cumulative_time.sec / 3600} hours.")
#     return cumulative_time