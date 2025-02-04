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
    "write_booking_file",
    "get_exoplanet_dict",
    "get_star_dict",
    "obsid_to_timerange",
    "get_source_exposure_time",
    "update_exposure_time"
]


from astropy.time import Time, TimeDelta
from typing import List, Tuple
import numpy as np
import operator
import os
import glob
import json
import logging
log = logging.getLogger("exoschedule")

from nenupy.schedule import Schedule, ReservedBlock


# ============================================================= #
# ---------------------- argument_parser ---------------------- #
# Minimum booking slot is 1h
MIN_DURATION = TimeDelta(3600, format="sec")

# No booking between 11-15 UTC except Saturday and Sunday
# to leave room for maintenance operations (and because the Sun
# is quite high and ppolluting).
MID_DAY_START_HOUR = 11
MID_DAY_STOP_HOUR = 15

KP_CODE = "LT02"
# ============================================================= #
# ============================================================= #

# ============================================================= #
# --------------------- bookings_from_vcr --------------------- #
def bookings_from_vcr(
        vcr_current_booking: str,
        start_time: Time,
        stop_time: Time,
        merge_slots: bool = True,
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
    min_dt : TimeDelta, optional
        Minimal time resolution of the booking time definition, by default TimeDelta(60, format="sec")

    Returns
    -------
    List[Tuple[Time, Time]]
        List of bookings in the shape [(start_0, stop_0), (start_1, stop_1), ...].
    """
    
    # Get the reserved slots currently scheduled within the VCR
    booked_slots = ReservedBlock.from_VCR(vcr_current_booking)
    if merge_slots:
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

        if (group_duration <= MIN_DURATION):
            log.info(f"Booking < min_duration ({MIN_DURATION.sec}s): {group_start.isot} - {group_stop.isot} --> skipped.")

        free_bookings.append( (group_start, group_stop) )

    return free_bookings


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
    file_name = os.path.join(output_path, f"{KP_CODE.lower()}_booking.csv")
    with open(file_name, "w") as wfile:
        for booking in booking_starts_stops:
            wfile.write(f"{booking[0]},{booking[1]},{KP_CODE},Nenupy Weekly Schedule\n")
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
def find_source_directories(source_name: str, data_path: str = "/databf/nenufar/LT02/") -> List[str]:
    directories = []
    year_paths = glob.glob(os.path.join(data_path, "*"))
    for year_path in sorted(year_paths):
        month_paths = glob.glob(os.path.join(year_path, "*"))
        for month_path in sorted(month_paths):
            observations = glob.glob(os.path.join(month_path, "*"))
            for observation in sorted(observations):
                if source_name.upper() not in observation:
                    continue
                directories.append(observation)
    return directories

# ============================================================= #
# ----------------- get_source_exposure_time ------------------ #
def get_source_exposure_time(source_name: str, data_path: str = "/databf/nenufar/LT02/") -> TimeDelta:
    """Get the total amount of exposure time spent on source_name.
    This will look into all repo within data_path and search for any directory containing source_name.

    Parameters
    ----------
    source_name : str
        The name of the source as it should appear in the observaiton repositories
    data_path : str, optional
        The path where all observations are put together, it is better to search the nenufar repo (containing only statistical data) as it lists everything, by default "/databf/nenufar/LT02/"

    Returns
    -------
    TimeDelta
        The exposure time on source_name
    """
    exposure_time = TimeDelta(0, format="sec")
    log.debug(f"Looking out for data repositories involving '{source_name}'...")
    observations = find_source_directories(source_name=source_name, data_path=data_path)
    for observation in observations:
        start, stop = obsid_to_timerange(observation)
        exposure_time += stop - start
    return exposure_time


# ============================================================= #
# ------------------- update_exposure_time -------------------- #
def update_exposure_time(source_dict: dict, data_path: str = "/databf/nenufar/LT02/", save: bool = False) -> dict:
    """Update the exposure time of every source in source_dict with observations found in datapath.
    For every source get_source_exposure_time() is called: the data repository names are parsed to figure out the exposure time from the obsid.

    Parameters
    ----------
    source_dict : dict
        The source dictionnary (such as returned by get_exoplanet_dict() for instance)
    data_path : str
        The path where all the data repositories are, by default "/databf/nenufar/LT02/"
    save : bool, optional
        If True the JSON files are updated with the new exposure times, by default False

    Returns
    -------
    dict
        The source dictionnary with the exposure time updated
    """

    for src in source_dict:
        dt = get_source_exposure_time(src, data_path=data_path)
        source_dict["exposure_time_hours"] = dt.sec / 3600

    if save:
        # Find out which of the json file this dict corresponds with
        this_directory = os.path.dirname(__file__)
        if src in list(get_star_dict().keys()):            
            json_file = os.path.join(this_directory, "stars.json")
        else:
            json_file = os.path.join(this_directory, "exoplanets.json")
        
        log.info(f"Updating '{json_file}'...")
        
        with open(json_file, "w") as wfile:
            json.dump(source_dict, wfile)
    
    return source_dict

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