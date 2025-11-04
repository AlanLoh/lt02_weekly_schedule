#! /usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Alan Loh"
__copyright__ = "Copyright 2025, exoschedule"
__credits__ = ["Alan Loh"]
__maintainer__ = "Alan"
__email__ = "alan.loh@obspm.fr"
__status__ = "Production"
__all__ = [
    "build_observation_blocks",
    "build_imaging_observation_blocks",
    "initialize_schedule",
    "sort_blocks_by_priorities",
    "book_observations",
    "plot_schedule",
    "is_night",
    "write_parset_file"
]


from nenupy.schedule import (
    Schedule,
    ObsBlock,
    ESTarget,
    Constraints,
    ElevationCnst,
    NightTimeCnst
)
from nenupy.observation import ParsetUser
from nenupy.astro.target import FixedTarget
from nenupy.schedule.open_time import _sort_night_time

from exoschedule import (
    SCHEDULE_DT,
    TYPE_PRIORITY_ORDER,
    TYPE_FACTOR,
    EXPOSURE_FACTOR,
    LAST_OBSERVE_FACTOR,
    DISABLE_POINTING_CORRECTION_DECLINATION,
    OBSERVATION_CONFIGURATION,
    CONFIG_EXCEPTION,
    MIN_OBS_DURATION,
    MIN_DURATION,
    A_TEAM_SOURCES,
    CALIBRATION_DURATION
)

from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, Angle
from astropy.table import Table
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patheffects as PathEffects
import numpy as np
import os
import copy
from typing import Tuple, List
import logging
log = logging.getLogger("exoschedule")


# ============================================================= #
# ----------------- build_observation_blocks ------------------ #
def build_observation_blocks(source_dict: dict, min_duration: TimeDelta = MIN_OBS_DURATION) -> ObsBlock:
    src_obs = sum([
        ObsBlock(
            name=src,
            program="LT02",
            target=ESTarget(SkyCoord(source_dict[src]["ra"], source_dict[src]["dec"], unit=(u.hourangle, u.deg))),
            # duration=TimeDelta(source_dict[src]["minimal_obs_time_minutes"] * 60, format="sec"),
            duration=min_duration,
            max_extended_duration=TimeDelta(source_dict[src]["nominal_obs_time_minutes"] * 60, format="sec"),
            constraints=Constraints(
                ElevationCnst(elevationMin=source_dict[src]["minimal_elevation"] * u.deg, scale_elevation=True, weight=1),
                # MeridianTransitCnst(weight=1)
            )
        ) for src in source_dict if source_dict[src]["main_target"] and (source_dict[src]["todo"] in ["observation"])
    ])
    return src_obs

def build_imaging_observation_blocks(source_dict: dict, source_list: List[str], min_duration: TimeDelta = MIN_OBS_DURATION) -> ObsBlock:
    source_names = list(source_dict.keys())
    for source in source_list:
        if source not in source_names:
            raise ValueError(f"Source {source} not known. Check that it could be found in exoplanet or star JSON files.")
    src_obs = sum([
        ObsBlock(
            name=src,
            program="LT02",
            target=ESTarget(SkyCoord(source_dict[src]["ra"], source_dict[src]["dec"], unit=(u.hourangle, u.deg))),
            duration=min_duration,
            max_extended_duration=TimeDelta(source_dict[src]["nominal_obs_time_minutes"] * 60, format="sec"),
            constraints=Constraints(
                ElevationCnst(elevationMin=source_dict[src]["minimal_elevation"] * u.deg, scale_elevation=True, weight=1),
                NightTimeCnst(scale_elevation=False)
            )
        ) for src in source_list if source_dict[src]["main_target"] and (source_dict[src]["todo"] in ["observation"])
    ])
    return src_obs

# ============================================================= #
# -------------------- initialize_schedule -------------------- #
def initialize_schedule(booking_starts_stops: List[Tuple[Time, Time]], dt: TimeDelta = SCHEDULE_DT) -> Schedule:
    start = booking_starts_stops[0][0]
    stop = booking_starts_stops[-1][1]
    
    log.info(f"Creating a Schedule from {start.isot} to {stop.isot}...")
    
    schedule =  Schedule(
        time_min=start,
        time_max=stop,
        dt=dt
    )

    block_starts = Time([blk[0].isot for blk in booking_starts_stops])
    block_stops = Time([blk[1].isot for blk in booking_starts_stops])

    schedule.set_free_slots(
        start_times=block_starts,
        stop_times=block_stops
    )

    return schedule

# ============================================================= #
# ----------------- sort_blocks_by_priorities ----------------- #
def sort_blocks_by_priorities_v0(source_dict: dict, schedule: Schedule) -> np.array:
    """Given a schedule containing observation_blocks, sort them out in order of booking priority.
    The priority is firt given to the source type as defined in TYPE_PRIORITY_ORDER.
    Then, the exposure time already spent on each source is looked for in source_dict.
    The sources are further sorted by increasing exposure time, hence prioritizing observations of lesser exposed targets.

    Parameters
    ----------
    source_dict : dict
        Source dictionnary (may be a concatenation of exoplanets and stars: {**stars, **exoplanets})
    schedule : Schedule
        The Schedule with the observation blocks already inserted.

    Returns
    -------
    np.array
        The sorted indices
    """

    source_names = np.array([blk.name for blk in schedule.observation_blocks])
    source_types = np.array([source_dict[name]["source_type"] for name in source_names if source_dict[name]["main_target"]])
    exposure = np.array([source_dict[name]["exposure_time_hours"] for name in source_names if source_dict[name]["main_target"]])
    
    # Sort by source type    
    # Create a list of index lists by source type, in order of TYPE_PRIORITY_ORDER
    list_idx_by_type = [
        np.argwhere(source_types == TYPE_PRIORITY_ORDER[i])[:, 0]
        for i in range(len(TYPE_PRIORITY_ORDER))
    ]
    # For each of those list, order by exposure time
    for i, type_idx in enumerate(list_idx_by_type):
        type_exposures = exposure[type_idx]
        exposure_sort = np.argsort(type_exposures)
        list_idx_by_type[i] = type_idx[exposure_sort]

    # Return the concatenated 1D index array
    return np.concatenate(list_idx_by_type)


def sort_blocks_by_priorities(source_dict: dict, schedule: Schedule, priority_sources: List[str] = []) -> np.array:
    log.info(f"Sorting observation blocks using type_factor={TYPE_FACTOR}, exposure_factor={EXPOSURE_FACTOR}, last_observation_factor={LAST_OBSERVE_FACTOR}...")

    source_names = np.array([blk.name for blk in schedule.observation_blocks])

    source_types = np.array([source_dict[name]["source_type"] for name in source_names if source_dict[name]["main_target"]])
    source_last_observation = Time([source_dict[name]["last_observation"] for name in source_names if source_dict[name]["main_target"]], format="isot")
    exposure = np.array([source_dict[name]["exposure_time_hours"] for name in source_names if source_dict[name]["main_target"]])

    type_index = np.array([TYPE_PRIORITY_ORDER.index(src_type) for src_type in source_types], dtype=float)
    type_index = 1 - type_index / type_index.max() # 1=exoplanet-l

    lastobs_index = Time.now() - source_last_observation
    lastobs_index /= lastobs_index.max() # 1=longest since last observed

    exposure_index = 1 - exposure / exposure.max() # 1=0h exposure, 0=max hour exposure

    # If a source is listed as a priority target, artificially increase its score
    for priority_source in priority_sources:
        if priority_source not in source_names:
            log.warning(f"Source {priority_source} not known. Check that it could be found in exoplanet or star JSON files.")
            continue
        source_index = np.argwhere(source_names == priority_source)
        type_index[source_index] += 1000
        exposure_index[source_index] += 1000
        lastobs_index[source_index] += 1000
        log.info(f"Observation of {priority_source} defined as a priority.")

    sorted_indices = np.argsort(TYPE_FACTOR * type_index + EXPOSURE_FACTOR * exposure_index + LAST_OBSERVE_FACTOR * lastobs_index)

    return sorted_indices[::-1] # return the lowest sorted


# ============================================================= #
# --------------------- book_observations --------------------- #
def book_observations(booking_starts_stops: List[Tuple[Time, Time]], source_dict: dict, observation_blocks: ObsBlock, priority_sources: List[str] = [], release_duration: bool = True) -> Schedule:

    log.info("Booking the schedule...")

    schedule = initialize_schedule(booking_starts_stops)
    schedule.insert(observation_blocks)
    
    sort_indices = sort_blocks_by_priorities(
        source_dict=source_dict,
        schedule=schedule,
        priority_sources=priority_sources
    )

    schedule.book(
        very_strict=False,
        reset_booking=True,
        sort_by_indices=sort_indices,
        minimal_score=0.9
    )

    # At the second pass, allow for less duration observations
    if release_duration:
        for obs_block in schedule.observation_blocks:
            if not obs_block.isBooked:
                obs_block.duration = MIN_DURATION

    schedule.book(
        very_strict=False,
        reset_booking=False,
        sort_by_availability=True,
        minimal_score=0.8
    )

    schedule.extend_scheduled_observations(
        relative_score_threshold=0.5,
    )

    results = schedule.export()
    log.info(f"Schedule proposed:\n{results}")

    return schedule


# ============================================================= #
# ----------------------- plot_schedule ----------------------- #
def plot_schedule(schedule: Schedule, source_dict: dict, output_path: str, file_prefix: str = "") -> None:
    """Plot the schedule

    Parameters
    ----------
    schedule : Schedule
        Schedule with booked slots
    output_path : str
        Path where the png will be stored
    file_prefix : str, optional
        Prefix of the name so that the resulting file is "<file_prefix>schedule_plot.png", by default "".
    """

    result = schedule.export()

    plt.close("all")

    dt_day = TimeDelta(1, format="jd")
    schedule_start = Time(f"{schedule.time_min.isot.split('T')[0]}T00:00:00")
    schedule_stop = schedule.time_max if schedule.time_max.isot.split("T")[1].startswith("00:00:00") else Time(f"{(schedule.time_max + dt_day).isot.split('T')[0]}T00:00:00")
    n_days = int(np.ceil((schedule_stop - schedule_start) / dt_day))
    days = np.array([(schedule.time_min + i * dt_day).isot.split("T")[0] for i in range(n_days)])

    fig, axs = plt.subplots(
        nrows=1,
        ncols=n_days,
        figsize=(2 * n_days, 15),
        sharey=True
    )
    fig.subplots_adjust(wspace=0, top=0.95)
    if n_days == 1:
        axs = np.array([axs])

    text = []

    # Show the source types
    type_cmap = matplotlib.colormaps["Dark2"]
    type_colors = type_cmap(np.linspace(0, 1, len(TYPE_PRIORITY_ORDER)))
    type_col_dict = {}
    for i, (src_type, type_color) in enumerate(zip(TYPE_PRIORITY_ORDER, type_colors)):
        type_col_dict[src_type] = type_color
        text.append(
            fig.text(0.07 + 0.09 * (i + 1), 0.96, src_type.upper(), ha="center", va="center", color=type_color)
        )

    cmap = matplotlib.cm.get_cmap("RdYlGn")

    for day_id, day in enumerate(days):

        day_start = Time(day)
        day_end = day_start + TimeDelta(1, format="jd")

        time_mask = ((result["start"] < day_start) * (result["stop"] > day_end)) + \
                    ((result["start"] <= day_start) * (result["stop"] > day_start)) +\
                    ((result["start"] > day_start) * (result["stop"] < day_end)) +\
                    ((result["start"] < day_end) * (result["stop"] >= day_end))
        
        for block in result[time_mask]:
            start = block["start"]
            if start < day_start:
                start = day_start
            start.precision = 0
            stop = block["stop"]
            if stop >= day_end:
                stop = day_end
            stop.precision = 0

            ymin = start.datetime.hour + start.datetime.minute / 60
            ymax = stop.datetime.hour + stop.datetime.minute / 60
            if ymax == 0:
                ymax = 24
            xmin = 0 # day_id / slots
            xmax = 1 # (day_id + 1) / slots
            axs[day_id].axhspan(
                ymin, ymax,
                xmin, xmax,
                edgecolor="black",
                facecolor="white"#cmap(block["score"])
            )
            axs[day_id].axhspan(
                ymin, ymax,
                0.9, xmax,
                edgecolor="black",
                facecolor=cmap(block["score"])
            )
            axs[day_id].axhline(ymin, color="black", linewidth=1.5)
            axs[day_id].axhline(ymax, color="black", linewidth=1.5)

            
            pad = 0.1
            small_font = 6
            large_font = 8
            text_color = "black"

            # print source name
            text.append(
                axs[day_id].text(
                    (xmax + xmin) / 2,
                    (ymax + ymin) / 2,
                    "\n".join([block["name"][i: i + 12] for i in range(0, len(block["name"]), 12)] ),
                    va="center",
                    ha="center",
                    fontsize=large_font,
                    color=type_col_dict[source_dict[block["name"]]["source_type"]]#text_color
                )
            )

            # Print elevation min / max
            skycoord = SkyCoord(
                source_dict[block["name"]]["ra"],
                source_dict[block["name"]]["dec"],
                unit=(u.hourangle, u.deg)
            )
            target = FixedTarget(skycoord, time=Time([block["start"], block["start"] + (block["stop"] - block["start"]) / 2, block["stop"]]))
            text.append(
                axs[day_id].text(
                    xmin + 0.05,
                    ymin + pad,
                    f"{target.horizontal_coordinates[0].alt.deg:.0f}",
                    va="top",
                    ha="left",
                    rotation="vertical",
                    fontsize=small_font,
                    color=text_color
                )
            )
            text.append(
                axs[day_id].text(
                    xmin + 0.05,
                    (ymax + ymin) / 2,
                    f"{target.horizontal_coordinates[1].alt.deg:.0f}",
                    va="center",
                    ha="left",
                    rotation="vertical",
                    fontsize=small_font,
                    color=text_color
                )
            )
            text.append(
                axs[day_id].text(
                    xmin + 0.05,
                    ymax - pad,
                    f"{target.horizontal_coordinates[2].alt.deg:.0f}",
                    va="bottom",
                    ha="left",
                    rotation="vertical",
                    fontsize=small_font,
                    color=text_color
                )
            )

            # Print time min / max
            text.append(
                axs[day_id].text(
                    (xmax + xmin) / 2,
                    ymin + pad,
                    start.isot.split("T")[1],
                    va="top",
                    ha="center",
                    fontsize=small_font,
                    color=text_color
                )
            )
            text.append(
                axs[day_id].text(
                    (xmax + xmin) / 2,
                    ymax - pad,
                    stop.isot.split("T")[1],
                    va="bottom",
                    ha="center",
                    fontsize=small_font,
                    color=text_color
                )
            )
        
        # for t in text:
        #     t.set_path_effects([PathEffects.withStroke(linewidth=1.5, foreground="black")])

        axs[day_id].set_xticklabels([])
        axs[day_id].set_xlabel(day, rotation="horizontal", va="center", ha="center")

        axs[day_id].set_ylim(24, 0)
        axs[day_id].set_yticks(np.arange(0, 24, 4))
        axs[day_id].set_yticks(np.arange(0, 24, 1), minor=True)
        axs[day_id].tick_params(axis="y", which="both", direction="inout", length=5)
        axs[day_id].grid(axis="y", which="major", color="gray", alpha=0.8, linestyle="--")
        axs[day_id].grid(axis="y", which="minor", color="gray", alpha=0.5, linestyle=":")

    axs[0].set_ylabel("UTC hour")

    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=axs.ravel().tolist(), ticks=np.linspace(0, 1, 11), boundaries=np.arange(0, 1, .01), pad=0.03, aspect=50)
    cb.set_label("Score")

    if output_path != "":
        figname = os.path.join(output_path, f"{file_prefix}schedule_plot.png")
        fig.savefig(figname, dpi=300, bbox_inches="tight")


# ============================================================= #
# ------------------------- is_night -------------------------- #
def is_night(start: Time, stop: Time) -> bool:
    """Determine if an observation is by night or not.

    Parameters
    ----------
    start : Time
        Start time of the observation
    stop : Time
        Stop time of the observation

    Returns
    -------
    bool
        If the obsrvation is considered night-time or not
    """
    start_angle = Angle(start.isot.split('T')[1], unit="hourangle")
    stop_angle = Angle(stop.isot.split('T')[1], unit="hourangle")
    if stop_angle < start_angle:
        day1, night1 = _sort_night_time(
            start_hour=start_angle,
            stop_hour=Angle("24", unit="hourangle"),
            current_day=start
        )
        day2, night2 = _sort_night_time(
            start_hour=Angle("00", unit="hourangle"),
            stop_hour=stop_angle,
            current_day=start
        )
        day = day1 + day2
        night = night1 + night2
    else:
        day, night = _sort_night_time(
            start_hour=start_angle,
            stop_hour=stop_angle,
            current_day=start
        )

    if day > night:
        return False
    else:
        return True

# ============================================================= #
# ------------------------ find_a_team ------------------------ #
def find_a_team(time: Time) -> str:
    """Return the A team source name that is at the highest elevation at time.

    Parameters
    ----------
    time : Time
        Time of observation.

    Returns
    -------
    str
        Name of the A team source
    """
    id_max_elevation = np.argmax([FixedTarget.from_name(a_team, time=time).horizontal_coordinates.alt.deg for a_team in A_TEAM_SOURCES])
    return A_TEAM_SOURCES[id_max_elevation]

# ============================================================= #
# ----------------------- find_transit ------------------------ #
def find_transit(target: FixedTarget, time: Time) -> Time:
    transits = [
        target.previous_meridian_transit(time=time).isot,
        target.next_meridian_transit(time=time).isot
    ]
    transit_idx = np.argmin(np.abs(Time(transits) - time))
    transit = Time(transits[transit_idx])
    transit.precision = 0
    return transit

# ============================================================= #
# --------------------- write_parset_file --------------------- #
def write_parset_file(schedule_row: Table, source_dict: dict, output_path: str, file_prefix: str = "", add_imaging: bool = False) -> None:

    source_name = schedule_row["name"]
    time_min = schedule_row["start"]
    time_max = schedule_row["stop"]

    skycoord = SkyCoord(source_dict[source_name]["ra"], source_dict[source_name]["dec"], unit=(u.hourangle, u.deg))

    time_min.precision = 0
    time_max.precision = 0
    init_delta = TimeDelta(70, format="sec")
    end_dt = TimeDelta(60, format="sec")

    # Find the transit source time
    source = FixedTarget(skycoord)
    transit = find_transit(source, time_min)

    # Consider time to make a calibration observation in case of imaging
    if add_imaging:
        cal_src_name = find_a_team(time_min)
        # Check whether to observe it before or after
        cal_src = FixedTarget.from_name(cal_src_name, time=Time([time_min.isot, time_max.isot], format="isot"))
        alt_max_idx = np.argmax(cal_src.horizontal_coordinates.alt.deg)
        if alt_max_idx == 0:
            # Observation before
            cal_time_min = time_min
            time_min += CALIBRATION_DURATION + TimeDelta(60 * 2, format="sec")
            time_min.precision = 0
            cal_time_max = time_min
        else:
            # Observation after
            cal_time_max = time_max
            time_max -= CALIBRATION_DURATION + TimeDelta(60 * 2, format="sec")
            time_max.precision = 0
            cal_time_min = time_max
        cal_transit = find_transit(cal_src, cal_time_min)

    parset = ParsetUser()
    obs = parset.observation
    obs["contactName"] = "philippe.zarka"
    main_target_name = f"{source_name.split()[0].upper()}_TRACKING"
    obs["name"] = main_target_name
    if "exoplanet" in source_dict[source_name]["source_type"]:
        source_type = "EXOPLANET"
    elif "star" in source_dict[source_name]["source_type"]:
        source_type = "STAR"
    elif "both" in source_dict[source_name]["source_type"]:
        source_type = "STAR_PLANETS" # for ceti group
    else:
        raise Exception(f"Unkown source type {source_dict['source_type']}!")
    obs["title"] = source_type
    obs["contactEmail"] = "philippe.zarka@obspm.fr"
    obs["topic"] = "lt02_exoplanets_and_stars"

    output = parset.output
    output["hd_receivers"] = "[undysputed]"

    # The configuration is different with respect to night and day
    day_period = "night" if is_night(time_min, time_max) else "day"

    # No azel correction if dec > 60
    corazel = "default" if skycoord.dec.deg < DISABLE_POINTING_CORRECTION_DECLINATION else "disabled"

    # Observation configuration
    if source_type.lower() == "star_planets":
        source_type = "exoplanet" # Ceti group... it doesnt matter to change it after this point
    config = copy.copy(OBSERVATION_CONFIGURATION[source_type.lower()][day_period])
    # Modify the config by the exceptions
    try:
        config.update(CONFIG_EXCEPTION[source_name.upper()][day_period])
    except KeyError:
        pass  

    parset.add_analog_beam(
        target=main_target_name,
        trackingType="tracking",
        directionType="j2000",
        ra=f"'{skycoord.ra.to_string(unit='hour', sep=':')}'",
        dec=f"'{skycoord.dec.to_string(unit='degree', sep=':')}'",
        startTime=(time_min + init_delta).isot + "Z",
        transitDate=transit.isot + "Z",
        filterStart=config["filter"],
        optFrq=config["squint"],
        corAzel=corazel,
        duration=f"{int(np.round((time_max - time_min - init_delta - end_dt).sec))}s"
    )

    # If CETI group observation
    if source_dict[source_name]["secondary_targets"] is not None:
        for secondary_name in source_dict[source_name]["secondary_targets"]:
            secondary_source = source_dict[secondary_name]
            skycoord_scd = SkyCoord(secondary_source["ra"], secondary_source["dec"], unit=(u.hourangle, u.deg))
            source = FixedTarget(skycoord_scd)
            transits = [
                source.previous_meridian_transit(time=time_min).isot,
                source.next_meridian_transit(time=time_min).isot
            ]
            transit_idx = np.argmin(np.abs(Time(transits) - time_min))
            transit = Time(transits[transit_idx])
            transit.precision = 0
            parset.add_numerical_beam(
                anabeam_index=0,
                target=f"{secondary_name.upper()}_TRACKING",
                trackingType="tracking",
                directionType="j2000",
                ra=f"'{secondary_source['ra']}'",
                dec=f"'{secondary_source['dec']}'",
                startTime=(time_min + init_delta).isot + "Z",
                transitDate=transit.isot + "Z",
                duration=f"{int(np.round((time_max - time_min - init_delta - end_dt).sec))}s",
                toDo="dynamicspectrum",
                subbandList=f"[{config['frequency'][0]}..{config['frequency'][1]}]",
                parameters=f"{config['parameters']}",
                corAzel=corazel
            )

    else:
        # Main
        parset.add_numerical_beam(
            anabeam_index=0,
            target=main_target_name,
            trackingType="tracking",
            directionType="j2000",
            ra=f"'{skycoord.ra.to_string(unit='hour', sep=':')}'",
            dec=f"'{skycoord.dec.to_string(unit='degree', sep=':')}'",
            startTime=(time_min + init_delta).isot + "Z",
            transitDate=transit.isot + "Z",
            duration=f"{int(np.round((time_max - time_min - init_delta - end_dt).sec))}s",
            toDo="dynamicspectrum",
            subbandList=f"[{config['frequency'][0]}..{config['frequency'][1]}]",
            parameters=f"{config['parameters']}",
            corAzel=corazel
        )

        # Add off beams
        off_beams = source_dict[source_name]["off_beams"]["option_1"]
        for off_beam in off_beams:
            skycoord_off = SkyCoord(off_beam["ra"], off_beam["dec"], unit=(u.hourangle, u.deg))
            source = FixedTarget(skycoord_off)
            transits = [
                source.previous_meridian_transit(time=time_min).isot,
                source.next_meridian_transit(time=time_min).isot
            ]
            transit_idx = np.argmin(np.abs(Time(transits) - time_min))
            transit = Time(transits[transit_idx])
            transit.precision = 0
            parset.add_numerical_beam(
                anabeam_index=0,
                target="J2000_TRACKING",
                trackingType="tracking",
                directionType="j2000",
                ra=f"'{off_beam['ra']}'",
                dec=f"'{off_beam['dec']}'",
                startTime=(time_min + init_delta).isot + "Z",
                transitDate=transit.isot + "Z",
                duration=f"{int(np.round((time_max - time_min - init_delta - end_dt).sec))}s",
                toDo="dynamicspectrum",
                subbandList=f"[{config['frequency'][0]}..{config['frequency'][1]}]",
                parameters=f"{config['parameters']}",
                corAzel=corazel
            )

    if add_imaging:
        imaging_freq = "[30.8-78.3]"
        imaging_params = "env_file=env_dp3-5.2.sh avg_timestep=1 avg_freqstep=30 startchan=2 nchan=60 compress=true flag_strategy=nenufar64c1s.lua sws=[158-221,222-281,282-341,342-401] stat_pols=[snr_xx,snr_yy,rfipercentage_xx]"
        cal_target = cal_src_name.replace(' ', "_").upper()

        parset.add_phase_center(
            target=main_target_name,
            subbandFrq=imaging_freq,
            useParentPointing=True,
            parameters=imaging_params
        )
        output["nri_receivers"] = True

        # Create the calibration parset
        cal_parset = ParsetUser()
        cal_obs = cal_parset.observation
        cal_obs["name"] = f"{cal_target}_TRACKING"
        cal_obs["contactName"] = "philippe.zarka"
        cal_obs["title"] = source_type
        cal_obs["contactEmail"] = "philippe.zarka@obspm.fr"
        cal_obs["topic"] = "lt02_exoplanets_and_stars"
        cal_output = cal_parset.output
        cal_output["nri_receivers"] = True
        cal_parset.add_analog_beam(
            target=f"{cal_target}_TRACKING",
            trackingType="tracking",
            directionType="j2000",
            ra=f"'{cal_src.coordinates.ra.to_string(unit='hour', sep=':')}'",
            dec=f"'{cal_src.coordinates.dec.to_string(unit='degree', sep=':')}'",
            startTime=(cal_time_min + init_delta).isot + "Z",
            transitDate=cal_transit.isot + "Z",
            filterStart=config["filter"],
            optFrq=config["squint"],
            corAzel=corazel,
            duration=f"{int(np.round((cal_time_max - cal_time_min - init_delta - end_dt).sec))}s"
        )
        cal_parset.add_phase_center(
            target=f"{cal_target}_TRACKING",
            subbandFrq=imaging_freq,
            useParentPointing=True,
            parameters=imaging_params
        )
        tmin_str = cal_time_min.isot.replace(":", "").replace("-", "")
        tmax_str = cal_time_max.isot.replace(":", "").replace("-", "")
        cal_parset_name = os.path.join(
            output_path,
            f"{file_prefix}{tmin_str}_{tmax_str}_{cal_target}_calibration.parset_user"
        )
        cal_parset.write( cal_parset_name )
        log.info(f"Parset '{cal_parset_name}' written.")
    

    parset.output["hd_bitMode"] = 16

    tmin_str = time_min.isot.replace(":", "").replace("-", "")
    tmax_str = time_max.isot.replace(":", "").replace("-", "")
    parset_name = os.path.join(
        output_path,
        f"{file_prefix}{tmin_str}_{tmax_str}_{main_target_name}.parset_user"
    )

    parset.write( parset_name )

    log.info(f"Parset '{parset_name}' written.")
