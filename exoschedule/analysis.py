#! /usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Alan Loh"
__copyright__ = "Copyright 2025, exoschedule"
__credits__ = ["Alan Loh"]
__maintainer__ = "Alan"
__email__ = "alan.loh@obspm.fr"
__status__ = "Production"
__all__ = [
    "source_statistics"
]

from exoschedule.functions import (
    find_directories,
    obsid_to_timerange,
    get_source_exposure_time,
    get_source_directories,
    get_piggyback_directories,
    is_imaging_observation
)

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from typing import Tuple, List
import logging
log = logging.getLogger("exoschedule")


# ============================================================= #
# --------------------- source_statistics --------------------- #
def source_statistics(sources: List[str], data_path: str, time_min: Time = None, time_max: Time = None, only_imaging: bool = True) -> dict:

    # Gather all repositories
    obs_directories = find_directories(
        data_path
    )

    # Filter by time range and by obs type if required
    obs_directories_filtered = []
    for obs_dir in obs_directories:
        start, stop = obsid_to_timerange(obs_dir)
        if time_min is not None:
            if stop <= time_min:
                continue
        if time_max is not None:
            if start >= time_max:
                continue
        if only_imaging:
            if not is_imaging_observation(obs_dir):
                continue
        obs_directories_filtered.append(obs_dir)

    # Populate the result dictionnary
    results = {}
    for source in sources:
        results[source] = {}
        
        results[source]["exposure_time"] = get_source_exposure_time(
            source_name=source,
            observation_directories=obs_directories
        )

        results[source]["n_observations"] = len(
            get_source_directories(
                source_name=source,
                observation_directories=obs_directories
            )
        )

        piggyback_directories = get_piggyback_directories(
            source_name=source,
            search_radius=8 * u.deg, # ~NenuFAR MA HWHM at 50 MHz
            observation_directories=obs_directories
        )
        results[source]["exposure_time_piggyback"] = get_source_exposure_time(
            source_name=None, # if None, compute the full exposure of the selected directories
            observation_directories=piggyback_directories
        )
