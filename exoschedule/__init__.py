__author__ = "Alan Loh"
__copyright__ = "Copyright 2025, exoschedule"
__credits__ = ["Alan Loh"]
__license__ = "MIT"
__version__ = "1.1.1"
__maintainer__ = "Alan Loh"
__email__ = "alan.loh@obspm.fr"


import logging
import sys
from astropy.time import TimeDelta

# ============================================================= #
# ------------------- Logging configuration ------------------- #
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ============================================================= #
# --------------------- Global Variables ---------------------- #
# Minimum booking slot
MIN_DURATION = TimeDelta(3600, format="sec")
# MIN_OBS_DURATION bypasses any 'minimal_obs_time_minutes' key in JSON files.
# Observations will preferentially be booked at this duration.
# In a second pass, if ther eis time left, they will be booked at MIN_DURATION
MIN_OBS_DURATION = TimeDelta(3600 * 2, format="sec")

# Time resolution of the Schedule
SCHEDULE_DT = TimeDelta(1 * 60, format="sec")

# No booking between 11-15 UTC except Saturday and Sunday
# to leave room for maintenance operations (and because the Sun
# is quite high and ppolluting).
MID_DAY_START_HOUR = 11
MID_DAY_STOP_HOUR = 15

KP_CODE = "LT02"

REMOTE_SERVER = "nancep7"

DATA_REPOSITORIES = (
    "/databf/nenufar/ES02",
    "/databf/nenufar/LT02"
)

# Disable the pointing correction for sources at declination > 60 deg
DISABLE_POINTING_CORRECTION_DECLINATION = 60

OBSERVATION_CONFIGURATION = {
    "exoplanet": {
        "day": {
            "filter": 3,
            "frequency": [154, 345], # 30.1-67.4 MHz
            "squint": 40,
            "parameters": "tf: df=3.05 dt=21.0 hamm"
        },
        "night": {
            "filter": 2,
            "frequency": [108, 299], # 21.1-58.4 MHz
            "squint": 40,
            "parameters": "tf: df=3.05 dt=21.0 hamm"
        }
    },
    "star": {
        "day": {
            "filter": 2,
            "frequency": [61, 444], # 12-86.8 MHz
            "squint": 60,
            "parameters": "tf: df=3.05 dt=42.0 hamm"
        },
        "night": {
            "filter": 1,
            "frequency": [61, 444], # 12-86.8 MHz
            "squint": 60,
            "parameters": "tf: df=3.05 dt=42.0 hamm"
        }
    }
}

CONFIG_EXCEPTION = {
    "CR_DRA": {
        "day": { 
            "frequency": [187, 378] # 36.5-73.8 MHz
        }
    },
    "TAU_BOO": {
        "night": {
            "frequency": [77, 268], # 15.1-52.4 MHz
            "filter": 1
        }
    }
} 

# The booking will be made in this order of priority (from left to right):
# "both" is for the special Ceti group case
TYPE_PRIORITY_ORDER = ["exoplanet-l", "exoplanet-p", "star-l", "exoplanet", "star", "both"]

# Tune how the source priority are computed (see schedule.sort_blocks_by_priorities)
TYPE_FACTOR = 3
EXPOSURE_FACTOR = 2
LAST_OBSERVE_FACTOR = 1

SPECIAL_CASES = {
    "LP212-62": {
        "include": ["LP_212-62"],
        "exclude": []
    },
    "TAU_BOOTIS": {
        "include": ["TAU_BOO", "TAU_BOOTIS_JOINT_CAMPAIGN_TAU_BOOTIS"],
        "exclude": ["TAU_BOO_CALIBRATOR_B0809+74", "TAU_BOOTIS_JOINT_CAMPAIGN_CALIBRATOR_B0809+74"]
    },
    "B0809+74_CALIBRATOR": {
        "include": ["TAU_BOO_CALIBRATOR_B0809+74", "TAU_BOOTIS_JOINT_CAMPAIGN_CALIBRATOR_B0809+74"],
        "exclude": []
    },
    "WOLF_28 (VANMAANEN)": {
        "include": ["WOLF_28", "VAN_MAANEN_STAR"],
        "exclude": []
    },
    "V830_TAU": {
        "include": ["V830TAU"],
        "exclude": []
    },
    "WASP-49": {
        "include": ["WASP_49"],
        "exclude": []
    },
    "V0612_CYG": {
        "include": ["V612_CYG"],
        "exclude": []
    },
    "T_CRB": {
        "include": ["T_CRB_TRACKING2"],
        "exclude": []
    },
    "GJ_436 (LALANDE_21185)": {
        "include": ["GJ_436"],
        "exclude": []
    },
    "LTT_3780 (TOI-732)": {
        "include": ["LTT_3780"],
        "exclude": []
    },
    "GJ_806 (TOI-4481)": {
        "include": ["GJ_806"],
        "exclude": []
    } 
}
# ============================================================= #
# ============================================================= #
