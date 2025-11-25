__author__ = "Alan Loh"
__copyright__ = "Copyright 2025, exoschedule"
__credits__ = ["Alan Loh"]
__license__ = "MIT"
__version__ = "1.2.0"
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

# A team sources used for calibration purposes during imaging observations
A_TEAM_SOURCES = ["Cas A", "Cyg A"]
CALIBRATION_DURATION = TimeDelta(15 * 60, format="sec")
IMAGING_MIN_DURATION = TimeDelta(2 * 3600, format="sec")

IMAGING_SOURCES = [
    "ROSS_248",
    "HD_189733",
    "GJ_687",
    "UPS_AND",
    "MCC_351",
    "LP_355-051",
    "YZ_CMI",
    "G_41-8",
    "NLTT_20670",
    "AD_LEO",
    "GJ_1134",
    "BK_CrB",
    "HD_239960A",
    "EV_LAC",
    "EQ_PEG"
]
# Mail 12 novembre 2025 17:20
# Les cibles Carmenes:
# Nom			RA			Dec
# MCC_351	00:08:27.3	+17:25:27.5
# LP_355-051	03:17:45.2	+25:15:06.4
# YZ_CMI		07:44:40.1	+03:33:08.9
# G_41-8		08:56:19.6	+12:39:49.8
# NLTT_20670	08:58:56.3	+08:28:26.0
# AD_LEO		10:19:36.2	+19:52:12.1
# GJ_1134		10:41:37.9	+37:36:39.2
# BK_CrB		15:36:50.3	+37:34:49.4
# HD_239960A	22:27:59.6	+57:41:42.1
# EV_LAC		22:46:49.6	+44:20:02.3
# EQ_PEG		23:31:52.1	+19:56:14

# Les cibles détectées par Xiang avec NenuFAR:
# UPS AND
# GJ 687 (same stellar rotation phase as Cyril’s detection)
# Ross 248 (brightest - 15 sigma detection - definitely resolved structure in burst)
# HAT-P-36 (two bursts with reversed polarimetry; suspicious due to terrible data quality - bad weather)
# 2MASS J04555897+2140007 (Resolved structure in burst)
# LP 298-45 (Potentially resolved structure)
# LSPM J1937+1747
# 2MASS J16071984+4801293
# 2MASS J22152477+3506576
# 2MASS J01220158+4943248
# 2MASS J09381115+4913012
# SDSS J101353.96+500359.5
# Unknown source 16:17:06.579 +68:55:51.168 (no M dwarf in pencil beam)

# Les cibles observés avec FAST:
# WISEPJ1122	11:22:54.70	+25:50:21.5
# BDR_J1750	17:50:00.00	+38:09:00.0
# Eps_Eri		03:32:55.84	-09:27:29.7
# GJ1151		11:50:57.72	+48:22:38.5

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
