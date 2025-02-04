__author__ = "Alan Loh"
__copyright__ = "Copyright 2025, exoschedule"
__credits__ = ["Alan Loh"]
__license__ = "MIT"
__version__ = "0.0.2"
__maintainer__ = "Alan Loh"
__email__ = "alan.loh@obspm.fr"

import logging
import sys

# ============================================================= #
# ------------------- Logging configuration ------------------- #
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)