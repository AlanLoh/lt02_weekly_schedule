# lt02_weekly_schedule
 
## Installation on nancep

### Python virtual environment

The first step is to create a virtual environment in which all the python modules will be installed cleanly (e.g. `<name_of_virtual_env> = alan_python312`). In your home directory in nancep:
```
python3 -m venv .pyenv/versions/<name_of_virtual_env>
```

This environment can be activated deactivated using:
```
source .pyenv/versions/<name_of_virtual_env>/bin/activate
deactivate
``` 
It may be useful to create an alias in your bashrc-equivalent:
```
alias activate_my_env="source .pyenv/versions/<name_of_virtual_env>/bin/activate"
```

### Installing the codes

Once in your environment, install the required python modules (others may be required):
```
pip install nenupy --upgrade
pip install matplotlib
```

Install the LT02 scheduling package:
```
pip install https://github.com/AlanLoh/lt02_weekly_schedule/tarball/master --upgrade
```
This installation needs to be re-done whenever you want to apply a modification done on github.

## Running the code

Activate the virtual environment.

### Download the current VCR booking

VCR -> Schedule -> "Current Booking.csv"

### Run

From anywhere

```
exoschedule -t0 2025-04-02T00:00:00 -t1 2025-04-08T00:00:00 -n 2 -b <path_to_vcr_booking.csv> -o <path_to_store_the_results> -u
```

### Options
Here are all the available options:
* `-t0` (or `--time_min`): Start time (UTC, ISOT format). E.g. `-t0 2025-04-02T00:00:00`
* `-t1` (or `--time_max`): Stop time (UTC, ISOT format). If not user-defined, will be set to t0 + 7 days.
* `-u` (or `--update_data`): Update the exposure time in the JSON files (default option). 
* `-nu` (or `--no-update_data`): Do not update the exposure time in the JSON files (either select `-u` or `-nu`)
* `-b` (or `--vcr_booking_csv`): CSV booking file downloaded from the VCR.
* `-o` (or `--output_path`): Where every result / parset will be written.
* `-n` (or `--max_repetitions`): Maximum number of observations per target (default=2).
* `-s` (or `--priority_source`): Set these sources as priority targets (e.g. `-s "GJ_687" "TOI-3884"`)
