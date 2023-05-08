#!/bin/sh
python3 Optimisation.py -e 3100 -s uni200 -z no_export
python3 Optimisation.py -e 3100 -s uni100 -z no_export
python3 Optimisation.py -e 3100 -s uni50 -z no_export
python3 Optimisation.py -e 3100 -s uni0 -z no_export
python3 Optimisation.py -e 3200 -s uni200 -z no_export
python3 Optimisation.py -e 3200 -s uni100 -z no_export
python3 Optimisation.py -e 3200 -s uni50 -z no_export
python3 Optimisation.py -e 3200 -s uni0 -z no_export
python3 Optimisation.py -e 3100 -s uni200 -z export
python3 Optimisation.py -e 3100 -s uni100 -z export
python3 Optimisation.py -e 3100 -s uni50 -z export
python3 Optimisation.py -e 3100 -s uni0 -z export
python3 Optimisation.py -e 3200 -s uni200 -z export
python3 Optimisation.py -e 3200 -s uni100 -z export
python3 Optimisation.py -e 3200 -s uni50 -z export
python3 Optimisation.py -e 3200 -s uni0 -z export
