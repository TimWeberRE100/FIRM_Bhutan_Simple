#!/bin/sh
python3 Optimisation.py -e 23 -s existing
python3 Optimisation.py -e 26 -s existing
python3 Optimisation.py -e 220 -s existing
python3 Optimisation.py -e 23 -s construction
python3 Optimisation.py -e 26 -s construction
python3 Optimisation.py -e 220 -s construction
python3 Optimisation.py -e 23 -s all
python3 Optimisation.py -e 26 -s all
python3 Optimisation.py -e 220 -s all