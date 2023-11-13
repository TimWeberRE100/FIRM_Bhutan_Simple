#!/bin/sh
python3 Optimisation.py -e 3 -s existing -y no_import -p 10 -i 1000
python3 Optimisation.py -e 6 -s existing -y no_import -p 10 -i 1000
python3 Optimisation.py -e 10 -s existing -y no_import -p 10 -i 1000
python3 Optimisation.py -e 15 -s existing -y no_import -p 10 -i 1000
python3 Optimisation.py -e 3 -s construction -y no_import -p 10 -i 1000
python3 Optimisation.py -e 6 -s construction -y no_import -p 10 -i 1000
python3 Optimisation.py -e 10 -s construction -y no_import -p 10 -i 1000
python3 Optimisation.py -e 15 -s construction -y no_import -p 10 -i 1000
python3 Optimisation.py -e 3 -s all -y no_import -p 10 -i 1000
python3 Optimisation.py -e 6 -s all -y no_import -p 10 -i 1000
python3 Optimisation.py -e 10 -s all -y no_import -p 10 -i 1000
python3 Optimisation.py -e 15 -s all -y no_import -p 10 -i 1000
