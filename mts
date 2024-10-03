#!/bin/bash

source /SNS/NOM/shared/Dev/mantid_total_scattering/.venv/bin/activate
python /SNS/NOM/shared/Dev/mantid/build/bin/AddPythonPath.py > /dev/null 2>&1

cwd=${PWD}

cd /SNS/NOM/shared/Dev/mantid_total_scattering/
python setup.py develop

cd ${cwd}

mantidtotalscattering $1
