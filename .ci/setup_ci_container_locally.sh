#!/usr/bin/env bash

# Used to setup the Travis-CI environment in a docker container
# Usage:
#   docker run -it -v $(pwd):/mts ubuntu:xenial
#   cd /mts
#   . .ci/setup_ci_container_locally.sh


CONDA=2.7
PKG="mantid-total-scattering-python-wrapper"

apt update -y 
apt install wget bzip2 vim freeglut3-dev libglu1-mesa -y

MINICONDA_URL="https://repo.continuum.io/miniconda"
MINICONDA_FILE="Miniconda${CONDA:0:1}-latest-Linux-x86_64.sh"
wget "${MINICONDA_URL}/${MINICONDA_FILE}"
bash ${MINICONDA_FILE} -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

conda config --set always_yes yes --set changeps1 no --set anaconda_upload no;
conda config --add channels conda-forge --add channels mantid;
conda update  -q conda;
conda info -a;

conda create -q -n ${PKG} python=${CONDA} flake8 conda-build=3.17 conda-verify anaconda-client pytest;
source activate ${PKG};

# Build recipe and install 
PKG_PATH="./conda.recipe/${PKG}";
conda build ${PKG_PATH};
export PKG_FILE=$(conda build ${PKG_PATH} --output);
conda install ${PKG_FILE};

cd $HOME
conda install -c mantid/label/nightly mantid-framework=4 python=${CONDA};
which python
python -V
python setup.py test
