| DevOps | Main Release (w/ Mantid Framework) | Dev Release (w/o Mantid Framework) |
|--------|------------------------------------|------------------------------------|
| [![Build Status](https://travis-ci.org/marshallmcdonnell/mantid_total_scattering.svg?branch=master)](https://travis-ci.org/marshallmcdonnell/mantid_total_scattering) | [![Anaconda-Server Badge](https://anaconda.org/marshallmcdonnell/mantid-total-scattering/badges/version.svg)](https://anaconda.org/marshallmcdonnell/mantid-total-scattering) | [![Anaconda-Server Badge](https://anaconda.org/marshallmcdonnell/mantid-total-scattering-python-wrapper/badges/version.svg)](https://anaconda.org/marshallmcdonnell/mantid-total-scattering-python-wrapper) |
| [![codecov](https://codecov.io/gh/marshallmcdonnell/mantid_total_scattering/branch/master/graph/badge.svg)](https://codecov.io/gh/marshallmcdonnell/mantid_total_scattering) | [![Anaconda-Server Badge](https://anaconda.org/marshallmcdonnell/mantid-total-scattering/badges/platforms.svg)](https://anaconda.org/marshallmcdonnell/mantid-total-scattering) | [![Anaconda-Server Badge](https://anaconda.org/marshallmcdonnell/mantid-total-scattering-python-wrapper/badges/platforms.svg)](https://anaconda.org/marshallmcdonnell/mantid-total-scattering-python-wrapper) |
| |  | [![PyPI version](https://badge.fury.io/py/mantid-total-scattering.svg)](https://badge.fury.io/py/mantid-total-scattering) |


Total Scattering Data Reduction using Mantid Framework
-----------------------------------------------------------

This project is trying to implement total scattering data reduction for neutron time-of-flight diffractometers using the algorithms currently available in the [Mantid framework](https://github.com/mantidproject/mantid)


This entails taking raw neutron counts from detectors in the diffraction experiment and turning them into the reciprocal-space structure factor patterns, F(Q) or S(Q), and applying a Fourier Transform to real-space to give the pair distribution fuction, PDF.

This is the future backend for the [ADDIE project](https://github.com/neutrons/addie) and hopes to support multiple diffractometers performing total scattering measurements.


Structure factor S(Q) -> Pair Distribution Function G(r)
-----------------------------------------------------------
![alt text](https://raw.githubusercontent.com/marshallmcdonnell/mantid_total_scattering/master/images/sofq_to_gofr.png)


# Installation

## Mantid Framework Included

### Anaconda (Recommended)

**Setup** 

Add channels with dependencies, create a conda environment with `python_version` set to either `2.7.14` or `3.6`, and activate the environment

```
conda config --add channels conda-forge --add channels mantid --add channels mantid/label/nightly
conda create -n mantidts_env python=${python_version}
source activate mantidts_env
```

**Install (or Update)**

```
conda install -c marshallmcdonnell mantid-total-scattering
```

Go here for how to delete an [environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#removing-an-environment) or use:

```
conda remove --name mantidts_env --all
```

**Notes**

If you have an error (see below for example) related to the `libGL` library, you may not have it installed for the Mantid Framework to work. See instructions [here](https://github.com/mantidproject/conda-recipes/#gl-and-glu-libs) for installing the necessary libraries for different OS

Example error:
`ImportError: First import of "._api" failed with "libGL.so.1: cannot open shared object file...`

If you have an error that another version of Mantid is installed on the machine and being imported via `PYTHONPATH`, you can use the following as a workaround for CLI tool:

```
PYTHONPATH="" mantidtotalscattering
```

# Usage (CLI reduction tool)

To launch the total scattering script, complete the input JSON file (found in `examples` directory), and run:

```bash
mantidtotalscattering examples/sns/nomad_simple.json
```

If you need to specify the path to Mantid build, use:
```
MANTIDPATH=/path/to/mantid/build/bin PATH=$MANTIDPATH:$PATH PYTHONPATH=$MANTIDPATH:$PATH mantidtotalscattering <json input>
```

## Mantid Framework Not Included (for development)

This is mainly for development if you want to use a local development build of Mantid Framework instead of one included.

## PyPI (Recommended)

**Install**

```
pip install mantid-total-scattering
```

## Anaconda 

**Setup**

Add channels with dependencies, create a conda environment with `python_version` set to either `2.7` or `3.6`, and activate the environment

```
conda config --add channels conda-forge
conda create -n mantidts_env python=${python_version}
conda activate mantidts_env
```

**Install**

```
conda install -c marshallmcdonnell mantid-total-scattering-python-wrapper
```

# Development

Clone the repository to a local directory

```bash
git clone https://github.com/marshallmcdonnell/mantid_total_scattering.git
cd mantid_total_scattering
```

## Pipenv

To setup the development environment with [pipenv](https://pipenv.readthedocs.io):

1. Install `pipenv`: `pip install --user pipenv` or read [this](https://pipenv.readthedocs.io/en/latest/install/)
2. Setup virtualenv with dependencies and `mantid_total_scattering` installed: `pipenv install -e .`
3. Activate the environment to run interatively: `pipenv shell`

NOTE: On Step 3, if you get something like
"Shell for UNKNOWN_VIRTUAL_ENVIRONMENT already activated.", the shell is already running from install.
Usually, do a `deactivate` and then repeat Step 3.

# Tests
To build and run the tests via [pytest](https://docs.pytest.org), use:
```bash
/path/to/mantid/build/bin/mantidpython setup.py test
```

To build and run tests via [Docker](https://docs.docker.com/), use:

```bash
docker build -t unit-test-env -f .ci/Dockerfile.nightly_ubuntu16.04 . && docker run -t unit-test-env /bin/bash -c "mantidpython -m pytest"
```


