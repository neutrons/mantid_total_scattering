| DevOps |
|--------|
| [![Build Status](https://travis-ci.org/marshallmcdonnell/mantid_total_scattering.svg?branch=master)](https://travis-ci.org/marshallmcdonnell/mantid_total_scattering) |
| [![codecov](https://codecov.io/gh/marshallmcdonnell/mantid_total_scattering/branch/master/graph/badge.svg)](https://codecov.io/gh/marshallmcdonnell/mantid_total_scattering) |


Total Scattering Data Reduction using Mantid Framework
-----------------------------------------------------------

This project is trying to implement total scattering data reduction for neutron time-of-flight diffractometers using the algorithms currently available in the [Mantid framework](https://github.com/mantidproject/mantid)


This entails taking raw neutron counts from detectors in the diffraction experiment and turning them into the reciprocal-space structure factor patterns, F(Q) or S(Q), and applying a Fourier Transform to real-space to give the pair distribution fuction, PDF.

This is the future backend for the [ADDIE project](https://github.com/neutrons/addie) and hopes to support multiple diffractometers performing total scattering measurements.


Structure factor S(Q) -> Pair Distribution Function G(r)
-----------------------------------------------------------
![alt text](https://raw.githubusercontent.com/marshallmcdonnell/mantid_total_scattering/master/images/sofq_to_gofr.png)


## Getting started for development

Clone the repository to a local directory

```bash
git clone https://github.com/marshallmcdonnell/mantid_total_scattering.git
cd mantid_total_scattering
```

### Setup using Pipenv

To setup the development environment with [pipenv](https://pipenv.readthedocs.io):

1. Install `pipenv`: `pip install --user pipenv` or read [this](https://pipenv.readthedocs.io/en/latest/install/)
2. Setup virtualenv with dependencies and `mantid_total_scattering` installed: `pipenv install -e .`
3. Activate the environment to run interatively: `pipenv shell`

NOTE: On Step 3, if you get something like
"Shell for UNKNOWN_VIRTUAL_ENVIRONMENT already activated.", the shell is already running from install.
Usually, do a `deactivate` and then repeat Step 3.

### Use CLI reduction tool

To launch the total scattering script, complete the input JSON file (found in `examples` directory), and run:

```bash
mantidtotalscattering examples/sns/nomad_simple.json
```

If you need to specify the path to Mantid build, use:
```
MANTIDPATH=/path/to/mantid/build/bin PATH=$MANTIDPATH:$PATH PYTHONPATH=$MANTIDPATH:$PATH mantidtotalscattering <json input>
```

## Running the tests
To build and run the tests via [pytest](https://docs.pytest.org), use:
```bash
/path/to/mantid/build/bin/mantidpython setup.py test
```

To build and run tests via [Docker](https://docs.docker.com/), use:

```bash
docker build -t unit-test-env -f .ci/Dockerfile.nightly_ubuntu16.04 . && docker run -t unit-test-env /bin/bash -c "mantidpython -m pytest"
```


