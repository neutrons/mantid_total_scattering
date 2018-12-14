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


## Getting started

To launch the total scattering script, complete the input JSON file, and run:

    python mantid_total_scattering.py input.json
    
## Running the tests
To build and run tests via [Docker](https://docs.docker.com/), use:

```bash
docker build -t unit-test-env -f docker/Dockerfile.test . && docker run -t unit-test-env "pytest"
```


