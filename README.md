Total Scattering Data Reduction using Mantid Framework
-----------------------------------------------------------

This project is trying to implement total scattering datareduction for neutron time-of-flight diffractometers using the algorithms currently available in the Mantid framework (https://github.com/mantidproject/mantid)


This entails taking raw neutron counts from detectors in the diffraction experiment and turning them into the reciprocal-space structure factor, F(Q) or S(Q), and applying a Fourier Transform to real-space to give the pair distribution fuction, PDF.


Structure factor S(Q) -> Pair Distribution Function G(r)
-----------------------------------------------------------
![alt text](https://raw.githubusercontent.com/marshallmcdonnell/mantid_total_scattering/master/images/sofq_to_gofr.png)



Usage
-------

To launch the total scattering script, complete the input JSON file, and run:

    python mantid_total_scattering.py input.json


[![Build Status](https://travis-ci.org/marshallmcdonnell/mantid_total_scattering.svg?branch=master)](https://travis-ci.org/marshallmcdonnell/mantid_total_scattering)
