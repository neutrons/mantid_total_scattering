| Health | Main Release (w/ Mantid Framework) | Dev Release (w/o Mantid Framework) |
|--------|------------------------------------|------------------------------------|
|[![Build Status](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fneutrons%2Fmantid_total_scattering%2Fbadge%3Fref%3Dmain&style=plastic)](https://actions-badge.atrox.dev/neutrons/mantid_total_scattering/goto?ref=main)  | [![Anaconda-Server Badge](https://anaconda.org/neutrons/mantid-total-scattering/badges/version.svg)](https://anaconda.org/neutrons/mantid-total-scattering) | [![Anaconda-Server Badge](https://anaconda.org/neutrons/mantid-total-scattering-python-wrapper/badges/version.svg)](https://anaconda.org/neutrons/mantid-total-scattering-python-wrapper) |
| [![codecov](https://codecov.io/gh/neutrons/mantid_total_scattering/branch/main/graph/badge.svg)](https://codecov.io/gh/neutrons/mantid_total_scattering) | [![Anaconda-Server Badge](https://anaconda.org/neutrons/mantid-total-scattering/badges/platforms.svg)](https://anaconda.org/neutrons/mantid-total-scattering) | [![Anaconda-Server Badge](https://anaconda.org/neutrons/mantid-total-scattering-python-wrapper/badges/platforms.svg)](https://anaconda.org/neutrons/mantid-total-scattering-python-wrapper) |
| |  | [![PyPI version](https://badge.fury.io/py/mantid-total-scattering.svg)](https://badge.fury.io/py/mantid-total-scattering) |

Total Scattering Data Reduction using Mantid Framework
-----------------------------------------------------------

This project is trying to implement total scattering data reduction for neutron time-of-flight diffractometers using the algorithms currently available in the [Mantid framework](https://github.com/mantidproject/mantid)

This entails taking raw neutron counts from detectors in the diffraction experiment and turning them into the reciprocal-space structure factor patterns, F(Q) or S(Q), and applying a Fourier Transform to real-space to give the pair distribution fuction, PDF.

This is the future backend for the [ADDIE project](https://github.com/neutrons/addie) and hopes to support multiple diffractometers performing total scattering measurements.

Structure factor S(Q) -> Pair Distribution Function G(r)
-----------------------------------------------------------
![alt text](https://raw.githubusercontent.com/neutrons/mantid_total_scattering/next/images/sofq_to_gofr.png)

Running `mantidtotalscattering` will generate the total scattering data in
reciprocal space saved as NeXus file, if one is using the multiple bank mode.
Another alternative mode is the single bank mode which will merge spectra from
all detectors into a single pattern. The single bank mode will be taken in the
auto reduction implementation since otherwise manual efforts are needed to merge
the different banks data. When using the single bank mode, no Bragg data will be
saved. With the multiple bank mode, both the total scattering data and the Bragg
diffraction data will be saved. For the Bragg data output, there are two
available formats -- the GSS format and the XYE format. The GSS data file could
be loaded into the ADDIE interface (see the link above) for visualization. Here,
it is noteworthy that if one is loading in the GSS data in Mantid, one has to
rebin the loaded in workspace first (since the output GSS data is ragged, i.e.,
different banks are not sharing the common x-axis, for the purpose of removing
the non-sense data in each bank), followed by running the Mantid algorithm
`ConvertToPointData`.

Installation
===========================================================

Pre-requiste: Installing `pixi`
---------------------------------

Please follow the instructions provided at https://pixi.sh/latest/installation/

Mantid Framework Included
-----------------------------------------------------------

### Options

#### 1. Setup from Source

The following will clone the repo and instruct pixi to download and setup the dependencies in a sandboxed environment.

```bash
git clone https://github.com/neutrons/mantid_total_scattering.git
cd mantid_total_scattering
pixi install
```

To update
```bash
pixi update
```

> Cheat sheet for running local tasks with pixi,

```bash
pixi shell  # activate the pixi shell, something analogous to `conda activate`
mantidtotalscattering input.json  # the command is only available when pixi shell is active
pytest  # run tests
flake8 total_scattering/  # format source codes with `flake8`
exit  # exit the pixi shell
pixi run mantidtotalscattering input.json  # run application task without activating pixi shell
pixi run pytest
pixi run flake8 total_scattering/
pixi run test  # specifically defined test task
pixi run link  # specifically defined link task
pixi shell -e local-mantid  # activate the pixi shell with the `local-mantid` environment
```

#### Setup from Anaconda

```bash
conda install -c neutrons mantid-total-scattering
```

To update
```bash
conda update -c neutrons mantid-total-scattering
```

#### Notes

If you have an error (see below for example) related to the `libGL` library, you may not have it installed for the Mantid Framework to work. See instructions [here](https://github.com/mantidproject/conda-recipes/#gl-and-glu-libs) for installing the necessary libraries for different OS

Example error:
`ImportError: First import of "._api" failed with "libGL.so.1: cannot open shared object file...`


Usage (CLI reduction tool)
===========================================================

To launch the total scattering script, complete the input JSON file (found in `examples` directory), and run:

#### From Source
```bash
pixi run mantidtotalscattering examples/sns/nomad_simple.json
```

#### From Conda Install
```bash
mantidtotalscattering examples/sns/nomad_simple.json
```

Development
===========================================================

Please follow the directions provided [here[(https://developer.mantidproject.org/GettingStarted/GettingStartedCondaLinux.html#gettingstartedcondalinux)] to build a local version of mantid.
Once built, you must update some values in your `pyproject.toml`:
 1. `MANTID_BUILD_DIR`: The path to your local mantid's build directory.
 2. `MANTID_SRC_DIR`: The path to your local mantid's repository root.

You may run `mantidtotalscattering` as normal via:

```bash
pixi run mantidtotalscattering-local examples/sns/nomad_simple.json
```

Or you may launch `mantidworkbench` with the following command:

```bash
pixi run workbench-local
```

This will launch `mantidworkbench` in a seperate sandboxed pixi environement.
Then we can load the following script into `mantidworkbench` and execute it,

```python
import json
from total_scattering.reduction import TotalScatteringReduction
from total_scattering.reduction import validateConfig

with open("FULL_PATH_TO_INPUT_JSON_FILE", "r") as handle:
    config = json.load(handle)

# validate the config
validateConfig(config)

# Run total scattering reduction
TotalScatteringReduction(config)
```

where we need to replace `FULL_PATH_TO_INPUT_JSON_FILE` with the full path to
our input json file.

Or, we can use the following script inside a general version of `MantidWorkbench`,

```python
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

mts_path = os.path.join(
    "/gpfs/neutronsfs/instruments/NOM/shared/Dev/mantid_total_scattering"
)
sys.path.append(mts_path)
mts_path = os.path.join(
    "/gpfs/neutronsfs/instruments/NOM/shared/Dev",
    "mantid_total_scattering/.venv/lib/python3.10/site-packages"
)
sys.path.append(mts_path)
from total_scattering.reduction.total_scattering_reduction import \
    TotalScatteringReduction  # noqa: E407
from total_scattering.reduction.validator import validateConfig  # noqa: E407
from total_scattering import __version__ as mts_version  # noqa: E407
from total_scattering.reduction.normalizations import Material  # noqa: E407

import json

with open("/SNS/PG3/IPTS-34378/shared/autoreduce/MTSRed/Input/PG3_59120.json", "r") as f:
    config = json.load(f)

sofq = TotalScatteringReduction(config)
```


Tests
===========================================================
To build and run the tests via [pytest](https://docs.pytest.org), use:
```bash
pixi run test
```
**N. B.** This is assuming that the pixi environment mentioned
above is installed.

Tagging a New Version
===========================================================


```bash
git branch --track main origin/main  #  create a local main branch set to follow remote main
git checkout main
git fetch -p -t  # fetch all changes from the remote repo
git rebase -v origin/main  # sync with remote main branch
git merge --ff-only origin/qa  # merge in all of the changes in branch next
git tag v1.0.7  # create the tag in the format that versioneer has been configured
git push origin v1.0.7  # push the tag to remote to kick off the deploy step
```

The steps above will tag a new main release (the production release) and will
kick off the pipeline action to push the release to the `neutrons` Anaconda channel. 


If we want to tag a release candidate, follow the steps below.

```bash
git branch --track qa origin/qa  #  create a local main branch set to follow remote main
git checkout qa
git fetch -p -t  # fetch all changes from the remote repo
git rebase -v origin/qa  # sync with remote main branch
git merge --ff-only origin/next  # merge in all of the changes in branch next
git tag v1.0.8rc1  # create the tag in the format that versioneer has been configured
git push origin v1.0.8rc1  # push the tag to remote to kick off the deploy step
```

This will kick off the pipeline action to push the release candidate to the `neutrons` Anaconda channel.

Releases will be deployed on the Analysis Cluster via a seperate deployment pipeline.
