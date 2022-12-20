from __future__ import absolute_import, division, print_function

import json
from total_scattering.reduction import TotalScatteringReduction
from total_scattering.reduction import validateConfig
from total_scattering import __version__ as mts_version
from total_scattering import __date__ as mts_date
from total_scattering import __commit__ as mts_commit
import sys

info_print = '''\n=======================MantidTotalScattering=======================

Version      => {0:s}
Time Stamp   => {1:s}
Commit SHA   => {2:s}

=======================MantidTotalScattering=======================

Go to https://github.com/neutrons/mantid_total_scattering/commits
Select a branch and check out the commit history with reference
to the SHA value printed above.

=======================MantidTotalScattering=======================

Here follows is listed the correspondence between the repo branch
and the deployed versions on SNS analysis.

=======================MantidTotalScattering=======================

|-------------------|-----------------------------|---------------|-------------|
| Deployed Version  | Command                     | GitHub Branch | Release Tag |
|-------------------|-----------------------------|---------------|-------------|
| Production        | mantidtotalscattering       | main          | vx.x.x      |
| Release Candidate | mantidtotalscattering --qa  | qa            | vx.x.xrcx   |
| Development       | mantidtotalscattering --dev | next          | N/A         |
|-------------------|-----------------------------|---------------|-------------|
'''


def main(config=None):

    version = mts_version
    up_date = mts_date
    up_commit = mts_commit

    # Read in JSON if not provided to main()
    if not config:
        if len(sys.argv) == 1:
            print(info_print.format(version, up_date, up_commit))
            return
        elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
            print(info_print.format(version, up_date, up_commit))
            return

        import argparse

        parser = argparse.ArgumentParser(
            description="Absolute normalization PDF generation"
        )
        parser.add_argument("json", help="Input json file")
        options = parser.parse_args()
        print("loading config from '%s'" % options.json)
        with open(options.json, "r") as handle:
            config = json.load(handle)

    # validate the config
    validateConfig(config)

    # Run total scattering reduction
    TotalScatteringReduction(config)


if __name__ == "__main__":
    main()
