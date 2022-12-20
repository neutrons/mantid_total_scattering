from __future__ import absolute_import, division, print_function

import json
from total_scattering.reduction import TotalScatteringReduction
from total_scattering.reduction import validateConfig
from total_scattering import __version__ as mts_version
from total_scattering import __date__ as mts_date
from total_scattering import __commit__ as mts_commit
import sys


def main(config=None):

    version = mts_version
    up_date = mts_date
    up_commit = mts_commit

    # Read in JSON if not provided to main()
    if not config:
        if len(sys.argv) == 1:
            print("=======================MantidTotalScattering=======================")
            print("")
            print("Version      => " + version)
            print("Time Stamp   => " + mts_date)
            print("Commit SHA   => " + mts_commit)
            print("")
            print("=======================MantidTotalScattering=======================")
            print("")
            print("Go to https://github.com/neutrons/mantid_total_scattering/commits")
            print("Select a branch and check out the commit history with reference")
            print("to the SHA value printed above.")
            print("")
            print("=======================MantidTotalScattering=======================")
            print("")
            print("Here follows is listed the correspondence between the repo branch")
            print("and the deployed versions on SNS analysis.")
            print("")
            print("=======================MantidTotalScattering=======================")
            print("")
            print("|-------------------|-----------------------------|---------------|-------------|")
            print("| Deployed Version  | Command                     | GitHub Branch | Release Tag |")
            print("|-------------------|-----------------------------|---------------|-------------|")
            print("| Production        | mantidtotalscattering       | main          | vx.x.x      |")
            print("| Release Candidate | mantidtotalscattering --qa  | qa            | vx.x.xrcx   |")
            print("| Development       | mantidtotalscattering --dev | next          | N/A         |")
            print("|-------------------|-----------------------------|---------------|-------------|")
            return
        elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
            print("=======================MantidTotalScattering=======================")
            print("")
            print("Version      => " + version)
            print("Time Stamp   => " + mts_date)
            print("Commit SHA   => " + mts_commit)
            print("")
            print("=======================MantidTotalScattering=======================")
            print("")
            print("Go to https://github.com/neutrons/mantid_total_scattering/commits")
            print("Select a branch and check out the commit history with reference")
            print("to the SHA value printed above.")
            print("")
            print("=======================MantidTotalScattering=======================")
            print("")
            print("Here follows is listed the correspondence between the repo branch")
            print("and the deployed versions on SNS analysis.")
            print("")
            print("=======================MantidTotalScattering=======================")
            print("")
            print("|-------------------|-----------------------------|---------------|-------------|")
            print("| Deployed Version  | Command                     | GitHub Branch | Release Tag |")
            print("|-------------------|-----------------------------|---------------|-------------|")
            print("| Production        | mantidtotalscattering       | main          | vx.x.x      |")
            print("| Release Candidate | mantidtotalscattering --qa  | qa            | vx.x.xrcx   |")
            print("| Development       | mantidtotalscattering --dev | next          | N/A         |")
            print("|-------------------|-----------------------------|---------------|-------------|")
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
