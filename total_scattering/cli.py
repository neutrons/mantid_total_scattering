from __future__ import absolute_import, division, print_function

import json
from total_scattering.reduction import TotalScatteringReduction
from total_scattering.reduction import validateConfig
import sys


def main(config=None):

    version = "1.0.3.12012022.1"

    # Read in JSON if not provided to main()
    if not config:
        if len(sys.argv) == 1:
            print("mantidtotalscattering version => " + version)
            return
        elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
            print("mantidtotalscattering version => " + version)
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
