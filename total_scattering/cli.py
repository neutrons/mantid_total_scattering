from __future__ import (absolute_import, division, print_function)

import json
from total_scattering.reduction import TotalScatteringReduction


def main(config=None):

    # Read in JSON if not provided to main()
    if not config:
        import argparse
        parser = argparse.ArgumentParser(
            description="Absolute normalization PDF generation")
        parser.add_argument('json', help='Input json file')
        options = parser.parse_args()
        print("loading config from '%s'" % options.json)
        with open(options.json, 'r') as handle:
            config = json.load(handle)

    # Run total scattering reduction
    TotalScatteringReduction(config)


if __name__ == "__main__":
    main()
