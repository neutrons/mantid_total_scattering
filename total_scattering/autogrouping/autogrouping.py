import json
import numpy as np

from mantid.dataobjects import \
    EventWorkspace, \
    MaskWorkspace
from mantid.simpleapi import \
    Load


def similarity_matrix_crosscorr(wksp: EventWorkspace, mask_wksp: MaskWorkspace):
    return


def similarity_matrix_degelder(wksp: EventWorkspace, mask_wksp: MaskWorkspace):
    return


def euclidean_distance(wksp: EventWorkspace, mask_wksp: MaskWorkspace,
                       threshold: float):
    return


def get_key(key, config):
    '''Returns the value of key in config. Raises error if not found.'''
    value = config.get(key, None)
    if value is None:
        raise RuntimeError("Expected '{}' to be set in inputs".format(key))
    return value


def get_grouping_method(grouping):
    '''Make sure grouping is one of the supported options'''
    method = grouping.split("_")
    clusterings = ["DBSCAN", "KMEANS"]
    methods = ["DG", "CC", "ED"]

    if len(method) == 2 and method[0] in clusterings and method[1] in methods:
        return method
    else:
        raise ValueError("Invalid grouping method")


def Autogrouping(config):

    diamond_file = get_key("DiamondFile", config)
    masking_file = get_key("MaskFile", config)

    grouping_method = get_key("GroupingMethod", config)
    method = get_grouping_method(grouping_method)

    num_groups = get_key("NumberOutputGroups", config)
    threshold = get_key("Threshold", config)

    output_file = get_key("OutputGroupingFile", config)
    output_mask = get_key("OutputMaskFile", config)

    wksp = Load(Filename=diamond_file)
    mask_wksp = Load(Filename=masking_file)

    if method[1] == "DG":
        # Compute the deGelder similarity matrix
        grouping = similarity_matrix_degelder(wksp, mask_wksp)
    elif method[1] == "CC":
        # Compute the cross-correlation
        grouping = similarity_matrix_crosscorr(wksp, mask_wksp)
    elif method[1] == "ED":
        # Compute the euclidean distance
        grouping = euclidean_distance(wksp, mask_wksp, threshold)

    # TODO:
    # Pass result to clustering algorithm
    # Convert to grouping workspace
    # Output grouping file

    return


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
    Autogrouping(config)


if __name__ == "__main__":
    main()
