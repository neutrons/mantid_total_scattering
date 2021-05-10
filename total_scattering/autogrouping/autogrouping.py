import json
import numpy as np

from Calibration.tofpd import diagnostics
from mantid.dataobjects import \
    EventWorkspace, \
    MaskWorkspace
from mantid.simpleapi import \
    ConvertUnits, \
    FitPeaks, \
    LoadEventAndCompress, \
    Load, \
    Rebin

# Diamond peak positions to perform multiple peak fitting on
DIAMOND_PEAKS = (0.8920, 1.0758, 1.2615)


def similarity_matrix_crosscorr(wksp: EventWorkspace, mask_wksp: MaskWorkspace):
    return


def similarity_matrix_degelder(wksp: EventWorkspace, mask_wksp: MaskWorkspace):
    return


def euclidean_distance(wksp: EventWorkspace, mask_wksp: MaskWorkspace,
                       threshold: float):
    # Create the peak windows from the diamond peak positions
    peakwindows = diagnostics.get_peakwindows(np.asarray(DIAMOND_PEAKS))

    # Convert from TOF to dspacing
    wksp = ConvertUnits(InputWorkspace=wksp, Target="dSpacing", EMode="Elastic")

    # Perform multiple peak fitting
    output = FitPeaks(InputWorkspace=wksp,
                      PeakFunction="Bk2BkExpConvPV",
                      RawPeakParameters=True,
                      HighBackground=False,
                      ConstrainPeakPositions=False,
                      MinimumPeakHeight=3,
                      PeakCenters=np.asarray(DIAMOND_PEAKS),
                      FitWindowBoundaryList=peakwindows,
                      FittedPeaksWorkspace='fitted',
                      OutputPeakParametersWorkspace='parameters',
                      OutputParameterFitErrorsWorkspace='fiterrors')

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

    wksp = LoadEventAndCompress(Filename=diamond_file, FilterBadPulses=0)
    mask_wksp = Load(Filename=masking_file)

    # Rebin (in TOF)
    # TODO: Double check if this is needed, and if this should only be done for ED case
    wksp = Rebin(InputWorkspace=wksp, Params=(300, -0.001, 16666.7))

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
