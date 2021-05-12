import json
import numpy as np

from Calibration.tofpd import diagnostics
from mantid.dataobjects import \
    EventWorkspace, \
    MaskWorkspace, \
    TableWorkspace
from mantid.simpleapi import \
    ConvertUnits, \
    CreateEmptyTableWorkspace, \
    FitPeaks, \
    LoadEventAndCompress, \
    Load, \
    Rebin, \
    mtd

# Diamond peak positions to perform multiple peak fitting on
DIAMOND_PEAKS = (0.8920, 1.0758, 1.2615)


def similarity_matrix_crosscorr(wksp: EventWorkspace, mask_wksp: MaskWorkspace):
    return


def similarity_matrix_degelder(wksp: EventWorkspace, mask_wksp: MaskWorkspace):
    return


def peakfitting(wksp: EventWorkspace, peaks, param_names, param_values):
    # Create the peak windows from the diamond peak positions
    peakwindows = diagnostics.get_peakwindows(np.asarray(peaks))

    # Convert from TOF to dspacing
    wksp = ConvertUnits(InputWorkspace=wksp, Target="dSpacing", EMode="Elastic")

    # Perform multiple peak fitting
    output = FitPeaks(InputWorkspace=wksp,
                      PeakFunction="PseudoVoigt",
                      PeakParameterNames=param_names,
                      PeakParameterValues=param_values,
                      RawPeakParameters=True,
                      HighBackground=False,
                      ConstrainPeakPositions=False,
                      MinimumPeakHeight=3,
                      PeakCenters=np.asarray(DIAMOND_PEAKS),
                      FitWindowBoundaryList=peakwindows,
                      FittedPeaksWorkspace='fitted',
                      OutputPeakParametersWorkspace='parameters',
                      OutputParameterFitErrorsWorkspace='fiterrors')

    return 'parameters', 'fiterrors'


def gather_fitparameters(paramws: TableWorkspace, errorws: TableWorkspace, mask,
                         peaks, threshold: float):
    '''Generate an array of peak fitting parameters over all spectra from FitPeaks results'''

    paramws = mtd[str(paramws)]
    errorws = mtd[str(errorws)]

    # The fit function parameters to gather
    cols = ["Mixing", "Intensity", "PeakCentre", "FWHM"]
    nprops = len(cols)*len(peaks)

    wsindex = paramws.column("wsindex")
    wsindex_unique = np.unique(wsindex)
    n = len(wsindex_unique)

    result = np.ndarray(shape=(n, nprops))

    for i in range(n):
        index = int(np.searchsorted(wsindex, wsindex_unique[i]))

        # Get the table rows corresponding to each peak index for this ws index
        nskip = 0
        for j in range(len(peaks)):
            row = paramws.row(index+j)

            # Skip bad fitting parameters based on fiterror values
            skip = False
            nzero = 0
            for k in range(len(cols)):
                p = errorws.row(index + j)[cols[k]]
                if np.isnan(p):
                    skip = True
                    break
                if p <= 0.0:
                    nzero += 1
            if skip or nzero == len(cols):
                nskip += 1
                continue

            for k in range(len(cols)):
                result[i, j*len(cols)+k] = row[cols[k]]

    return result


def gathered_parameters_to_tablewksp(wsname: str, result: np.ndarray, peaks):
    '''Converts the result from gather_fitparameters() to a Mantid TableWorkspace'''
    cols = ["Mixing", "Intensity", "PeakCentre", "FWHM"]
    if result.shape != (len(result), len(peaks) * len(cols)):
        raise ValueError("input array does not match expected size of {}"
                         .format((len(result), len(peaks) * len(cols))))

    tab = CreateEmptyTableWorkspace(OutputWorkspace=wsname)
    tab.addColumn("double", "wsindex")
    for peak in peaks:
        for col in cols:
            tab.addColumn("double", "{}-{}".format(col, peak))
    for i in range(len(result)):
        tab.addRow(np.insert(result[i], 0, i))
    return tab


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
    mask_wksp = None
    if masking_file:
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
        params, fiterrors = peakfitting(wksp, DIAMOND_PEAKS, "Mixing", "0.6")
        clustering_input = gather_fitparameters(params, fiterrors, None, DIAMOND_PEAKS, 1.0)

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
