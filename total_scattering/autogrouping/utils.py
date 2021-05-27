import numpy as np

from Calibration.tofpd import diagnostics
from mantid.dataobjects import \
    EventWorkspace, \
    TableWorkspace
from mantid.simpleapi import \
    ConvertUnits, \
    CreateEmptyTableWorkspace, \
    FitPeaks, \
    mtd
from total_scattering.reduction.total_scattering_reduction import compress_ints
from total_scattering.autogrouping.similarity import similarity_metric


def is_badfit(row, colnames, thresholds=dict(), max_chi=None):
    '''Determine if a single peakindex row for a workspace
    index has a bad fit result'''
    isbad = False
    nzero = 0

    # Filter by chi^2 if enabled
    if max_chi is not None and row["chi2"] > max_chi:
        return True

    for k in range(len(colnames)):
        p = row[colnames[k]]
        if np.isnan(p):
            isbad = True
            break
        if p <= 0.0:
            nzero += 1

        # Restrict fit against parameter thresholds, if any
        for threshold in thresholds:
            if colnames[k] == threshold:
                if p < thresholds[threshold][0] or \
                        p > thresholds[threshold][1]:
                    isbad = True
                    break
        if isbad:
            break
    # Count as "bad" if parameter errors are nan, or all 0
    if nzero == len(colnames):
        isbad = True
    return isbad


def get_badfitcount(index, ws, peaks, colnames, thresholds=dict(),
                    max_chi=None):
    '''Returns the number of bad fits in for a given workspace index
    (including all peakindex rows)'''
    # Note: index is the global table row (i.e, if ws
    # index = 3000, 3 peaks = 9000)
    nbad = 0
    for j in range(len(peaks)):
        if is_badfit(ws.row(index + j), colnames, thresholds, max_chi):
            nbad += 1
    return nbad


def get_goodfits(paramws, peaks, colnames, thresholds=dict(),
                 max_chi=None):
    '''
    Return a list of ws indices containing good fit parameters to include
    '''
    fitlist = []
    masklist = []

    wsindex = np.asarray(paramws.column("wsindex"))
    wsindex_unique = np.unique(wsindex)
    n = len(wsindex_unique)
    for i in range(n):
        ind = int(np.searchsorted(wsindex, wsindex_unique[i]))
        # Add ws index to list if fit result for at least one peak is "good"
        if get_badfitcount(ind, paramws, peaks, colnames, thresholds,
                           max_chi) != len(peaks):
            fitlist.append(i)
        else:
            # If all fits were bad, then we want to mask this index
            masklist.append(i)

    return fitlist, masklist


def peakfitting(wksp: EventWorkspace, peaks: np.ndarray, **fitpeaks_args):
    '''
    Converts wksp to dspacing and runs FitPeaks, returns parameter and
    parameter error fit workspaces.
    :param wksp: Input event workspace to convert and run FitPeaks on
    :param peaks: An array containing diamond peaks to fit
    :param fitpeaks_args: Dictionary of options to pass to FitPeaks
    '''
    # Create the peak windows from the diamond peak positions
    peakwindows = diagnostics.get_peakwindows(peaks)

    # Convert from TOF to dspacing
    wksp = ConvertUnits(InputWorkspace=wksp, Target="dSpacing",
                        EMode="Elastic")

    # Perform multiple peak fitting
    FitPeaks(InputWorkspace=wksp,
             RawPeakParameters=True,
             PeakCenters=peaks,
             FitWindowBoundaryList=peakwindows,
             FittedPeaksWorkspace='fitted',
             OutputWorkspace='output',
             OutputPeakParametersWorkspace='parameters',
             OutputParameterFitErrorsWorkspace='fiterrors',
             **fitpeaks_args)

    return 'parameters', 'fiterrors'


def gather_fitparameters(paramws: TableWorkspace, cols, mask,
                         peaks, thresholds=dict(), max_chi=None):
    '''Generate an array of peak fitting parameters over all spectra from
    FitPeaks results'''

    paramws = mtd[str(paramws)]

    # The fit function parameters to gather
    nprops = len(cols) * len(peaks) + 1

    wsindex = paramws.column("wsindex")
    wsindex_unique = np.unique(wsindex)

    fitlist, mask = get_goodfits(paramws, peaks, cols, thresholds, max_chi)
    n = len(fitlist)

    print("Number of good fits: {}".format(len(fitlist)))
    print("Number of bad fits (to be masked): {}".format(len(mask)))
    print("Mask list (wsindices): {}".format(compress_ints(mask)))

    frac = int(n / 10)  # get a tenth of total spectra for progress updates

    result = np.zeros(shape=(n, nprops))

    print("Gathering fit parameters...")
    for i in range(n):
        # Get mapping to ws index
        ws_index = fitlist[i]
        index = int(np.searchsorted(wsindex, wsindex_unique[ws_index]))

        # Get the table rows corresponding to each peak index for this ws index
        result[i, 0] = ws_index
        for j in range(len(peaks)):
            row = paramws.row(index + j)

            # Skip bad fitting parameters based on fiterror values
            if is_badfit(row, cols, thresholds, max_chi):
                continue

            for k in range(len(cols)):
                result[i, j * len(cols) + k + 1] = row[cols[k]]

        if i % frac == 0:
            print("-- {}/{} spectra".format(i, n))

    return result, mask


def gathered_parameters_to_tablewksp(wsname: str, result: np.ndarray,
                                     peaks, cols):
    '''Converts the result from gather_fitparameters() to a Mantid
    TableWorkspace'''
    if result.shape != (len(result), len(peaks) * len(cols) + 1):
        raise ValueError("input array does not match expected size of {}"
                         .format((len(result), len(peaks) * len(cols) + 1)))

    tab = CreateEmptyTableWorkspace(OutputWorkspace=wsname)
    tab.addColumn("double", "wsindex")
    for peak in peaks:
        for col in cols:
            tab.addColumn("double", "{}-{}".format(col, peak))
    for i in range(len(result)):
        tab.addRow(result[i])
    return tab


def display_parameter_stats(parameters: np.ndarray, peaks, cols):
    '''Prints out statistics about fit parameters for each peak'''
    for i in range(len(peaks)):
        for j in range(len(cols)):
            param = parameters[..., i * len(cols) + j]
            print("Param {}-{}:".format(cols[j], peaks[i]))
            print("   MAX = {}".format(np.max(param)))
            print("   MIN = {}".format(np.min(param)))
            print("   AVG = {}".format(np.mean(param)))
            print("   STD = {}".format(np.std(param)))
            print("")


def similarity_matrix_degelder(wksp):
    '''
    Generate the similarity matrix using the deGelder function.
    Returns an nxn matrix where n=number of spectra in wksp, where
    each entry of the matrix contains the deGelder similarity of that
     ith pixel to the jth pixel.
    Since this matrix is symmetric, this calculates the upper
    triangle of the matrix

    Note: for large input workspaces, this method can take quite a long time
    (~100 min for n=10000)
    '''
    print("Computing deGelder similarity matrix...")

    sm = similarity_metric()
    n = wksp.getNumberHistograms()
    y = wksp.extractY()

    frac = int(n / 10)  # get a tenth of total spectra for progress updates

    result = np.zeros(shape=(n, n))
    for i in range(n):
        for j in range(i, n):
            result[i][j] = sm.de_gelder_similarity(y[i], y[j])

        if i % frac == 0:
            print("-- {}/{} spectra".format(i, n))

    # Zero out any nans since this causes problems with clustering
    result[np.isnan(result)] = 0.0
    return result


def similarity_matrix_crosscorr(wksp):
    print("Computing pointwise crosscorr similarity matrix...")

    sm = similarity_metric()
    n = wksp.getNumberHistograms()
    y = wksp.extractY()

    frac = int(n / 10)  # get a tenth of total spectra for progress updates

    result = np.zeros(shape=(n, n))
    for i in range(n):
        for j in range(i, n):
            result[i][j] = sm.pointwise_squared_difference_similarity(
                y[i], y[j])

        if i % frac == 0:
            print("-- {}/{} spectra".format(i, n))

    # Zero out any nans since this causes problems with clustering
    result[np.isnan(result)] = 0.0
    return result


def get_detector_mask(wksp, ws_mask):
    '''
    Converts the mask of ws indices to their corresponding detector ids
    :param wksp: Donor workspace containing detector mapping
    :param ws_mask: List of workspace indices to mask
    :return: List of detector IDs to mask
    '''
    mask = []
    for index in ws_mask:
        # Get the detector id from this index
        det_id = wksp.getDetector(index).getID()
        mask.append(det_id)
    return mask
