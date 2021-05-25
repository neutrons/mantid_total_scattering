import json
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, DBSCAN

from Calibration.tofpd import diagnostics
from mantid.dataobjects import \
    EventWorkspace, \
    MaskWorkspace, \
    TableWorkspace
from mantid.simpleapi import \
    ConvertUnits, \
    CreateEmptyTableWorkspace, \
    CreateCacheFilename, \
    CreateGroupingWorkspace, \
    FitPeaks, \
    LoadEventAndCompress, \
    Load, \
    MaskSpectra, \
    Rebin, \
    RemoveMaskedSpectra, \
    SaveDetectorsGrouping, \
    SaveNexusProcessed, \
    mtd
from total_scattering.autogrouping.similarity import similarity_metric
from total_scattering.reduction.total_scattering_reduction import compress_ints


def is_badfit(row, colnames, thresholds=dict()):
    '''Determine if a single peakindex row for a workspace index has a bad fit result'''
    isbad = False
    nzero = 0
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
                if p < thresholds[threshold][0] or p > thresholds[threshold][1]:
                    isbad = True
                    break
        if isbad:
            break
    # Count as "bad" if parameter errors are nan, or all 0
    if nzero == len(colnames):
        isbad = True
    return isbad


def get_badfitcount(index, ws, peaks, colnames, thresholds=dict()):
    '''Returns the number of bad fits in for a given workspace index (including all peakindex rows)'''
    # Note: index is the global table row (i.e, if ws index = 3000, 3 peaks = 9000)
    nbad = 0
    for j in range(len(peaks)):
        if is_badfit(ws.row(index + j), colnames, thresholds):
            nbad += 1
    return nbad


def get_goodfits(paramws, peaks, colnames, thresholds=dict()):
    '''Return a list of ws indices containing good fit parameters to include'''
    fitlist = []
    masklist = []

    wsindex = np.asarray(paramws.column("wsindex"))
    wsindex_unique = np.unique(wsindex)
    n = len(wsindex_unique)
    for i in range(n):
        ind = int(np.searchsorted(wsindex, wsindex_unique[i]))
        # Add the ws index to list if fit result for at least one peak is "good"
        if get_badfitcount(ind, paramws, peaks, colnames, thresholds) != len(peaks):
            fitlist.append(i)
        else:
            # If all fits were bad, then we want to mask this index
            masklist.append(i)

    return fitlist, masklist


def get_cachename(prefix, cache_dir, props):
    filename, sig = CreateCacheFilename(OtherProperties=props, Prefix=prefix, CacheDir=cache_dir)
    return filename, sig


def peakfitting(wksp: EventWorkspace, peaks: np.ndarray, **fitpeaks_args):
    '''
    Converts wksp to dspacing and runs FitPeaks, returns parameter and parameter error fit workspaces.
    :param wksp: Input event workspace to convert and run FitPeaks on
    :param peaks: An array containing diamond peaks to fit
    :param fitpeaks_args: Dictionary of options to pass to FitPeaks
    '''
    # Create the peak windows from the diamond peak positions
    peakwindows = diagnostics.get_peakwindows(peaks)

    # Convert from TOF to dspacing
    wksp = ConvertUnits(InputWorkspace=wksp, Target="dSpacing", EMode="Elastic")

    # Perform multiple peak fitting
    output = FitPeaks(InputWorkspace=wksp,
                      RawPeakParameters=True,
                      PeakCenters=peaks,
                      FitWindowBoundaryList=peakwindows,
                      FittedPeaksWorkspace='fitted',
                      OutputPeakParametersWorkspace='parameters',
                      OutputParameterFitErrorsWorkspace='fiterrors',
                      **fitpeaks_args)

    return 'parameters', 'fiterrors'


def gather_fitparameters(paramws: TableWorkspace, cols, mask,
                         peaks, thresholds=dict()):
    '''Generate an array of peak fitting parameters over all spectra from FitPeaks results'''

    paramws = mtd[str(paramws)]

    # The fit function parameters to gather
    nprops = len(cols) * len(peaks) + 1

    wsindex = paramws.column("wsindex")
    wsindex_unique = np.unique(wsindex)

    fitlist, mask = get_goodfits(paramws, peaks, cols, thresholds)
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
            if is_badfit(row, cols, thresholds):
                continue

            for k in range(len(cols)):
                result[i, j * len(cols) + k + 1] = row[cols[k]]

        if i % frac == 0:
            print("-- {}/{} spectra".format(i, n))

    return result, mask


def gathered_parameters_to_tablewksp(wsname: str, result: np.ndarray, peaks, cols):
    '''Converts the result from gather_fitparameters() to a Mantid TableWorkspace'''
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


def plot_features(data, peaks, colnames, labels=None, centroids=None):
    n = len(colnames)
    npeaks = len(peaks)
    for p in range(npeaks):
        fig, ax = plt.subplots(n, n)

        fig.suptitle("peak {}".format(peaks[p]))

        # Plot each peak feature against its other features
        for i in range(n):
            for j in range(n):
                if labels is not None:
                    ax[i, j].scatter(data[..., p * n + i], data[..., p * n + j], c=labels, alpha=0.5)
                else:
                    ax[i, j].scatter(data[..., p * n + i], data[..., p * n + j])

                if centroids is not None:
                    ax[i, j].scatter(centroids[..., p * n + i], centroids[..., p * n + j],
                                     marker="X", c=np.unique(labels), s=200, alpha=0.75,
                                     edgecolors="red")
                ax[i, j].set_xlabel(colnames[i])
                ax[i, j].set_ylabel(colnames[j])

    return fig, ax


def similarity_matrix_degelder(wksp):
    '''
    Generate the similarity matrix using the deGelder function.
    Returns an nxn matrix where n=number of spectra in wksp, where
    each entry of the matrix contains the deGelder similarity of that ith pixel
    to the jth pixel.
    Since this matrix is symmetric, this calculates the upper triangle of the matrix

    Note: for large input workspaces, this method can take quite a long time (~100 min for n=10000)
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
            result[i][j] = sm.pointwise_squared_difference_similarity(y[i], y[j])

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
    diamond_file = get_key("DiamondFile", config)  # .nxs.h5 or .nxs input file
    masking_file = get_key("MaskFile", config)

    grouping_method = get_key("GroupingMethod", config)
    method = get_grouping_method(grouping_method)

    num_groups = int(get_key("NumberOutputGroups", config))  # k parameter for KMeans clustering
    epsilon = float(get_key("DBSCANEpsilon", config))  # eps parameter for DBSCAN clustering

    fitparams = get_key("FittingFunctionParameters", config).split(",")
    fitpeaks_args = get_key("FitPeaksArgs", config)
    diamond_peaks = np.asarray(get_key("DiamondPeaks", config).strip().split(","), dtype=float)
    thresholds = get_key("ParameterThresholds", config)
    for threshold in thresholds:
        thresholds[threshold] = tuple(float(x) for x in thresholds[threshold].strip("()").split(","))

    cache_dir = ""
    if "CacheDir" in config:
        cache_dir = os.path.abspath(get_key("CacheDir", config))

    output_file = os.path.abspath(get_key("OutputGroupingFile", config))
    output_mask = os.path.abspath(get_key("OutputMaskFile", config))

    output_fittable = ""
    if "OutputFitParamFile" in config:
        output_fittable = os.path.abspath(get_key("OutputFitParamFile", config))

    if diamond_file.endswith(".nxs.h5"):
        wksp = LoadEventAndCompress(Filename=diamond_file, FilterBadPulses=0)
    else:
        wksp = Load(Filename=diamond_file)

    mask = None
    new_mask = None
    if masking_file:
        # Check whether to load masking file with mantid or to an array with numpy
        if masking_file.endswith(".txt") or masking_file.endswith(".out"):
            mask = np.loadtxt(masking_file, dtype=int)
            print("Loaded detector mask: {}".format(mask))
            # mask numbers are detector ids - so convert them to workspace indices
            mask_ws = wksp.getIndicesFromDetectorIDs(mask.tolist())
            print("Mask converted to workspace indices: {}".format(mask_ws))
            # compress down to ranges
            mask_ws = compress_ints(mask_ws)
            print("Compressed mask of wsindex: {}".format(mask_ws))
            # mask spectra in workspace
            wksp = MaskSpectra(InputWorkspace=wksp, InputWorkspaceIndexType="WorkspaceIndex",
                               InputWorkspaceIndexSet=mask_ws)
            wksp = RemoveMaskedSpectra(InputWorkspace=wksp)

    # Rebin (in TOF)
    wksp = Rebin(InputWorkspace=wksp, Params=(300, -0.001, 16666.7))

    wsindex = list(range(wksp.getNumberHistograms()))
    clustering_input = None
    use_cache = False
    if method[1] == "DG":
        # Compute the deGelder similarity matrix
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms())]
            prefix = diamond_file.rstrip(".h5").rstrip(".nxs").split("/")[-1] + "_DG"
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            print("Checking for cachefile '{}'".format(cachefile))
            if os.path.exists(cachefile):
                # Load cached similarity matrix
                use_cache = True
                print("Found cached deGelder similarity matrix, loading from '{}'".format(cachefile))
                clustering_input = np.loadtxt(cachefile)

        if not use_cache:
            clustering_input = similarity_matrix_degelder(wksp)
            if cache_dir:
                # Save matrix if using caching
                print("Caching deGelder similarity matrix to '{}'".format(cachefile))
                np.savetxt(cachefile, clustering_input)
    elif method[1] == "CC":
        # Compute the cross-correlation matrix using a pointwise squared difference method
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms())]
            prefix = diamond_file.rstrip(".h5").rstrip(".nxs").split("/")[-1] + "_CC"
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            print("Checking for cachefile '{}'".format(cachefile))
            if os.path.exists(cachefile):
                # Load cached similarity matrix
                use_cache = True
                print("Found cached pointwise crosscorr similarity matrix, loading from '{}'".format(cachefile))
                clustering_input = np.loadtxt(cachefile)

        if not use_cache:
            clustering_input = similarity_matrix_crosscorr(wksp)
            if cache_dir:
                # Save matrix if using caching
                print("Caching pointwise crosscorr similarity matrix to '{}'".format(cachefile))
                np.savetxt(cachefile, clustering_input)
    elif method[1] == "ED":
        # Compute the euclidean distance
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms()),
                     "fitpeaksargs={}".format(fitpeaks_args),
                     "peaks={}".format(diamond_peaks)]

            prefix = diamond_file.rstrip(".h5").rstrip(".nxs").split("/")[-1]
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            print("Cachefile = {}".format(cachefile))
            if os.path.exists(cachefile):
                # Load cachefile
                print("Found FitPeaks cached result, loading '{}'".format(cachefile))
                use_cache = True
                Load(Filename=cachefile, OutputWorkspace="parameters")
                params = mtd["parameters"]

        # Perform peak fitting if we aren't caching, OR we are but the cache file doesn't exist yet
        if not use_cache:
            params, fiterrors = peakfitting(wksp, diamond_peaks, **fitpeaks_args)
            # Save result parameter wksp to cache file if we want to use caching
            if cache_dir:
                SaveNexusProcessed(InputWorkspace=params, Filename=cachefile)

        use_cache = False
        if cache_dir:
            cachefile = cachefile.replace(".nxs", ".txt")
            if os.path.exists(cachefile):
                # Load clustering cachefile
                print("Found gathered params cache, loading '{}'".format(cachefile))
                use_cache = True
                clustering_input = np.loadtxt(cachefile)
                new_mask = np.loadtxt(cachefile.replace(".txt", "_mask.txt"), dtype=int).tolist()

        if not use_cache:
            clustering_input, new_mask = gather_fitparameters(params, fitparams, None, diamond_peaks, thresholds)
            if cache_dir:
                np.savetxt(cachefile, clustering_input)
                np.savetxt(cachefile.replace(".txt", "_mask.txt"), new_mask, fmt='%i')

        if output_fittable:
            print("Exporting fit parameter table to '{}'".format(output_fittable))
            tablews = gathered_parameters_to_tablewksp("table", clustering_input, diamond_peaks, fitparams)
            SaveNexusProcessed(Filename=output_fittable, InputWorkspace=tablews)

        new_mask = get_detector_mask(wksp, new_mask)
        print("New mask (detector IDs) = {}".format(compress_ints(new_mask)))

        # Extract the first column to keep wsindex for later
        wsindex = clustering_input[..., 0]
        clustering_input = clustering_input[..., 1:]

        display_parameter_stats(clustering_input, diamond_peaks, fitparams)

        #fig, ax = plot_features(clustering_input, diamond_peaks, fitparams)

    model = None
    centroids = None
    if method[0] == "KMEANS":
        model = KMeans(n_clusters=num_groups).fit(clustering_input)
        centroids = model.cluster_centers_
        print("KMeans centroids: {}".format(centroids))
        print("KMeans inertia: {}".format(model.inertia_))
    elif method[0] == "DBSCAN":
        model = DBSCAN(eps=epsilon).fit(clustering_input)
    else:
        raise ValueError("Invalid grouping method '{}'. Must be KMEANS or DBSCAN".format(method[0]))

    labels = model.labels_
    unique_labels = np.unique(labels)
    skip_grouping = False  # flag to skip generate grouping file

    print("Labels: {}".format(labels))
    print("Unique labels (clusters) = {} ({})".format(len(unique_labels), unique_labels))
    if len(unique_labels) == 1 and unique_labels[0] == -1:
        # Print warning that the only labels found were noisy
        skip_grouping = True
        print("NOTE: only noisy clusters were found.. skipping grouping generation!")

    if method[1] == "ED":
        fig, ax = plot_features(clustering_input, diamond_peaks, fitparams, labels, centroids)

    fig, ax = plt.subplots()
    for i in range(len(wsindex)):
        wsindex[i] = wksp.getDetector(i).getID()
    ax.scatter(wsindex, labels+1, c=labels)
    ax.set_title("{}: {}".format(grouping_method, diamond_file.split("/")[-1]))
    ax.set_xlabel("pixel")
    ax.set_ylabel("cluster/grouping")

    # Export mask to file
    if new_mask is not None:
        print("Saving new mask to '{}'".format(output_mask))
        # new mask is the union of original + newly masked items
        if mask is not None:  # check in case no masking file was given
            new_mask = np.union1d(mask, new_mask)
        np.savetxt(output_mask, new_mask, fmt='%10i')

    # Export grouping based on clustering result. Use the input workspace as the donor
    if not skip_grouping:
        print("Generating grouping file '{}'".format(output_file))
        CreateGroupingWorkspace(InputWorkspace=wksp, OutputWorkspace="grouping")
        grouping = mtd["grouping"]
        for i in range(len(labels)):
            det_id = wksp.getDetector(i).getID()
            # Skip if label is -1 (DBSCAN labels these as noise)
            if labels[i] == -1:
                continue
            grouping.setY(det_id, [int(labels[i]+1)])  # Shift by 1 since group 0 is unused

        SaveDetectorsGrouping(InputWorkspace=grouping, OutputFile=output_file)

    plt.show(block=True)

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
