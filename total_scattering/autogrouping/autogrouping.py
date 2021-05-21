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
    FitPeaks, \
    LoadEventAndCompress, \
    Load, \
    Rebin, \
    SaveNexusProcessed, \
    mtd
from total_scattering.autogrouping.similarity import similarity_metric


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
    skiplist = []

    wsindex = np.asarray(paramws.column("wsindex"))
    wsindex_unique = np.unique(wsindex)
    n = len(wsindex_unique)
    for i in range(n):
        ind = int(np.searchsorted(wsindex, wsindex_unique[i]))
        # Add the ws index to list if fit result for ALL peaks are bad
        if get_badfitcount(ind, paramws, peaks, colnames, thresholds) != len(peaks):
            skiplist.append(i)

    return skiplist


def get_cachename(prefix, cache_dir, props):
    filename, sig = CreateCacheFilename(OtherProperties=props, Prefix=prefix, CacheDir=cache_dir)
    return filename, sig


def peakfitting(wksp: EventWorkspace, peaks: np.ndarray, fitfunction, param_names, param_values):
    # Create the peak windows from the diamond peak positions
    peakwindows = diagnostics.get_peakwindows(peaks)

    # Convert from TOF to dspacing
    wksp = ConvertUnits(InputWorkspace=wksp, Target="dSpacing", EMode="Elastic")

    # Perform multiple peak fitting
    output = FitPeaks(InputWorkspace=wksp,
                      PeakFunction=fitfunction,
                      PeakParameterNames=param_names,
                      PeakParameterValues=param_values,
                      RawPeakParameters=True,
                      HighBackground=False,
                      ConstrainPeakPositions=False,
                      MinimumPeakHeight=3,
                      PeakCenters=peaks,
                      FitWindowBoundaryList=peakwindows,
                      FittedPeaksWorkspace='fitted',
                      OutputPeakParametersWorkspace='parameters',
                      OutputParameterFitErrorsWorkspace='fiterrors')

    return 'parameters', 'fiterrors'


def gather_fitparameters(paramws: TableWorkspace, cols, mask,
                         peaks, thresholds=dict()):
    '''Generate an array of peak fitting parameters over all spectra from FitPeaks results'''

    paramws = mtd[str(paramws)]

    # The fit function parameters to gather
    nprops = len(cols) * len(peaks) + 1

    wsindex = paramws.column("wsindex")
    wsindex_unique = np.unique(wsindex)

    fitlist = get_goodfits(paramws, peaks, cols, thresholds)
    n = len(fitlist)

    result = np.ndarray(shape=(n, nprops))

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

    return result


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
            param = parameters[..., i * len(peaks) + j + 1]
            print("Param {}-{}:".format(cols[j], peaks[i]))
            print("   MAX = {}".format(np.max(param)))
            print("   MIN = {}".format(np.min(param)))
            print("   AVG = {}".format(np.mean(param)))
            print("   STD = {}".format(np.std(param)))
            print("")


def plot_features(data, peaks, colnames, labels=None):
    n = len(colnames)
    # n = 2
    fig, ax = plt.subplots(n, n)

    # Plot each peak feature against its other features
    for i in range(n):
        for j in range(n):
            if labels is not None:
                ax[i, j].scatter(data[..., i], data[..., j], c=labels, alpha=0.5)
            else:
                ax[i, j].scatter(data[..., i], data[..., j])
            ax[i, j].set_xlabel(colnames[i])
            ax[i, j].set_ylabel(colnames[j])

    return fig, ax


def similarity_matrix_degelder(wksp, mask_wksp):
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

    frac = int(n/10)  # get a tenth of total spectra for progress updates

    result = np.zeros(shape=(n, n))
    for i in range(n):
        for j in range(i, n):
            result[i][j] = sm.de_gelder_similarity(y[i], y[j])

        if i % frac == 0:
            print("-- {}/{} spectra".format(i, n))

    # Zero out any nans since this causes problems with clustering
    result[np.isnan(result)] = 0.0
    return result


def similarity_matrix_crosscorr(wksp, mask_wksp):
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

    fitfunction = get_key("FittingFunction", config)
    fitparams = get_key("FittingFunctionParameters", config).split(",")

    fitpeaks_names = get_key("InitialParameterNames", config)
    fitpeaks_values = get_key("InitialParameterValues", config)
    diamond_peaks = np.asarray(get_key("DiamondPeaks", config).strip().split(","), dtype=float)
    thresholds = get_key("ParameterThresholds", config)
    for threshold in thresholds:
        thresholds[threshold] = tuple(float(x) for x in thresholds[threshold].strip("()").split(","))

    cache_dir = ""
    if "CacheDir" in config:
        cache_dir = get_key("CacheDir", config)

    output_file = get_key("OutputGroupingFile", config)
    output_mask = get_key("OutputMaskFile", config)

    if diamond_file.endswith(".nxs.h5"):
        wksp = LoadEventAndCompress(Filename=diamond_file, FilterBadPulses=0)
    else:
        wksp = Load(Filename=diamond_file)
    mask_wksp = None
    if masking_file:
        mask_wksp = Load(Filename=masking_file)

    # Rebin (in TOF)
    # TODO: Double check if this is needed, and if this should only be done for ED case
    wksp = Rebin(InputWorkspace=wksp, Params=(300, -0.001, 16666.7))

    wsindex = range(wksp.getNumberHistograms())
    clustering_input = None
    use_cache = False
    if method[1] == "DG":
        # Compute the deGelder similarity matrix
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms())]
            prefix = diamond_file.rstrip(".nxs").rstrip(".h5").split("/")[-1] + "_DG"
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            print("Checking for cachefile '{}'".format(cachefile))
            if os.path.exists(cachefile):
                # Load cached similarity matrix
                use_cache = True
                print("Found cached deGelder similarity matrix, loading from '{}'".format(cachefile))
                clustering_input = np.loadtxt(cachefile)

        if not use_cache:
            clustering_input = similarity_matrix_degelder(wksp, mask_wksp)
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
            prefix = diamond_file.rstrip(".nxs").rstrip(".h5").split("/")[-1] + "_CC"
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            print("Checking for cachefile '{}'".format(cachefile))
            if os.path.exists(cachefile):
                # Load cached similarity matrix
                use_cache = True
                print("Found cached pointwise crosscorr similarity matrix, loading from '{}'".format(cachefile))
                clustering_input = np.loadtxt(cachefile)

        if not use_cache:
            clustering_input = similarity_matrix_crosscorr(wksp, mask_wksp)
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
                     "function={}".format(fitfunction),
                     "param_names={}".format(fitpeaks_names),
                     "param_vals={}".format(fitpeaks_values),
                     "peaks={}".format(diamond_peaks)]

            prefix = diamond_file.rstrip(".nxs").rstrip(".h5").split("/")[-1]
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
            params, fiterrors = peakfitting(wksp, diamond_peaks, fitfunction, fitpeaks_names, fitpeaks_values)
            # Save result parameter wksp to cache file if we want to use caching
            if cache_dir:
                SaveNexusProcessed(InputWorkspace=params, Filename=cachefile)

        clustering_input = gather_fitparameters(params, fitparams, None, diamond_peaks, thresholds)

        # Extract the first column to keep wsindex for later
        wsindex = clustering_input[..., 0]
        clustering_input = clustering_input[..., 1:]

        display_parameter_stats(clustering_input, diamond_peaks, fitparams)

        fig, ax = plot_features(clustering_input, diamond_peaks, fitparams)

    model = None
    if method[0] == "KMEANS":
        model = KMeans(n_clusters=num_groups).fit(clustering_input)
        centroids = model.cluster_centers_
        print("KMeans centroids: {}".format(centroids))
        print("KMeans inertia: {}".format(model.inertia_))
    elif method[0] == "DBSCAN":
        model = DBSCAN(eps=1.0).fit(clustering_input)
    else:
        raise ValueError("Invalid grouping method '{}'. Must be KMEANS or DBSCAN".format(method[0]))

    labels = model.labels_
    print("Labels: {}".format(labels))
    print("Unique labels = {}".format(np.unique(labels)))

    if method[1] == "ED":
        fig, ax = plot_features(clustering_input, diamond_peaks, fitparams, labels)

    fig, ax = plt.subplots()
    ax.scatter(wsindex, labels, c=labels)
    ax.set_xlabel("pixel")
    ax.set_ylabel("cluster/grouping")

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
