import json
import numpy as np
import os

import matplotlib.pyplot as plt
from mantid.simpleapi import \
    ConvertUnits, \
    CreateCacheFilename, \
    CreateGroupingWorkspace, \
    LoadEventNexus, \
    Load, \
    MaskSpectra, \
    Rebin, \
    RemoveMaskedSpectra, \
    SaveDetectorsGrouping, \
    SaveNexusProcessed, \
    mtd
from sklearn.cluster import KMeans, DBSCAN
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

from total_scattering.autogrouping.utils import \
    peakfitting, \
    gather_fitparameters, \
    gathered_parameters_to_tablewksp, \
    display_parameter_stats, \
    similarity_matrix_degelder, \
    similarity_matrix_crosscorr, \
    get_detector_mask
from total_scattering.reduction.total_scattering_reduction import compress_ints


def plot_features(data, peaks, colnames, labels=None, centroids=None):
    """
    Creates a grid plot of the fitted peak parameters against all other params
    for each peak
    :param data: array of parameters for each peak, see gather_fitparameters()
    :param peaks: array of diamond peak positions
    :param colnames: list of parameter names
    :param labels: optional labels array from clustering used to color data
    :param centroids: optional centroids from clustering algorithm
    :return:
    """
    n = len(colnames)
    npeaks = len(peaks)
    for p in range(npeaks):
        fig, ax = plt.subplots(n, n)

        fig.suptitle("peak {}".format(peaks[p]))

        # Plot each peak feature against its other features
        for i in range(n):
            for j in range(n):
                if labels is not None:
                    ax[i, j].scatter(data[..., p * n + i],
                                     data[..., p * n + j],
                                     c=labels, alpha=0.5)
                else:
                    ax[i, j].scatter(data[..., p * n + i],
                                     data[..., p * n + j])

                if centroids is not None:
                    ax[i, j].scatter(centroids[..., p * n + i],
                                     centroids[..., p * n + j],
                                     marker="X", c=np.unique(labels),
                                     s=200, alpha=0.75,
                                     edgecolors="red")
                ax[i, j].set_xlabel(colnames[i])
                ax[i, j].set_ylabel(colnames[j])

    return fig, ax


def get_cachename(prefix, cache_dir, props):
    filename, sig = CreateCacheFilename(
        OtherProperties=props, Prefix=prefix, CacheDir=cache_dir)
    return filename, sig


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


def autogrouping(config):
    diamond_file = get_key("DiamondFile", config)  # .nxs.h5 or .nxs file
    masking_file = get_key("MaskFile", config)

    grouping_method = get_key("GroupingMethod", config)
    method = get_grouping_method(grouping_method)

    # k parameter for KMeans clustering
    num_groups = int(get_key("NumberOutputGroups", config))
    # eps parameter for DBSCAN clustering
    epsilon = float(get_key("DBSCANEpsilon", config))
    # whether to use scikit standard scaler for ED method
    standard_scaling = get_key("StandardScaling", config)

    fitparams = get_key("FittingFunctionParameters", config).split(",")
    fitpeaks_args = get_key("FitPeaksArgs", config)
    diamond_peaks = np.asarray(
        get_key("DiamondPeaks", config).strip().split(","), dtype=float)
    if "DSpaceBin" in config:
        dspace_bin = get_key("DSpaceBin", config)
    else:
        dspace_bin = "0.001"
    thresholds = get_key("ParameterThresholds", config)
    for threshold in thresholds:
        thresholds[threshold] = tuple(
            float(x) for x in thresholds[threshold].strip("()").split(","))

    max_chi = get_key("FilterByChi2", config)
    if max_chi["Enable"]:
        max_chi = max_chi["Value"]
    else:
        max_chi = None

    wks_index_range = None
    if "WorkspaceIndexRange" in config:
        wks_index_range = get_key("WorkspaceIndexRange", config)

    cache_dir = ""
    if "CacheDir" in config:
        cache_dir = os.path.abspath(get_key("CacheDir", config))

    plots = ""
    if "Plots" in config:
        plots = get_key("Plots", config)

    output_file = os.path.abspath(get_key("OutputGroupingFile", config))
    output_mask = os.path.abspath(get_key("OutputMaskFile", config))

    output_fittable = ""
    if "OutputFitParamFile" in config:
        output_fittable = os.path.abspath(get_key("OutputFitParamFile",
                                                  config))

    if diamond_file.endswith(".nxs.h5"):
        wksp = LoadEventNexus(Filename=diamond_file,
                              CompressTolerance=0.01,
                              SpectrumList=wks_index_range)
    else:
        wksp = Load(Filename=diamond_file,
                    SpectrumList=wks_index_range)

    mask = None
    new_mask = None
    if masking_file:
        if masking_file.endswith(".txt") or masking_file.endswith(".out"):
            mask = np.loadtxt(masking_file, dtype=int)
            print("Loaded detector mask: {}".format(mask))
            # mask numbers are detector ids - so convert them to ws indices
            mask_ws = wksp.getIndicesFromDetectorIDs(mask.tolist())
            if len(mask_ws) == 0:
                # Skip masking if no ws indices are masked
                print("No detector masks applied to workspace indices, "
                      "skipping masking")
            else:
                # compress down to ranges
                mask_ws = compress_ints(mask_ws)
                print("Compressed mask of wsindex: {}".format(mask_ws))
                # mask spectra in workspace
                wksp = MaskSpectra(InputWorkspace=wksp,
                                   InputWorkspaceIndexType="WorkspaceIndex",
                                   InputWorkspaceIndexSet=mask_ws)
                wksp = RemoveMaskedSpectra(InputWorkspace=wksp)

    # Convert from TOF to dspacing
    wksp = ConvertUnits(InputWorkspace=wksp, Target="dSpacing",
                        EMode="Elastic")

    # Rebin (in d-space)
    dspace_bin = float(dspace_bin)
    wksp = Rebin(InputWorkspace=wksp, Params=(0.1, dspace_bin, 3.0))

    wsindex = list(range(wksp.getNumberHistograms()))
    clustering_input = None
    use_cache = False
    if method[1] == "DG":
        # Compute the deGelder similarity matrix
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms()),
                     "wksprange={}".format(wks_index_range)]
            prefix = diamond_file.rstrip(".h5").rstrip(
                ".nxs").split("/")[-1] + "_DG"
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            print("Checking for cachefile '{}'".format(cachefile))
            if os.path.exists(cachefile):
                # Load cached similarity matrix
                use_cache = True
                print(
                    "Found cached deGelder similarity matrix, "
                    "loading from '{}'".format(cachefile))
                clustering_input = np.loadtxt(cachefile)

        if not use_cache:
            clustering_input = 1 - similarity_matrix_degelder(wksp)
            if cache_dir:
                # Save matrix if using caching
                print("Caching deGelder similarity matrix to "
                      "'{}'".format(cachefile))
                np.savetxt(cachefile, clustering_input)
    elif method[1] == "CC":
        # Compute the cross-correlation matrix using a
        # pointwise squared difference method
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms()),
                     "wksprange={}".format(wks_index_range)]
            prefix = diamond_file.rstrip(".h5").rstrip(
                ".nxs").split("/")[-1] + "_CC"
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            print("Checking for cachefile '{}'".format(cachefile))
            if os.path.exists(cachefile):
                # Load cached similarity matrix
                use_cache = True
                print("Found cached pointwise crosscorr similarity matrix, "
                      "loading from '{}'".format(cachefile))
                clustering_input = np.loadtxt(cachefile)

        if not use_cache:
            clustering_input = 1 - similarity_matrix_crosscorr(wksp)
            if cache_dir:
                # Save matrix if using caching
                print(
                    "Caching pointwise crosscorr similarity matrix "
                    "to '{}'".format(cachefile))
                np.savetxt(cachefile, clustering_input)
    elif method[1] == "ED":
        # Compute the euclidean distance
        if cache_dir:
            props = ["filename={}".format(diamond_file),
                     "mask={}".format(masking_file),
                     "nhisto={}".format(wksp.getNumberHistograms()),
                     "wksprange={}".format(wks_index_range),
                     "fitpeaksargs={}".format(fitpeaks_args),
                     "peaks={}".format(diamond_peaks)]

            prefix = diamond_file.rstrip(".h5").rstrip(".nxs").split("/")[-1]
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            print("Cachefile = {}".format(cachefile))
            if os.path.exists(cachefile):
                # Load cachefile
                print("Found FitPeaks cached result, "
                      "loading '{}'".format(cachefile))
                use_cache = True
                Load(Filename=cachefile, OutputWorkspace="parameters")
                params = mtd["parameters"]

        # Perform peak fitting if we aren't caching, OR we are but the
        # cache file doesn't exist yet
        if not use_cache:
            params, fiterrors = peakfitting(
                wksp, diamond_peaks, **fitpeaks_args)
            # Save result parameter wksp to cache file if
            # we want to use caching
            if cache_dir:
                SaveNexusProcessed(InputWorkspace=params, Filename=cachefile)

        # After peak fitting, assemble clustering input to an array of
        # fit parameters for each peak
        use_cache = False
        if cache_dir:
            # cache properties change slightly since these are ones
            # that the fit peaks result do not depend upon
            # these are cached separately since this might take awhile
            # to generate depending on workspace size
            props.append("chifiltering={}".format(max_chi))
            props.append("thresholds={}".format(thresholds))
            cachefile, sig = get_cachename(prefix, cache_dir, props)
            cachefile = cachefile.replace(".nxs", ".txt")
            if os.path.exists(cachefile):
                # Load clustering cachefile
                print("Found gathered params cache, "
                      "loading '{}'".format(cachefile))
                use_cache = True
                clustering_input = np.loadtxt(cachefile)
                # ensure masking cache is present
                if not os.path.exists(cachefile.replace(".txt", "_mask.txt")):
                    print("could not find cached masked file, regenerating..")
                    use_cache = False
                else:
                    new_mask = np.loadtxt(cachefile.replace(
                        ".txt", "_mask.txt"), dtype=int).tolist()

        if not use_cache:
            # if we are not using caching, or we are but haven't
            # saved a cached result yet
            clustering_input, new_mask = \
                gather_fitparameters(params, fitparams, None, diamond_peaks,
                                     thresholds, max_chi)
            if cache_dir:
                np.savetxt(cachefile, clustering_input)
                np.savetxt(cachefile.replace(
                    ".txt", "_mask.txt"), new_mask, fmt='%i')

        if output_fittable:
            # convert the clustering input array to a TableWorkspace
            # so it can be used in Mantid
            print("Exporting fit parameter table to "
                  "'{}'".format(output_fittable))
            tablews = gathered_parameters_to_tablewksp(
                "table", clustering_input, diamond_peaks, fitparams)
            SaveNexusProcessed(Filename=output_fittable,
                               InputWorkspace=tablews)

        # convert workspace index mask to detector index mask
        new_mask = get_detector_mask(wksp, new_mask)
        print("New mask (detector IDs) = {}"
              "".format(compress_ints(new_mask)))

        # Extract the first column to keep wsindex for later
        wsindex = clustering_input[..., 0]
        clustering_input = clustering_input[..., 1:]

        if standard_scaling:
            clustering_input = \
                StandardScaler().fit_transform(clustering_input)

        display_parameter_stats(clustering_input, diamond_peaks,
                                fitparams)

    model = None
    centroids = None
    if method[0] == "KMEANS":
        model = KMeans(n_clusters=num_groups).fit(clustering_input)
        centroids = model.cluster_centers_
        print("KMeans iterations ran: {}".format(model.n_iter_))
        print("KMeans centroids: {}".format(centroids))
        print("KMeans inertia: {}".format(model.inertia_))

        if plots and plots["KMeans_Elbow"]:
            distortion = []
            for i in range(2, 12):
                kmean_tmp = KMeans(n_clusters=i).fit(clustering_input)
                distortion.append(kmean_tmp.inertia_)
            fig, ax = plt.subplots()
            ax.plot(range(2, 12), distortion, marker='o')
            ax.set_xlabel("number of clusters (k)")
            ax.set_ylabel("model inertia")
        if plots and plots["KMeans_Silhouette"]:
            scores = []
            for i in range(2, 12):
                kmean_tmp = KMeans(n_clusters=i).fit(clustering_input)
                scores.append(silhouette_score(
                    clustering_input, kmean_tmp.labels_))
            fig, ax = plt.subplots()
            ax.plot(range(2, 12), scores, marker='o')
            ax.set_xlabel("number of clusters (k)")
            ax.set_ylabel("silhouette score")

    elif method[0] == "DBSCAN":
        model = DBSCAN(eps=epsilon)
        if method[1] != "ED":
            model.set_params(metric="precomputed")
        model.fit(clustering_input)
    else:
        raise ValueError(
            "Invalid grouping method '{}'. Must be KMEANS or "
            "DBSCAN".format(method[0]))

    labels = model.labels_
    unique_labels = np.unique(labels)
    skip_grouping = False  # flag to skip generate grouping file

    print("Labels: {}".format(labels))
    print("Unique labels (clusters) = {} ({})".format(
        len(unique_labels), unique_labels))
    if len(unique_labels) == 1 and unique_labels[0] == -1:
        # Print warning that the only labels found were noisy
        skip_grouping = True
        print("NOTE: only noisy clusters were found.. skipping "
              "grouping generation!")

    if method[1] == "ED":
        if plots and plots["ED_Features"]:
            fig, ax = plot_features(
                clustering_input, diamond_peaks, fitparams,
                labels, centroids)

    if plots and plots["Grouping"]:
        fig, ax = plt.subplots()
        wsindex_init = []
        for item in wsindex:
            wsindex_init.append(wksp.getDetector(int(item)).getID())
        ax.scatter(wsindex_init, labels + 1, c=labels)
        ax.set_title("{}: {}".format(
            grouping_method, diamond_file.split("/")[-1]))
        ax.set_xlabel("pixel")
        ax.set_ylabel("cluster/grouping")

    if plots and plots["PCA"]:
        pca = PCA(n_components=2)
        pca_res = pca.fit_transform(clustering_input)
        print("PCA variation of components: {}".format(
            pca.explained_variance_ratio_))

        fig, ax = plt.subplots()
        ax.set_xlabel("PCA 1")
        ax.set_ylabel("PCA 2")
        ax.scatter(pca_res[:, 0], pca_res[:, 1], c=labels, alpha=0.5)
        if method[0] == "KMEANS":
            centroids_pca = pca.transform(centroids)
            ax.scatter(centroids_pca[:, 0], centroids_pca[:, 1],
                       marker='X', s=200, alpha=0.75,
                       edgecolors='red', c=np.unique(labels))

    # Export mask to file
    if new_mask is not None:
        print("Saving new mask to '{}'".format(output_mask))
        # new mask is the union of original + newly masked items
        if mask is not None:  # check in case no masking file was given
            new_mask = np.union1d(mask, new_mask)
        np.savetxt(output_mask, new_mask, fmt='%10i')

    # Export grouping based on clustering result. Use the input
    # workspace as the donor
    if not skip_grouping:
        print("Generating grouping file '{}'".format(output_file))
        instr_name = wksp.getInstrument().getName()
        print("Instrument - {}".format(instr_name))
        CreateGroupingWorkspace(InputWorkspace=wksp,
                                OutputWorkspace="grouping")
        grouping = mtd["grouping"]
        num_monitors = int(np.sum(wksp.detectorInfo().detectorIDs() < 0))
        for i in range(len(labels)):
            det_id = wksp.getDetector(int(wsindex[i])).getID()
            det_id = wksp.detectorInfo().indexOf(det_id) - num_monitors
            # Skip if label is -1 (DBSCAN labels these as noise)
            if labels[i] == -1:
                continue
            # Shift by 1 since group 0 is unused
            grouping.setY(det_id, [int(labels[i] + 1)])

        SaveDetectorsGrouping(InputWorkspace=grouping,
                              OutputFile=output_file)

    if plots:
        # block only if we have plotted something
        nplot = sum(plot for plot in plots.values())
        if nplot > 0:
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
    autogrouping(config)


if __name__ == "__main__":
    main()
