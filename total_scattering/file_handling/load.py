import logging
from mantid import mtd
from mantid.simpleapi import \
    AlignAndFocusPowderFromFiles, \
    ApplyDiffCal, \
    ConvertUnits, \
    DeleteWorkspace, \
    DiffractionFocussing, \
    Divide, \
    Load, \
    LoadDetectorsGroupingFile, \
    LoadDiffCal, \
    MaskDetectors, \
    MultipleScatteringCorrection, \
    NormaliseByCurrent, \
    PDDetermineCharacterizations, \
    PDLoadCharacterizations, \
    PropertyManagerDataService, \
    SetSample, \
    Rebin, \
    RebinToWorkspace, \
    CreateGroupingWorkspace, \
    SaveDetectorsGrouping, \
    GroupDetectors, \
    LoadNexus, \
    SaveNexusProcessed, \
    CloneWorkspace, \
    Plus, \
    MaskDetectors
from mantid.utils import absorptioncorrutils
from sklearn.cluster import KMeans
import numpy as np
import os
import re
import hashlib
import base64

_shared_shape_keys = ["Shape", "Height", "Center"]
required_shape_keys = {
    "FlatPlate": _shared_shape_keys + ["Width", "Thick", "Angle"],
    "Cylinder": _shared_shape_keys + ["Radius"],
    "HollowCylinder": _shared_shape_keys + ["InnerRadius", "OuterRadius"]
}


def load(ws_name, input_files, group_wksp,
         facility=None, instr_name=None, ipts=None, group_num=None,
         geometry=None, chemical_formula=None, mass_density=None,
         absorption_wksp='', out_group_dict=None,
         qparams='0.01,0.001,40.0',
         wlparams='0.06,0.0001,2.98',
         auto_red=False,
         group_all_file=None,
         sam_files=None,
         re_cache=False,
         cache_dir=None,
         **align_and_focus_args):
    '''Routine for loading workspace'''

    # Figure out the list of run number
    run_list = list()
    potential_cache_bn = list()
    for item in input_files.split(","):
        if "/" in item or "\\" in item:
            file_n = os.path.basename(item)
            reg_test = re.search('.*_[0-9]+\\.nxs\\.h5', file_n)
            if reg_test is None:
                run_list.append("-1")
                hash_obj = hashlib.sha256(file_n.encode())
                hash_str = hash_obj.hexdigest()
                vt = hash_str.encode()
                s_hash_str = base64.urlsafe_b64encode(vt).decode()[:12]
                potential_cache_bn.append(s_hash_str)
            else:
                run_list.append(reg_test.group(0).split(".")[0].split("_")[1])
                potential_cache_bn.append(None)
        else:
            run_list.append(item.split("_")[1])
            potential_cache_bn.append(None)

    str_tmp = sam_files.split(',')[0]
    if "/" in str_tmp or "\\" in str_tmp:
        sfile_n = os.path.basename(sam_files.split(',')[0])
        hash_obj = hashlib.sha256(sfile_n.encode())
        hash_str = hash_obj.hexdigest()
        vt = hash_str.encode()
        sf_cfn_part = base64.urlsafe_b64encode(vt).decode()[:12]
    else:
        sf_cfn_part = str_tmp

    if ipts is not None:
        if cache_dir is None:
            cache_dir = os.path.join("/" + facility,
                                     instr_name,
                                     ipts,
                                     "shared",
                                     "autoreduce",
                                     "cache")
        else:
            cache_dir = os.path.abspath(cache_dir)
        os.makedirs(cache_dir, exist_ok=True)

    # For autoreduction, we want to fill the low-q part with the asymptotic
    # value according to the paper by D. Keen,
    # J. Appl. Cryst. (2001). 34, 172-177
    # to guarantee a robust Fourier transform automatically. In such a
    # situation, we have to change the binning parameter to extend to the
    # very low-q region and later on (at the end of the main program) we will
    # fill the low-q part (below the original Qmin as specified in the `qparams`
    # parameter) with the asymptotic value. Here, the temporary Qmin value is
    # set to the provided Q-space interval.
    if auto_red:
        qminmin = qparams.split(",")[1]
        qmaxmax = qparams.split(",")[2]
        qparams_use = f"{qminmin},{qminmin},{qmaxmax}"
    else:
        qparams_use = qparams

    if group_wksp is None:
        # Given the current implementation mechanism, there are two cases when
        # the input `group_wksp` will be `None`,
        #   - no absorption calculation was ever initiated
        #       > In this case, we can cache and load the sum
        #   - manual grouping of detectors was specified
        #       > In this case, no sense for the summed caching

        # Figure out a unique name for the summed cache file
        # The codes here were originating from GPT3.5-turbo API
        hash_obj = hashlib.sha256(str(run_list).encode())
        hash_str = hash_obj.hexdigest()
        short_hash_str = base64.urlsafe_b64encode(hash_str.encode()).decode()[:12]

        cache_sf_bn = f"{instr_name}_mts_summed_"
        if auto_red:
            cache_sf_bn += f"{short_hash_str}_sgb"
        else:
            cache_sf_bn += f"{short_hash_str}"
        cache_sf_bn += f"_{ws_name}_{sam_files.split(',')[0]}.nxs"
        if ipts is not None:
            cache_sf_fn = os.path.join(
                cache_dir,
                cache_sf_bn
            )

        cond_1 = ipts is not None
        cond_2 = os.path.isfile(cache_sf_fn)
        cond_3 = absorption_wksp == ''
        if cond_1 and cond_2 and cond_3:
            LoadNexus(Filename=cache_sf_fn, OutputWorkspace=ws_name)
        else:
            if auto_red:
                align_and_focus_args["GroupFilename"] = group_all_file

            if "GroupFilename" in align_and_focus_args:
                proc_group = "proc_group"
                LoadDetectorsGroupingFile(
                    InputFile=align_and_focus_args["GroupFilename"],
                    OutputWorkspace="proc_group"
                )
            else:
                proc_group = None

            out_group = align_and_focus_args["GroupingWorkspace"]

            for run_i, run in enumerate(run_list):
                if run == "-1":
                    cache_f_bn = potential_cache_bn[run_i]
                else:
                    cache_f_bn = f"{instr_name}_{run}"

                if auto_red:
                    cache_f_bn += f"_mts_no_subg_sgb_{ws_name}"
                    cache_f_bn += f"_{sf_cfn_part}.nxs"
                else:
                    cache_f_bn += f"_mts_no_subg_{ws_name}"
                    cache_f_bn += f"_{sf_cfn_part}.nxs"

                if ipts is not None:
                    cache_f_fn = os.path.join(
                        cache_dir,
                        cache_f_bn
                    )
                if ipts is not None and os.path.isfile(cache_f_fn):
                    cache_f_exist = True
                else:
                    cache_f_exist = False

                if cache_f_exist and absorption_wksp == '':
                    wksp_tmp = "wksp_tmp_qrb"
                    wksp_tmp = LoadNexus(
                        OutputWorkspace=wksp_tmp,
                        Filename=cache_f_fn
                    )
                else:
                    wksp_tmp = "wksp_tmp"

                    tmin_tmp = align_and_focus_args["TMin"]
                    tmax_tmp = align_and_focus_args["TMax"]
                    params = f"{tmin_tmp}, -0.0008, {tmax_tmp}"

                    if absorption_wksp == '':
                        proc_group_in = proc_group
                    else:
                        proc_group_in = None

                    align_focus_mts(
                        wksp_tmp,
                        instr_name,
                        input_files.split(",")[run_i],
                        align_and_focus_args["CalFilename"],
                        params,
                        group_wksp_in=proc_group_in,
                        pres_events=align_and_focus_args["PreserveEvents"]
                    )
                    Rebin(
                        InputWorkspace=wksp_tmp,
                        OutputWorkspace="wksp_tmp_qrb",
                        Params=qparams_use,
                        PreserveEvents=align_and_focus_args["PreserveEvents"]
                    )
                    if ipts is not None:
                        SaveNexusProcessed(
                            InputWorkspace="wksp_tmp_qrb",
                            Filename=cache_f_fn,
                            Title=f"{run}_cached_no_abs",
                            WorkspaceIndexList=range(
                                mtd[wksp_tmp].getNumberHistograms()
                            )
                        )

                # Accumulate individual files
                if run_i == 0:
                    CloneWorkspace(InputWorkspace="wksp_tmp_qrb",
                                   OutputWorkspace=ws_name + "_tmp")
                else:
                    Plus(LHSWorkspace="wksp_tmp_qrb",
                         RHSWorkspace=wksp_tmp,
                         OutputWorkspace=ws_name + "_tmp")
                    DeleteWorkspace(Workspace="wksp_tmp_qrb")

            if absorption_wksp != '':
                ConvertUnits(
                    InputWorkspace=absorption_wksp,
                    OutputWorkspace=absorption_wksp,
                    Target="MomentumTransfer",
                    EMode="Elastic")
                Rebin(InputWorkspace=absorption_wksp,
                      OutputWorkspace="absorption_wksp_rb",
                      Params=qparams_use,
                      PreserveEvents=align_and_focus_args["PreserveEvents"])
                Divide(LHSWorkspace=mtd[ws_name + "_tmp"],
                       RHSWorkspace="absorption_wksp_rb",
                       OutputWorkspace=ws_name,
                       AllowDifferentNumberSpectra=True)

                DeleteWorkspace(Workspace=absorption_wksp)
                DeleteWorkspace(Workspace="absorption_wksp_rb")

            DiffractionFocussing(
                InputWorkspace=ws_name + "_tmp",
                OutputWorkspace=ws_name,
                GroupingWorkspace=out_group
            )

            DeleteWorkspace(Workspace=wksp_tmp)
            DeleteWorkspace(Workspace=ws_name + "_tmp")

            if ipts is not None and absorption_wksp == '':
                SaveNexusProcessed(ws_name, cache_sf_fn, Title="cache_summed")
    else:
        for run_i, run in enumerate(run_list):
            if run == "-1":
                cache_f_bn = potential_cache_bn[run_i]
            else:
                cache_f_bn = f"{instr_name}_{run}_mts_subg"

            cache_f_bn += f"_{ws_name}"
            cache_f_bn += f"_{sf_cfn_part}.nxs"
            if ipts is not None:
                cache_f_fn = os.path.join(
                    cache_dir,
                    cache_f_bn
                )
            if os.path.isfile(cache_f_fn):
                cache_f_exist = True
            else:
                cache_f_exist = False

            # `group_num == 0` means no regeneration of grouping
            # was ever initialized and only in such cases will we
            # move ahead to load in the existing cache file.
            if cache_f_exist and group_num == 0:
                wksp_tmp = "wksp_tmp_qrb"
                LoadNexus(OutputWorkspace=wksp_tmp, Filename=cache_f_fn)
            else:
                wksp_tmp = "wksp_tmp"

                tmin_tmp = align_and_focus_args["TMin"]
                tmax_tmp = align_and_focus_args["TMax"]
                params = f"{tmin_tmp}, -0.0008, {tmax_tmp}"
                align_focus_mts(
                    wksp_tmp,
                    instr_name,
                    input_files.split(",")[run_i],
                    align_and_focus_args["CalFilename"],
                    params,
                    group_wksp_in=group_wksp,
                    pres_events=align_and_focus_args["PreserveEvents"]
                )
                Rebin(InputWorkspace=wksp_tmp,
                      OutputWorkspace="wksp_tmp_qrb",
                      Params=qparams_use,
                      PreserveEvents=align_and_focus_args["PreserveEvents"])
                SaveNexusProcessed(
                    InputWorkspace="wksp_tmp_qrb",
                    Filename=cache_f_fn,
                    Title=f"{run}_cached",
                    WorkspaceIndexList=range(
                        mtd[wksp_tmp].getNumberHistograms()))

            # Accumulate individual files
            if run_i == 0:
                CloneWorkspace(InputWorkspace="wksp_tmp_qrb",
                               OutputWorkspace=ws_name)
            else:
                Plus(LHSWorkspace=ws_name,
                     RHSWorkspace="wksp_tmp_qrb",
                     OutputWorkspace=ws_name)
                DeleteWorkspace(Workspace="wksp_tmp_qrb")

        if absorption_wksp != '':
            ConvertUnits(
                InputWorkspace=absorption_wksp,
                OutputWorkspace=absorption_wksp,
                Target="MomentumTransfer",
                EMode="Elastic")
            Rebin(InputWorkspace=absorption_wksp,
                  OutputWorkspace="absorption_wksp_rb",
                  Params=qparams_use,
                  PreserveEvents=align_and_focus_args["PreserveEvents"])
            Divide(LHSWorkspace=mtd[ws_name],
                   RHSWorkspace="absorption_wksp_rb",
                   OutputWorkspace=ws_name,
                   AllowDifferentNumberSpectra=True)

        if auto_red:
            min_min = 0
            max_max = mtd[ws_name].getNumberHistograms() - 1
            g_pattern = f"{min_min}-{max_max}"
        else:
            spec_map = dict()
            for i in range(mtd[ws_name].getNumberHistograms()):
                spec_id = mtd[ws_name].getSpectrum(i).getSpectrumNo()
                spec_map[spec_id] = i
            g_pattern = ""

            for _, item in out_group_dict.items():
                min_tmp = item[0]
                max_tmp = item[1]
                list_tmp = [num for num in spec_map.keys() if num >= min_tmp]
                min_min = spec_map[min(list_tmp)]
                list_tmp = [num for num in spec_map.keys() if num <= max_tmp]
                max_max = spec_map[max(list_tmp)]

                g_pattern += f"{min_min}-{max_max},"

            g_pattern = g_pattern[:-1]

        GroupDetectors(InputWorkspace=ws_name,
                       OutputWorkspace=ws_name,
                       GroupingPattern=g_pattern)

    NormaliseByCurrent(
        InputWorkspace=ws_name,
        OutputWorkspace=ws_name,
        RecalculatePCharge=True)

    if group_wksp is None:
        Rebin(
            InputWorkspace=ws_name,
            OutputWorkspace=ws_name,
            Params=qparams_use,
            PreserveEvents=align_and_focus_args["PreserveEvents"]
        )

    if geometry and chemical_formula and mass_density:
        set_sample(ws_name, geometry, chemical_formula, mass_density)

    return ws_name


def align_focus_mts(out_wksp,
                    instr_name,
                    file_name,
                    cal_file_name,
                    tof_bin_params,
                    group_wksp_in=None,
                    pres_events=True):
    """The MantidTotalScattering internal version of the align and focus
    algorithm. Simple enough but does the job just as what it should do.
    """
    wksp_proc = Load(file_name)

    LoadDiffCal(
        InstrumentName=instr_name,
        Filename=cal_file_name,
        WorkspaceName="calib_wksp",
        TofMin=float(tof_bin_params.split(",")[0])
    )

    ApplyDiffCal(
        InstrumentWorkspace="wksp_proc",
        CalibrationWorkspace="calib_wksp_cal"
    )

    MaskDetectors(
        Workspace="wksp_proc",
        MaskedWorkspace="calib_wksp_mask"
    )

    Rebin(
        InputWorkspace="wksp_proc",
        OutputWorkspace="wksp_proc_rebin",
        Params=tof_bin_params,
        PreserveEvents=True
    )

    # When the DIFA term is non zero in the calibration table,
    # the `ConvertUnits` algorithm cannot convert from TOF to
    # Q properly. As a workaround, we can convert to d first,
    # then to Q.
    ConvertUnits(
        InputWorkspace="wksp_proc_rebin",
        OutputWorkspace="wksp_proc_rebin_d",
        Target="dSpacing"
    )

    if group_wksp_in is not None:
        DiffractionFocussing(
            InputWorkspace="wksp_proc_rebin_d",
            OutputWorkspace="wksp_proc_focus",
            GroupingWorkspace=group_wksp_in
        )
        wksp_to_convert = "wksp_proc_focus"
    else:
        wksp_to_convert = "wksp_proc_rebin_d"

    ConvertUnits(
        InputWorkspace=wksp_to_convert,
        OutputWorkspace=out_wksp,
        Target="MomentumTransfer"
    )

    DeleteWorkspace(Workspace="wksp_proc")
    DeleteWorkspace(Workspace="wksp_proc_rebin")
    DeleteWorkspace(Workspace=wksp_to_convert)

    return


def set_sample(ws_name, geometry=None, chemical_formula=None,
               mass_density=None):
    '''Sets sample'''
    if 'Center' not in geometry:
        geometry.update({'Center': [0., 0., 0., ]})
    if "Shape" not in geometry:
        geometry.update({'Shape': 'Cylinder'})
    geometry = configure_geometry(geometry)
    SetSample(
        InputWorkspace=ws_name,
        Geometry=geometry,
        Material={
            'ChemicalFormula': chemical_formula,
            'SampleMassDensity': mass_density})


def configure_geometry(geo):
    '''Configure geometry'''
    new_geo = dict()
    shape = geo['Shape'].lower().replace(" ", "")
    if shape == 'cylinder':
        new_geo = add_required_shape_keys(geo, "Cylinder")

    if shape == 'hollowcylinder':
        new_geo = add_required_shape_keys(geo, "HollowCylinder")
        if 'Radius' in geo:
            new_geo['OuterRadius'] = geo['Radius']
        if 'Radius2' in geo:
            new_geo['InnerRadius'] = geo['Radius2']

    if shape == 'flatplate':
        new_geo = add_required_shape_keys(geo, "FlatPlate")
    return new_geo


def add_required_shape_keys(mydict, shape):
    '''Add the required 'Shape' key'''
    new_dict = dict()
    for key in required_shape_keys[shape]:
        if key not in mydict:
            new_dict[key] = None
        else:
            new_dict[key] = mydict[key]
    new_dict['Shape'] = shape
    return new_dict


def mask_generator(group_wksp, abs_ref_file):
    '''Generate a mask list given the input file containing
    the list of detectors that we don't want to mask, together with
    a group workspace. We don't need much information from the
    provided grouping workspace beyond all those detector IDs. In
    fact, this input workspace can be any workspace containing the
    instrument geometry.
    '''
    with open(abs_ref_file, "r") as f:
        ref_list = [int(line.strip()) for line in f.readlines()]

    num_dets = group_wksp.getNumberHistograms()
    all_det_ids = list()
    for i in range(num_dets):
        all_det_ids.append(group_wksp.getDetector(i).getID())

    return [i for i in all_det_ids if i not in ref_list]


def abs_grouping(sam_abs_ws,
                 con_abs_ws,
                 num_clusters,
                 group_out_file,
                 group_ref_det_out_file,
                 subgroup_index_file):

    '''The routine for clustering absorption spectra. The clustering
    will be performed for each of the physical banks separately as this will
    make the implementation at the data reduction stage relatively more
    straightforward in terms of cooperating with the grouped absorption spectra
    scenario.

    :param sam_abs_ws: sample abs workspace containing spectra for all pixels
    :type sam_abs_ws: str
    :param con_abs_ws: container abs workspace containing spectra for all pixels
    :type con_abs_ws: str
    :param num_clusters: number of clusters in total
    :type num_clusters: int
    :param group_out_file: output for the generated group per abs similarity
    :type group_out_file: str
    :param group_ref_det_out_file: reference spectrum for each generated group
    :type group_ref_det_out_file: str
    :param subgroup_index_file: file for the belonging of abs group to banks
    :type subgroup_index_file: str
    '''

    # Rebin the spectra in prep for the clustering
    if isinstance(sam_abs_ws, str):
        xmin = (int(mtd[sam_abs_ws].readX(0)[0] / 0.05) + 1.) * 0.05
    else:
        xmin = (int(sam_abs_ws.readX(0)[0] / 0.05) + 1.) * 0.05
    Rebin(InputWorkspace=sam_abs_ws,
          OutputWorkspace="abs_pbp",
          Params="{0:4.2F},0.05,3.0".format(xmin))
    
    # Figure out the number of spectra and monitors
    num_spectra = mtd['abs_pbp'].getNumberHistograms()
    num_monitors = int(np.sum(mtd['abs_pbp'].detectorInfo().detectorIDs() < 0))
    
    (cw_out) = CreateGroupingWorkspace(InputWorkspace="abs_pbp",
                                       GroupDetectorsBy="Group")
    output_group = cw_out[0]
    spectra_count = cw_out[1]
    group_count = cw_out[2]

    # Figure out the number of sub groups in each physical bank to meet the
    # overall number of groups to identify. Then we will save the index of the
    # sub groups and their belongings to the physical banks into a file for
    # later use at the reduction stage.
    sub_num_clusters = int(num_clusters / group_count)
    with open(subgroup_index_file, "w") as f:
        f.write("# Start & Stop index of subgroups\n")
        for i in range(group_count):
            start_index = sub_num_clusters * i + 1
            stop_index = sub_num_clusters * (i + 1)
            f.write("{0:<10s}{1:<10d}{2:<10d}\n".format("Bank-" + str(i),
                                                        start_index,
                                                        stop_index))

    # Clustering the absorption spectra using the KMEANS algorithm
    print("[Info] Clustering the absorption spectra for all detectors...")
    all_spectra_id = dict()
    all_spectra_all = list()
    labels = list()
    for j in range(group_count):
        all_spectra = list()
        for i in range(num_spectra):
            # Here we only pick up those spectra belonging to a certain physical
            # bank in the loop.
            if output_group.readY(i) == float(j + 1):
                y_tmp = mtd['abs_pbp'].readY(i)
                all_spectra.append(y_tmp)
                # Save the correspondence between the spectra index and the
                # detector ID which will be used later for creating the grouping
                # workspace.
                all_spectra_id[i] = mtd['abs_pbp'].getDetector(i).getID()
        all_spectra_all.extend(all_spectra)

        X = np.array(all_spectra)
        clustering = KMeans(n_clusters=sub_num_clusters, n_init='auto').fit(X)
    
        # For each physical bank, the identified clustering label will be
        # starting from 0. We need to connect continuously the labelling of the
        # subgroups belonging to each of the physical banks one after another.
        labels_tmp = [x + j * sub_num_clusters for x in clustering.labels_]
        labels.extend(labels_tmp)
        msg_tmp = "[Info] Done with the clustering of absorption "
        msg_tmp += f"spectra for bank-{j + 1}."
        print(msg_tmp)

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    
    # Save the identified group based on the absorption spectra similarity to
    # grouping file.
    CreateGroupingWorkspace(InputWorkspace="abs_pbp",
                            OutputWorkspace="grouping")
    grouping = mtd["grouping"]
    for i, label in enumerate(labels):
        if label == -1:
            continue

        # Given the detector ID, figure out the spectra index
        # It should be noticed that the index of a certain detector ID in the
        # returned list from `detectorInfo()` method takes the monitors into
        # account, but `setY()` method takes the spectrum index as input which
        # does not take into account of the monitors. So, we need the following
        # operations to remove number of monitors from the obtained spectra
        # index.
        det_id = all_spectra_id[i]
        det_id = mtd['abs_pbp'].detectorInfo().indexOf(det_id) - num_monitors
        grouping.setY(det_id, [int(label + 1)])
    
    SaveDetectorsGrouping(InputWorkspace=grouping,
                          OutputFile=group_out_file)

    # Group the absorption spectra within each of the identified groups
    # following the average scheme and save the averaged absorption spectra to
    # an external file.
    GroupDetectors(InputWorkspace=sam_abs_ws,
                   OutputWorkspace="sam_abs_grouped", Behaviour="Average",
                   CopyGroupingFromWorkspace=grouping)
    sam_abs_grouped = mtd['sam_abs_grouped']
    if con_abs_ws != "":
        GroupDetectors(InputWorkspace=con_abs_ws,
                       OutputWorkspace="con_abs_grouped",
                       Behaviour="Average",
                       CopyGroupingFromWorkspace=grouping)
        con_abs_grouped = mtd['con_abs_grouped']
    else:
        con_abs_grouped = ""

    # Here we want to check the similarity of the spectra within each of the
    # identified group and see how different they are.
    print("[Info] Check the similarity of absorption spectra within groups.")
    ref_id = list()
    for i in range(group_count * sub_num_clusters):
        ref_id.append(list(labels).index(i))
    with open(group_ref_det_out_file, "w") as f:
        for i, item in enumerate(ref_id):
            if i == len(ref_id) - 1:
                f.write("{0:d}".format(all_spectra_id[item]))
            else:
                f.write("{0:d}\n".format(all_spectra_id[item]))
    
    perct_max_tmp = list()
    for i, group in enumerate(labels):
        perct_tmp = list()
        for j in range(len(all_spectra_all[i])):
            top = abs(all_spectra_all[i][j] - all_spectra_all[ref_id[group]][j])
            bottom = all_spectra_all[ref_id[group]][j]
            perct_tmp.append(top / bottom * 100.)
        perct_max_tmp.append(max(perct_tmp))
    perct_max = max(perct_max_tmp)
    
    print(f"[Info] # of clusters: {n_clusters_}")
    print(f"[Info] # of noise: {n_noise_}")
    print("[Info] Maximum percentage difference: {0:6.4F}%".format(perct_max))

    return sam_abs_grouped, con_abs_grouped, grouping


def create_absorption_wksp(filename, abs_method, geometry, material,
                           environment=None, props=None,
                           characterization_files=None,
                           ms_method=None,
                           elementsize=1.0,  # mm
                           con_elementsize=1.0,  # mm
                           group_wksp_in=None,
                           num_groups=0,
                           group_out_file=None,
                           group_ref_det_out_file=None,
                           sg_index_f=None,
                           **align_and_focus_args):
    group_wksp_out = group_wksp_in

    '''Create absorption workspace'''
    if abs_method is None:
        return '', '', ''

    abs_input = Load(filename, MetaDataOnly=True)

    # If no run characterization properties given, load any provided files
    if not props and characterization_files:
        msg = "No props were given, but determining from characterization files"
        print(msg)

        charfile = characterization_files
        # Reduce to a string if multiple files were provided
        if isinstance(charfile, list):
            charfile = ','.join(characterization_files)
        charTable = PDLoadCharacterizations(Filename=charfile)
        chars = charTable[0]

        # Create the properties for the absorption workspace
        # NOTE
        # WaveLengthLogNames used here will be the standard default one in the
        # future, however let's keep them in until Mantid_v6.3 comes out
        PDDetermineCharacterizations(
            InputWorkspace=abs_input,
            Characterizations=chars,
            ReductionProperties="__absreductionprops",
            WaveLengthLogNames="LambdaRequest,lambda,skf12.lambda,"
                               "BL1B:Det:TH:BL:Lambda,freq"
            )
        props = PropertyManagerDataService.retrieve("__absreductionprops")

    # If neither run characterization properties or files, guess from input
    if not (props and characterization_files):
        msg = ("No props or characterizations were given, "
               "determining props from input file")
        print(msg)
        # NOTE
        # WaveLengthLogNames used here will be the standard default one in the
        # future, however let's keep them in until Mantid_v6.3 comes out
        PDDetermineCharacterizations(
            InputWorkspace=abs_input,
            ReductionProperties="__absreductionprops",
            WaveLengthLogNames="LambdaRequest,lambda,skf12.lambda,"
                               "BL1B:Det:TH:BL:Lambda,freq"
            )
        props = PropertyManagerDataService.retrieve("__absreductionprops")

        # Default to wavelength from JSON input / align and focus args
        if "AlignAndFocusArgs" in align_and_focus_args:
            input_wl = align_and_focus_args["AlignAndFocusArgs"]
            if "TMin" and "TMax" in input_wl:
                props["tof_min"] = input_wl["TMin"]
                props["tof_max"] = input_wl["TMax"]

        # But set wavelength max from logs if not set in JSON or elsewhere
        else:
            wl_lognames = [
                "LambdaRequest",
                "lambda",
                "skf12.lambda",
                "BL1B:Det:TH:BL:Lambda",
                "frequency"]

            for logname_wl in wl_lognames:
                run = abs_input.run()
                is_max_wavelength_zero = props["wavelength_max"].value == 0
                if logname_wl in run and is_max_wavelength_zero:
                    props["wavelength_max"] = run[logname_wl].lastValue()

    # NOTE: We have two options from this point forward.
    #       As of 10-04-2021, use option 2 to bypass the automated caching

    # Option 1: use top level API from absorptioncorrutils for easy caching
    # abs_s, abs_c = absorptioncorrutils.calculate_absorption_correction(
    #                   filename,
    #                   abs_method,
    #                   props,
    #                   sample_formula=material['ChemicalFormula'],
    #                   mass_density=material['SampleMassDensity'],
    #                   cache_dir=align_and_focus_args["CacheDir"],
    #                   ms_method=ms_method,
    # )

    # Option 2 (Original method)
    # Use low level API from absorptioncorrutils to bypass the automated
    # caching
    # 1. Setup the donor workspace for absorption correction
    re_gen_group = num_groups > 0
    try:
        if isinstance(filename, str):
            list_filenames = filename.split(",")
            filename = list_filenames[0]

        cond1 = material["ChemicalFormula"] == "V"
        cond2 = material["ChemicalFormula"] == "V1"
        find_env = not (cond1 or cond2)

        donor_ws = absorptioncorrutils.create_absorption_input(
            filename,
            props,
            material=material,
            geometry=geometry,
            environment=environment,
            find_environment=find_env)

        if not (group_wksp_in is None or re_gen_group):
            mask_list = mask_generator(group_wksp_in, group_ref_det_out_file)
            MaskDetectors(Workspace=donor_ws, DetectorList=mask_list)

    except RuntimeError as e:
        msg = "Could not create absorption correction donor workspace: {}"
        raise RuntimeError(msg.format(e))

    # 2. calculate the absorption workspace (first order absorption) without
    #    calling to cache
    if "CacheDir" in align_and_focus_args:
        if "BL11A:CS:ITEMS:HeightInContainer" not in mtd[donor_ws].run():
            h_val = geometry["Height"]
            mtd[donor_ws].run()["BL11A:CS:ITEMS:HeightInContainer"] = h_val
        abs_s, abs_c = absorptioncorrutils.calc_absorption_corr_using_wksp(
            donor_ws,
            abs_method,
            element_size=elementsize,
            cache_dirs=align_and_focus_args["CacheDir"]
        )
    else:
        abs_s, abs_c = absorptioncorrutils.calc_absorption_corr_using_wksp(
            donor_ws,
            abs_method,
            element_size=elementsize
        )

    if not (group_wksp_in is None or re_gen_group):
        if abs_s != "":
            MaskDetectors(Workspace=abs_s, DetectorList=mask_list)
        if abs_c != "":
            MaskDetectors(Workspace=abs_c, DetectorList=mask_list)

    # 3. Convert to effective absorption correction workspace if multiple
    # scattering correction is requested
    # NOTE:
    #   Multiple scattering and absorption correction are using the same
    #   element size when discretizing the volume.
    if ms_method is not None:
        try:
            donor_ws_tmp = mtd[donor_ws].getItem(0)
        except AttributeError:
            donor_ws_tmp = mtd[donor_ws]
        MultipleScatteringCorrection(
            InputWorkspace=donor_ws_tmp,
            ElementSize=elementsize,
            ContainerElementSize=con_elementsize,
            method=ms_method,
            OutputWorkspace="ms_tmp"
        )
        if ms_method == "SampleOnly":
            ms_sampleOnly = mtd["ms_tmp_sampleOnly"]
            ms_sampleOnly = 1 - ms_sampleOnly
            # abs_s now point to the effective absorption correction
            # A = A / (1 - ms_s)
            Divide(
                LHSWorkspace=abs_s,  # str
                RHSWorkspace=ms_sampleOnly,  # workspace
                OutputWorkspace=abs_s,  # str
                )
            # nothing need to be done for container
            mtd.remove("ms_tmp_sampleOnly")
        elif ms_method == "SampleAndContainer":
            ms_sampleAndContainer = mtd["ms_tmp_sampleAndContainer"]
            ms_sampleAndContainer = 1 - ms_sampleAndContainer
            Divide(
                LHSWorkspace=abs_s,  # str
                RHSWorkspace=ms_sampleAndContainer,  # workspace
                OutputWorkspace=abs_s,  # str
            )
            mtd.remove("ms_tmp_sampleAndContainer")
            ms_containerOnly = mtd["ms_tmp_containerOnly"]
            ms_containerOnly = 1 - ms_containerOnly
            Divide(
                LHSWorkspace=abs_c,  # str
                RHSWorkspace=ms_containerOnly,  # workspace
                OutputWorkspace=abs_c,  # str
            )
            mtd.remove("ms_tmp_containerOnly")
        else:
            logging.warning(
                f"multiple scattering correction {ms_method}"
                "is performed independent from absorption correction."
                )

    if not (group_wksp_in is None or re_gen_group):
        if abs_s != "":
            GroupDetectors(InputWorkspace=abs_s,
                           OutputWorkspace=abs_s, Behaviour="Sum",
                           CopyGroupingFromWorkspace=group_wksp_in)
        if abs_c != "":
            GroupDetectors(InputWorkspace=abs_c,
                           OutputWorkspace=abs_c,
                           Behaviour="Sum",
                           CopyGroupingFromWorkspace=group_wksp_in)

    if re_gen_group:
        abs_s_g, abs_c_g, group_wksp_out = abs_grouping(abs_s,
                                                        abs_c,
                                                        num_groups,
                                                        group_out_file,
                                                        group_ref_det_out_file,
                                                        sg_index_f)
        return abs_s_g, abs_c_g, group_wksp_out

    return abs_s, abs_c, group_wksp_out
