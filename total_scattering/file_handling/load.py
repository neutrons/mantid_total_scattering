import logging
from mantid import mtd
from mantid.simpleapi import \
    AlignAndFocusPowderFromFiles, \
    ConvertUnits, \
    Divide, \
    Load, \
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
         **align_and_focus_args):
    '''Routine for loading workspace'''

    # Figure out the list of run number
    run_list = list()
    for item in input_files.split(","):
        if "/" in item or "\\" in item:
            file_n = os.path.basename(item)
            reg_test = re.search('.*_[0-9]+\\.nxs\\.h5', file_n)
            if reg_test is None:
                run_list.append("-1")
            else:
                run_list.append(reg_test.group(0).split(".")[0].split("_")[1])
        else:
            run_list.append(item.split("_")[1])

    if group_wksp is None:
        # Given the current implementation mechanism, if the input
        # `group_wksp` is None, that means there was no absorption
        # calculation involved in the preparation stage (i.e., before
        # calling current load method). In such a situation, we go
        # ahead to perform the align and focus in the old way by
        # calling the `AlignAndFocusPowderFromFiles` algorithm,
        # followed by saving the summed file to cache. For sure,
        # if such a cache file already exists, we load it in.

        # Figure out a unique name for the summed cache file
        # The codes here were originating from GPT3.5-turbo API
        hash_obj = hashlib.sha256(str(run_list).encode())
        hash_str = hash_obj.hexdigest()
        short_hash_str = base64.urlsafe_b64encode(hash_str.encode()).decode()[:12]
        cache_sf_bn = f"{instr_name}_mts_summed_{short_hash_str}.nxs"
        cache_sf_fn = os.path.join("/" + facility,
                                   instr_name,
                                   ipts,
                                   "shared",
                                   "autoreduce",
                                   cache_sf_bn)

        if os.path.isfile(cache_sf_fn):
            LoadNexus(Filename=cache_sf_fn, OutputWorkspace=ws_name)
        else:
            for run_i, run in enumerate(run_list):
                if run == "-1":
                    cache_f_exist = False
                else:
                    cache_f_bn = f"{instr_name}_{run}_mts_reduced_no_subg.nxs"
                    cache_f_fn = os.path.join("/" + facility,
                                              instr_name,
                                              ipts,
                                              "shared",
                                              "autoreduce",
                                              cache_f_bn)
                    if os.path.isfile(cache_f_fn):
                        cache_f_exist = True
                    else:
                        cache_f_exist = False

                if cache_f_exist:
                    wksp_tmp = LoadNexus(Filename=cache_f_fn)
                else:
                    wksp_tmp = "wksp_tmp"
                    AlignAndFocusPowderFromFiles(
                        OutputWorkspace=wksp_tmp,
                        Filename=input_files.split(",")[run_i],
                        AbsorptionWorkspace=absorption_wksp,
                        **align_and_focus_args)
                    ConvertUnits(
                        InputWorkspace=wksp_tmp,
                        OutputWorkspace=wksp_tmp,
                        Target="MomentumTransfer",
                        EMode="Elastic")
                    Rebin(
                        InputWorkspace=wksp_tmp,
                        OutputWorkspace=wksp_tmp,
                        Params=qparams)
                    SaveNexusProcessed(
                        InputWorkspace=wksp_tmp,
                        Filename=cache_f_fn,
                        Title=f"{run}_cached_no_abs",
                        WorkspaceIndexList=range(
                            mtd[wksp_tmp].getNumberHistograms()))

                # Accumulate individual files
                if run_i == 0:
                    CloneWorkspace(InputWorkspace=wksp_tmp,
                                   OutputWorkspace=ws_name)
                else:
                    Plus(LHSWorkspace=ws_name,
                         RHSWorkspace=wksp_tmp,
                         OutputWorkspace=ws_name)

            SaveNexusProcessed(ws_name, cache_sf_fn, Title="cache_summed")
    else:
        for run_i, run in enumerate(run_list):
            if run == "-1":
                cache_f_exist = False
            else:
                cache_f_bn = f"{instr_name}_{run}_mts_reduced.nxs"
                cache_f_fn = os.path.join("/" + facility,
                                          instr_name,
                                          ipts,
                                          "shared",
                                          "autoreduce",
                                          cache_f_bn)
                if os.path.isfile(cache_f_fn):
                    cache_f_exist = True
                else:
                    cache_f_exist = False

            # `group_num == 0` means no regeneration of grouping
            # was ever initialized and only in such cases will we
            # move ahead to load in the existing cache file.
            if cache_f_exist and group_num == 0:
                wksp_tmp = LoadNexus(Filename=cache_f_fn)
            else:
                wksp_tmp = "wksp_tmp"
                AlignAndFocusPowderFromFiles(
                    OutputWorkspace=wksp_tmp,
                    Filename=input_files.split(",")[run_i],
                    GroupingWorkspace=group_wksp,
                    **align_and_focus_args)
                tmin_tmp = align_and_focus_args["TMin"]
                tmax_tmp = align_and_focus_args["TMax"]
                params = f"{tmin_tmp}, -0.0008, {tmax_tmp}"
                Rebin(InputWorkspace=wksp_tmp,
                      OutputWorkspace=wksp_tmp,
                      Params=params)
                ConvertUnits(
                    InputWorkspace=wksp_tmp,
                    OutputWorkspace=wksp_tmp,
                    Target="Wavelength",
                    EMode="Elastic")
                SaveNexusProcessed(
                    InputWorkspace=wksp_tmp,
                    Filename=cache_f_fn,
                    Title=f"{run}_cached",
                    WorkspaceIndexList=range(
                        mtd[wksp_tmp].getNumberHistograms()))

            # Accumulate individual files
            if run_i == 0:
                CloneWorkspace(InputWorkspace=wksp_tmp,
                               OutputWorkspace=ws_name)
            else:
                Plus(LHSWorkspace=ws_name,
                     RHSWorkspace=wksp_tmp,
                     OutputWorkspace=ws_name)

        if absorption_wksp != '':
            RebinToWorkspace(
                WorkspaceToRebin=absorption_wksp,
                WorkspaceToMatch=ws_name,
                OutputWorkspace=absorption_wksp)
            Divide(LHSWorkspace=mtd[ws_name],
                   RHSWorkspace=absorption_wksp,
                   OutputWorkspace=ws_name,
                   AllowDifferentNumberSpectra=True)

    if group_wksp is not None:
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

        ConvertUnits(
            InputWorkspace=ws_name,
            OutputWorkspace=ws_name,
            Target="MomentumTransfer",
            EMode="Elastic")
        Rebin(
            InputWorkspace=ws_name,
            OutputWorkspace=ws_name,
            Params=qparams)
        GroupDetectors(InputWorkspace=ws_name,
                       OutputWorkspace=ws_name,
                       GroupingPattern=g_pattern)

    NormaliseByCurrent(
        InputWorkspace=ws_name,
        OutputWorkspace=ws_name,
        RecalculatePCharge=True)

    if geometry and chemical_formula and mass_density:
        set_sample(ws_name, geometry, chemical_formula, mass_density)

    return ws_name


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
    Rebin(InputWorkspace=sam_abs_ws,
          OutputWorkspace="abs_pbp",
          Params="0.1,0.05,3.0")

    num_spectra = mtd['abs_pbp'].getNumberHistograms()
    num_monitors = int(np.sum(mtd['abs_pbp'].detectorInfo().detectorIDs() < 0))

    output_group, spectra_count, group_count = CreateGroupingWorkspace(
        InputWorkspace="abs_pbp",
        GroupDetectorsBy="Group")
    sub_num_clusters = int(num_clusters / group_count)

    with open(subgroup_index_file, "w") as f:
        f.write("# Start & Stop index of subgroups\n")
        for i in range(group_count):
            start_index = sub_num_clusters * i + 1
            stop_index = sub_num_clusters * (i + 1)
            f.write("{0:<10s}{1:<10d}{2:<10d}\n".format("Bank-" + str(i),
                                                        start_index,
                                                        stop_index))

    print("[Info] Clustering the absorption spectra for all detectors...")
    all_spectra_id = dict()
    all_spectra_all = list()
    labels = list()
    for j in range(group_count):
        all_spectra = list()
        for i in range(num_monitors, num_spectra):
            if output_group.readY(i) == float(j + 1):
                y_tmp = mtd['abs_pbp'].readY(i)
                all_spectra.append(y_tmp)
                all_spectra_id[i] = mtd['abs_pbp'].getDetector(i).getID()
        all_spectra_all.extend(all_spectra)

        X = np.array(all_spectra)
        clustering = KMeans(n_clusters=sub_num_clusters, n_init=10).fit(X)

        labels_tmp = [x + j * sub_num_clusters for x in clustering.labels_]
        labels.extend(labels_tmp)
        print(f"[Info] Done with the clustering of absorption spectra for bank-{j + 1}.")

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    CreateGroupingWorkspace(InputWorkspace="abs_pbp",
                            OutputWorkspace="grouping")
    grouping = mtd["grouping"]
    for i, label in enumerate(labels):
        if label == -1:
            continue
        det_id = all_spectra_id[i + num_monitors]
        det_id = mtd['abs_pbp'].detectorInfo().indexOf(det_id) - num_monitors
        grouping.setY(det_id, [int(label + 1)])

    SaveDetectorsGrouping(InputWorkspace=grouping,
                          OutputFile=group_out_file)
    GroupDetectors(InputWorkspace=sam_abs_ws,
                   OutputWorkspace="sam_abs_grouped", Behaviour="Average",
                   CopyGroupingFromWorkspace=grouping)
    sam_abs_grouped = mtd['sam_abs_grouped']
    if con_abs_ws != "":
        GroupDetectors(InputWorkspace=con_abs_ws,
                       OutputWorkspace="con_abs_grouped", Behaviour="Average",
                       CopyGroupingFromWorkspace=grouping)
        con_abs_grouped = mtd['con_abs_grouped']
    else:
        con_abs_grouped = ""

    print("[Info] Check the similarity of absorption spectra within groups.")
    ref_id = list()
    for i in range(group_count * sub_num_clusters):
        ref_id.append(list(labels).index(i) + 2)
    with open(group_ref_det_out_file, "w") as f:
        for i, item in enumerate(ref_id):
            if i == len(ref_id) - 1:
                f.write("{0:d}".format(item))
            else:
                f.write("{0:d}\n".format(item))

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
    print("[Info] Maximum percentage difference: {0:3.1F}".format(perct_max))

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
    abs_s, abs_c = absorptioncorrutils.calc_absorption_corr_using_wksp(
            donor_ws,
            abs_method,
            element_size=elementsize)

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
