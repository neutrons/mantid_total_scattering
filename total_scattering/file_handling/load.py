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
    SetSample
from mantid.utils import absorptioncorrutils

_shared_shape_keys = ["Shape", "Height", "Center"]
required_shape_keys = {
    "FlatPlate": _shared_shape_keys + ["Width", "Thick", "Angle"],
    "Cylinder": _shared_shape_keys + ["Radius"],
    "HollowCylinder": _shared_shape_keys + ["InnerRadius", "OuterRadius"]
}


def load(ws_name, input_files,
         geometry=None, chemical_formula=None, mass_density=None,
         absorption_wksp='', **align_and_focus_args):
    '''Load workspace'''
    AlignAndFocusPowderFromFiles(
        OutputWorkspace=ws_name,
        Filename=input_files,
        AbsorptionWorkspace=absorption_wksp,
        **align_and_focus_args)
    NormaliseByCurrent(
        InputWorkspace=ws_name,
        OutputWorkspace=ws_name,
        RecalculatePCharge=True)
    if geometry and chemical_formula and mass_density:
        set_sample(ws_name, geometry, chemical_formula, mass_density)

    ConvertUnits(
        InputWorkspace=ws_name,
        OutputWorkspace=ws_name,
        Target="MomentumTransfer",
        EMode="Elastic")
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


def create_absorption_wksp(filename, abs_method, geometry, material,
                           environment=None, props=None,
                           characterization_files=None,
                           ms_method=None,
                           elementsize=1.0, # mm
                           **align_and_focus_args):
    '''Create absorption workspace'''
    if abs_method is None:
        return '', ''

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
    try:
        if isinstance(filename, str):
            list_filenames = filename.split(",")
            filename = list_filenames[0]

        find_environment = not material["ChemicalFormula"] == "V"

        donor_ws = absorptioncorrutils.create_absorption_input(
            filename,
            props,
            material=material,
            geometry=geometry,
            environment=environment,
            find_environment=find_environment)

    except RuntimeError as e:
        msg = "Could not create absorption correction donor workspace: {}"
        raise RuntimeError(msg.format(e))

    # 2. calculate the absorption workspace (first order absorption) without
    #    calling to cache
    abs_s, abs_c = absorptioncorrutils.calc_absorption_corr_using_wksp(
            donor_ws,
            abs_method,
            element_size=elementsize)

    # 3. Convert to effective absorption correction workspace if multiple
    # scattering correction is requested
    # NOTE:
    #   Multiple scattering and absorption correction are using the same
    #   element size when discretizing the volume.
    if ms_method is not None:
        MultipleScatteringCorrection(
            InputWorkspace=donor_ws,
            ElementSize=elementsize,
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

    return abs_s, abs_c
