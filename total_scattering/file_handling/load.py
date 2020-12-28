from mantid.simpleapi import \
    AlignAndFocusPowderFromFiles, Load, NormaliseByCurrent, PDDetermineCharacterizations, \
    PDLoadCharacterizations, PropertyManagerDataService, SetSample, ConvertUnits
from mantid.utils import AbsorptionCorrUtils


required_shape_keys = {
    "FlatPlate": ["Shape", "Width", "Height", "Thick", "Center", "Angle"],
    "Cylinder": ["Shape", "Height", "Radius", "Center"],
    "HollowCylinder": ["Shape", "Height", "InnerRadius", "OuterRadius", "Center"]
}


def load(
        ws_name,
        input_files,
        geometry=None,
        chemical_formula=None,
        mass_density=None,
        absorption_wksp='',
        **align_and_focus_args):
    AlignAndFocusPowderFromFiles(OutputWorkspace=ws_name,
                                 Filename=input_files,
                                 AbsorptionWorkspace=absorption_wksp,
                                 **align_and_focus_args)
    NormaliseByCurrent(InputWorkspace=ws_name,
                       OutputWorkspace=ws_name,
                       RecalculatePCharge=True)
    if geometry and chemical_formula and mass_density:
        set_sample(ws_name, geometry, chemical_formula, mass_density)

    ConvertUnits(InputWorkspace=ws_name,
                 OutputWorkspace=ws_name,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    return ws_name


def set_sample(ws_name, geometry=None, chemical_formula=None, mass_density=None):
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
    new_dict = dict()
    for key in required_shape_keys[shape]:
        if key not in mydict:
            new_dict[key] = None
        else:
            new_dict[key] = mydict[key]
    new_dict['Shape'] = shape
    return new_dict


def create_absorption_wksp(filename, abs_method, geometry, material,
                           environment=None, props=None, characterization_files=None):
    if abs_method is None:
        return '', ''

    # Check against supported absorption correction methods so we can error out early if needed
    valid_methods = ["SampleOnly", "SampleAndContainer", "FullPaalmanPings"]
    if abs_method not in valid_methods:
        raise RuntimeError("Unrecognized absorption correction method '{}'".format(abs_method))

    abs_input = Load(filename, MetaDataOnly=True)

    # If no run characterization properties were given, load them if any files were provided
    if not props and characterization_files:
        print("No props were given, but determining from characterization files")
        charfile = characterization_files
        # Reduce to a string if multiple files were provided
        if isinstance(charfile, list):
            charfile = ','.join(characterization_files)
        charTable = PDLoadCharacterizations(Filename=charfile)
        chars = charTable[0]

        # Create the properties for the absorption workspace
        #  Note: Should BackRun, NormRun, and NormBackRun be specified here?
        PDDetermineCharacterizations(InputWorkspace=abs_input,
                                     Characterizations=chars,
                                     ReductionProperties="__absreductionprops",
                                     WaveLengthLogNames="LambdaRequest,lambda,skf12.lambda,"
                                                        "BL1B:Det:TH:BL:Lambda,freq")
        props = PropertyManagerDataService.retrieve("__absreductionprops")

    # If neither run characterization properties or char files were given, try to guess from input
    if not (props and characterization_files):
        print("No props or characterizations were given, determining props from input file")
        PDDetermineCharacterizations(InputWorkspace=abs_input,
                                     ReductionProperties="__absreductionprops",
                                     WaveLengthLogNames="LambdaRequest,lambda,skf12.lambda,"
                                                        "BL1B:Det:TH:BL:Lambda,freq")
        props = PropertyManagerDataService.retrieve("__absreductionprops")

        # Set wavelength max from logs if not set
        wl_lognames = ["LambdaRequest", "lambda", "skf12.lambda", "BL1B:Det:TH:BL:Lambda", "freq"]
        for logname_wl in wl_lognames:
            if logname_wl in abs_input.run() and props["wavelength_max"].value == 0:
                props["wavelength_max"] = abs_input.run()[logname_wl].lastValue()

    # Setup the donor workspace for absorption correction
    try:
        donor_ws = AbsorptionCorrUtils.create_absorption_input(filename, props, material=material,
                                                               geometry=geometry,
                                                               environment=environment)
    except RuntimeError as e:
        raise RuntimeError("Could not create absorption correction donor workspace: {}".format(e))

    abs_s, abs_c = AbsorptionCorrUtils.calc_absorption_corr_using_wksp(donor_ws, abs_method)

    return abs_s, abs_c
