from mantid.simpleapi import \
    AlignAndFocusPowderFromFiles, NormaliseByCurrent, SetSample, ConvertUnits


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
        **align_and_focus_args):
    AlignAndFocusPowderFromFiles(OutputWorkspace=ws_name,
                                 Filename=input_files,
                                 Absorption=None,
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
