from mantid.simpleapi import AlignAndFocusPowderFromFiles, NormaliseByCurrent, SetSample, ConvertUnits


def load(ws_name, input_files, geometry=None, material=None, mass_density=None, **align_and_focus_args):
    AlignAndFocusPowderFromFiles(OutputWorkspace=ws_name,
                                 Filename=input_files,
                                 Absorption=None,
                                 **align_and_focus_args)
    NormaliseByCurrent(InputWorkspace=ws_name,
                       OutputWorkspace=ws_name,
                       RecalculatePCharge=True)
    if geometry and material and mass_density:
        set_sample(ws_name, geometry, material, mass_density)

    ConvertUnits(InputWorkspace=ws_name,
                 OutputWorkspace=ws_name,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    return ws_name


def set_sample(ws_name, geometry, material, mass_density):
    new_geometry = dict()
    for k, v in geometry.items():
        key = str(k)
        v = str(v)
        new_geometry[key] = v
    geometry = new_geometry
    geometry.update({'Center': [0., 0., 0., ]})
    if "Shape" not in geometry:
        geometry.update({'Shape': 'Cylinder'})
    material = str(material)
    SetSample(
        InputWorkspace=ws_name,
        Geometry=geometry,
        Material={
            'ChemicalFormula': material,
            'SampleMassDensity': mass_density})
