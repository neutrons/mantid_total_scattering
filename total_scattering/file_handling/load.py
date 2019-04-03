from mantid.simpleapi import \
    AlignAndFocusPowderFromFiles, NormaliseByCurrent, SetSample, ConvertUnits


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
    SetSample(
        InputWorkspace=ws_name,
        Geometry=geometry,
        Material={
            'ChemicalFormula': chemical_formula,
            'SampleMassDensity': mass_density})
