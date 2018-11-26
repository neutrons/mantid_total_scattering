from mantid.simpleapi import AlignAndFocusPowderFromFiles, NormaliseByCurrent, SetSample, ConvertUnits


def load(ws_name, input_files, Geometry=None, ChemicalFormula=None, MassDensity=None, **align_and_focus_args):
    AlignAndFocusPowderFromFiles(OutputWorkspace=ws_name,
                                 Filename=input_files,
                                 Absorption=None,
                                 **align_and_focus_args)
    NormaliseByCurrent(InputWorkspace=ws_name,
                       OutputWorkspace=ws_name,
                       RecalculatePCharge=True)
    if Geometry and ChemicalFormula and MassDensity:
        set_sample(ws_name, Geometry, ChemicalFormula, MassDensity)

    ConvertUnits(InputWorkspace=ws_name,
                 OutputWorkspace=ws_name,
                 Target="MomentumTransfer",
                 EMode="Elastic")
    return ws_name


def set_sample(ws_name, Geometry=None, ChemicalFormula=None, MassDensity=None):
    if 'Center' not in Geometry:
        Geometry.update({'Center': [0., 0., 0., ]})
    if "Shape" not in Geometry:
        Geometry.update({'Shape': 'Cylinder'})
    SetSample(
        InputWorkspace=ws_name,
        Geometry=Geometry,
        Material={
            'ChemicalFormula': ChemicalFormula,
            'SampleMassDensity': MassDensity})
