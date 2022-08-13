# 3rd-party imports
from mantid.api import mtd
from mantid.dataobjects import Workspace2D
from mantid.simpleapi import (
    CreateWorkspace, DeleteWorkspace,
    RenameWorkspace, Fit)

# standard imports
from typing import Tuple, Union


class Material:
    r"""Thin wrapper of mantid.kernel.Material, along with some
     convenient properties"""

    def __init__(self, input_workspace: Union[str, Workspace2D]):
        material = mtd[str(input_workspace)].sample().getMaterial()
        if material.name() is None or len(material.name().strip()) == 0:
            raise RuntimeError('Sample material was not set')
        self.__dict__['_material'] = material

    def __getattr__(self, item):
        if item not in self.__dict__:
            return getattr(self._material, item)

    @property
    def bcoh_avg_sqrd(self):
        return self.cohScatterLength()**2 / 100.

    @property
    def btot_sqrd_avg(self):
        return self.totalScatterLengthSqrd() / 100.

    @property
    def laue_monotonic_diffuse_scat(self):
        return self.btot_sqrd_avg / self.bcoh_avg_sqrd


def to_absolute_scale(
        input_workspace: Union[str, Workspace2D],
        output_workspace: str) -> Workspace2D:
    r"""
    :param input_workspace: structure factor
    :param output_workspace: structure factor
    """
    m = Material(input_workspace)
    s_of_q = mtd[str(input_workspace)]
    _normalized = (1. / m.bcoh_avg_sqrd) * s_of_q +\
        (1 - m.laue_monotonic_diffuse_scat)
    RenameWorkspace(InputWorkspace=_normalized,
                    OutputWorkspace=output_workspace)
    return mtd[output_workspace]


def calculate_and_apply_fitted_levels(
        input_workspace: Union[str, Workspace2D],
        q_ranges: dict) -> Tuple[Workspace2D, dict]:
    r"""Fits a straight line to each bank and region specified in q_ranges.

    :param input_workspace: workspace containing S(Q) with spectra for each
    bank, in MomentumTransfer units
    :param q_ranges: dictionary of bank number with tuple range of Q fitting
    :return: fitted offset and slope for all banks
    """
    input_workspace = mtd[str(input_workspace)]

    # Check units (and for now, throw error if not in MomentumTransfer)
    if input_workspace.getAxis(0).getUnit().unitID() != "MomentumTransfer":
        raise RuntimeError(
            "S(Q) workspace was expected to be in MomentumTransfer, "
            "but has units '{}'".format(
                input_workspace.getAxis(0).getUnit().unitID()))

    # Perform the horizontal fitting to each bank in q_ranges
    offset = dict()
    slope = dict()
    for key, value in q_ranges.items():
        # get corresponding workspace index number in S(Q) from the bank number
        ws_index = key - 1

        start_index = input_workspace.yIndexOfX(value[0], ws_index) + 1
        end_index = input_workspace.yIndexOfX(value[1], ws_index) + 1
        # Extract the bank between the fit regions
        bank_x = input_workspace.readX(ws_index)[start_index:end_index]
        bank_y = input_workspace.readY(ws_index)[start_index:end_index]

        tmp_wks = CreateWorkspace(bank_x, bank_y)
        Fit("name=UserFunction, Formula=a+b*x, a=1, b=0",
            tmp_wks,
            Output='tmp_wks_fitted')
        offset_tmp = mtd['tmp_wks_fitted_Parameters'].row(0)['Value']
        slope_tmp = mtd['tmp_wks_fitted_Parameters'].row(1)['Value']

        offset[key] = offset_tmp
        slope[key] = slope_tmp

        DeleteWorkspace(tmp_wks)
        DeleteWorkspace('tmp_wks_fitted_Parameters')
        DeleteWorkspace('tmp_wks_fitted_Workspaces')
        DeleteWorkspace('tmp_wks_fitted_NormalisedCovarianceMatrix')

    return offset, slope


def to_f_of_q(
        input_workspace: Union[str, Workspace2D],
        output_workspace: str) -> Workspace2D:
    r"""
    :param input_workspace:
    :param output_workspace:
    """
    m = Material(input_workspace)
    s_of_q = mtd[str(input_workspace)]
    _f_of_q = m.bcoh_avg_sqrd * (s_of_q - 1)
    RenameWorkspace(InputWorkspace=_f_of_q,
                    OutputWorkspace=output_workspace)
    return mtd[output_workspace]
