# 3rd-party imports
from mantid.api import mtd
from mantid.dataobjects import Workspace2D
from mantid.simpleapi import (
    CloneWorkspace, CreateWorkspace, DeleteWorkspace, MatchSpectra,
    RenameWorkspace)

# standard imports
import numpy as np
from typing import Optional, Tuple, Union


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
        q_ranges: dict,
        output_workspace: Optional[str] = None) -> Tuple[Workspace2D, dict]:
    r"""Fits a horizontal line to each bank and region specified in q_ranges.
    Scales the full bank data in input_workspace by the offset (fitted level)
    and returns a new workspace with the scaled data (s_q_norm)

    :param input_workspace: workspace containing S(Q) with spectra for each
    bank, in MomentumTransfer units
    :param q_ranges: dictionary of bank number with tuple range of Q fitting
    :param output_workspace: name of output workspace. If `None`, the name is
    that of the input_workspace plus the `_scaled` suffix.
    :return: Clone of input_workspace scaled by fitted level for banks
    specified in q_ranges, and a dict of banks with bad fitted levels (<= 0)
    """
    if output_workspace is None:
        output_workspace = str(input_workspace) + '_scaled'

    input_workspace = mtd[str(input_workspace)]
    bad_levels = dict()

    # Check units (and for now, throw error if not in MomentumTransfer)
    if input_workspace.getAxis(0).getUnit().unitID() != "MomentumTransfer":
        raise RuntimeError(
            "S(Q) workspace was expected to be in MomentumTransfer, "
            "but has units '{}'".format(
                input_workspace.getAxis(0).getUnit().unitID()))

    # Setup the output workspace. This gets returned as-is if q_ranges is empty
    s_q_norm = CloneWorkspace(InputWorkspace=input_workspace,
                              OutputWorkspace=str(output_workspace))

    # Perform the horizontal fitting to each bank in q_ranges
    for key, value in q_ranges.items():
        # get corresponding workspace index number in S(Q) from the bank number
        ws_index = key - 1

        start_index = input_workspace.yIndexOfX(value[0], ws_index) + 1
        end_index = input_workspace.yIndexOfX(value[1], ws_index) + 1
        # Extract the bank between the fit regions
        bank_x = input_workspace.readX(ws_index)[start_index:end_index]
        bank_y = input_workspace.readY(ws_index)[start_index:end_index]
        bank_e = input_workspace.readE(ws_index)[start_index:end_index]

        # Create reference spectrum which is flat line at 1.0
        ref_y = np.zeros(bank_x.shape) + 1.0

        # Create a workspace with 2 spectra (bank data and reference spectrum)
        # as input to MatchSpectra so that the bank data can be fit to the flat
        # line at 1.0
        CreateWorkspace(np.concatenate((bank_x, bank_x)),
                        np.concatenate((ref_y, bank_y)),
                        np.concatenate((ref_y, bank_e)),
                        NSpec=2,
                        OutputWorkspace="__tmp_bank_fit")
        match_result = MatchSpectra(InputWorkspace="__tmp_bank_fit",
                                    OutputWorkspace="__tmp_bank_fit",
                                    ReferenceSpectrum=1,
                                    CalculateOffset=True,
                                    CalculateScale=False)
        offset = match_result.Offset[1]
        scale = match_result.Scale[1]
        # Scale away from 1.0 means the baseline is tilting and therefore we
        # may have some uncorrected effect (e.g., hydrogen background).
        if abs(scale - 1.0) > 0.5:
            bad_levels[key] = scale
        else:
            # scale full bank data by offset
            s_q_norm.setY(ws_index, s_q_norm.dataY(ws_index) * (1.0 / (1.0 - offset)))

        DeleteWorkspace("__tmp_bank_fit")

    return s_q_norm, bad_levels


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
