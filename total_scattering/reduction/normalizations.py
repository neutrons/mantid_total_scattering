# 3rd-party imports
from mantid.api import mtd
from mantid.dataobjects import Workspace2D
from mantid.simpleapi import Divide, RenameWorkspace

# standard imports
from typing import Union


class Material:
    r"""Wrapper of mantid.kernel.Material with some convenient properties"""

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
        return self.cohScatterLength()**2

    @property
    def btot_sqrd_avg(self):
        return self.totalScatterLengthSqrd()

    @property
    def laue_monotonic_diffuse_scat(self):
        return self.btot_sqrd_avg / self.bcoh_avg_sqrd


def to_absolute_scale(
        input_workspace: Union[str, Workspace2D],
        output_workspace: str) -> None:
    r"""
    :param input_workspace: structure factor
    :param output_workspace: structure factor
    """
    m = Material(input_workspace)
    s_of_q = mtd[str(input_workspace)]
    normalized = (1. / m.bcoh_avg_sqrd) * s_of_q - \
        m.laue_monotonic_diffuse_scat + 1.
    RenameWorkspace(InputWorkspace=normalized,
                    OutputWorkspace=output_workspace)


# TODO #354 (This mock returns 1.0 everywhere)
def calculate_fitted_levels(
        input_workspace: Union[str, Workspace2D],
        q_ranges: dict,
        output_workspace: str = 'fitted_levels') -> Workspace2D:
    r"""
    :param input_workspace:
    :param output_workspace:
    """
    return Divide(LHSWorkspace=input_workspace,
                  RHSWorkspace=input_workspace,
                  OutputWorkspace=output_workspace)


def to_f_of_q(
        input_workspace: Union[str, Workspace2D],
        output_workspace: str) -> None:
    r"""
    :param input_workspace:
    :param output_workspace:
    """
    m = Material(input_workspace)
    s_of_q = mtd[str(input_workspace)]
    f_of_q = m.bcoh_avg_sqrd * (s_of_q - 1)
    RenameWorkspace(InputWorkspace=f_of_q,
                    OutputWorkspace=output_workspace)
