import os
from mantid import mtd
from mantid.api import IEventWorkspace
from mantid.simpleapi import \
    CloneWorkspace, Rebin, ConvertToDistribution, DiffractionFocussing, \
    SaveNexusProcessed, SaveAscii, DeleteWorkspace


def save_banks(InputWorkspace, Filename, Title, OutputDir,
               Binning=None, GroupingWorkspace=None):
    """
    Saves input workspace to processed NeXus file in specified
    output directory with optional rebinning and grouping
    (to coarsen) the output in a bank-by-bank manner. Mainly
    wraps Mantid `SaveNexusProcessed` algorithm.

    :param InputWorkspace: Mantid workspace to save out
    :type InputWorkspace: MatrixWorkspace
    :param Filename: Filename to save output
    :type Filename: str
    :param Title: A title to describe the saved workspace
    :type Title: str
    :param OutputDir: Output directory to save the processed NeXus file
    :type OutputDir: path str
    :param Binning: Optional rebinning of event workspace.
                    See `Rebin` in Mantid for options
    :type Binning: dbl list
    :param GroupingWorkspace: A workspace with grouping
                              information for the output spectra
    :type GroupWorkspace: GroupWorkspace
    """

    # Make a local clone
    CloneWorkspace(InputWorkspace=InputWorkspace, OutputWorkspace="__tmp")
    tmp_wksp = mtd["__tmp"]

    # Rebin if requested
    if Binning:
        tmp_wksp = Rebin(InputWorkspace=tmp_wksp,
                         Params=Binning,
                         PreserveEvents=True)

    # Convert to distributions to remove bin width dependence
    yunit = tmp_wksp.YUnit()
    if yunit == "Counts":
        try:
            ConvertToDistribution(tmp_wksp)
        except BaseException:
            pass

    # Output to desired level of grouping
    isEventWksp = isinstance(tmp_wksp, IEventWorkspace)
    if isEventWksp and GroupingWorkspace and yunit == "Counts":
        tmp_wksp = DiffractionFocussing(InputWorkspace=tmp_wksp,
                                        GroupingWorkspace=GroupingWorkspace,
                                        PreserveEvents=False)

    # Save out wksp to file
    filename = os.path.join(os.path.abspath(OutputDir), Filename)
    SaveNexusProcessed(
        InputWorkspace=tmp_wksp,
        Filename=filename,
        Title=Title,
        Append=True,
        PreserveEvents=False,
        WorkspaceIndexList=range(
            tmp_wksp.getNumberHistograms()))
    DeleteWorkspace(tmp_wksp)


def save_file(ws, filename, header=None):
    """
    Small wrapper to Mantid `SaveAscii` algorithm to add a header lines.

    :param ws: Mantid workspace to save out
    :type ws: MatrixWorkspace
    :param filename: Filename to save output
    :type filename: str
    :param headers: A list of header comments
    :type headers: list
    """
    with open(filename, 'w') as f:
        if header:
            for line in header:
                f.write('# %s \n' % line)
    SaveAscii(
        InputWorkspace=ws,
        Filename=filename,
        Separator='Space',
        ColumnHeader=False,
        AppendToFile=True)
