#!/usr/bin/env python

from pdf_plotter.io.nexus_load     import ExperimentFileInput
from pdf_plotter.ui.control_panel import ControlPanel, ControlPanelView, ControlPanelHandler


# Use the ControlPanel to View the Measurement
cp = ControlPanel(experiment_file=ExperimentFileInput())
cp.configure_traits(view=ControlPanelView, handler=ControlPanelHandler)

