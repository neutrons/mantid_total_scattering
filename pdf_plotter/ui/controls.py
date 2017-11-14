#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function)

import numpy as np

# Traits
from traits.api \
    import HasTraits, Instance, List, Property, Any, Bool, \
    property_depends_on, cached_property, on_trait_change

from traitsui.api \
    import InstanceEditor, \
    View, VSplit, UItem

# Local
from pdf_plotter.ui.models \
    import Dataset, CorrectedDatasets, Experiment

from pdf_plotter.ui.nodes.base_node \
    import NodeControls, NodeButtons

from pdf_plotter.ui.editors \
    import ExperimentTreeEditor

# -----------------------------------------------------------#
# Controls Model

ControlsView = View(
    VSplit(
        # Experiment Tree
        UItem(
            name='experiment',
            editor=ExperimentTreeEditor,
            resizable=True,
            show_label=False,
            height=0.7,
        ),
        UItem('node_controls',
              editor=InstanceEditor(),
              style='custom',
              resizable=True,
              height=0.2,
              ),
        UItem('node_buttons',
              editor=InstanceEditor(),
              style='custom',
              resizable=True,
              height=0.1,
              show_label=False,
              ),
    ),
)


class Controls(HasTraits):
    # View
    view = ControlsView

    # -------------------------------------------------------#
    # Traits

    # Passed in measurement
    experiment = Instance(Experiment, ())

    # Flag to help with unnecessary trait notification on initialization
    intialized = Bool(False)

    # Controls for selected node type
    node_controls = Instance(NodeControls)

    # Buttons for the selected node type
    node_buttons = Instance(NodeButtons)

    # The currently selected node
    selected = Any

    # Current selected dataset
    current_dataset = Dataset

    # Node controls used for different types of TreeNodes
    # Cached plots we keep on plot
    cached_plots = List

    # Cached properties for the loaded experiment
    datasets = Property(depends_on='experiment')
    xlist = Property(depends_on='datasets')
    ylist = Property(depends_on='datasets')

    exp_xmin = Property(depends_on='xlist')
    exp_xmax = Property(depends_on='xlist')

    # -------------------------------------------------------#
    # Utilities

    # Adds nodes to the Tree View in Controls
    def add_plot_to_node(self, dataset, parents):
        if dataset is None or parents is None:
            return

        # Get the pointer to the right Measurements and CorrectedDatasets
        for m in self.experiment.measurements:
            if m == parents['measurement']:
                measurement = m

        # Create the 'Other' CorrectedDatasets Node if it does not exist
        if 'Other' not in [m.title for m in measurement.corrected_datasets]:
            other = CorrectedDatasets(datasets=[dataset], title='Other')
            measurement.corrected_datasets.append(other)
        else:
            other = [m for m in measurement.corrected_datasets
                     if m.title == 'Other']
            if len(other) != 1:
                print("WARNING: More than 1 'Other' CorrectedDatsets...")
            other = other[0]
            other.datasets.append(dataset)

    # -------------------------------------------------------#
    # Dynamic
    @on_trait_change('node_controls.dataset_selected')
    def get_current_dataset(self):
        if isinstance(self.selected, Dataset):
            self.current_dataset = self.selected
        elif isinstance(self.selected, CorrectedDatasets):
            if self.node_controls.dataset_selected:
                self.current_dataset = self.node_controls.dataset_selected[0]
            else:
                self.current_dataset = self.selected.datasets[0]


    # Extracts Datasets models that are stored in the Experiment model
    @property_depends_on('experiment')
    def _get_datasets(self):
        datasets = list()
        for measurement in self.experiment.measurements:
            for corrected_dataset in measurement.corrected_datasets:
                for dataset in corrected_dataset.datasets:
                    # Strip +- inf
                    y = dataset.y
                    y[y == np.inf] = 0.0
                    y[y == -np.inf] = 0.0

                    datasets.append(dataset)
        return datasets

    # Get a list of all x values for all datasets
    @cached_property
    def _get_xlist(self):
        xlist = []
        for d in self.datasets:
            xlist = np.append(xlist, d.x, axis=None)
        return xlist

    @cached_property
    def _get_exp_xmin(self):
        return min(self.xlist)

    @cached_property
    def _get_exp_xmax(self):
        return max(self.xlist)

    @on_trait_change("selected")
    def updated_selected_in_node_controls(self):
        self.node_controls.selected = self.selected
