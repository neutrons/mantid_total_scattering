#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function)

import numpy as np

# Traits
from traits.api \
    import HasTraits, Instance, Property, CFloat, \
    Dict, on_trait_change

from traitsui.api \
    import View, Item, UItem, StatusItem, HSplit, VSplit,\
    InstanceEditor, Handler

# Matplotlib
from matplotlib import cm
from matplotlib.figure import Figure


# Local
from pdf_plotter.utils.mpl_utilities \
    import ZoomOnWheel, DraggableLegend, MPLFigureEditor

from pdf_plotter.ui.models \
    import Dataset, CorrectedDatasets

from pdf_plotter.ui.controls \
    import Controls

from pdf_plotter.ui.nodes.dataset \
    import DatasetNodeControls, DatasetNodeButtons

from pdf_plotter.ui.nodes.corrected_datasets \
    import CorrectedDatasetsNodeControls, CorrectedDatasetsNodeButtons

from pdf_plotter.io.nexus_load import ExperimentFileInput

# -----------------------------------------------------------#
# Figure Model

SofqPlotView = View(
    Item('figure',
         editor=MPLFigureEditor(),
         show_label=False,
         resizable=True,
         springy=True
         ),
)


class SofqPlot(HasTraits):
    # View
    view = SofqPlotView

    # Cached xlims
    xlims = None

    # Cached ylims
    ylims = None

    # Figure to display selected dataset and axes for figure
    figure = Instance(Figure, ())


# -----------------------------------------------------------#
# Main Control Panel View

ControlPanelView = View(
    HSplit(
        UItem('sofq_plot', width=0.8, style='custom', editor=InstanceEditor()),
        VSplit(
            UItem('experiment_file',
                  height=0.1,
                  style='custom',
                  editor=InstanceEditor()
                  ),
            UItem('controls',
                  height=0.9,
                  style='custom',
                  editor=InstanceEditor()
                  ),
        ),
    ),
    resizable=True,
    statusbar=[StatusItem(name='load_status')]
)

# -----------------------------------------------------------#
# Main Control Panel Controller


class ControlPanelHandler(Handler):
    def get_parents(self, info, node):
        # Grab selected Dataset
        selected = info.object.selected

        # Get the TreeEditor for the Experiment
        controls_editors = info.ui.get_editors("controls")
        experiment_editors = list()
        for editor in controls_editors:
            experiment_editors.extend(editor._ui.get_editors("experiment"))
        experiment_editor = experiment_editors[0]  # just grab first

        # Get the parents
        corrected_dataset = experiment_editor.get_parent(selected)
        measurement = experiment_editor.get_parent(corrected_dataset)
        experiment = experiment_editor.get_parent(measurement)

        parents = {'corrected_dataset': corrected_dataset,
                   'measurement': measurement,
                   'experiment': experiment}

        return parents

    def object_button_pressed_changed(self, info):

        # Map of button's name to the function it calls
        name2func = {'cache_plot': self.cache_plot,
                     'cache_plots': self.cache_plots,
                     'clear_cache': self.clear_cache,
                     }

        #
        if info.initialized:
            if info.object.controls.node_buttons.button_event:
                name = info.object.controls.node_buttons.button_name
                button_func = name2func[name]
                button_func(info)

    def cache_plot(self, info):
        dataset = info.object.current_dataset

        # Get info for selected Dataset (=a) and create new Dataset (=b)
        a = dataset
        shift = info.object.controls.node_controls.shift_factor
        scale = info.object.controls.node_controls.scale_factor
        b = Dataset(x=a.x, 
                    y=scale * a.y + shift, 
                    xmin_filter=a.xmin_filter,
                    xmax_filter=a.xmax_filter,
                    title=a.title)

        # Apply x-range filter
        b.x, b.y = info.object.controls.node_controls.filter_xrange(b.x, b.y, b)

        # If we have modified Dataset 'a', change title of 'b' for
        # differences in...
        tmp_title = str(b.title)

        try:
            # Shift
            if shift != 0.0:
                b.title += " shift: {0:>5.2f}".format(shift)

            # Scale
            if scale != 1.0:
                b.title += " scale: {0:>5.2f}".format(scale)

            # Xmin
            if min(b.x) != min(a.x):
                b.title += " xmin: {0:.2f}".format(min(b.x))

            # Xmax
            if max(b.x) != max(a.x):
                b.title += " xmax: {0:.2f}".format(max(b.x))

            # Check if title changed.
            # If so, add as a different Node in 'Other' Measurement
            if tmp_title != b.title:
                parents = self.get_parents(info, b)
                info.object.controls.add_plot_to_node(dataset=b,
                                                      parents=parents)

            # Add 'b' to cached plots
            info.object.controls.cached_plots.append(b)

            # Add to plot and refresh
            axes = info.object.get_axes()
            line = axes.plot(b.x, b.y, label=b.title)
            info.object.line2dataset[line[0]] = b
            info.object._setupColorMap(info.object.controls.cached_plots)
            info.object.plot_cached()
            info.object.plot_modification_to_selected_change()

        except ValueError:
            pass

    def cache_plots(self, info):
        datasets = info.object.controls.selected.datasets
        axes = info.object.get_axes()
        for dataset in datasets:
            x, y = info.object.controls.node_controls.filter_xrange(dataset.x, dataset.y, dataset)
            line = axes.plot(x, y, label=dataset.title)
            info.object.line2dataset[line[0]] = dataset
            info.object.controls.cached_plots.append(Dataset(x=x,y=y,title=dataset.title))

        info.object._setupColorMap(info.object.controls.cached_plots)
        info.object.plot_cached()
        info.object.plot_modification_to_selected_change()

    def clear_cache(self, info):
        info.object.controls.cached_plots = []
        axes = info.object.get_axes()
        axes.cla()
        info.object.plot_selected_node_change()

    def object_selected_changed(self, info):
        if not info.initialized:
            return

        # Keep selected color map
        selected_cmap = info.object.controls.node_controls.selected_cmap
        selected_node = info.object.controls.selected

        # Get Experiment xmin and xmax for min and max slider limits
        exp_xmin = info.object.controls.exp_xmin
        exp_xmax = info.object.controls.exp_xmax

        # Keep x,y lock on axes
        freeze_xlims = info.object.controls.node_controls.freeze_xlims
        freeze_ylims = info.object.controls.node_controls.freeze_ylims

        # Inform that we initializating node_controls to turn off trait notifications
        info.object.controls.node_controls.initialized = False

        if isinstance(info.object.selected, Dataset):
            if freeze_xlims:
                xmin = info.object.controls.node_controls.xmin
                xmax = info.object.controls.node_controls.xmax

            else:
                xmin = info.object.selected.xmin_filter
                xmax = info.object.selected.xmax_filter

            info.object.controls.node_controls = DatasetNodeControls(
                selected=selected_node,
                xmin=xmin,
                xmin_min=exp_xmin,
                xmin_max=exp_xmax,
                xmax=xmax,
                xmax_min=exp_xmin,
                xmax_max=exp_xmax,
                selected_cmap=selected_cmap,
                freeze_xlims=freeze_xlims,
                freeze_ylims=freeze_ylims,
            )

            info.object.controls.node_buttons = DatasetNodeButtons()

        elif isinstance(info.object.selected, CorrectedDatasets):
            if freeze_xlims:
                xmin = info.object.controls.node_controls.xmin
                xmax = info.object.controls.node_controls.xmax

            else:
                xmin = info.object.selected.datasets[0].xmin_filter
                xmax = info.object.selected.datasets[0].xmax_filter

            datasets = [ (dataset, dataset.title) for dataset in info.object.selected.datasets ]

            info.object.controls.node_controls = CorrectedDatasetsNodeControls(
                selected=selected_node,
                datasets=datasets,
                xmin=xmin,
                xmax=xmax,
                selected_cmap=selected_cmap,
                freeze_xlims=freeze_xlims,
                freeze_ylims=freeze_ylims,
            )

            info.object.controls.node_buttons = CorrectedDatasetsNodeButtons()

        # Inform that we are done with initialization to turn back on trait notifications
        info.object.controls.node_controls.initialized = True

# -----------------------------------------------------------#
# Main Control Panel


class ControlPanel(HasTraits):

    # -------------------------------------------------------#
    # Traits

    experiment_file = Instance(ExperimentFileInput)

    # S(Q) Plot
    sofq_plot = Instance(SofqPlot, ())

    # Controls for adjusting plots
    controls = Instance(Controls)

    # Status
    load_status = Property(depends_on='experiment_file.load_status')

    # Selected node
    selected = Property(depends_on='controls.selected')

    # Selected dataset
    current_dataset = Property(depends_on='selected')

    # Button event fired
    button_pressed = Property(depends_on='controls.node_buttons.button_event')

    # Current colors of contents that is setup by _setupColorMap
    _colors = list()

    # Stuff for plots
    cache_start_index = 0
    dataset = Instance(Dataset)
    corrected_datasets = Instance(CorrectedDatasets)
    line2dataset = Dict

    # Plot's overall x, y min and max
    plot_xmin = CFloat
    plot_xmax = CFloat
    plot_ymin = CFloat
    plot_ymax = CFloat

    # -------------------------------------------------------#
    # Utilities

    def _get_selected(self):
        return self.controls.selected

    def _get_current_dataset(self):
        return self.controls.current_dataset

    # Updates the displayed load status based on notifications from Controls
    def _get_load_status(self):
        return ("%s" % self.experiment_file.load_status)

    # Updates the buttons fired from Controls
    def _get_button_pressed(self):
        return self.controls.node_buttons

    # Initialize the ColorMap selected for the list of plots
    def _setupColorMap(self, plot_list, num_lines=8, reset_cycler=True):
        # If there are more lines than
        if num_lines < len(plot_list):
            num_lines = len(plot_list)

        # Reset the cycler so we don't get the same colors in order of the
        # ColorMap
        if reset_cycler:
            self.sofq_plot.figure.gca().set_prop_cycle(None)

        # Get the currently selected ColorMap from the Controls
        myCMap = cm.get_cmap(
            name=self.controls.node_controls.selected_cmap_contents)

        # Create a value for each line that corresponds to the color in the
        # ColorMap, ranged from 0->1
        cm_subsection = np.linspace(0.0, 1.0, num_lines)

        # Pull the colors from the ColorMap selected using the values created
        # in cm_subsection
        self._colors = [myCMap(x) for x in cm_subsection]

    # Sets the limits for the plots in the figure (cached and selected)
    def _get_limits_on_plot(self, x, y, not_selected=None):

        xlist = list()
        ylist = list()

        # Get current selected / input x values
        xlist = np.append(xlist, x)
        ylist = np.append(ylist, y)

        # Get cached plot x values
        if len(self.controls.cached_plots) > 0:
            for a in self.controls.cached_plots:

                xlist = np.append(xlist, a.x, axis=None)
                ylist = np.append(ylist, a.y, axis=None)

        # Get all datasets not selected but plotted
        if not_selected:
            for a in not_selected:
                x, y = self.controls.node_controls.filter_xrange(a.x, a.y, a)
                xlist = np.append(xlist, x, axis=None)
                ylist = np.append(ylist, y, axis=None)

        # Get xmin and xmax specified by the Controls
        xlist = np.append(xlist, self.controls.node_controls.xmin)
        xlist = np.append(xlist, self.controls.node_controls.xmax)

        # Convert infinities and NaN to 0.0
        ylist[ylist == np.inf] = 0.0
        ylist[ylist == -np.inf] = 0.0
        ylist[np.isnan(ylist)] = 0.0

        # Get plot x, y min and max
        if len(xlist) > 0:
            self.plot_xmin = min(xlist)
            self.plot_xmax = max(xlist)

        if len(ylist) > 0:
            self.plot_ymin = min(ylist)
            self.plot_ymax = max(ylist)

    # Get the axes we are plotting on
    def get_axes(self):
        # If no Axes for the figure, 1) Make the figure zoomable by mouse wheel
        # and 2) initialize an Axes
        if len(self.sofq_plot.figure.axes) == 0:
            self.sofq_plot.figure.pan_zoom = ZoomOnWheel(self.sofq_plot.figure)
            axes = self.sofq_plot.figure.add_subplot(111)

        # Get only Axes currently available
        else:
            axes = self.sofq_plot.figure.axes[0]

        return axes

    # Clear the plot
    def clear_plot(self):
        # Get the Axes
        axes = self.get_axes()

        # Save axes if we are freezing the view
        if self.controls.node_controls.freeze_xlims:
            self.sofq_plot.xlims = axes.get_xlim()
        if self.controls.node_controls.freeze_ylims:
            self.sofq_plot.ylims = axes.get_ylim()

        axes.cla()

        # Clear the line -> dataset dictionary
        self.line2dataset = dict()

    # Get the datasets we want to plot
    def get_datasets_to_plot(self):
        # Get the Datasets
        if isinstance(self.controls.selected, Dataset):
            datasets = [ self.controls.selected ]

        elif isinstance(self.controls.selected, CorrectedDatasets):
            datasets = self.controls.selected.datasets

        else:
            datasets = list()

        return datasets

    def apply_selected_mask(self, datasets):
        selected_dataset = self.controls.current_dataset

        selected = None
        not_selected = list()
        for dataset in datasets:
            if dataset == selected_dataset:
                selected = dataset
            else:
                not_selected.append(dataset)

        return selected, not_selected

    def plot_all_datasets(self,datasets):
        # Get the Axes
        axes = self.get_axes()

        for i, dataset in enumerate(datasets):
            x, y = self.controls.node_controls.filter_xrange(dataset.x, dataset.y, dataset)
            line = axes.plot(x, y, label=dataset.title)
            self.line2dataset[line[0]] = dataset

    # Plots the X, Y data (cached and selected) on the given axes
    def plot_selected(self, x, y, dataset):
        # Get axes to plot on
        axes = self.get_axes()

        # Add first line to plot if none already
        if not axes.lines:
            line = axes.plot(x, y, label=dataset.title)
            self.line2dataset[line[0]] = dataset

        # Get the corresponding  line to the dataset selected
        for line in axes.lines:
            if dataset == self.line2dataset[line]:
                l = line

        # Plot xy
        l.set_data(x, y)
        #l.set_color('b')
        l.set_marker('o')
        l.set_markersize(4.0)
        #l.set_linestyle('--')
        l.set_label(dataset.title)

    # Add the cached lines back to the plot (style taken care of in
    # plot_cached)
    def add_cached(self):
        axes = self.get_axes()
        for cached_plot in self.controls.cached_plots:
            line = axes.plot(cached_plot.x, cached_plot.y)
            self.line2dataset[line[0]] = cached_plot

    # Loop over lines in the Axes for the cached plots and modify
    # accordingly
    def plot_cached(self):
        axes = self.get_axes()
        for i, (cached_plot, color) in enumerate(
                zip(self.controls.cached_plots, self._colors)):
            l = axes.lines[i + self.cache_start_index]
            l.set_data(cached_plot.x, cached_plot.y)
            l.set_color(color)
            l.set_marker(None)
            l.set_linestyle('-')
            l.set_label(cached_plot.title)

    # Apply adjustments to the currently selected dataset
    def adjust_selected(self, selected, not_selected=None):
        axes = self.get_axes()
        scale = self.controls.node_controls.scale_factor
        shift = self.controls.node_controls.shift_factor

        # Set the x-range filter
        selected.xmin_filter = self.controls.node_controls.xmin
        selected.xmax_filter = self.controls.node_controls.xmax

        # Get x and y. Apply scale and shift to y
        x = selected.x
        y = scale * (selected.y) + shift

        # Apply x-range filter
        x, y = self.controls.node_controls.filter_xrange(x, y, selected)

        # Get the X, Y limits from all plots (selected + cached)
        self._get_limits_on_plot(x, y, not_selected)

        # Set the limits
        if not self.controls.node_controls.freeze_xlims \
           and self.plot_xmin <= self.plot_xmax:
            axes.set_xlim(self.plot_xmin, self.plot_xmax)

        if not self.controls.node_controls.freeze_ylims \
           and  self.plot_ymin <= self.plot_ymax:
            axes.set_ylim(self.plot_ymin, self.plot_ymax)

        return x, y, selected.title

    # After updating the figure with all changes, redraw the plot
    def redraw_canvas(self):
        axes = self.get_axes()

        # Make the Legend: draggable, wheel-adjusted font, right-click hidable,
        leg = axes.legend()
        leg = DraggableLegend(leg)

        # Re-draw the Canvas of the Figure with our changes
        canvas = self.sofq_plot.figure.canvas
        if canvas is not None:
            canvas.draw()

    # -------------------------------------------------------#
    # Dynamic

    # If the experiment changes, initialize the controls if needed
    @on_trait_change('experiment_file.experiment')
    def experiment_updated(self):
        experiment = self.experiment_file.experiment
        if experiment:
            if self.controls:
                self.controls.experiment = experiment
            else:
                self.controls = Controls(
                    experiment=experiment,
                    node_controls=DatasetNodeControls(),
                    node_buttons=DatasetNodeButtons())

    # Upon loading the experiment file, watch for updates from the multi-thread load
    @on_trait_change('experiment_file.load_status')
    def update_status(self):
        self.status = self.experiment_file.load_status

    # If the 'lock x-axis' radio button is pressed, save the current x-axis limits
    @on_trait_change('controls.node_controls.freeze_xlims')
    def cache_xlims(self):
        axes = self.get_axes()
        self.sofq_plot.xlims = axes.get_xlim()

    # If the 'lock y-axis' radio button is pressed, save the current y-axis limits
    @on_trait_change('controls.node_controls.freeze_ylims')
    def cache_ylims(self):
        axes = self.get_axes()
        self.sofq_plot.ylims = axes.get_ylim()

    # Re-plot when selecting another node or if we change the ColorMap
    @on_trait_change('controls:selected,controls:node_controls:[selected_cmap,dataset_selected]')
    def plot_selected_node_change(self):
        try:

            self.clear_plot()

            # Force an update in the selected node
            self.controls.updated_selected_in_node_controls()
            self.controls.get_current_dataset()

            # Update node_controls xmin xmax filters
            self.controls.node_controls.initialized = False

            # Update xmin filter
            self.controls.node_controls.xmin = self.current_dataset.xmin_filter
            self.controls.node_controls.xmin_min = self.current_dataset.xmin_filter
            self.controls.node_controls.xmin_max = self.current_dataset.xmax_filter

            # Update xmax filter
            self.controls.node_controls.xmax = self.current_dataset.xmax_filter
            self.controls.node_controls.xmax_min = self.current_dataset.xmin_filter
            self.controls.node_controls.xmax_max = self.current_dataset.xmax_filter

            # Get a list of datasets we will plot
            datasets = self.get_datasets_to_plot()

            # Go ahead and initially plot all the datasets
            self.plot_all_datasets(datasets)

            # Separate the list into selected and not selected datasets
            selected, not_selected = self.apply_selected_mask(datasets)

            # Adjust/modify the currently selected dataset
            x, y, title = self.adjust_selected(selected, not_selected)

            self.controls.node_controls.initialized = True

            # Plot the selected dataset after modification
            self.plot_selected(x, y, selected)

             # Add cached lines back to plot
            self.cache_start_index = len(datasets)
            self.add_cached()

            # Setup the ColorMap for the cache plots
            self._setupColorMap(self.controls.cached_plots)

            # Plot the cached lines
            self.plot_cached()

            # Redraw the canvas of the figure
            self.redraw_canvas()


        except AttributeError:
            pass

    # Re-plot when adjusting the currently selected node
    @on_trait_change('controls:node_controls:[scale_factor,shift_factor],'
                     'controls.node_controls.[xmin,xmax]')
    def plot_modification_to_selected_change(self):

        try:
            if not self.controls.node_controls.initialized:
                return
            datasets = self.get_datasets_to_plot()

            selected, not_selected = self.apply_selected_mask(datasets)

            x, y, title = self.adjust_selected(selected, not_selected)

            self.plot_selected(x, y, selected)

            # Redraw the canvas of the figure
            self.redraw_canvas()

        except AttributeError:
            pass


if __name__ == "__main__":
    cp = ControlPanel(experiment_file=ExperimentFileInput())
    cp.configure_traits(
        view=ControlPanelView,
        handler=ControlPanelHandler)
