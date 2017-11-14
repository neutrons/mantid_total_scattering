import logging
import math
import weakref

import numpy

from pyface.qt import QtGui, QtCore
from traits.etsconfig.api \
    import ETSConfig
ETSConfig.toolkit = "qt4"

if ETSConfig.toolkit == "qt4":
    from traitsui.qt4.editor import Editor
    from traitsui.qt4.basic_editor_factory import BasicEditorFactory
elif ETSConfig.toolkit == "wx":
    from traitsui.wx.editor import Editor
    from traitsui.qt4.basic_editor_factory import BasicEditorFactory

import matplotlib.pyplot as _plt
from matplotlib.backends.backend_qt4agg \
    import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg \
    import NavigationToolbar2QT


# -----------------------------------------------------------#
# Matplotlib w/ Qt4 classes for TraitsUI Editor

class _MPLFigureEditor(Editor):

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas
        frame = QtGui.QWidget()
        mpl_canvas = FigureCanvas(self.value)
        mpl_canvas.setParent(frame)
        mpl_toolbar = NavigationToolbar2QT(mpl_canvas, frame)

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(mpl_canvas)
        vbox.addWidget(mpl_toolbar)
        frame.setLayout(vbox)

        return frame


class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor


# -----------------------------------------------------------#
# Matplotlib Interactions

class MplInteraction(object):
    """Base class for class providing interaction to a matplotlib Figure."""

    def __init__(self, figure):
        """Initializer

        :param Figure figure: The matplotlib figure to attach the behavior to.
        """
        self._fig_ref = weakref.ref(figure)
        self._cids = []

    def __del__(self):
        self.disconnect()

    def _add_connection(self, event_name, callback):
        """Called to add a connection to an event of the figure

        :param str event_name: The matplotlib event name to connect to.
        :param callback: The callback to register to this event.
        """
        cid = self.figure.canvas.mpl_connect(event_name, callback)
        self._cids.append(cid)

    def disconnect(self):
        """Disconnect interaction from Figure."""
        if self._fig_ref is not None:
            figure = self._fig_ref()
            if figure is not None:
                for cid in self._cids:
                    figure.canvas.mpl_disconnect(cid)
            self._fig_ref = None

    @property
    def figure(self):
        """The Figure this interaction is connected to or
        None if not connected."""
        return self._fig_ref() if self._fig_ref is not None else None

    def _axes_to_update(self, event):
        """Returns two sets of Axes to update according to event.

        Takes care of multiple axes and shared axes.

        :param MouseEvent event: Matplotlib event to consider
        :return: Axes for which to update xlimits and ylimits
        :rtype: 2-tuple of set (xaxes, yaxes)

        """
        x_axes, y_axes = set(), set()

        # Go through all axes to enable zoom for multiple axes subplots
        for ax in self.figure.axes:
            if ax.contains(event)[0]:
                # For twin x axes, makes sure the zoom is applied once
                shared_x_axes = set(ax.get_shared_x_axes().get_siblings(ax))
                if x_axes.isdisjoint(shared_x_axes):
                    x_axes.add(ax)

                # For twin y axes, makes sure the zoom is applied once
                shared_y_axes = set(ax.get_shared_y_axes().get_siblings(ax))
                if y_axes.isdisjoint(shared_y_axes):
                    y_axes.add(ax)

        return x_axes, y_axes

    def _draw(self):
        """Conveninent method to redraw the figure"""
        self.figure.canvas.draw()


class ZoomOnWheel(MplInteraction):
    """Class providing zoom on wheel interaction to a matplotlib Figure.

    Supports subplots, twin Axes and log scales.
    """

    def __init__(self, figure=None, scale_factor=1.1):
        """Initializer

        :param Figure figure: The matplotlib figure to attach the behavior to.
        :param float scale_factor: The scale factor to apply on wheel event.
        """
        super(ZoomOnWheel, self).__init__(figure)
        self._add_connection('scroll_event', self._on_mouse_wheel)

        self.scale_factor = scale_factor

    @staticmethod
    def _zoom_range(begin, end, center, scale_factor, scale):
        """Compute a 1D range zoomed around center.

        :param float begin: The begin bound of the range.
        :param float end: The end bound of the range.
        :param float center: The center of the zoom (i.e., invariant point)
        :param float scale_factor: The scale factor to apply.
        :param str scale: The scale of the axis
        :return: The zoomed range (min, max)
        """
        if begin < end:
            min_, max_ = begin, end
        else:
            min_, max_ = end, begin

        if scale == 'linear':
            old_min, old_max = min_, max_
        elif scale == 'log':
            old_min = numpy.log10(min_ if min_ > 0. else numpy.nextafter(0, 1))
            center = numpy.log10(
                center if center > 0. else numpy.nextafter(0, 1))
            old_max = numpy.log10(max_) if max_ > 0. else 0.
        else:
            logging.warning(
                'Zoom on wheel not implemented for scale "%s"' % scale)
            return begin, end

        offset = (center - old_min) / (old_max - old_min)
        range_ = (old_max - old_min) / scale_factor
        new_min = center - offset * range_
        new_max = center + (1. - offset) * range_

        if scale == 'log':
            try:
                new_min, new_max = 10. ** float(new_min), 10. ** float(new_max)
            except OverflowError:  # Limit case
                new_min, new_max = min_, max_
            if new_min <= 0. or new_max <= 0.:  # Limit case
                new_min, new_max = min_, max_

        if begin < end:
            return new_min, new_max
        else:
            return new_max, new_min

    def _on_mouse_wheel(self, event):
        # Skip if we are in the legend space
        eventInLegend, t = \
            self.figure.gca().get_legend().legendPatch.contains(event)
        if eventInLegend:
            return

        if event.step > 0:
            scale_factor = self.scale_factor
        else:
            scale_factor = 1. / self.scale_factor

        # Go through all axes to enable zoom for multiple axes subplots
        x_axes, y_axes = self._axes_to_update(event)

        for ax in x_axes:
            transform = ax.transData.inverted()
            xdata, ydata = transform.transform_point((event.x, event.y))

            xlim = ax.get_xlim()
            xlim = self._zoom_range(xlim[0], xlim[1],
                                    xdata, scale_factor,
                                    ax.get_xscale())
            ax.set_xlim(xlim)

        for ax in y_axes:
            ylim = ax.get_ylim()
            ylim = self._zoom_range(ylim[0], ylim[1],
                                    ydata, scale_factor,
                                    ax.get_yscale())
            ax.set_ylim(ylim)

        if x_axes or y_axes:
            self._draw()


class PanAndZoom(ZoomOnWheel):
    """Class providing pan & zoom interaction to a matplotlib Figure.

    Left button for pan, right button for zoom area and zoom on wheel.

    Support subplots, twin Axes and log scales.
    """

    def __init__(self, figure=None, scale_factor=1.1):
        """Initializer

        :param Figure figure: The matplotlib figure to attach the behavior to.
        :param float scale_factor: The scale factor to apply on wheel event.
        """
        super(PanAndZoom, self).__init__(figure, scale_factor)
        self._add_connection('button_press_event', self._on_mouse_press)
        self._add_connection('button_release_event', self._on_mouse_release)
        self._add_connection('motion_notify_event', self._on_mouse_motion)

        self._pressed_button = None  # To store active button
        self._axes = None  # To store x and y axes concerned by interaction
        self._event = None  # To store reference event during interaction

    @staticmethod
    def _pan_update_limits(ax, axis_id, event, last_event):
        """Compute limits with applied pan."""
        assert axis_id in (0, 1)
        if axis_id == 0:
            lim = ax.get_xlim()
            scale = ax.get_xscale()
        else:
            lim = ax.get_ylim()
            scale = ax.get_yscale()

        pixel_to_data = ax.transData.inverted()
        data = pixel_to_data.transform_point((event.x, event.y))
        last_data = pixel_to_data.transform_point((last_event.x, last_event.y))

        if scale == 'linear':
            delta = data[axis_id] - last_data[axis_id]
            new_lim = lim[0] - delta, lim[1] - delta
        elif scale == 'log':
            try:
                delta = math.log10(data[axis_id]) - \
                    math.log10(last_data[axis_id])
                new_lim = [pow(10., (math.log10(lim[0]) - delta)),
                           pow(10., (math.log10(lim[1]) - delta))]
            except (ValueError, OverflowError):
                new_lim = lim  # Keep previous limits
        else:
            logging.warning('Pan not implemented for scale "%s"' % scale)
            new_lim = lim
        return new_lim

    def _pan(self, event):
        if event.name == 'button_press_event':  # begin pan
            self._event = event

        elif event.name == 'button_release_event':  # end pan
            self._event = None

        elif event.name == 'motion_notify_event':  # pan
            if self._event is None:
                return

            if event.x != self._event.x:
                for ax in self._axes[0]:
                    xlim = self._pan_update_limits(ax, 0, event, self._event)
                    ax.set_xlim(xlim)

            if event.y != self._event.y:
                for ax in self._axes[1]:
                    ylim = self._pan_update_limits(ax, 1, event, self._event)
                    ax.set_ylim(ylim)

            if event.x != self._event.x or event.y != self._event.y:
                self._draw()

            self._event = event

    def _zoom_area(self, event):
        if event.name == 'button_press_event':  # begin drag
            self._event = event
            self._patch = _plt.Rectangle(
                xy=(event.xdata, event.ydata), width=0, height=0,
                fill=False, linewidth=1., linestyle='solid', color='black')
            self._event.inaxes.add_patch(self._patch)

        elif event.name == 'button_release_event':  # end drag
            self._patch.remove()
            del self._patch

            if (abs(event.x - self._event.x) < 3 or
                    abs(event.y - self._event.y) < 3):
                return  # No zoom when points are too close

            x_axes, y_axes = self._axes

            for ax in x_axes:
                pixel_to_data = ax.transData.inverted()
                begin_pt = pixel_to_data.transform_point((event.x, event.y))
                end_pt = pixel_to_data.transform_point(
                    (self._event.x, self._event.y))

                min_ = min(begin_pt[0], end_pt[0])
                max_ = max(begin_pt[0], end_pt[0])
                if not ax.xaxis_inverted():
                    ax.set_xlim(min_, max_)
                else:
                    ax.set_xlim(max_, min_)

            for ax in y_axes:
                pixel_to_data = ax.transData.inverted()
                begin_pt = pixel_to_data.transform_point((event.x, event.y))
                end_pt = pixel_to_data.transform_point(
                    (self._event.x, self._event.y))

                min_ = min(begin_pt[1], end_pt[1])
                max_ = max(begin_pt[1], end_pt[1])
                if not ax.yaxis_inverted():
                    ax.set_ylim(min_, max_)
                else:
                    ax.set_ylim(max_, min_)

            self._event = None

        elif event.name == 'motion_notify_event':  # drag
            if self._event is None:
                return

            if event.inaxes != self._event.inaxes:
                return  # Ignore event outside plot

            self._patch.set_width(event.xdata - self._event.xdata)
            self._patch.set_height(event.ydata - self._event.ydata)

        self._draw()

    def _on_mouse_press(self, event):
        if self._pressed_button is not None:
            return  # Discard event if a button is already pressed

        if event.button in (1, 3):  # Start
            x_axes, y_axes = self._axes_to_update(event)
            if x_axes or y_axes:
                self._axes = x_axes, y_axes
                self._pressed_button = event.button

                if self._pressed_button == 1:  # pan
                    self._pan(event)
                elif self._pressed_button == 3:  # zoom area
                    self._zoom_area(event)

    def _on_mouse_release(self, event):
        if self._pressed_button == event.button:
            if self._pressed_button == 1:  # pan
                self._pan(event)
            elif self._pressed_button == 3:  # zoom area
                self._zoom_area(event)
            self._pressed_button = None

    def _on_mouse_motion(self, event):
        if self._pressed_button == 1:  # pan
            self._pan(event)
        elif self._pressed_button == 3:  # zoom area
            self._zoom_area(event)


def figure_pz(*args, **kwargs):
    """matplotlib.pyplot.figure with pan and zoom interaction"""
    fig = _plt.figure(*args, **kwargs)
    fig.pan_zoom = PanAndZoom(fig)
    return fig


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    fig = figure_pz()
    # Alternative:
    # fig = plt.figure()
    # pan_zoom = PanAndZoom(fig)

    nrow, ncol = 2, 3

    ax1 = fig.add_subplot(nrow, ncol, 1)
    ax1.set_title('basic')
    ax1.plot((1, 2, 3))

    ax2 = fig.add_subplot(nrow, ncol, 2)
    ax2.set_title('log + twinx')
    ax2.set_yscale('log')
    ax2.plot((1, 2, 1))

    ax2bis = ax2.twinx()
    ax2bis.plot((3, 2, 1), color='red')

    ax3 = fig.add_subplot(nrow, ncol, 3)
    ax3.set_title('inverted y axis')
    ax3.plot((1, 2, 3))
    lim = ax3.get_ylim()
    ax3.set_ylim(lim[1], lim[0])

    ax4 = fig.add_subplot(nrow, ncol, 4)
    ax4.set_title('keep ratio')
    ax4.axis('equal')
    ax4.imshow(numpy.arange(100).reshape(10, 10))

    ax5 = fig.add_subplot(nrow, ncol, 5)
    ax5.set_xlabel('symlog scale + twiny')
    ax5.set_xscale('symlog')
    ax5.plot((1, 2, 3))
    ax5bis = ax5.twiny()
    ax5bis.plot((3, 2, 1), color='red')

    # The following is taken from:
    # http://matplotlib.org/examples/axes_grid/demo_curvelinear_grid.html
    from mpl_toolkits.axisartist import Subplot
    from mpl_toolkits.axisartist.grid_helper_curvelinear import \
        GridHelperCurveLinear

    def tr(x, y):  # source (data) to target (rectilinear plot) coordinates
        x, y = numpy.asarray(x), numpy.asarray(y)
        return x + 0.2 * y, y - x

    def inv_tr(x, y):
        x, y = numpy.asarray(x), numpy.asarray(y)
        return x - 0.2 * y, y + x

    grid_helper = GridHelperCurveLinear((tr, inv_tr))

    ax6 = Subplot(fig, nrow, ncol, 6, grid_helper=grid_helper)
    fig.add_subplot(ax6)
    ax6.set_title('non-ortho axes')

    xx, yy = tr([3, 6], [5.0, 10.])
    ax6.plot(xx, yy)

    ax6.set_aspect(1.)
    ax6.set_xlim(0, 10.)
    ax6.set_ylim(0, 10.)

    ax6.axis["t"] = ax6.new_floating_axis(0, 3.)
    ax6.axis["t2"] = ax6.new_floating_axis(1, 7.)
    ax6.grid(True)

    plt.show()


class DraggableLegend:
    def __init__(self, legend):
        self.legend = legend
        self.gotLegend = False
        legend.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        legend.figure.canvas.mpl_connect('pick_event', self.on_picker)
        legend.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        legend.set_picker(self.my_legend_picker)

    # ----------------------------------------------------#
    # Connected event handlers

    def on_motion(self, event):
        if self.gotLegend:
            dx = event.x - self.mouse_x
            dy = event.y - self.mouse_y
            loc_in_canvas = self.legend_x + dx, self.legend_y + dy
            loc_in_norm_axes = \
                self.legend.parent.transAxes.inverted().transform_point(loc_in_canvas)
            self.legend._loc = tuple(loc_in_norm_axes)
            self.legend.figure.canvas.draw()

    def my_legend_picker(self, legend, event):
        return self.legend.legendPatch.contains(event)

    def on_picker(self, event):
        if event.artist == self.legend:
            # left-click
            if event.mouseevent.button == 1:
                self._move_legend(event)

            # mouse button pressed
            if event.mouseevent.button == 2:
                pass

            # right-click
            if event.mouseevent.button == 3:
                self._hideLegend()

            # mouse up
            if event.mouseevent.button == 'up':
                self._scaleUpLegendFont()

            # mouse down
            if event.mouseevent.button == 'down':
                self._scaleDownLegendFont()

    def on_release(self, event):
        if self.gotLegend:
            self.gotLegend = False

    # ----------------------------------------------------#
    # Utility functions

    def _move_legend(self, event):
        bbox = self.legend.get_window_extent()
        self.mouse_x = event.mouseevent.x
        self.mouse_y = event.mouseevent.y
        self.legend_x = bbox.xmin
        self.legend_y = bbox.ymin
        self.gotLegend = 1

    def _scaleUpLegendFont(self, size_step=4):
        size = self.legend.get_texts()[0].get_fontsize()
        size += size_step
        _plt.setp(self.legend.get_texts(),
                  fontsize=size)  # legend 'list' fontsize
        self.legend.figure.canvas.draw()

    def _scaleDownLegendFont(self, size_step=4):
        size = self.legend.get_texts()[0].get_fontsize()
        size -= size_step
        _plt.setp(self.legend.get_texts(),
                  fontsize=size)  # legend 'list' fontsize
        self.legend.figure.canvas.draw()

    def _hideLegend(self):
        if self.legend.get_visible():
            self.legend.set_visible(False)
        else:
            self.legend.set_visible(True)
        self.legend.figure.canvas.draw()
