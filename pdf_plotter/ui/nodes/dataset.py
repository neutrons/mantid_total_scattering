#!/usr/bin/env python

from traits.api \
    import Button, CFloat

from traitsui.api \
    import View, VGroup, VSplit, HGroup, HSplit, Item, \
    RangeEditor, CheckListEditor

from pdf_plotter.ui.nodes.base_node \
    import NodeButtonHandler, NodeButtons, NodeControls


# -----------------------------------------------------------#
# Dataset Node Buttons

class DatasetNodeButtonHandler(NodeButtonHandler):
    def object_cache_button_changed(self, info):
        info.object.button_name = 'cache_plot'
        self.trigger_button_event(info)

    def object_clear_cache_button_changed(self, info):
        info.object.button_name = 'clear_cache'
        self.trigger_button_event(info)


class DatasetNodeButtons(NodeButtons):
    cache_button = Button
    clear_cache_button = Button

    traits_view = View(
        HGroup(
            Item('cache_button',
                 label="Cache Plot",
                 show_label=False,
                 ),
            Item('clear_cache_button',
                 label="Clear Cache",
                 show_label=False,
                 ),
        ),
        handler=DatasetNodeButtonHandler(),
    )


# -----------------------------------------------------------#
# Dataset Node Controls

class DatasetNodeControls(NodeControls):



    traits_view = View(
        VGroup(

            # Scale group
            HSplit(
                Item('scale_min', width=0.1, label='Min'),
                Item(
                    'scale_factor',
                    editor=RangeEditor(
                        mode='slider',
                        low_name='scale_min',
                        high_name='scale_max',
                        format='%4.2f',
                    ),
                    width=0.8,
                    show_label=False,
                ),
                Item('scale_max', width=0.1, label='Max'),
                show_border=True,
                label='Scale',
            ),

            # Shift group
            HSplit(
                Item('shift_min', width=0.1, label='Min'),
                Item(
                    'shift_factor',
                    editor=RangeEditor(
                        mode='slider',
                        low_name='shift_min',
                        high_name='shift_max',
                        format='%4.2f',
                    ),
                    width=0.8,
                    show_label=False,
                ),
                Item('shift_max', width=0.1, label='Max'),
                show_border=True,
                label='Shift',
            ),

            # X range
            VSplit(
                HGroup(
                    Item('freeze_xlims',
                         label="X axis",
                    ),
                    Item('freeze_ylims',
                         label="Y axis",
                    ),
                    show_border=True,
                    label='Lock',
                ),

                # Xmin
                HSplit(
                    Item('xmin_min',
                         width=0.1,
                         editor=RangeEditor(auto_set=False,mode='text'),
                         label='Min',
                         ),
                    Item('xmin',
                         editor=RangeEditor(
                             mode='slider',
                             low_name='xmin_min',
                             high_name='xmin_max',
                             format='%4.2f',
                         ),
                         width=0.8,
                         show_label=False,
                         ),
                    Item('xmin_max',
                         width=0.1,
                         editor=RangeEditor(auto_set=False,mode='text'),
                         label='Max',
                         ),
                    label='Xmin',
                ),

                # Xmax
                HSplit(
                    Item('xmax_min',
                         width=0.1,
                         editor=RangeEditor(auto_set=False,mode='text'),
                         label='Min',
                         ),
                    Item('xmax',
                         editor=RangeEditor(
                             mode='slider',
                             low_name='xmax_min',
                             high_name='xmax_max',
                             format='%4.2f',
                         ),
                         width=0.8,
                         show_label=False,
                         ),
                    Item('xmax_max',
                         width=0.1,
                         editor=RangeEditor(auto_set=False,mode='text'),
                         label='Max',
                         ),
                    label='Xmax',
                ),
                show_border=True,
                label='X-range',
            ),

            # Color map
            HSplit(
                Item('selected_cmap',
                     editor=CheckListEditor(name='cmap_list'),
                     show_label=False,
                     ),
                show_border=True,
                label='ColorMap',
            ),
        ),
    )
