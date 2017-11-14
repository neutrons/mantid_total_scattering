#!/usr/bin/env python

from traits.api \
    import Button, Any, List, Dict, on_trait_change

from traitsui.api \
    import View, HGroup, HSplit, VGroup, VSplit, Item, \
    CheckListEditor, TextEditor, RangeEditor

from pdf_plotter.ui.nodes.base_node \
    import NodeButtonHandler, NodeButtons, NodeControls

# -----------------------------------------------------------#
# CorrectedDatasets Node Buttons


class CorrectedDatasetsNodeButtonHandler(NodeButtonHandler):
    def object_cache_button_changed(self, info):
        info.object.button_name = 'cache_plots'
        self.trigger_button_event(info)

    def object_clear_cache_button_changed(self, info):
        info.object.button_name = 'clear_cache'
        self.trigger_button_event(info)


class CorrectedDatasetsNodeButtons(NodeButtons):
    cache_button = Button
    clear_cache_button = Button

    traits_view = View(
        HGroup(
            Item('cache_button',
                 label="Cache Plots",
                 show_label=False,
                 ),
            Item('clear_cache_button',
                 label="Clear Cache",
                 show_label=False,
                 ),
        ),
        handler=CorrectedDatasetsNodeButtonHandler(),
    )


# -----------------------------------------------------------#
# CorrectedDatasets Node Controls

class CorrectedDatasetsNodeControls(NodeControls):

    # List of CorrectedDatasets Node's datasets
    datasets   = List

    # Selected Dataset
    dataset_selected = Any

    traits_view = View(

        VGroup(
            # List of Datasets
            HSplit(
                Item('dataset_selected',
                     editor=CheckListEditor(name='datasets'),
                     show_label=False,
                     ),
                show_border=True,
                label='Datasets',
            ),

            # Color map for cached plots
            HSplit(
                Item('selected_cmap',
                     editor=CheckListEditor(name='cmap_list'),
                     show_label=False,
                     ),
                show_border=True,
                label='CacheColorMap',
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
                         editor=TextEditor(auto_set=False,),
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
                         editor=TextEditor(auto_set=False,),
                         label='Max',
                         ),
                    label='Xmin',
                ),

                # Xmax
                HSplit(
                    Item('xmax_min',
                         width=0.1,
                         editor=TextEditor(auto_set=False,),
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
                         editor=TextEditor(auto_set=False,),
                         label='Max',
                         ),
                    label='Xmax',
                ),
                show_border=True,
                label='X-range',
            ),
        ),
    )
