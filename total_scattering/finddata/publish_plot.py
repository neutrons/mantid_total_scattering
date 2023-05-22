#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)
import json
import logging
import os

try:
    from postprocessing.Configuration import Configuration, CONFIG_FILE, CONFIG_FILE_ALTERNATE
except ImportError:
    # Local dev only, mocking Configuration as needed
    CONFIG_FILE = '/etc/autoreduce/post_processing.conf'
    CONFIG_FILE_ALTERNATE = 'post_processing.conf'

    class Configuration(object):
        """
        Read and process configuration file and provide an easy
        way to hold the various options for a client. This is a
        heavily abridged version of what is found in postprocessing.
        """
        def __init__(self, config_file):
            if os.access(config_file, os.R_OK) is False:
                raise RuntimeError("Configuration file doesn't exist or is not readable: %s" % config_file)
            with open(config_file, 'r') as cfg:
                json_encoded = cfg.read()
            config = json.loads(json_encoded)

            # Keep a record of which config file we are using
            self.config_file = config_file

            # plot publishing
            self.publish_url = config.get('publish_url_template', '')
            self.publisher_username = config.get('publisher_username', '')
            self.publisher_password = config.get('publisher_password', '')
            self.publisher_certificate = config.get('publisher_certificate', '')


def _determine_config_file(config_file):
    # put together the list of all choices
    choices = [config_file, CONFIG_FILE, CONFIG_FILE_ALTERNATE]

    # filter out bad choices
    choices = [name for name in choices
               if name is not None]
    choices = [name for name in choices
               if len(name) > 0]
    choices = [name for name in choices
               if os.access(name, os.R_OK)]

    # first one is a winner
    if len(choices) > 0:
        return choices[0]
    else:
        return None


def read_configuration(config_file=None):
    """
    Returns a new configuration object for a given
    configuration file
    @param config_file: configuration file to process
    """
    config_file = _determine_config_file(config_file)
    if config_file is None:
        raise RuntimeError('Failed to find Configuration file')

    logging.info('Loading configuration \'%s\'' % config_file)
    return Configuration(config_file)


def _loadDiv(filename):
    if not os.path.exists(filename):
        raise RuntimeError('\'%s\' does not exist' % filename)
    print('loading \'%s\'' % filename)
    with open(filename, 'r') as handle:
        div = handle.read()
    return div


def _getURL(url_template, instrument, run_number):
    import string
    url_template = string.Template(url_template)
    url = url_template.substitute(instrument=instrument,
                                  run_number=str(run_number))
    return url


def publish_plot(instrument, run_number, files, config=None):
    # read the configuration if one isn't provided
    if config is None:
        config = read_configuration()
    # verify that it has an attribute that matters
    try:
        config.publish_url
    except AttributeError:  # assume that it is a filename
        config = read_configuration(config)

    run_number = str(run_number)
    url = _getURL(config.publish_url, instrument, run_number)
    print('posting to \'%s\'' % url)

    # these next 2 lines are explicity bad - and doesn't seem
    # to do ANYTHING
    # https://urllib3.readthedocs.org/en/latest/security.html
    import urllib3
    urllib3.disable_warnings()

    import requests
    if config.publisher_certificate:
        request = requests.post(url, data={'username': config.publisher_username,
                                           'password': config.publisher_password},
                                files=files, cert=config.publisher_certificate)
    else:
        request = requests.post(url, data={'username': config.publisher_username,
                                           'password': config.publisher_password},
                                files=files, verify=False)
    return request


def plot1d(run_number, data_list, data_names=None, x_title='', y_title='',
           x_log=False, y_log=False, instrument='', show_dx=True, title='', publish=True):
    """
        Produce a 1D plot
        @param data_list: list of traces [ [x1, y1], [x2, y2], ...]
        @param data_names: name for each trace, for the legend
    """
    from plotly.offline import plot
    import plotly.graph_objs as go

    # Create traces
    if not isinstance(data_list, list):
        raise RuntimeError("plot1d: data_list parameter is expected to be a list")

    # Catch the case where the list is in the format [x y]
    data = []
    show_legend = False
    if len(data_list) == 2 and not isinstance(data_list[0], list):
        label = ''
        if isinstance(data_names, list) and len(data_names) == 1:
            label = data_names[0]
            show_legend = True
        data = [go.Scatter(name=label, x=data_list[0], y=data_list[1])]
    else:
        for i in range(len(data_list)):
            label = ''
            if isinstance(data_names, list) and len(data_names) == len(data_list):
                label = data_names[i]
                show_legend = True
            err_x = {}
            err_y = {}
            if len(data_list[i]) >= 3:
                err_y = dict(type='data', array=data_list[i][2], visible=True)
            if len(data_list[i]) >= 4:
                err_x = dict(type='data', array=data_list[i][3], visible=True)
                if show_dx is False:
                    err_x['thickness'] = 0
            data.append(go.Scatter(name=label, x=data_list[i][0], y=data_list[i][1],
                                   error_x=err_x, error_y=err_y))

    x_layout = dict(title=x_title, zeroline=False, exponentformat="power",
                    showexponent="all", showgrid=True,
                    showline=True, mirror="all", ticks="inside")
    if x_log:
        x_layout['type'] = 'log'
    y_layout = dict(title=y_title, zeroline=False, exponentformat="power",
                    showexponent="all", showgrid=True,
                    showline=True, mirror="all", ticks="inside")
    if y_log:
        y_layout['type'] = 'log'

    layout = go.Layout(
        showlegend=show_legend,
        autosize=True,
        width=600,
        height=400,
        margin=dict(t=40, b=40, l=80, r=40),  # noqa: E741
        hovermode='closest',
        bargap=0,
        xaxis=x_layout,
        yaxis=y_layout,
        title=title
    )

    fig = go.Figure(data=data, layout=layout)
    plot_div = plot(fig, output_type='div', include_plotlyjs=False, show_link=False)
    if publish:
        try:
            return publish_plot(instrument, run_number, files={'file': plot_div})
        except:  # noqa: E722
            logging.error("Publish plot failed: %s", sys.exc_value)
            return None
    else:
        return plot_div


def plot_heatmap(run_number, x, y, z, x_title='', y_title='', surface=False,
                 x_log=False, y_log=False, instrument='', title='', publish=True):
    """
        Produce a 2D plot
    """
    from plotly.offline import plot
    import plotly.graph_objs as go

    x_layout = dict(title=x_title, zeroline=False, exponentformat="power",
                    showexponent="all", showgrid=True,
                    showline=True, mirror="all", ticks="inside")
    if x_log:
        x_layout['type'] = 'log'

    y_layout = dict(title=y_title, zeroline=False, exponentformat="power",
                    showexponent="all", showgrid=True,
                    showline=True, mirror="all", ticks="inside")
    if y_log:
        y_layout['type'] = 'log'

    layout = go.Layout(
        showlegend=False,
        autosize=True,
        width=600,
        height=500,
        margin=dict(t=40, b=40, l=80, r=40),  # noqa: E741
        hovermode='closest',
        bargap=0,
        xaxis=x_layout,
        yaxis=y_layout,
        title=title
    )

    colorscale = [[0, "rgb(0,0,131)"], [0.125, "rgb(0,60,170)"], [0.375, "rgb(5,255,255)"],
                  [0.625, "rgb(255,255,0)"], [0.875, "rgb(250,0,0)"], [1, "rgb(128,0,0)"]]
    # plot_type = 'surface' if surface else 'heatmap'
    trace = go.Heatmap(z=z, x=x, y=y, autocolorscale=False,  # type=plot_type,
                       hoverinfo="none", colorscale=colorscale)
    fig = go.Figure(data=[trace], layout=layout)
    plot_div = plot(fig, output_type='div', include_plotlyjs=False, show_link=False)

    # The following would remove the hover options, which are not accessible through python
    # https://github.com/plotly/plotly.js/blob/master/src/components/modebar/buttons.js
    # plot_div = plot_div.replace('modeBarButtonsToRemove:[]',
    #                             'modeBarButtonsToRemove:["hoverClosestCartesian",
    #                                                      "hoverCompareCartesian"]')

    if publish:
        try:
            return publish_plot(instrument, run_number, files={'file': plot_div})
        except:  # noqa: E722
            logging.error("Publish plot failed: %s", sys.exc_value)
            return None
    else:
        return plot_div


if __name__ == '__main__':
    import sys
    div = _loadDiv(sys.argv[1])

    # run information is generated from the filename
    name = os.path.split(sys.argv[1])[-1]
    (instr, runnumber) = name.split('_')[:2]

    config = read_configuration('post_processing.conf')
    request = publish_plot(instr, runnumber, {'file': div}, config)
    print('request returned', request.status_code)
