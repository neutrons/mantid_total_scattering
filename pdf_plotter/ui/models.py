#!/usr/bin/env python

from traits.api \
    import HasTraits, Float, Array, Str, List, Dict

# -----------------------------------------------------------#
# Models


class Dataset(HasTraits):
    x = Array
    y = Array
    title = Str
    info = Dict

    xmin_filter = Float
    xmax_filter = Float


class CorrectedDatasets(HasTraits):
    datasets = List(Dataset)
    title = Str
    info = Dict


class Measurement(HasTraits):
    corrected_datasets = List(CorrectedDatasets)
    title = Str


class Experiment(HasTraits):
    measurements = List(Measurement)
    title = Str
    info = Dict
