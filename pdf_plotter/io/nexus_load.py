from __future__ import (absolute_import, division, print_function)

import os
import collections
import threading
import h5py
import time
import numpy as np

# Traits
from traits.api \
    import HasTraits, Instance, Str, Button

from traitsui.api \
    import View, Item

from traitsui.file_dialog \
    import open_file, FileInfo

# Local
from pdf_plotter.ui.models \
    import Dataset, CorrectedDatasets, Measurement, Experiment

from pdf_plotter.utils.instrument_filters \
    import NomadFilters

# -----------------------------------------------------------#
# Map of instrument name to the filters it can use 
# based on Dataset title

instrument2filters = { 'NOM' : NomadFilters(), }

# -----------------------------------------------------------#
# Measurement-type to workspace-title-startswith Map
#   NOTE: Need to search *Background 1st, hence the order

mtype2title = collections.OrderedDict()
mtype2title['Sample'] = 'sample'
mtype2title['Container Background'] = 'container_background'
mtype2title['Container'] = 'container'
mtype2title['Vanadium Background'] = 'vanadium_background'
mtype2title['Vanadium'] = 'vanadium'
mtype2title['Correction'] = 'correction'

# workspace-title startswith to Measurement-type Map
title2mtype = collections.OrderedDict()
for k, v in mtype2title.iteritems():
    title2mtype[v] = k


# List of measurement types in order of the Controlsview Tree View
mtype_list = ['Sample',
              'Container',
              'Container Background',
              'Vanadium',
              'Vanadium Background',
              'Correction']

# -----------------------------------------------------------#
# Filters dictionary

#NOM_title2filters = { 'Bank: {0:.2f}'.format(theta)

# -----------------------------------------------------------#
# Experiment File Input view

ExperimentFileInputView = View(
    Item('load_button',
         show_label=False,
         ),
)


# -----------------------------------------------------------#
# Thread to handle loading in Experiment files


class NexusFileThread(threading.Thread):
    def __init__(self, f):
        self.f = f
        self.nxresult = None
        threading.Thread.__init__(self)

    
    # Extract input file
    def extract_input_file(self):
        print("No input extraction yet")
        self.instrument = 'NOM'

    # Multithreaded extraction of Datasets from Nexus file
    def extract_datasets_nexus(self):
        nx = self.nxresult
        wksps = [nx[wksp] for wksp in nx
                 if wksp.startswith("mantid_workspace")]

        t_list = list()
        for i, wksp in enumerate(wksps):
            t = DatasetThread(wksp)
            t.corrected_datasets = self.corrected_datasets
            t.instrument = self.instrument
            t_list.append(t)
            t.start()

        for t in t_list:
            t.join()

    # Main thread opens and extracts Nexus and then launchs threads to extract
    # Datasets
    def run(self):
        self.update_status("Loading...")
        self.nxresult = h5py.File(self.f, 'r')
        self.update_status("Done Loading!")
        self.update_status("Extracting input file...")
        self.extract_input_file()
        self.update_status("Extracted input file...")
        self.update_status('Extracting Datasets...')
        self.extract_datasets_nexus()
        self.update_status("Done Extracting!")
        return

# -----------------------------------------------------------#
# Thread to handle extracting Datasets


class DatasetThread(threading.Thread):
    def __init__(self, wksp):
        self.wksp = wksp
        threading.Thread.__init__(self)

    def sort_lists(self, sorter=None, sortee=None):
        if sorter is None or sortee is None:
            return
        sorter_result = [x for x, y in sorted(
            zip(sorter, sortee), key=lambda pair: pair[0])]
        sortee_result = [y for x, y in sorted(
            zip(sorter, sortee), key=lambda pair: pair[0])]
        return sorter_result, sortee_result

    def get_measurement_type(self, title):
        for key in title2mtype:
            if title.startswith(key):
                return title2mtype[key]
        return 'Other'

    def get_filters(self, title):
        instrument_filters = instrument2filters[self.instrument]
        filters = dict()
        for ifilter in instrument_filters.list_of_filters:
            if ifilter.title == title:
                filters['xmin'] = ifilter.xmin
                filters['xmax'] = ifilter.xmax
        return filters
        

    def run(self):
        wksp = dict(self.wksp)

        # Get title and detector group info (L1s, Thetas, and Phis)
        title = str(wksp['title'].value[0])
        groups = wksp['instrument']['detector']['detector_positions'].value

        # Extract detector group info
        L1 = [float(l1) for (l1, theta, phi) in groups]
        Theta = [float(theta) for (l1, theta, phi) in groups]
        Phi = [float(phi) for (l1, theta, phi) in groups]

        # Get detector group data
        x = np.array(wksp['workspace']['axis1'].value)
        err_groups = wksp['workspace']['errors'].value
        y_groups = wksp['workspace']['values'].value   # 0==error, 1==values

        # Re-sort based on Theta degrees
        tmp, err_groups = self.sort_lists(sorter=Theta, sortee=err_groups)
        tmp, y_groups = self.sort_lists(sorter=Theta, sortee=y_groups)
        tmp, L1 = self.sort_lists(sorter=Theta, sortee=L1)
        tmp, Phi = self.sort_lists(sorter=Theta, sortee=Phi)

        # Theta must be sorted last
        tmp, Theta = self.sort_lists(sorter=Theta, sortee=Theta)

        dataset_list = list()
        for i, (l1, theta, phi, y, err) in enumerate(
                zip(L1, Theta, Phi, y_groups, err_groups)):
            info_dict = {
                'L1': l1,
                'Theta': theta,
                'Phi': Phi,
                'yerr': err}
            dataset_title = "Bank: {0:.2f}".format(theta)

            my_x = list()
            if len(x) > len(y):
                diff = len(x) - len(y)
                my_x = x[diff:]

            xmin = min(x)
            xmax = max(x)

            filters = self.get_filters(dataset_title)
            if 'xmin' in filters:
                xmin = filters['xmin']
            if 'xmax' in filters:
                xmax = filters['xmax']

            dataset_list.append(
                Dataset(x=my_x, y=y,
                        title=dataset_title,
                        xmin_filter = xmin,
                        xmax_filter = xmax,
                        info=info_dict,
                        )
            )

        # Get Measurement-type based on title
        measurement_type = self.get_measurement_type(title)
        info_dict = {'measurement_type': measurement_type}
        self.corrected_datasets[title] = CorrectedDatasets(
            datasets=dataset_list,
            title=title,
            info=info_dict
        )
        return


# -----------------------------------------------------------------#
# Thread to handle creating Measurement from CoorectedDatasets
class MeasurementThread(threading.Thread):
    def __init__(self, measurement_type):
        self.measurement_type = measurement_type
        threading.Thread.__init__(self)

    def get_measurement_of_type(self, my_type):
        cd_list = [cd for title, cd in self.corrected_datasets.items()
                   if cd.info['measurement_type'] == my_type]
        self.measurements[my_type] = Measurement(corrected_datasets=cd_list,
                                                 title=my_type)

    def get_other_measurement(self):
        cd_list = [cd for title, cd in self.corrected_datasets.items()
                   if cd.info['measurement_type'] not in mtype_list]

        if len(cd_list) == 0:
            return

        self.measurements['Other'] = Measurement(corrected_datasets=cd_list,
                                                 title='Other')

    def run(self):
        if self.measurement_type == 'Other':
            self.get_other_measurement()
        else:
            self.get_measurement_of_type(self.measurement_type)

# -----------------------------------------------------------------#
# Thread to handle creating the Experiment from CorrectedDatasets


class ExperimentThread(threading.Thread):
    def __init__(self):
        self.experiment = None
        self.measurements = dict()
        threading.Thread.__init__(self)

    # Multithreaded extraction of Datasets from Nexus file
    def create_measurements(self):
        my_list = mtype_list + ['Other']
        t_list = list()
        for i, measurement_type in enumerate(my_list):
            t = MeasurementThread(measurement_type)
            t.corrected_datasets = self.corrected_datasets
            t.measurements = self.measurements
            t_list.append(t)
            t.start()

        for t in t_list:
            t.join()

    def create_experiments(self):
        measurement_list = list()
        for next_title_in_list in mtype_list:
            for title, measurement in self.measurements.iteritems():
                if title == next_title_in_list:
                    measurement_list.extend([measurement])
        name, ext = os.path.splitext(self.filename)
        name = os.path.basename(name)
        self.experiment = Experiment(measurements=measurement_list, title=name)

    # Main thread opens and extracts Nexus and then launchs threads to extract
    # Datasets
    def run(self):
        self.update_status("Creating Measurements...")
        self.create_measurements()
        self.update_status("Done with Measurements!")
        self.update_status("Creating Experiment...")
        self.create_experiments()
        self.update_status("Done loading Experiment!")
        self.update_experiment(self.experiment)
        return


# -----------------------------------------------------------#
# Experiment File Input Model

class ExperimentFileInput(HasTraits):
    # View
    view = ExperimentFileInputView

    # Load button
    load_button = Button("Load Experiment...")

    # NeXus filename that is loaded in
    filename = Str

    # Thread to handle loading in the file
    file_thread = Instance(NexusFileThread)

    # Thread to handle putting together the Experiment from CorrectedDatasets
    experiment_thread = Instance(ExperimentThread)

    # status
    load_status = Str('Load in a file.')

    # Dataset dictionary: key=title, value={CorrectedDatasets, tag}
    corrected_datasets = dict()

    # Returned Experiment
    experiment = Instance(Experiment)

    # Get update from thread
    def update_status(self, status):
        self.load_status = status

    # Get experiment update from thread
    def update_experiment(self, experiment):
        self.experiment = experiment

    # Handle the user clicking the Load button
    def _load_button_changed(self):
        f = open_file(file_name=os.getcwd(),
                      extensions=[FileInfo()],
                      filter=['*.nxs', '*.dat'])
        self.filename = f
        if f != '':
            name, ext = os.path.splitext(f)
            if ext == '.nxs':
                self.load_and_extract_nexus(f)
                while self.file_thread.isAlive():
                    time.sleep(1)
                self.form_experiment()
            elif ext == '.dat':
                self.load_and_extract_dat_file(f)

    # Parse the Experiment NeXus file
    def load_and_extract_nexus(self, f):
        self.file_thread = NexusFileThread(f)
        self.file_thread.update_status = self.update_status
        self.file_thread.corrected_datasets = self.corrected_datasets
        self.file_thread.start()

    def form_experiment(self):
        self.experiment_thread = ExperimentThread()
        self.experiment_thread.update_status = self.update_status
        self.experiment_thread.update_experiment = self.update_experiment
        self.experiment_thread.corrected_datasets = self.corrected_datasets
        self.experiment_thread.filename = self.filename
        self.experiment_thread.experiment = self.experiment
        self.experiment_thread.start()
