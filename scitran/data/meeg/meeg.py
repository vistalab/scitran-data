"""
scitran.data.meeg
=================

Data format for M/EEG data using mne-python.

"""

import logging
import tempfile
import zipfile
import warnings
import os
from os import path as op
from datetime import datetime, date
import shutil

from mne.io import read_raw_fif

from .. import data

log = logging.getLogger(__name__)  # root logger already configured

# see data.py for expected project properties
project_properties = data.project_properties

# add additional session properties, which should be added as attributes of
# the Reader object
_session_properties = {
    'subject': {
        'type': 'object',
        'properties': {
            'firstname': {
                'field': 'subj_firstname',
                'title': 'First Name',
                'type': 'string',
            },
            'lastname': {
                'field': 'subj_lastname',
                'title': 'Last Name',
                'type': 'string',
            },
            'dob': {
                'field': 'subj_dob',
                'title': 'Date of Birth',
                'type': 'string',
                'format': 'date',  # i.e., datetime object
            },
            'sex': {
                'field': 'subj_sex',
                'title': 'Sex',
                'type': 'string',
                'enum': ['male', 'female'],
            },
            'hand': {
                'field': 'subj_hand',
                'title': 'Handedness',
                'type': 'string',
                'enum': ['right', 'left'],
            },
        },
    },
}
session_properties = data.dict_merge(data.session_properties,
                                     _session_properties)

_acquisition_properties = {  # add custom acquisition properties
}
acquisition_properties = data.dict_merge(data.acquisition_properties,
                                         _acquisition_properties)


class MEEGError(data.DataError):
    pass


class MEEGReader(data.Reader):

    """
    Parameters
    ----------
    path : str
        Path to input file.
    load_data : boolean
        Indicate if a reader should attempt to immediately load all data.
        Default False.
    """
    project_properties = project_properties
    session_properties = session_properties
    acquisition_properties = acquisition_properties

    domain = u'meeg'
    filetype = u'meeg'
    state = ['orig']

    def __init__(self, path, load_data=False, timezone=None):
        super(MEEGReader, self).__init__(path, load_data, timezone)

        #
        # Process the incoming data
        #

        self._temp_dir = tempfile.mkdtemp()
        os.mkdir(op.join(self._temp_dir, 'reap'))
        try:
            with zipfile.ZipFile(self.filepath, 'r') as zip_file:
                zip_fnames = [op.join('reap', op.basename(fname))
                              for fname in zip_file.namelist()]
                fnames = [zip_file.extract(fname, self._temp_dir)
                          for fname in zip_fnames if fname.endswith('.fif')]
        except Exception as e:
            raise MEEGError(e)

        # load information and optionally data from the files
        with warnings.catch_warnings(record=True):
            self._raws = [read_raw_fif(fname, allow_maxshield=True,
                                       preload=load_data)
                          for fname in fnames]
        info = self._raws[0].info
        subject_info = info['subject_info']
        hand_dict = {1: 'right', 2: 'left'}
        sex_dict = {1: 'male', 2: 'female'}

        #
        # Parameters required by NIMS
        #

        # pick a unique filename
        meas_date = datetime.fromtimestamp(info['meas_date'][0])
        fname = meas_date.strftime('%Y_%m_%d_%H_%M_%S')
        self.filename = fname
        self.group_id = info['experimenter']  # XXX always "neuromag", !useful
        self.project_name = info['proj_name']
        self.session_id = meas_date.strftime('%Y%m%d')
        self.acquisition = info['description']
        self.subject = subject_info['his_id']

        #
        # Additional session properties
        #

        self.subj_firstname = subject_info['first_name']
        self.subj_lastname = subject_info['last_name']
        self.subj_dob = \
            datetime.combine(date(*subject_info['birthday']),
                             datetime.min.time())
        self.subj_hand = hand_dict[subject_info['hand']]
        self.subj_sex = sex_dict[subject_info['sex']]

    def __del__(self):
        shutil.rmtree(self._temp_dir)

    def load_data(self):
        super(MEEGReader, self).load_data()
        for raw in self._raws:
            raw.preload_data()

    @property
    def nims_group_id(self):
        return self.group_id

    @property
    def nims_project(self):
        return self.project_name

    @property
    def nims_session_id(self):
        return self.session_id

    @property
    def nims_session_label(self):
        return self.session_id

    @property
    def nims_session_subject(self):
        return self.session_subject

    @property
    def nims_acquisition_id(self):
        return self.acquisition

    @property
    def nims_acquisition_label(self):
        return self.acquisition

    @property
    def nims_acquisition_description(self):
        return self.acquisition

    @property
    def nims_file_name(self):
        return self.filename

    @property
    def nims_file_kinds(self):
        return ['FIF']

    # the following are all handled by the super class Reader

    @property
    def nims_metadata_status(self):
        return super(MEEGReader, self).nims_metadata_status

    @property
    def nims_file_ext(self):
        return super(MEEGReader, self).nims_file_ext

    @property
    def nims_file_domain(self):
        return super(MEEGReader, self).nims_file_domain

    @property
    def nims_file_type(self):
        return super(MEEGReader, self).nims_file_type

    @property
    def nims_file_state(self):
        return super(MEEGReader, self).nims_file_state

    @property
    def nims_timestamp(self):
        return super(MEEGReader, self).nims_timestamp

    @property
    def nims_timezone(self):
        return super(MEEGReader, self).nims_timezone
