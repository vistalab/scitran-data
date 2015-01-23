# @author:  Gunnar Schaefer

"""
nimsdata.medimg.nimsgephysio
============================

Parse and identify GE MR Physio files.

nimsdata.medimg.nimsbehavior will need to be paired with a writer that is capable of outputting
to text or csv.  Currently, there is no such writer.

"""

import bson
import json
import logging
import tarfile

from .. import data

log = logging.getLogger(__name__)


class GEPhysioError(data.DataError):
    pass


class GEPhysio(data.Reader):

    """Parse and identify GE Physio data."""


    domain = u'mr'
    filetype = u'gephysio'
    state = ['orig']

    def __init__(self, path, load_data=False):
        super(GEPhysio, self).__init__(path, load_data)
        with tarfile.open(path) as archive:
            for ti in archive:
                if not ti.isreg():
                    continue
                try:
                    self._hdr = json.loads(archive.extractfile(ti).read(), object_hook=bson.json_util.object_hook)['header']
                except (KeyError, ValueError):
                    pass
                else:
                    break
            else:
                raise GEPhysioError('no header found')
        self.group = self._hdr['group']
        self.project = self._hdr['project']
        self.session = self._hdr['session']
        self.acquisition = self._hdr['acquisition']
        self.timestamp = self._hdr['timestamp']
        self.metadata_status = None

    def load_data(self):
        log.debug('nimsgephysio.load_data() has not been implemented yet. sorry!')

    @property
    def nims_metadata_status(self):
        return self.metadata_status

    @property
    def nims_group_id(self):
        return self.group

    @property
    def nims_project(self):
        return self.project

    @property
    def nims_session_id(self):
        return self.session

    @property
    def nims_session_label(self):
        return None

    @property
    def nims_session_subject(self):
        return None

    @property
    def nims_acquisition_id(self):
        return self.acquisition

    @property
    def nims_acquisition_label(self):
        return None

    @property
    def nims_acquisition_description(self):
        return None

    @property
    def nims_file_name(self):
        return self.acquisition + '_' + self.filetype

    @property
    def nims_file_ext(self):
        return '.tgz'

    @property
    def nims_file_domain(self):
        return self.domain

    @property
    def nims_file_type(self):
        return self.filetype

    @property
    def nims_file_kinds(self):
        return ['resp', 'ecg']

    @property
    def nims_file_state(self):
        return self.state

    @property
    def nims_timestamp(self): # FIXME: should return UTC time and timezone
        return self.timestamp.replace(tzinfo=bson.tz_util.FixedOffset(-7*60, 'pacific')) #FIXME: use pytz

    @property
    def nims_timezone(self):
        return None
