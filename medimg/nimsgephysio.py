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

from .. import nimsdata

log = logging.getLogger(__name__)


class NIMSGEPhysioError(nimsdata.NIMSDataError):
    pass


class NIMSGEPhysio(nimsdata.NIMSReader):

    """Parse and identify GE Physio data."""


    domain = u'mr'
    filetype = u'gephysio'
    state = ['orig']

    def __init__(self, path, load_data=False):
        super(NIMSGEPhysio, self).__init__(path, load_data)
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
                raise NIMSGEPhysioError('no header found')
        self.group = self._hdr['group']
        self.experiment = self._hdr['experiment']
        self.session = self._hdr['session']
        self.epoch = self._hdr['epoch']
        self.timestamp = self._hdr['timestamp']
        self.nims_metadata_status = None

    def load_data(self):
        log.debug('nimsgephysio.load_data() has not been implemented yet. sorry!')

    @property
    def nims_group_id(self):
        return self.group

    @property
    def nims_experiment(self):
        return self.experiment

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
    def nims_epoch_id(self):
        return self.epoch

    @property
    def nims_epoch_label(self):
        return None

    @property
    def nims_epoch_description(self):
        return None

    @property
    def nims_file_name(self):
        return self.epoch + '_' + self.filetype

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
