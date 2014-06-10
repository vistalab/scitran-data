# @author:  Gunnar Schaefer

"""
nimsdata.nimsgephysio
=====================

Parse and identify GE MR Physio files.

"""

import bson
import json
import logging
import tarfile

import nimsdata

log = logging.getLogger(__name__)


class NIMSGEPhysioError(nimsdata.NIMSDataError):
    pass


class NIMSGEPhysio(nimsdata.NIMSReader):

    filetype = 'gephysio'

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
    def nims_epoch_id(self):
        return self.epoch

    @property
    def nims_type(self):
        return ('original', 'physio', self.filetype)

    @property
    def nims_filename(self):
        return self.epoch + '_' + self.filetype

    @property
    def nims_file_ext(self):
        return '.tgz'

    @property
    def nims_timestamp(self): # FIXME: should return UTC time and timezone
        return self.timestamp.replace(tzinfo=bson.tz_util.FixedOffset(-7*60, 'pacific')) #FIXME: use pytz

    @property
    def nims_timezone(self):
        return None

    @property
    def nims_session_label(self):
        return None

    @property
    def nims_session_subject(self):
        return None

    @property
    def nims_session_type(self):
        return None

    @property
    def nims_epoch_label(self):
        return None

    @property
    def nims_epoch_description(self):
        return None

    @property
    def nims_epoch_type(self):
        return None
