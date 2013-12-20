# @author:  Gunnar Schaefer
#           Bob Dougherty

import abc
import json
import datetime
import bson.json_util


class NIMSDataError(Exception):
    pass


class NIMSData(object):

    __metaclass__ = abc.ABCMeta

    parse_priority = 0

    _session_properties = {
            'exam': {
                'attribute': 'exam_no',
                'title': 'Exam Number',
                'type': 'integer',
            },
            'patient_id': {
                'attribute': 'patient_id',
                'title': 'Patient ID',
                'type': 'string',
            },
            'firstname': {
                'attribute': 'subj_firstname',
                'title': 'First Name',
                'type': 'string',
            },
            'lastname': {
                'attribute': 'subj_lastname',
                'title': 'Last Name',
                'type': 'string',
            },
            'dob': {
                'attribute': 'subj_dob',
                'title': 'Date of Birth',
                'type': 'string',
                'format': 'date-time',
            },
            'sex': {
                'attribute': 'subj_sex',
                'title': 'Sex',
                'type': 'string',
                'enum': ['male', 'female'],
            },
    }
    session_properties = _session_properties

    _epoch_properties = {
            'timestamp': {
                'attribute': 'timestamp',
                'title': 'Timestamp',
                'format': 'date-time',
            },
            'series': {
                'attribute': 'series_no',
                'title': 'Series',
                'type': 'integer',
            },
            'acquisition': {
                'attribute': 'acq_no',
                'title': 'Acquisition',
                'type': 'integer',
            },
            'description': {
                'attribute': 'series_desc',
                'title': 'Description',
                'type': 'string',
                'maxLength': 64,
            },
    }
    epoch_properties = _epoch_properties

    _file_properties = {
            'datakind': {
                'attribute': 'datakind',
                'type': 'string',
            },
            'datatype': {
                'attribute': 'datatype',
                'type': 'string',
            },
            'filetype': {
                'attribute': 'filetype',
                'type': 'string',
            },
    }
    file_properties = _file_properties

    @staticmethod
    def parse(filepath):
        for sc in subclasses:
            try:
                dataset = sc(filepath)
            except NIMSDataError:
                dataset = None
            else:
                break
        else:
            raise NIMSDataError('%s could not be parsed' % filepath)
        return dataset

    @abc.abstractmethod
    def __init__(self):
        self.default_subj_code = None

    @property
    def session_uid(self):
        return self.exam_uid.replace('.', '_')

    @property
    def epoch_uid(self):
        return self.series_uid.replace('.', '_') + '_' + str(self.acq_no)

    @property
    def canonical_filename(self):
        return '%s_%s_%s_%s' % (self.exam_uid.replace('.', '_'), self.series_no, self.acq_no, self.filetype)

    def get_json_metadata(self, tgt_cls=None):
        tgt_cls = tgt_cls or self.__class__
        field_names  = ['exam_uid', 'series_uid']
        field_names += [attrs['attribute'] for attrs in tgt_cls.session_properties.itervalues()]
        field_names += [attrs['attribute'] for attrs in tgt_cls.epoch_properties.itervalues()]
        return json.dumps(dict(((fn, getattr(self, fn, None)) for fn in field_names)), default=bson.json_util.default)

    def set_json_metadata(self, json_metadata):
        metadata = json.loads(json_metadata, object_hook=bson.json_util.object_hook)
        for attribute, value in metadata.iteritems():
            if isinstance(value, datetime.datetime):
                value = value.replace(tzinfo=None)
            setattr(self, attribute, value)
