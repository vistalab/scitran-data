#!/usr/bin/env python
# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S Hahn

import os
import abc
import logging

log = logging.getLogger('nimsdata')


class NIMSDataError(Exception):
    pass


class NIMSData(object):

    __metaclass__ = abc.ABCMeta

    parse_priority = 0

    _session_properties = {
            'name': {
                'attribute': 'nims_session_name',
                'title': 'Name',
                'type': 'string',
                'maxLength': 32,
            },
            'subject': {
                'attribute': 'nims_session_subject',
                'title': 'Subject',
                'type': 'string',
                'maxLength': 16,
            },
            'datatype': {
                'attribute': 'nims_session_type',
                'title': 'Type',
                'type': 'array',
            },
            'timestamp': {
                'title': 'Timestamp',
                'type': 'string',
                'format': 'date-time',
            },
            'timezone': {
                'title': 'Timezone',
                'type': 'string',
            },
    }
    session_properties = _session_properties

    _epoch_properties = {
            'name': {
                'attribute': 'nims_epoch_name',
                'title': 'Name',
                'type': 'string',
                'maxLength': 32,
            },
            'description': {
                'attribute': 'nims_epoch_description',
                'title': 'Description',
                'type': 'string',
                'maxLength': 64,
            },
            'datatype': {
                'attribute': 'nims_epoch_type',
                'title': 'Type',
                'type': 'string',
                'maxLength': 32,
            },
            'timestamp': {
                'attribute': 'nims_timestamp',
                'title': 'Timestamp',
                'type': 'string',
                'format': 'date-time',
            },
            'timezone': {
                'attribute': 'nims_timezone',
                'title': 'Timezone',
                'type': 'string',
            },
    }
    epoch_properties = _epoch_properties

    @abc.abstractmethod
    def __init__(self):
        self.data = None
        self.metadata = type('NIMSMetadata', (object, ), {})()
        self.status = 'metadata pending'                        # dict?
        self.status = {'metadata': 'pending'}
        pass

    @abc.abstractproperty
    def nims_group(self):
        pass

    @abc.abstractproperty
    def nims_experiment(self):
        pass

    @abc.abstractproperty
    def nims_session(self):
        pass

    @abc.abstractproperty
    def nims_epoch(self):
        pass

    @abc.abstractproperty
    def nims_type(self):
        pass

    @abc.abstractproperty
    def nims_filename(self):
        pass

    @abc.abstractproperty
    def nims_timestamp(self):
        pass

    @abc.abstractproperty
    def nims_timezone(self):
        pass

    @property
    def nims_session_name(self):
        pass

    @property
    def nims_session_subject(self):
        pass

    @property
    def nims_session_type(self):
        pass

    @property
    def nims_epoch_name(self):
        pass

    @property
    def nims_epoch_description(self):
        pass

    @property
    def nims_epoch_type(self):
        pass

    @abc.abstractmethod
    def load_data(self, archive):
        pass

    @abc.abstractmethod
    def convert(self, outbase):
        pass

    def getnestedattr(self, attr):
        for name in attr.split('.'):
            self = getattr(self, name)
        return self


if __name__ == '__main__':
    import sys
    import argparse

    sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../..')))

    import nimsdata
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='directory of dicoms to convert')
    parser.add_argument('outbase', nargs='?', help='basename for output files (default: input)')
    parser.add_argument('-p', '--parser', help='parser must be one of %s' % str(nimsdata.PARSERS.keys()), choices=nimsdata.PARSERS.keys())
    parser.add_argument('-i', '--ignore_json', help='do not use json metadata in output metadata', dest='ignore_json', action='store_true', default=False)
    parser.add_argument('-v', '--verbose', help='enable verbose logging', dest='verbose', action='store_true', default=False)
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    if not os.path.exists(args.input):
        log.error('could not find input file %s' % args.input)
    else:
        # TODO: add output-exists-check, and add option to overwrite

        # TODO: improve name mangling
        outbase = args.outbase or os.path.basename(os.path.splitext(args.input.rstrip('/'))[0])
        nimsdata.convert(nimsdata.parse(args.input, args.ignore_json, args.parser), outbase)
