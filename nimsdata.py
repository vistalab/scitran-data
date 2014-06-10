#!/usr/bin/env python
# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S. Hahn

"""
nimsdata.nimsdata
=================

provides base classes for NIMS related exceptions and an abstract implementation of
a NIMSData object. nimsdata.nimsdata module provides base classes for Exceptions,
Readers and Writers.  Also provides the nimsdata command line interface.

The NIMSData object is meant to be subclassed to create abstract classes for data domains
(e.g. MRData) and concrete classes for domain specific data type (e.g. Dicom).

In some cases, such as with MR data, there are functions and definitions that
are used by multiple data types, such as MR Dicoms, MR Nifti, and MR PFiles,
all may want to manipulate metadata in the same way.  In these cases, it may be
more prudent to add an abstract base class for your data domain, and several
seperate data type specific classes.

Adding a new reader or writer by subclassing NIMSReader.  A reader must implement
an `__init__` that parses enough information to sort a file, and a `load_data` method
that converts the data to an intermediate format compatible with a writer.
Adding a new domain specific type: **write this out**

A writer must implement a `write` method that accepts an intermediate format, and saves it
to a specific filetype.

Readers and Writers of the same group/datadomain should share an intermediate format that
is produced by a reader, and can be used by a Writer.

Some general rules about Reader and writer subclasses:
    - readers input should be a tgz of files of interest, and a json file
    - readers should read files, collect metadata, and define nims properties
    - readers should return the data as close to 'native' as possible
    - `__init__` should be able to quickly collect info for sorting
    - `load_data` should handle the rest

"""

import os
import abc
import copy
import logging
import traceback

log = logging.getLogger(__name__)


def dict_merge(a, b):
    """
    Recursively merge two dictionaries.

    If a key exists in both dictionary_A and dictionary_B, then the resulting
    dictionary will contain the value from dictionary_B.

    Parameters
    ----------
    a : dict
    b : dict

    Returns
    -------
    result : dict

    """
    result = copy.deepcopy(a)
    for k, v in b.iteritems():
        if k in result and isinstance(result[k], dict):
                result[k] = dict_merge(result[k], v)
        else:
            result[k] = copy.deepcopy(v)
    return result


class NIMSDataError(Exception):

    """
    Base class for NIMS exceptions.

    Parameters
    ----------
    message : str
        Message to report along with the exception.
    log_level : logging level, i.e. logging.ERROR
        What logging level to report the error with,
        must be a valid logging level class, such as
        `logging.DEBUG` or `logging.ERROR`.

    """

    def __init__(self, message, log_level=None):
        """instantiate nimsdataerror exception class."""
        super(NIMSDataError, self).__init__(message)
        if log_level is not None:
            message = '%s\n%s' % (message, traceback.format_exc())
            log.log(log_level, message)


class NIMSReader(object):

    """
    Abstract base class that provides Standard interfaces, properties and methods that are required for data readers.

    Cannot be instantiated.

    Parameters
    ----------
    path : str
        filepath as string
    load_data : boolean
        Indicate if parse should return with all data loaded.  Default is False.

    """

    __metaclass__ = abc.ABCMeta

    _experiment_properties = {
            'gid': {
                'field': 'nims_group_id',
                'type': 'string',
            },
            'group': {
                'title': 'Group',
                'type': 'string',
                # FIXME: dynamically add enum of all possible groups
            },
            'name': {
                'field': 'nims_experiment',
                'title': 'Name',
                'type': 'string',
                'maxLength': 32,
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
    experiment_properties = _experiment_properties

    _session_properties = {
            'uid': {
                'field': 'nims_session_id',
                'type': 'string',
            },
            'label': {
                'field': 'nims_session_label',
                'title': 'Label',
                'type': 'string',
                'maxLength': 32,
            },
            'subject': {
                'title': 'Subject',
                'type': 'object',
                'properties': {
                    'code': {
                        'field': 'nims_session_subject',
                        'title': 'Code',
                        'type': 'string',
                        'maxLength': 16,
                    },
                },
            },
            'datatype': {
                'field': 'nims_session_type',
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
            'uid': {
                'field': 'nims_epoch_id',
                'type': 'string',
            },
            'label': {
                'field': 'nims_epoch_label',
                'title': 'Label',
                'type': 'string',
                'maxLength': 32,
            },
            'description': {
                'field': 'nims_epoch_description',
                'title': 'Description',
                'type': 'string',
                'maxLength': 64,
            },
            'datatype': {
                'field': 'nims_epoch_type',
                'title': 'Datatype',
                'type': 'string',
                'maxLength': 32,
            },
            'timestamp': {
                'field': 'nims_timestamp',
                'title': 'Timestamp',
                'type': 'string',
                'format': 'date-time',
            },
            'timezone': {
                'field': 'nims_timezone',
                'title': 'Timezone',
                'type': 'string',
            },
    }
    epoch_properties = _epoch_properties

    def _schema_init(self, schema):
        for k, v in schema.iteritems():
            if isinstance(v, dict):
                    self._schema_init(v)
            elif k == 'field' and not v.startswith('nims_'):
                setattr(self, v, None)

    @abc.abstractmethod
    def __init__(self, path, load_data=False):
        self._schema_init(self.experiment_properties)
        self._schema_init(self.session_properties)
        self._schema_init(self.epoch_properties)
        self.filepath = os.path.abspath(path)
        self.data = None
        self.metadata_status = 'empty'

    @abc.abstractmethod
    def load_data(self):
        pass

    @abc.abstractproperty
    def nims_group_id(self):
        pass

    @abc.abstractproperty
    def nims_experiment(self):
        pass

    @abc.abstractproperty
    def nims_session_id(self):
        pass

    @abc.abstractproperty
    def nims_epoch_id(self):
        pass

    @abc.abstractproperty
    def nims_filename(self):
        pass

    @abc.abstractproperty
    def nims_file_ext(self):
        pass

    @abc.abstractproperty
    def nims_timestamp(self):
        pass

    @abc.abstractproperty
    def nims_timezone(self):
        pass

    @abc.abstractproperty
    def nims_session_label(self):
        pass

    @abc.abstractproperty
    def nims_session_subject(self):
        pass

    @abc.abstractproperty
    def nims_session_type(self):
        pass

    @abc.abstractproperty
    def nims_epoch_label(self):
        pass

    @abc.abstractproperty
    def nims_epoch_description(self):
        pass

    @abc.abstractproperty
    def nims_epoch_type(self):
        pass

    def __str__(self):
        # TODO: tighten this up.
        properties = []
        for prop, value in vars(self).iteritems():
            if prop in ['data', 'imagedata']:
                continue
            if prop.startswith('_'):
                continue
            if value is None:
                value = ''
            properties.append((prop, value))

        properties.sort()

        for prop in dir(self):
            value = None
            if prop in ['data', 'imagedata']:
                continue
            if not prop.startswith('nims'):
                continue
            value = getattr(self, prop)
            if value == None:
                value = ''
            properties.append((prop, str(value)))

        return '\n'.join(['%-30s: %s' % (p, v) for p, v in properties])


class abstractclassmethod(classmethod):

    """
    Decorator for abstract classmethods.

    Similar to abstractmethod.

    .. code-block:: python

        class C(metaclass=ABCMeta):
            @abstractclassmethod
            def my_abstract_classmethod(cls, ...):
                ...

    'abstractclassmethod' is deprecated. Use 'classmethod' with
    'abstractmethod' instead.

    """

    __isabstractmethod__ = True

    def __init__(self, callable):
        callable.__isabstractmethod__ = True
        super(abstractclassmethod, self).__init__(callable)


class NIMSWriter(object):

    """
    Abstract base class that provides standard interfaces, properties and methods for all writers.

    Cannot be instantiated.

    """

    __metaclass__ = abc.ABCMeta

    @abstractclassmethod
    def write(cls, metadata, imagedata, outbase, **kwargs):
        """
        Do stuff and then write the file.

        abstract implementaiton doesn't REALLY need to know about filepath.  Subclasses
        may way to build upon write, such as performing voxel reordering of MR data.

        Parameters
        ----------
        metadata : object
            fully loaded instance of a NIMSReader
        imagedata : dict
            dictionary of np.arrays.  primary dataset exists in imagedata['']
        outbase : str
            output name prefix
        """
        if metadata is None or imagedata is None:
            raise NIMSDataError('no metadata/data?')


if __name__ == '__main__':
    import sys
    import argparse
    import datetime             # mainly for debug messages
    log = logging.getLogger('nimsdata')

    sys.path.insert(0, os.path.abspath(os.path.join(__file__, '../..')))
    import nimsdata

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='file to convert')
    parser.add_argument('outbase', nargs='?', help='basename for output files (default: input)')
    parser.add_argument('-p', '--parser', help='parser to use', choices=nimsdata.READERS.keys())
    parser.add_argument('-i', '--ignore_json', help='do not use json metadata in output metadata', action='store_true', default=False)
    parser.add_argument('-w', '--writer', help='write to use', choices=nimsdata.WRITERS.keys())
    parser.add_argument('-v', '--verbose', help='enable verbose logging', dest='verbose', action='store_true', default=False)
    parser.add_argument(      '--parser_kwarg', action='append', help='keyword arguments to pass directly to the parser')
    parser.add_argument(      '--writer_kwarg', action='append', help='keyword arguments to pass directly to the writer')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if not os.path.exists(args.input):
        raise NIMSDataError('could not find input file %s' % args.input)

    outbase = args.outbase or os.path.basename(os.path.splitext(args.input.rstrip('/'))[0])

    p_kwargs = {}
    if args.parser_kwarg:
        for item in args.parser_kwarg:
            kw, val = item.split('=')
            p_kwargs[kw] = val
    log.debug(p_kwargs)

    w_kwargs = {}
    if args.writer_kwarg:
        for item in args.writer_kwarg:
            kw, val = item.split('=')
            w_kwargs[kw] = val
    log.debug(w_kwargs)

    log.debug('%20s - %s' % ('start parse', str(datetime.datetime.now())))
    ds = nimsdata.parse(args.input, load_data=True, ignore_json=args.ignore_json, filetype=args.parser, **p_kwargs)
    log.debug('%20s - %s' % ('done parse', str(datetime.datetime.now())))

    if not ds:
        raise NIMSDataError('%s could not be parsed' % args.input)
    if ds.data is None:
        raise NIMSDataError('%s has no data' % args.input)

    log.debug('%20s - %s' % ('start write', str(datetime.datetime.now())))
    nimsdata.write(ds, ds.data, outbase, filetype=args.writer, **w_kwargs)
    log.debug('%20s - %s' % ('done write', str(datetime.datetime.now())))
