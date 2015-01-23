#!/usr/bin/env python
# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S. Hahn

"""
nimsdata.nimsdata
=================

provides base classes for NIMS related exceptions and an abstract implementation of a NIMSData
object. nimsdata.nimsdata module provides base classes for Exceptions, Readers and Writers.
Also provides the nimsdata command line interface. The NIMSData object is meant to be subclassed
to create abstract classes for data domains (e.g. MRData) and concrete classes for domain specific
data type (e.g. Dicom).

In some cases, such as with MR data, there are functions and definitions that are used by multiple
data types, such as MR Dicoms, MR Nifti, and MR PFiles, all may want to manipulate metadata in the
same way.  In these cases, it may be more prudent to add an abstract base class for your data
domain, and several seperate data type specific classes.

see :doc:`extending_nimsdata` for more information on creating a new reader
subclass.

"""

import os
import abc
import copy
import json
import logging
import tarfile
import warnings
import datetime
import traceback
import bson.json_util


log = logging.getLogger(__name__)
warnings.simplefilter('ignore', FutureWarning)

# note: readers/writers.json will never contain datetime, or objects that require bson.json_util
READERS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'readers.json')))
WRITERS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'writers.json')))
MODULES = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'modules.json')))


project_properties = {
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
            'field': 'nims_project',
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

session_properties = {
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
        'type': {
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

acquisition_properties = {
        'uid': {
            'field': 'nims_acquisition_id',
            'type': 'string',
        },
        'label': {
            'field': 'nims_acquisition_label',
            'title': 'Label',
            'type': 'string',
            'maxLength': 32,
        },
        'description': {
            'field': 'nims_acquisition_description',
            'title': 'Description',
            'type': 'string',
            'maxLength': 64,
        },
        'datatype': {
            'title': 'Datatype',
            'type': 'array',
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
        'filetype': {
            'field': 'nims_file_type',
            'title': 'Filetype',
            'type': 'string',
        },
        'state': {
            'field': 'nims_file_state',
            'title': 'State',
            'type': 'string',
        },
}


def module_by_type(domain_kind):
    """
    Retrieve a module based on tuple(domain, kind).

    Parameters
    ----------
    domain_kind : 2-item tuple
        Tuple of data domain and kind. kind can be None.

    Returns
    -------
    mod : module
        module that defines the properties for (domain, kind).

    Raises
    ------
    NIMSDataError
        Unable to import the module specified by (domain, kind) in modules.json.

    """
    domain, kind = domain_kind[0], domain_kind[1]
    mod_name = MODULES.get('%s.%s' % (domain, kind)) if (domain and kind) else MODULES.get(domain)
    try:
        mod = __import__(mod_name, globals())
    except (ImportError, TypeError):
        raise DataError('no module matches that domain/datatype %s' % mod_name, log_level=logging.ERROR)
    return mod


def dict_merge(a, b=None, in_place=False):
    """
    Recursively merge dictionaries.

    If a and b are both dictionaries, merge them. If a is a list of dicionaries, b must be None,
    then merge all the dictionaries.

    Parameters
    ----------
    a : dict or list
        dictionary, or list of dicionaries to merge.
    b : dict or None
        dictionary or None.  if a is a list, b must be None.
    in_place : bool [default False]
        False creates a deep copy of dict a. dict a is not modified.
        True merges dict a with dict b, altering dict a. dict a is modified.

    Returns
    -------
    merged : dict
        dictionary of dict b merged into dict a.

    Raises
    ------
    NIMSDataError
        inputs were of unexpected types.  Inputs must be two dictionaries, or one list.

    """
    if isinstance(a, dict) and isinstance(b, dict):
        result = copy.deepcopy(a) if not in_place else a
        for k, v in b.iteritems():
            if k in result and isinstance(result[k], dict):
                result[k] = dict_merge(result[k], v, in_place=True)
            else:
                result[k] = copy.deepcopy(v)
        return result
    elif isinstance(a, list) and b is None:
        result = copy.deepcopy(a[0]) if not in_place else a[0]
        for b in a[1:]:
            result = dict_merge(result, b, in_place=True)
    else:
        raise DataError('inputs must be two dictionaries, or a list of dicionaries',
                            log_level=logging.ERROR)
    return result


def acquisition_properties_by_type_list(type_list):
    """
    Assemble acquisition properties from list of types.

    Parameters
    ----------
    type_list : list
        list of two-item (domain, kind) tuples.

    Returns
    -------
    merged_dict : dict
        merged acquisition properties for items in type_list.

    """
    return dict_merge([module_by_type(t).acquisition_properties for t in type_list])


def session_properties_by_type_list(type_list):
    """
    Assemble session properties from list of types.

    Parameters
    ----------
    type_list : list
        list of two item (domain, kind) tuples.

    Returns
    -------
    merged_dict : dict
        merged session properties for items in type_list.

    """
    return dict_merge([module_by_type(t).session_properties for t in type_list])


def project_properties_by_type_list(type_list):
    """
    Assemble project properties from list of types.

    Parameters
    ----------
    type_list : list
        list of two-item (domain, kind) tuples.

    Returns
    -------
    merged_dict : dict
        merged project properties for items in type_list.

    """
    return dict_merge([module_by_type(t).project_properties for t in type_list])


# FIXME; make this more flexible to deal with various levels of nesting
def get_handler(name, handlerdict):
    """
    Retrieve the class object of reader or writer from its string label.

    Parameters
    ----------
    name : str
        reader or writer name
    handlerdict : dict
        name of which handler dict to search from, READERS or WRITERS.

    Returns
    -------
    handler : NIMSReader or NIMSWriter subclass
        Reader or Writer object.

    Raises
    ------
    NIMSDataError
        The specified handler does not exist, or cannot be imported

    """
    try:
        # is there a more flexible way to do this? get nested attr?
        container, mod, klass = map(str, handlerdict.get(name).split('.'))
        handler = getattr(getattr(__import__(container, globals(), fromlist=[mod]), mod), klass)
    except (ImportError, AttributeError):
        raise DataError('no handler', log_level=logging.ERROR)  # XXX FAIL! unexpected no handler
    return handler


def _parse_dataset(path, filetype, load_data, debug, **kwargs):
    """
    Implementation detail to reduce repeated code, called by parse().

    Main parse interface has two paths, one that inspects a json file
    for filetype, and one that does not.

    Parameters
    ----------
    path : str
        path to input file
    filetype : str
        filetype
    load_data : bool
        whether or not to load all data
    debug : bool
        whether or not to mask exceptions
    **kwargs : dict
        passes **kwargs to filetype specific parsers

    Returns
    -------
    ds : NIMSReader subclass
        readers NIMSReader subclass of that data that could be parsed.
    """
    try:
        parser = get_handler(filetype, READERS)
        ds = parser(path, load_data, **kwargs)
    except Exception as e:
        if debug:
            raise e
        else:
            raise DataError(e)
    else:
        ds.compressed = True  # v1 compat: inform scheduler to not compress
        log.debug('parse end: %s' % str(datetime.datetime.now()))

    if hasattr(ds, 'failure_reason') and ds.failure_reason:
        log.warning('error during parsing: %s' % str(ds.failure_reason))
    return ds


def parse(path, filetype=None, load_data=False, ignore_json=False, debug=False, **kwargs):
    """
    Parse the file at path with a filetype-specific parser.

    The parse function expects a tgz that contains input file(s) and a metadata.json.  The json
    should be the first item in the tar achive, and the json should declare the filetype, and any
    metadata that should be overwritten.

    This function will pass the input file to the appropriate handler,
    if available.

    =========== ======== ===================================================
    ignore_json filetype behavior
    =========== ======== ===================================================
    False       None     **DEFAULT** try to read json, use parser from json.
    False       'dicom'  try to read json, use 'dicom' parser.
    True        None     **INVALID**
    True        'dicom'  don't try to read json, use 'dicom' parser.
    =========== ======== ===================================================

    Parameters
    ----------
    path : str
        path to input file.
    filetype : str
        string name of parser to use.  no special search or inferences will be made if this option
        is specified.
    load_data : bool  [default False]
        attempt to load all metadata and data at parse.  like running parse, then load_data
    ignore_json : bool [default False]
        True, don't look for json file. therefore the parser cannot be read from the json, and a
        parser MUST be specfied if ignore_json is true. False, do look for json file.
    debug : bool [default False]
        developer option, False masks all exceptions as NIMSDataError.  debug=True does not
        mask exceptions.
    kwargs : dict
        keyword arguments passed to reader.

    Returns
    -------
    parser_class : obj
        NIMSReader populated with data and metadata attributes. this NIMSReader object can be
        passed to any compatible writer to complete the conversion process.

    Raises
    ------
    NIMSDataError
        input file is of unacceptable type, or problem with combination of input parameters.

    Examples
    --------
    convert dicom.tgz to dicoms.nii.gz

    .. code-block:: python

        import nimsdata
        ds = nimsdata.parse('dicoms.tgz', load_data=False, filetype='dicom')
        ds.load_data()
        nimsdata.write(ds, ds.data, outbase, filetype='nifti')

    TODO: convert all in-line code EXAMPLES to using doc-string code.

    """
    log.debug('parse start: %s' % str(datetime.datetime.now()))
    if not os.path.exists(path):
        raise DataError('input path %s not found' % path, log_level=logging.ERROR)

    if ignore_json:
        log.debug('ignore_json = True')
        if not filetype:
            raise DataError('filetype must be specified if ignore_json=True', log_level=logging.ERROR)
        else:
            ds = _parse_dataset(path, filetype, load_data, debug, **kwargs)
            return ds

    json_data = {}
    if os.path.isfile(path) and tarfile.is_tarfile(path):
        log.debug('inspecting tarfile %s for json' % path)
        with tarfile.open(path) as archive:
            # TODO: explicit name check for Metadata.json or metadata.json
            # if ti.name.lower is metadata.json
            for ti in archive:
                try:
                    json_data = json.loads(archive.extractfile(ti).read(), object_hook=bson.json_util.object_hook)
                except (TypeError, KeyError, AttributeError, ValueError):
                    pass
                else:
                    log.debug('json found, %s' % ti.name)
                    if not filetype:
                        filetype = json_data.get('filetype')
                        log.debug('filetype from json: %s' % filetype)
                    break
            else:
                raise DataError('expected json file that indicates filetype in tgz')
    elif os.path.isdir(path):  # TODO: implement 'break-at-first-readable' json read for directory
        raise DataError('directory input not implemented', log_level=logging.ERROR)
    elif os.path.isfile(path):  # TODO: infer type
        if path.endswith('.7.gz') or path.endswith('.7'):   # single P12345.7 or P12345.7.gz
            filetype = 'pfile'
        else:
            raise DataError('non tar-files not implemented', log_level=logging.ERROR)

    ds = _parse_dataset(path, filetype, load_data, debug, **kwargs)
    for key, value in json_data.get('overwrite', {}).iteritems():  # FIXME: handle NESTED information
        setattr(ds, key, value)
    return ds


def write(metadata, imagedata, outbase, filetype, **kwargs):
    """
    Write the metadata, imagedata into ouput file named outbase.

    kwargs can be passed to the writer specified by filetype.

    Parameters
    ----------
    metadata : dataset object
        dataset from nimsdata.parse.
    imagedata : dict
        dictionary of data with string labels as keys.
    outbase : string
        base of name to use, without any file extension.
    filetype : string
        string name of writer to use.
    kwargs : dict
        keyword arguments to be passed to writer.

    Returns
    -------
    output_list : list
        list of created output filepaths.

    Raises
    ------
    NIMSDataError
        filetype parameters was not specified.

    Examples
    --------
    convert dicom.tgz to dicoms.nii.gz

    .. code-block:: python

        import nimsdata
        ds = nimsdata.parse('dicoms.tgz', load_data=False, filetype='dicom')
        ds.load_data()
        nimsdata.write(ds, ds.data, outbase, filetype='nifti')

    """
    log.debug('write start: %s' % str(datetime.datetime.now()))
    if not filetype:
        raise DataError('filetype cannot be None')  # XXX FAIL! unexpected to get no filetype

    if metadata is None:
        raise DataError('metadata cannot be None')
    if imagedata is None:
        raise DataError('imagedata cannot be None')

    writer = get_handler(filetype, WRITERS)  # raises exception if no handler
    output_list = writer.write(metadata, imagedata, outbase, **kwargs)  # writer checks if data is present

    log.debug('write end: %s' % str(datetime.datetime.now()))
    log.debug('generated: %s' % str(output_list))
    return output_list


class DataError(Exception):

    """
    Base class for NIMS exceptions.

    Parameters
    ----------
    message : str
        Message to report along with the exception.
    log_level : logging level [default None]
        What logging level to report the error with, must be a valid logging level class, such as
        `logging.DEBUG` or `logging.ERROR`.

    """

    def __init__(self, message, log_level=None):
        """instantiate nimsdataerror exception class."""
        super(DataError, self).__init__(message)
        if log_level is not None:
            message = '%s\n%s' % (message, traceback.format_exc())
            log.log(log_level, message)


class Reader(object):

    """
    Abstract base class that provides interfaces, and functions for readers.

    Cannot be instantiated.  See :doc:`extending_nimsdata` for more information on subclassing
    NIMSReader to create a new data reader.

    Parameters
    ----------
    path : str
        filepath as string
    load_data : boolean [default False]
        Indicate if parse should return with all data loaded.

    """

    __metaclass__ = abc.ABCMeta

    project_properties = project_properties
    session_properties = session_properties
    acquisition_properties = acquisition_properties

    def _schema_init(self, schema):
        for k, v in schema.iteritems():
            if isinstance(v, dict):
                    self._schema_init(v)
            elif k == 'field' and not v.startswith('nims_'):
                setattr(self, v, None)

    @abc.abstractmethod
    def __init__(self, path, load_data=False):
        self._schema_init(self.project_properties)
        self._schema_init(self.session_properties)
        self._schema_init(self.acquisition_properties)
        self.filepath = os.path.abspath(path)
        self.data = None
        self.metadata_status = 'empty'

    @abc.abstractmethod
    def load_data(self):
        pass

    @abc.abstractproperty
    def nims_metadata_status(self):
        pass

    @abc.abstractproperty
    def nims_group_id(self):
        pass

    @abc.abstractproperty
    def nims_project(self):
        pass

    @abc.abstractproperty
    def nims_session_id(self):
        pass

    @abc.abstractproperty
    def nims_session_label(self):
        pass

    @abc.abstractproperty
    def nims_session_subject(self):
        pass

    @abc.abstractproperty
    def nims_acquisition_id(self):
        pass

    @abc.abstractproperty
    def nims_acquisition_label(self):
        pass

    @abc.abstractproperty
    def nims_acquisition_description(self):
        pass

    @abc.abstractproperty
    def nims_file_name(self):
        pass

    @abc.abstractproperty
    def nims_file_ext(self):
        pass

    @abc.abstractproperty
    def nims_file_domain(self):
        pass

    @abc.abstractproperty
    def nims_file_type(self):
        pass

    @abc.abstractproperty
    def nims_file_kinds(self):
        pass

    @abc.abstractproperty
    def nims_file_state(self):
        pass

    @abc.abstractproperty
    def nims_timestamp(self):
        pass

    @abc.abstractproperty
    def nims_timezone(self):
        pass

    def __str__(self):  # TODO: tighten this up. or just get rid of it completely
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
            if value is None:
                value = ''
            properties.append((prop, str(value)))

        return '\n'.join(['%-30s: %s' % (p, v) for p, v in properties])


class abstractclassmethod(classmethod):

    """
    Decorator for abstract classmethods.

    .. code-block:: python

        class C(object):
            __metaclass__ = abc.ABCMeta
            @abstractclassmethod
            def my_abstract_classmethod(cls, ...):
                ...

    """

    __isabstractmethod__ = True

    def __init__(self, callable):
        callable.__isabstractmethod__ = True
        super(abstractclassmethod, self).__init__(callable)


class Writer(object):

    """
    Abstract base class that provides interfaces and functions for writers.

    Cannot be instantiated.  See :doc:`extending_nimsdata` for more information on subclassing
    NIMSReader to create a new data writer.

    """

    __metaclass__ = abc.ABCMeta

    @abstractclassmethod
    def write(cls, metadata, imagedata, outbase, **kwargs):
        """
        Write metadata and imagedata to output file outbase.

        Abstract implementaiton doesn't REALLY need to know about filepath.  Subclasses may want to
        build upon write, such as performing voxel reordering of MR data.

        Parameters
        ----------
        metadata : object
            fully loaded instance of a NIMSReader.
        imagedata : dict
            dictionary of np.darrays. label suffix as keys, with np.darrays as values.
        outbase : str
            output name prefix.
        **kwargs :
            all remaining keyword arguments will be passed to the underlying file writer.

        Returns
        -------
        None : NoneType

        Raises
        ------
        NIMSDataError
            metadata or data is None.

        """
        if metadata is None or imagedata is None:
            raise DataError('no metadata/data?')  # XXX FAIL! unexpected to get no data


if __name__ == '__main__':
    import sys
    import argparse

    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import data

    log = logging.getLogger('data')

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='file to convert')
    parser.add_argument('outbase', nargs='?', help='basename for output files (default: input)')
    parser.add_argument('-p', '--parser', help='parser to use', choices=data.data.READERS.keys())
    parser.add_argument('-i', '--ignore_json', help='do not use json metadata in output metadata', action='store_true', default=False)
    parser.add_argument('-w', '--writer', help='write to use', choices=data.data.WRITERS.keys())
    parser.add_argument('-v', '--verbose', help='enable verbose logging', dest='verbose', action='store_true', default=False)
    parser.add_argument(      '--parser_kwarg', action='append', help='keyword arguments to pass directly to the parser')
    parser.add_argument(      '--writer_kwarg', action='append', help='keyword arguments to pass directly to the writer')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if not os.path.exists(args.input):
        raise DataError('could not find input file %s' % args.input)

    outbase = args.outbase or os.path.basename(os.path.splitext(args.input.rstrip('/'))[0])

    def cast_if_number(s):
        try:
            return float(s) if '.' in s else int(s)
        except ValueError:
            return s

    p_kwargs = {}
    if args.parser_kwarg:
        for item in args.parser_kwarg:
            kw, val = item.split('=')
            p_kwargs[kw] = cast_if_number(val)
    log.debug(p_kwargs)

    w_kwargs = {}
    if args.writer_kwarg:
        for item in args.writer_kwarg:
            kw, val = item.split('=')
            w_kwargs[kw] = cast_if_number(val)
    log.debug(w_kwargs)

    ds = data.parse(args.input, load_data=True, ignore_json=args.ignore_json, filetype=args.parser, **p_kwargs)

    if not ds:
        raise DataError('%s could not be parsed' % args.input)
    if ds.data is None:
        raise DataError('%s has no data' % args.input)

    data.write(ds, ds.data, outbase, filetype=args.writer, **w_kwargs)
