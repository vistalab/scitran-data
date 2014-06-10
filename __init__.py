# @author:  Gunnar Schaefer
#           Kevin S Hahn

"""
nimsdata
========

The nimsdata package provides two interfaces; one for reading from file, and the other for writing to file.

Readers and Writers are defined in external json files, readers.json and writers.json.  New readers and writers can be added to
a running system by defining it in the appropriate json file.

"""

import os
import json
import logging
import tarfile

import bson.json_util

import nimsdata

log = logging.getLogger(__name__)


dict_merge = nimsdata.dict_merge
NIMSDataError = nimsdata.NIMSDataError

READERS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'readers.json')))
WRITERS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'writers.json')))
# note: readers/writers.json will never contain datetime, or objects that require bson.json_util
# TODO: reader/writer loader


def get_handler(name, handlerdict):
    """
    Retrieve the class object of reader or writer from its string label.

    Parameters
    ----------
    name : str
        reader or writer name, by string
    handlerdict : dict
        name of which handler dict to search from, READERS or WRITERS

    Returns
    -------
    handler : NIMSReader subclass
        NIMSReader subclass capable of reading
    """
    handler = None
    try:
        module, klass = handlerdict.get(name).split('.')
        handler = getattr(__import__(module, globals()), klass)
    except AttributeError:
        raise NIMSDataError('no handler')
    return handler


# TODO: do we even need this anymore?
def getnattr(obj, attr):
    """
    Get nested attribute.

    getnattr(foo, 'bar.face')

    This function is not needed anymore. safe to delete?

    """
    if not attr:
        return
    for name in attr.split('.'):
        obj = getattr(obj, name)
    return obj


def parse(path, load_data=False, ignore_json=False, filetype=None, **kwargs):
    """
    Infer the type of file at path and parse with filetype-specific parser.

    The parse function expects a tgz of raw input files and a metadata.json file that
    declares the filetype, and any metadata that should be overwritten.  The parse
    function will pass the input file to the appropriate handler, if available.

    Parameters
    ----------
    path : str
        input path
    load_data : bool  [default False]
        attempt to load all metadata and data at parse.  like running parse, then load_data
    ignore_json : bool [default False]
        True, don't look for json file. therefore the parser cannot be read from the json, and
        a parser MUST be specfied if ignore_json is true
        False, do look for json file.
    filetype : str
        string name of parser to use.  no special search or inferences will be made if this
        option is specified. if the specified parser does not exist, an error will be raised.

    Returns
    -------
    parser_class : obj
        NIMSReader populated with data and metadata attributes. this NIMSReader object can be
        passed to any compatible writer to complete the conversion process

    Notes
    -----
    =========== ======== ==================================================
    ignore_json filetype  behavior
    =========== ======== ==================================================
    False       None     **DEFAULT** try to read json, use parser from json
    False       'dicom'  try to read json, use 'dicom' parser
    True        None     **INVALID**
    True        'dicom'  don't try to read json AT ALL, use 'dicom' parser
    =========== ======== ==================================================

    TODO: expand the notes table section to be MUCH more descriptive of the behavior
    of nimsdata.parse's multiple combinations of options.

    Examples
    --------
    convert dicom.tgz to dicoms.nii.gz

    .. code-block:: python

        import nimsdata
        ds = nimsdata.parse('dicoms.tgz', load_data=False, filetype='dicom')
        ds.load_data()
        nimsdata.write(ds, ds.data, outbase, filetype='nifti')

    TODO: convert all in-line code EXAMPLES to using doc-string code.  and then run test with
    doctests. this makes sure that the provided examples actually work. The only downside,
    is that provided code examples must be very explicit.

    """
    if not os.path.exists(path):
        raise NIMSDataError('input path %s not found' % path, log_level=logging.ERROR)

    if ignore_json:
        if not filetype:
            raise NIMSDataError('filetype must be specified if ignore_json=True', log_level=logging.ERROR)
        else:
            log.debug('ignoring json entirely')
            nimsparser = get_handler(filetype, READERS)
            return nimsparser(path, load_data, **kwargs)

    json_data = {}
    if tarfile.is_tarfile(path):
        log.debug('inspecting tgz for json')
        with tarfile.open(path) as archive:
            log.debug('"%s" is tarfile' % path)
            for ti in archive:
                if not ti.isreg():
                    continue
                try:
                    json_data = json.loads(archive.extractfile(ti).read(), object_hook=bson.json_util.object_hook)
                except Exception:
                    pass
                else:
                    log.debug('json found, %s' % ti.name)
                    if not filetype:
                        filetype = json_data.get('filetype')
                        logging.debug('filetype from json: %s' % filetype)
                    break
    elif os.path.isdir(path):
        # TODO: implement 'break-at-first-readable' json read for directory
        raise NIMSDataError('directory input not implemented', log_level=logging.ERROR)

    elif os.path.isfile(path):
        # TODO: infer type
        # single file, by definition, will not have a json with it
        # some P files may come in as tgz, with parser inside, OR as single P12345.7 or P12345.7.gz
        if path.endswith('.7.gz') or path.endswith('.7'):
            filetype = 'pfile'
        else:
            raise NIMSDataError('non tar-files not implemented', log_level=logging.ERROR)

    nimsparser = get_handler(filetype, READERS)
    ds = nimsparser(path, load_data, **kwargs)

    # FIXME: handle NESTED information
    for key, value in json_data.get('overwrite', {}).iteritems():
        setattr(ds, key, value)

    return ds


def write(metadata, imagedata, outbase, filetype, **kwargs):
    """
    Write the metadata, imagedata, into file outbase, with the specified filetype writer.

    Parameters
    ----------
    metadata : dataset object
        dataset from parse
    imagedata : dict
        dictionary of np.arrays with string labels as keys.  The primary dataset should
        be in imagedata[''], while secondary data, such as fieldmap, should be stored in
        imagedata['_fieldmap'].
    outbase : string
        base of name to use
    filetype : string
        string name of writer to use

    Returns
    -------
    output_list : list
        list of created output filepaths

    Examples
    --------
    convert dicom.tgz to dicoms.nii.gz

    .. code-block:: python

        import nimsdata
        ds = nimsdata.parse('dicoms.tgz', load_data=False, filetype='dicom')
        ds.load_data()
        nimsdata.write(ds, ds.data, outbase, filetype='nifti')

    """
    if not filetype:
        raise NIMSDataError('filetype cannot be None')

    nimswriter = get_handler(filetype, WRITERS)
    # think about this return value. what does the processor really need?
    output_list = nimswriter.write(metadata, imagedata, outbase, **kwargs)
    log.debug('generated: %s' % str(output_list))
    return output_list
