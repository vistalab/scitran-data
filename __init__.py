# @author:  Gunnar Schaefer
#           Kevin S Hahn

import os
import json
import logging
import tarfile
import cStringIO

import nimsdata
import nimsdicom
import nimsmontage

NIMSDataError = nimsdata.NIMSDataError

log = logging.getLogger('nimsdata')


PARSERS = {'gedicom': 'nimsdicomge.NIMSDicomGE',
           'siemensdicom': 'nimsdicomsiemens.NIMSDicomSiemens'}


def _open(path):
    """open an input path.

    Parameters
    ----------
    path : filepath
        directory path of input as string.  accepts a directory, or a tarfile, or a tarfile like object

    Returns
    -------
    returns open tarfile object
    """
    if isinstance(path, cStringIO.OutputType):
        log.debug('opening cstringio file object')
        path.flush()
        path.seek(0)
        archive = tarfile.open(fileobj=path, mode='r')
    elif os.path.isfile(path) and tarfile.is_tarfile(path):
        log.debug('opening tar archive')
        archive = tarfile.open(path)
    elif os.path.isdir(path):
        log.debug('creating fake archive from dir')
        f = cStringIO.StringIO()
        archive = tarfile.open(fileobj=f, mode='w')     # write w/o gz; faster, but more mem
        archive.add(path)                               # add directory
        archive.close()                                 # close archive, to re-open as READ
        f.flush()                                       # clear internal buffer
        f.seek(0)                                       # reset position, so iteration starts at beginning
        archive = tarfile.open(fileobj=f, mode='r')
    return archive


def parse(path, ignore_json=False, parser=None):
    """locate json file inside tarfile to identify parser class.

    Parameters
    ----------
    path : filepath
        directory path of input as a string. accepts a directory, or a tarfile

    ignore_json : True|False
        option to ignore metadata in json. True disables. default False.

    parser_id : string id of parser
        parser_id as a string, corresponds to parsers dict.

    Returns
    -------
    instance of NIMSData subclass, which has load_data() and convert() methods
    """
    json_data = None
    ds = None
    archive = _open(path)
    for ti in (ti for ti in archive if ti.isreg()):
        try:
            json_data = json.loads(archive.extractfile(ti).read())
        except Exception:
            pass
        else:
            parser = json_data.get('parser') or parser
            break
    log.debug('looking up parser')
    # TODO: parse should give USEFUL feedback for errors!! re think structuring and error handling
    try:
        module, klass = PARSERS.get(parser).split('.')               # Attribute error for attempting to split None
        nimsparser = getattr(__import__(module, globals()), klass)      # Attribute error for no attribute (read as: class not found)
    except AttributeError:
        raise AttributeError('no matching parser')
    else:
        log.debug('parsing with %s for %s' % (nimsparser, parser))
        ds = nimsparser(archive)                        # lots of action here
        ds.filepath = archive.name or archive.fileobj   # keep filepath for tarfiles, or object reference for file_obj (in case of directory)
        archive.close()
    if ds and json_data and not ignore_json:
        log.debug('updating dataset metadata with json metadata')
        try:
            [setattr(ds.metadata, key, json_data[key]) for key in json_data if key != 'parser']
        except TypeError:
            pass                                        # TypeError, json contains ONLY parser, nothing else
    return ds


def load_data(ds):
    """invoke NIMSData subclass specific load_data method

    Parameters
    ----------
    ds : dataset
        subclass of NIMSData, returned by parse()

    Returns
    -------
    subclass of NIMSData object, which has load_data() and convert() methods
    """
    archive = _open(ds.filepath)
    ds.load_data(archive)
    archive.close()
    return ds


def convert(ds, outbase):
    """invoke NIMSData subclass specific convert method.

    Parameters
    ----------
    ds : dataset
        subclass of NIMSdata

    outbase : outbase
        string of desired base of output name, do not include file extension

    Returns
    -------
    None
    """
    if ds.data is None:
        ds = load_data(ds)
    result = ds.convert(outbase)
    return result
