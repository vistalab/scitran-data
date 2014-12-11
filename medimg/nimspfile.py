#!/usr/bin/env python
#
# @author:  Gunnar Schaefer
#           Kevin S Hahn

"""
nimsdata.medimg.nimspfile
=========================

This module provides functions, classes and errors for fully minimally parsing
and reconstructing pfiles. Additional modules are required to enable
full parsing of pfiles, spiral reconstruction, and mux_epi reconstruction.

"""

import os
import bson
import glob
import gzip
import json
import time
import shlex
import struct
import logging
import tarfile
import datetime
import subprocess
import bson.json_util

import numpy as np

import medimg
import dcm.mr.ge
import dcm.mr.generic_mr

from .. import tempdir as tempfile

log = logging.getLogger(__name__)


def unpack_uid(uid):
    """
    Convert packed PFile UID to standard DICOM UID.

    Parameters
    ----------
    uid : str
        packed PFile UID as a string
    Returns
    -------
    uid : str
        unpacked PFile UID as string

    """
    return ''.join([str(i-1) if i < 11 else '.' for pair in [(ord(c) >> 4, ord(c) & 15) for c in uid] for i in pair if i > 0])


def is_gzip(filepath):
    """
    Convert packed PFile UID to standard DICOM UID.

    Parameters
    ----------
    uid : str
        packed PFile UID as a string
    Returns
    -------
    uid : str
        unpacked PFile UID as string

    """
    with open(filepath, 'rb') as fp:
        compressed = (fp.read(2) == '\x1f\x8b')
    return compressed


def get_version(filepath):
    """
    Determine the pfile version of the file at filepath.

    An NIMSPFileError exception will be raised if the file is not a valid PFile.

    Parameters
    ----------
    filepath : str
        filepath of file to check

    Returns
    -------
    version : str
        PFile version number of file at filepath

    Raises
    ------
    NIMSPFileError : Exception
        error if the file is not a valid PFile

    """
    fileobj = gzip.open(filepath, 'rb') if is_gzip(filepath) else open(filepath, 'rb')

    version_bytes = fileobj.read(4)
    fileobj.seek(34); logo = (struct.unpack("10s", fileobj.read(struct.calcsize("10s")))[0]).split('\0', 1)[0]
    if version_bytes == '\x00\x00\xc0A':
        version = 24
    elif version_bytes == 'V\x0e\xa0A':
        version = 23
    elif version_bytes == 'J\x0c\xa0A':
        version = 22
    elif version_bytes == '\x00\x000A':
        version = 12
    else:
        raise NIMSPFileError(fileobj.name + ' is not a valid PFile or of an unsupported version')
    if logo != 'GE_MED_NMR' and logo != 'INVALIDNMR':
        raise NIMSPFileError(fileobj.name + ' is not a valid PFile')
    fileobj.close()
    return version


class NIMSPFileError(medimg.MedImgError):
    pass


class NIMSPFile(medimg.MedImgReader):

    """
    Parse and load data from a pfile.

    This class reads the data and/or header from a pfile, runs k-space reconstruction.

    NIMSPFile object can handle several different input types
        - .tgz of directory containing Pfile, and supporting files such as ref.dat, vrgf.dat and tensor.dat.
        - a single pfile, either gz or uncompressed.

    tgz cannot be "full parsed".  setting full_parse=True, with an input tgz, will raise an exception.
    nims2 input tgz format
    Pfile.7, Pfile.7ref.dat, Pfile.7vrgf.dat Pfile.7ref

    .. code:: python

        import nimsdata
        ds = nimsdata.parse('pfile.tgz', filetype='pfile', load_data=True)
        if not ds.failure_reason:
            nimsdata.write(ds, ds.data, outbase='output_name', filetype='nifti')

    Some pfiles require calibration files from another scan.  This 'aux_file' can be provided
    during `__init__()`, or `load_data()`.

    .. code:: python

        import nimsdata
        ds = nimsdata.parse('muxarcepi_nocal.tgz', filetype='pfile', load_data=True, aux_file='muxarcepi_cal.tgz')
        if not ds.failure_reason:
            nimsdata.write(ds, ds.data, outbase='output_name', filetype='nifti')

    .. code:: python

        import nimsdata
        ds = nimsdata.parse('muxarcepi_nocal.tgz', filetype='pfile', load_data=False)
        ds.load_data(aux_file='muxarcepi_cal.tgz')
        if no ds.failure_reason:
            nimsdata.write(ds, ds.data, outbase='output_name', filetype='nifti')

    """

    domain = u'mr'
    filetype = u'pfile'
    parse_priority = 5
    state = ['orig']

    def __init__(self, filepath, load_data=False, full_parse=False, tempdir=None, aux_file=None, num_jobs=4, num_virtual_coils=16, notch_thresh=0, recon_type=None):
        """
        Read basic sorting information.

        There are a lot of parameters; most of the parameters only apply to mux_epi scans. The muxepi only
        parameters are num_jobs, num_virtual_coils, notch_thresh, recon_type and aux_file.

        Parameters
        ----------
        filepath : str
            path to pfile.7 or pfile.tgz
        load_data : bool [default False]
            load all data and run reconstruction
        full_parse : bool [default False]
            full parse the input file, only applies to pfile.7 inputs
        tempdir : str
            path prefix to use for temp directory
        num_jobs : int
            muxepi only, number of simultaneous jobs
        num_virtual_coils : int
            muxepi only, number of virtual coils
        notch_thresh : int
            muxepi only, number of virtual coils
        recon_type : NoneType or str
            muxepi only, if recon_type is 'sense', then run sense recon
        aux_file : None or str
            path to pfile.tgz that contains valid vrgf.dat and ref.dat files

        """
        super(NIMSPFile, self).__init__(filepath)       # sets self.filepath
        self.full_parsed = False                        # indicates if fully parsed
        self.dirpath = os.path.dirname(self.filepath)   # what contains the input file
        self.basename = os.path.basename(self.filepath)
        # TODO setting the file name and extension should be different for .7 and .7.tgz
        # if pfile_arc.tgz, file_name = pfile_arc, file_ext = .tgz
        # if P?????.7,   file_name = P?????, file_ext = .7
        self.file_name, self.file_ext = os.path.splitext(self.filepath)
        self.num_jobs = num_jobs
        self.num_vcoils = num_virtual_coils
        self.notch_thresh = notch_thresh
        self.recon_type = recon_type
        self.aux_file = aux_file
        self.tempdir = tempdir
        self.data = None


        log.debug('parsing %s' % filepath)
        if tarfile.is_tarfile(self.filepath):  # tgz; find json with a ['header'] section
            log.debug('tgz')
            with tarfile.open(self.filepath) as archive:
                for ti in archive:
                    if not ti.isreg():
                        continue
                    try:
                        _hdr = json.load(archive.extractfile(ti), object_hook=bson.json_util.object_hook)['header']
                    except ValueError as e:  # json file does not exist
                        log.debug('%s; not a json file' % e)
                    except KeyError as e:  # header section does not exist
                        log.debug('%s; header section does not exist' % e)
                    else:
                        log.debug('_min_parse_tgz')
                        self.exam_uid = _hdr.get('session')
                        self.acquisition_id = _hdr.get('acquisition')
                        self.timestamp = _hdr.get('timestamp')
                        self.group_name = _hdr.get('group')
                        self.project_name = _hdr.get('project')
                        self.metadata_status = 'pending'
                        break
                else:
                    raise NIMSPFileError('no json file with header section found. bailing', log_level=logging.WARNING)
        else:  # .7 or .7.gz, doing it old world style
            try:
                self.version = get_version(self.filepath)
                self._full_parse(self.filepath) if full_parse else self._min_parse(self.filepath)  # full_parse arg indicates run full_parse
            except Exception as e:
                raise NIMSPFileError('not a PFile? %s' % str(e))

        if load_data:
            self.load_data()

    def infer_psd_type(self):
        """
        Infer the psd type based on self.psd_type.

        Also makes any corrections to the psd_type to account for mis-named psds.

        Returns
        -------
        None : NoneType
            sets self.psd_type

        """
        dcm.mr.ge.infer_psd_type(self)
        if self.psd_type == 'epi' and int(self._hdr.rec.user6) > 0:  # XXX HACK check for misnamed mux scans
            self.psd_type = 'muxepi'
        log.debug('psd_name: %s, psd_type: %s' % (self.psd_name, self.psd_type))

    def infer_scan_type(self):
        """
        Infer the scan type based on the dataset attributes.

        Returns
        -------
        None : NoneType
            sets self.scan_type

        """
        dcm.mr.generic_mr.infer_scan_type(self)
        log.debug('scan_type: %s' % self.scan_type)

    def _min_parse(self, filepath=None):
        """
        Parse the minimum sorting information from a pfile.7.

        Does not work if input file is a tgz.  If NIMSPfile was init'd with a tgz input, the tgz can be
        unpacked into a temporary directory, and then this function can parse the unpacked pfile.

        Parameters
        ----------
        filepath : str
            path to a pfile.7.  Does not accept pfile.tgz.

        """
        filepath = filepath or self.filepath  # use filepath if provided, else fall back to self.filepath
        if tarfile.is_tarfile(filepath):
            raise NIMSPFileError('_min_parse() expects a .7 or .7.gz')
        log.debug('_min_parse of %s' % filepath)

        fileobj = gzip.open(filepath, 'rb') if is_gzip(self.filepath) else open(filepath, 'rb')

        fileobj.seek(16); self.scan_date = str(struct.unpack("10s", fileobj.read(struct.calcsize("10s")))[0])
        fileobj.seek(26); self.scan_time = str(struct.unpack("8s", fileobj.read(struct.calcsize("8s")))[0])
        fileobj.seek(64); self.num_timepoints = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
        fileobj.seek(70); self.num_echos = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
        fileobj.seek(216); self.rec_user0 = struct.unpack("f", fileobj.read(struct.calcsize("f")))[0]
        fileobj.seek(240); self.rec_user6 = struct.unpack("f", fileobj.read(struct.calcsize("f")))[0]
        fileobj.seek(244); self.rec_user7 = struct.unpack("f", fileobj.read(struct.calcsize("f")))[0]
        fileobj.seek(914); self.ileaves = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
        if self.version in [24, 23, 22]:
            fileobj.seek(143516); self.exam_no = str(struct.unpack("H", fileobj.read(struct.calcsize("H")))[0])
            fileobj.seek(145622); self.series_no = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
            fileobj.seek(145762); self.series_desc = (struct.unpack("65s", fileobj.read(struct.calcsize("65s")))[0]).split('\0', 1)[0]
            fileobj.seek(145875); self.series_uid = unpack_uid(struct.unpack("32s", fileobj.read(struct.calcsize("32s")))[0])
            fileobj.seek(148388); self.im_datetime = struct.unpack("i", fileobj.read(struct.calcsize("i")))[0]
            fileobj.seek(148396); self.tr = struct.unpack("i", fileobj.read(struct.calcsize("i")))[0] / 1e6
            fileobj.seek(148834); self.acq_no = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
            fileobj.seek(148972); self.psd_name = os.path.basename(struct.unpack("33s", fileobj.read(struct.calcsize("33s")))[0]).split('\0', 1)[0].lower()
        if self.version in [24, 23]:
            fileobj.seek(144248); self.exam_uid = unpack_uid(struct.unpack("32s", fileobj.read(struct.calcsize("32s")))[0])
            fileobj.seek(144409); self.patient_id = (struct.unpack("65s", fileobj.read(struct.calcsize("65s")))[0]).split('\0', 1)[0]
        if self.version == 22:
            fileobj.seek(144240); self.exam_uid = unpack_uid(struct.unpack("32s", fileobj.read(struct.calcsize("32s")))[0])
            fileobj.seek(144401); self.patient_id = (struct.unpack("65s", fileobj.read(struct.calcsize("65s")))[0]).split('\0', 1)[0]
        if self.version == 12:
            fileobj.seek(61576); self.exam_no = str(struct.unpack("H", fileobj.read(struct.calcsize("H")))[0])
            fileobj.seek(61966); self.exam_uid = unpack_uid(struct.unpack("32s", fileobj.read(struct.calcsize("32s")))[0])
            fileobj.seek(62127); self.patient_id = (struct.unpack("65s", fileobj.read(struct.calcsize("65s")))[0]).split('\0', 1)[0]
            fileobj.seek(62710); self.series_no = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
            fileobj.seek(62786); self.series_desc = (struct.unpack("65s", fileobj.read(struct.calcsize("65s")))[0]).split('\0', 1)[0]
            fileobj.seek(62899); self.series_uid = unpack_uid(struct.unpack("32s", fileobj.read(struct.calcsize("32s")))[0])
            fileobj.seek(65016); self.im_datetime = struct.unpack("i", fileobj.read(struct.calcsize("i")))[0]
            fileobj.seek(65024); self.tr = struct.unpack("i", fileobj.read(struct.calcsize("i")))[0] / 1e6
            fileobj.seek(65328); self.acq_no = struct.unpack("h", fileobj.read(struct.calcsize("h")))[0]
            fileobj.seek(65374); self.psd_name = os.path.basename(struct.unpack("33s", flleobj.read(struct.calcsize("33s")))[0]).split('\0', 1)[0].lower()

        if self.im_datetime > 0:
            self.timestamp = datetime.datetime.utcfromtimestamp(self.im_datetime)
        else:
            month, day, year = map(int, self.scan_date.split('\0', 1)[0].split('/'))
            hour, minute = map(int, self.scan_time.split('\0', 1)[0].split(':'))
            self.timestamp = datetime.datetime(year + 1900, month, day, hour, minute)  # GE's epoch begins in 1900

        self.infer_psd_type()
        if self.psd_type == 'spiral':
            self.num_timepoints = int(self.rec_user0)
        elif self.psd_type == 'basic':
            self.num_timepoints = (self.num_timepoints * self.num_echos - 6) / 2
        elif self.psd_type == 'muxepi':
            self.num_timepoints = self.num_timepoints + int(self.rec_user6) * self.ileaves * (int(self.rec_user7) - 1)
        self.prescribed_duration = self.num_timepoints * self.tr
        self.subj_code, self.group_name, self.project_name = medimg.parse_patient_id(self.patient_id, 'ex' + self.exam_no)
        self.metadata_status = 'pending'

    def _full_parse(self, filepath=None):
        """
        Fully parse the input pfile.

        Attempts to import pfile version specific parser from pfile submodule.  Full parse is
        not possible without access to the pfile submodule.

        Does not work if input file is a tgz.  If NIMSPfile was init'd with a tgz input, the tgz can be
        unpacked into a temporary directory, and then this function can parse the unpacked pfile.

        Parameters
        ----------
        filepath : str
            path to a pfile.7.  Does not accept pfile.tgz.

        """
        filepath = filepath or self.filepath
        if tarfile.is_tarfile(filepath):
            raise NIMSPFileError('_full_parse() expects a .7 or .7.gz')
        log.debug('_full_parse of %s' % filepath)

        try:
            pfile = getattr(__import__('pfile.pfile%d' % self.version, globals()), 'pfile%d' % self.version)
        except ImportError:
            raise ImportError('no pfile parser for v%d' % self.version)

        with gzip.open(filepath, 'rb') if is_gzip(filepath) else open(filepath, 'rb') as fileobj:
            self._hdr = pfile.POOL_HEADER(fileobj)
            if not self._hdr:
                raise NIMSPFileError('no pfile was read', log_level=logging.WARNING)
            self.data = None    # data always starts as None

            self.pfilename = 'P%05d' % self._hdr.rec.run_int

            self.exam_no = self._hdr.exam.ex_no
            self.exam_uid = unpack_uid(self._hdr.exam.study_uid)
            self.series_no = self._hdr.series.se_no
            self.series_desc = self._hdr.series.se_desc.split('\0', 1)[0]
            self.series_uid = unpack_uid(self._hdr.series.series_uid)
            self.acq_no = self._hdr.image.scanactno
            self.patient_id = self._hdr.exam.patidff.split('\0', 1)[0]
            self.subj_code, self.group_name, self.project_name = medimg.parse_patient_id(self.patient_id, 'ex' + str(self.exam_no))
            self.subj_firstname, self.subj_lastname = medimg.parse_patient_name(self._hdr.exam.patnameff.split('\0', 1)[0])
            self.subj_dob = medimg.parse_patient_dob(self._hdr.exam.dateofbirth.split('\0', 1)[0])
            self.subj_sex = ('male', 'female')[self._hdr.exam.patsex-1] if self._hdr.exam.patsex in [1, 2] else None

            self.psd_name = os.path.basename(self._hdr.image.psdname.partition('\x00')[0]).lower()
            # self.scan_type = self._hdr.image.psd_iname.split('\0', 1)[0]  # XXX is this needed, it gets overwritten by end of fullparse
            if self._hdr.image.im_datetime > 0:
                self.timestamp = datetime.datetime.utcfromtimestamp(self._hdr.image.im_datetime)
            else:   # HOShims don't have self._hdr.image.im_datetime
                month, day, year = map(int, self._hdr.rec.scan_date.split('\0', 1)[0].split('/'))
                hour, minute = map(int, self._hdr.rec.scan_time.split('\0', 1)[0].split(':'))
                self.timestamp = datetime.datetime(year + 1900, month, day, hour, minute)  # GE's epoch begins in 1900

            # expose study date, study time, acquisition date, acquisition time
            month, day, year = map(int, self._hdr.rec.scan_date.split('\0', 1)[0].split('/'))
            hour, minute = map(int, self._hdr.rec.scan_time.split('\0', 1)[0].split(':'))
            self.study_date = '%4d%02d%02d' % (year + 1900, month, day)
            self.study_time = '%02d%02d%02d' % (hour, minute, 0)
            self.study_datetime = self.study_date and self.study_time and datetime.datetime.strptime(self.study_date + self.study_time[:6], '%Y%m%d%H%M%S')

            if self._hdr.image.im_datetime > 0:
                self.acq_datetime = datetime.datetime.utcfromtimestamp(self._hdr.image.im_datetime)
                self.acq_date = datetime.datetime.strftime(self.acq_datetime, '%Y%m%d')
                self.acq_time = datetime.datetime.strftime(self.acq_datetime, '%H%M%S')
            else:
                self.acq_datetime = None
                self.acq_date = None
                self.acq_time = None
            self.ti = self._hdr.image.ti / 1e6
            self.te = self._hdr.image.te / 1e6
            self.tr = self._hdr.image.tr / 1e6  # tr in seconds
            self.flip_angle = float(self._hdr.image.mr_flip)
            self.pixel_bandwidth = self._hdr.rec.bw
            # Note: the freq/phase dir isn't meaningful for spiral trajectories.
            # GE numbers the dims 1,2, so freq_dir==1 is the first dim. We'll use
            # the convention where first dim = 0, second dim = 1, etc. for phase_encode.
            self.phase_encode = 1 if self._hdr.image.freq_dir == 1 else 0
            self.mt_offset_hz = self._hdr.image.offsetfreq
            self.num_slices = self._hdr.image.slquant
            self.num_averages = self._hdr.image.averages
            self.num_echos = self._hdr.rec.nechoes
            self.receive_coil_name = self._hdr.image.cname.split('\0', 1)[0]
            self.num_receivers = self._hdr.rec.dab[0].stop_rcv - self._hdr.rec.dab[0].start_rcv + 1
            self.operator = self._hdr.exam.operator_new.split('\0', 1)[0]
            self.protocol_name = self._hdr.series.prtcl.split('\0', 1)[0]
            self.scanner_name = self._hdr.exam.hospname.split('\0', 1)[0] + ' ' + self._hdr.exam.ex_sysid.split('\0', 1)[0]
            self.scanner_type = 'GE MEDICAL SYSTEMS DISCOVERY MR750'    # FIXME: don't hardcode
            self.acquisition_type = None                                # hope this doesn't break anything...
            self.size = [self._hdr.image.dim_X, self._hdr.image.dim_Y]  # imatrix_Y
            self.fov = [self._hdr.image.dfov, self._hdr.image.dfov_rect]
            self.num_bands = 1
            self.num_mux_cal_cycle = 0
            self.num_timepoints = self._hdr.rec.npasses
            # Some sequences (e.g., muxepi) acuire more timepoints that will be available in the resulting data file.
            # The following will indicate how many to expect in the final image.
            self.num_timepoints_available = self.num_timepoints
            self.deltaTE = 0.0
            self.scale_data = False
            # Compute the voxel size rather than use image.pixsize_X/Y
            self.mm_per_vox = [self.fov[0] / self.size[0], self.fov[1] / self.size[1], self._hdr.image.slthick + self._hdr.image.scanspacing]
            image_tlhc = np.array([self._hdr.image.tlhc_R, self._hdr.image.tlhc_A, self._hdr.image.tlhc_S])
            image_trhc = np.array([self._hdr.image.trhc_R, self._hdr.image.trhc_A, self._hdr.image.trhc_S])
            image_brhc = np.array([self._hdr.image.brhc_R, self._hdr.image.brhc_A, self._hdr.image.brhc_S])

            # psd-specific params get set here
            self.infer_psd_type()
            if self.psd_type == 'spiral':
                self.num_timepoints = int(self._hdr.rec.user0)    # not in self._hdr.rec.nframes for sprt
                self.deltaTE = self._hdr.rec.user15
                self.band_spacing = 0
                self.scale_data = True
                # spiral is always a square encode based on the frequency encode direction (size_x)
                # Atsushi also likes to round up to the next higher power of 2.
                # self.size_x = int(pow(2,ceil(log2(pf.size_x))))
                # The rec.im_size field seems to have the correct reconned image size, but
                # this isn't guaranteed to be correct, as Atsushi's recon does whatever it
                # damn well pleases. Maybe we could add a check to ninfer the image size,
                # assuming it's square?
                self.size_x = self.size_y = self._hdr.rec.im_size
                self.mm_per_vox_x = self.mm_per_vox_y = self.fov_x / self.size_x
            elif self.psd_type == 'basic':
                # first 6 are ref scans, so ignore those. Also, two acquired timepoints are used
                # to generate each reconned time point.
                self.num_timepoints = (self._hdr.rec.npasses * self._hdr.rec.nechoes - 6) / 2
                self.num_echos = 1
            elif self.psd_type == 'muxepi':
                self.num_bands = int(self._hdr.rec.user6)
                self.num_mux_cal_cycle = int(self._hdr.rec.user7)
                self.band_spacing_mm = self._hdr.rec.user8
                # When ARC is used with mux, the number of acquired TRs is greater than what's Rxed.
                # ARC calibration uses multi-shot, so the additional TRs = num_bands*(ileaves-1)*num_mux_cal_cycle
                self.num_timepoints = self._hdr.rec.npasses + self.num_bands * (self._hdr.rec.ileaves-1) * self.num_mux_cal_cycle
                # The actual number of images returned by the mux recon is npasses - num_calibration_passes + num_mux_cal_cycle
                self.num_timepoints_available = self._hdr.rec.npasses - self.num_bands * self.num_mux_cal_cycle + self.num_mux_cal_cycle
                # TODO: adjust the image.tlhc... fields to match the correct geometry.
            elif self.psd_type == 'mrs':
                self._hdr.image.scanspacing = 0.
                self.mm_per_vox = [self._hdr.rec.roileny, self._hdr.rec.roilenx, self._hdr.rec.roilenz]
                image_tlhc = np.array((-self._hdr.rec.roilocx - self.mm_per_vox[0]/2.,
                                        self._hdr.rec.roilocy + self.mm_per_vox[1]/2.,
                                        self._hdr.rec.roilocz - self.mm_per_vox[1]/2.))
                image_trhc = image_tlhc - [self.mm_per_vox[0], 0., 0.]
                image_brhc = image_trhc + [0., self.mm_per_vox[1], 0.]
            # Tread carefully! Most of the stuff down here depends on various fields being corrected in the
            # sequence-specific set of hacks just above. So, move things with care!

            # Note: the following is true for single-shot planar acquisitions (EPI and 1-shot spiral).
            # For multishot sequences, we need to multiply by the # of shots. And for non-planar aquisitions,
            # we'd need to multiply by the # of phase encodes (accounting for any acceleration factors).
            # Even for planar sequences, this will be wrong (under-estimate) in case of cardiac-gating.
            self.prescribed_duration = self.num_timepoints * self.tr
            self.total_num_slices = self.num_slices * self.num_timepoints
            # The actual duration can only be computed after the data are loaded. Settled for rx duration for now.
            self.duration = self.prescribed_duration
            self.effective_echo_spacing = self._hdr.image.effechospace / 1e6
            self.phase_encode_undersample = 1. / self._hdr.rec.ileaves
            # TODO: Set this correctly! (it's in the dicom at (0x0043, 0x1083))
            self.slice_encode_undersample = 1.          # FIXME
            self.acquisition_matrix_x, self.acquisition_matrix_y = [self._hdr.rec.rc_xres, self._hdr.rec.rc_yres]
            # TODO: it looks like the pfile now has a 'grad_data' field!
            # Diffusion params
            self.dwi_numdirs = self._hdr.rec.numdifdirs
            # You might think that the b-value for diffusion scans would be stored in self._hdr.image.b_value.
            # But alas, this is GE. Apparently, that var stores the b-value of the just the first image, which is
            # usually a non-dwi. So, we had to modify the PSD and stick the b-value into an rhuser CV. Sigh.
            # NOTE: pre-dv24, the bvalue was stored in rec.user22.
            self.dwi_bvalue = self._hdr.rec.user1 if self.version == 24 else self._hdr.rec.user22
            self.is_dwi = True if self.dwi_numdirs >= 6 else False
            # if bit 4 of rhtype(int16) is set, then fractional NEX (i.e., partial ky acquisition) was used.
            self.partial_ky = self._hdr.rec.scan_type & np.uint16(16) > 0
            # was pepolar used to flip the phase encode direction?
            self.phase_encode_direction = 1 if np.bitwise_and(self._hdr.rec.dacq_ctrl,4)==4 else 0
            self.caipi = self._hdr.rec.user13   # true: CAIPIRINHA-type acquisition; false: Direct aliasing of simultaneous slices.
            self.cap_blip_start = self._hdr.rec.user14   # Starting index of the kz blips. 0~(mux-1) correspond to -kmax~kmax.
            self.cap_blip_inc = self._hdr.rec.user15   # Increment of the kz blip index for adjacent acquired ky lines.
            self.mica = self._hdr.rec.user17   # MICA bit-reverse?
            self.slice_duration = self.tr / self.num_slices
            lr_diff = image_trhc - image_tlhc
            si_diff = image_trhc - image_brhc
            if not np.all(lr_diff == 0) and not np.all(si_diff == 0):
                row_cosines =  lr_diff / np.sqrt(lr_diff.dot(lr_diff))
                col_cosines = -si_diff / np.sqrt(si_diff.dot(si_diff))
            else:
                row_cosines = np.array([1., 0, 0])
                col_cosines = np.array([0, -1., 0])
            self.slice_order = dcm.mr.generic_mr.SLICE_ORDER_UNKNOWN
            # FIXME: check that this is correct.
            if self._hdr.series.se_sortorder == 0:
                self.slice_order = dcm.mr.generic_mr.SLICE_ORDER_SEQ_INC
            elif self._hdr.series.se_sortorder == 1:
                self.slice_order = dcm.mr.generic_mr.SLICE_ORDER_ALT_INC
            # header geometry is LPS, but we need RAS, so negate R and A.
            slice_norm = np.array([-self._hdr.image.norm_R, -self._hdr.image.norm_A, self._hdr.image.norm_S])

            # This is either the first slice tlhc (image_tlhc) or the last slice tlhc. How to decide?
            # And is it related to wheather I have to negate the slice_norm?
            # Tuned this empirically by comparing spiral and EPI data with the same Rx.
            # Everything seems reasonable, except the test for axial orientation (start_ras==S|I).
            # I have no idea why I need that! But the flipping only seems necessary for axials, not
            # coronals or the few obliques I've tested.
            # FIXME: haven't tested sagittals!
            if (self._hdr.series.start_ras in 'SI' and self._hdr.series.start_loc > self._hdr.series.end_loc):
                self.reverse_slice_order = True
                slice_fov = np.abs(self._hdr.series.start_loc - self._hdr.series.end_loc)
                image_position = image_tlhc - slice_norm * slice_fov
                # FIXME: since we are reversing the slice order here, should we change the slice_order field below?
            else:
                image_position = image_tlhc
                self.reverse_slice_order = False

            # not sure why the following is needed.
            # TODO: * test non-slice-reversed coronals-- do they also need l/r flip?
            #       * test sagitals-- do they need any flipping?
            if (self._hdr.series.start_ras in 'AP' and self._hdr.series.start_loc > self._hdr.series.end_loc):
                slice_norm = -slice_norm
                self.flip_lr = True
            else:
                self.flip_lr = False

            if self.num_bands > 1:
                image_position = image_position - slice_norm * self.band_spacing_mm * (self.num_bands - 1.0) / 2.0

            # origin = image_position * np.array([-1, -1, 1])
            # Fix the half-voxel offset. Apparently, the p-file convention specifies coords at the
            # corner of a voxel. But DICOM/NIFTI convention is the voxel center. So offset by a half-voxel.
            origin = image_position + (row_cosines+col_cosines)*(np.array(self.mm_per_vox)/2)
            # The DICOM standard defines these two unit vectors in an LPS coordinate frame, but we'll
            # need RAS (+x is right, +y is anterior, +z is superior) for NIFTI. So, we compute them
            # such that self.row_cosines points to the right and self.col_cosines points up.
            row_cosines[0:2] = -row_cosines[0:2]
            col_cosines[0:2] = -col_cosines[0:2]
            if self.is_dwi and self.dwi_bvalue == 0:
                log.warning('the data appear to be diffusion-weighted, but image.b_value is 0! Setting it to 10.')
                # Set it to something other than 0 so non-dwi's can be distinguised from dwi's
                self.dwi_bvalue = 10.
            # The bvals/bvecs will get set later
            self.bvecs, self.bvals = (None, None)
            self.image_rotation = dcm.mr.generic_mr.compute_rotation(row_cosines, col_cosines, slice_norm)
            self.qto_xyz = dcm.mr.generic_mr.build_affine(self.image_rotation, self.mm_per_vox, origin)
            self.infer_psd_type()
            self.infer_scan_type()
            log.debug((self.psd_name, self.psd_type, self.scan_type))
            if self.psd_type == 'muxepi' and self.num_mux_cal_cycle < 2:
                if self.aux_file:
                    log.warning('muxepi without own calibration, will checking aux_file %s.' % self.aux_file)
                else:
                    log.warning('muxepi without own calibration. please provide an aux_file to load_data fxn.')
            self.full_parsed = True
            self.metadata_statue = 'complete'

    @property
    def canonical_filename(self):
        """Return the pfile name, without .7."""
        return self.pfilename

    @property
    def priority(self):
        """Return priority, 1 if can recon, -1 if cannot."""
        return int(bool(self.recon_func)) * 2 - 1  # 1 = can recon, -1 = cannot

    def get_bvecs_bvals(self, dirpath):
        """
        Parse tensor data from tensor file.

        Parameters
        ----------
        dirpath : str
            path to directory that contains tensor file.  This is usually the same directory that
            contains the P?????.7 file.

        Returns
        -------
        None : NoneType
            Set dataset.bvecs and dataset.bvals if tensor file is found.

        """
        tensor_name = '%s.7_tensor.dat' % self.pfilename  # pfilename set during _full_parse
        tensor_path = os.path.join(dirpath, tensor_name)

        if not os.path.exists(tensor_path):
            log.warning('tensor file %s not found' % tensor_path)
        else:
            log.warning('tensor file %s found' % tensor_path)
            with open(tensor_path) as fp:
                try:
                    uid = fp.readline().rstrip()
                    ndirs = int('0' + fp.readline.rstrip())
                except:
                    fp.seek(0, 0)
                    uid = None
                    ndirs = int('0' + fp.readline().rstrip())
                bvecs = np.fromfile(fp, sep=' ')

            if uid and uid != self.series_uid:  # uid provided does not match
                raise NIMSPFileError('tensor file UID does not match PFile UID!')
            if (ndirs or None) != self.dwi_numdirs or self.dwi_numdirs != bvecs.size / 3.:
                log.warning('tensor file numdirs does not match PFile header numdirs!')
                self.bvecs = None
                self.bvals = None
            else:
                num_nondwi = self.num_timepoints_available - self.dwi_numdirs
                bvals = np.concatenate((np.zeros(num_nondwi, dtype=float), np.tile(self.dwi_bvalue, self.dwi_numdirs)))
                bvecs = np.hstack((np.zeros((3, num_nondwi), dtype=float), bvecs.reshape(self.dwi_numdirs, 3).T))
                self.bvecs, self.bvals = dcm.mr.generic_mr.adjust_bvecs(bvecs, bvals, self.scanner_type, self.image_rotation)

    @property
    def recon_func(self):
        """Property that returns a member function that can then be executed."""
        if self.psd_type == 'spiral':
            return self.recon_spirec
        elif self.psd_type == 'muxepi':
            return self.recon_muxepi
        elif self.psd_type == 'mrs':
            return self.recon_mrs
        elif self.psd_type == 'hoshim':
            return self.recon_hoshim
        elif self.psd_type == 'basic':
            return self.recon_basic
        else:
            return None

    def do_recon(self, filepath=None, tempdir=None):
        """
        Run recon_func on filepath in the specified tempdir.

        Parameters
        ----------
        filepath : str
            path to pfile.7. input file must be pfile.7
        tempdir : str
            path to as base for temporary directory

        """
        pfilepath = filepath or self.filepath
        pfiledir = os.path.dirname(pfilepath)

        if self.is_dwi:
            self.get_bvecs_bvals(pfiledir)
        if self.recon_func:
            try:
                self.recon_func(pfilepath, tempdir)
            except Exception as e:
                log.debug('an error occured: pixel data could not be loaded from %s' % (self.filepath))
                self.data = None
                self.failure_reason = e

        # common stuff that can occur after the recon_func has been run, and data
        # loaded into self.data, if a recon func was successfull, self.data will be a dict, instead of None
        if self.data:
            for k in self.data.iterkeys():
                if self.reverse_slice_order:
                    self.data[k] = self.data[k][:,:,::-1,]
                if self.flip_lr:
                    self.data[k] = self.data[k][::-1,:,:,]

    def load_data(self, num_jobs=None, num_virtual_coils=None, tempdir=None, aux_file=None):
        """
        Load the data and run the appropriate reconstruction.

        Load data always works on the __init__ filepath.  it will determine if the file is a tgz, or not, and
        take the appropriate action to fully parse and prepare to reconstruct.

        Some parameters are repeated from __init__, to allow resetting those parameters at the time of data load time.

        Parameters
        ----------
        num_jobs : int
            override the number of jobs to use that was set during __init__
        num_virtual_coils : int
            override the number of virtual coils that was set during __init__
        tempdir : str
            override the temporary directory that was set during __init__
        aux_file : list
            override the list of potential aux files that was set during __init__

        """
        self.num_jobs = num_jobs or self.num_jobs
        self.num_vcoils = num_virtual_coils or self.num_vcoils
        self.aux_file = aux_file or self.aux_file
        self.tempdir = tempdir or self.tempdir

        if tarfile.is_tarfile(self.filepath):
            log.debug('loading data from tgz %s' % self.filepath)
            with tempfile.TemporaryDirectory(dir=self.tempdir) as temp_dirpath:
                log.debug('now working in temp_dirpath=%s' % temp_dirpath)
                with tarfile.open(self.filepath) as archive:
                    archive.extractall(path=temp_dirpath)

                temp_datadir = os.path.join(temp_dirpath, os.listdir(temp_dirpath)[0])    # tgz always has subdir that contains data
                for f in os.listdir(temp_datadir):
                    fpath = os.path.join(temp_datadir, f)
                    try:
                        self.version = get_version(fpath)
                    except Exception:
                        pass
                    else:
                        self._full_parse(fpath)  # provide input to parse
                        break
                self.do_recon(fpath, self.tempdir)

            log.debug('closing tempdir %s' % self.tempdir)
        else:
            log.debug('loading data from .7 %s' % self.filepath)
            if not self.full_parsed:
                self._full_parse()      # parse original input
            self.do_recon(self.filepath, self.tempdir)

    def load_imagedata_from_file(self, filepath):
        """
        Load raw image data from a file and do sanity checking on metadata values.

        Parameters
        ----------
        filepath : str
            path to *.mat, such as sl_001.mat

        Returns
        -------
        imagedata: np.array
            TODO: more details about np.array format?

        """
        # TODO confirm that the voxel reordering is necessary
        import scipy.io
        mat = scipy.io.loadmat(filepath)
        if 'd' in mat:
            sz = mat['d_size'].flatten().astype(int)
            slice_locs = mat['sl_loc'].flatten().astype(int) - 1
            imagedata = np.zeros(sz, mat['d'].dtype)
            raw = np.atleast_3d(mat['d'])
            if len(slice_locs)<raw.shape[2]:
                slice_locs = range(raw.shape[2])
                log.warning('Slice_locs is too short. Assuming slice_locs=[0,1,...,nslices]')
            imagedata[:,:,slice_locs,...] = raw[::-1,...]
        elif 'MIP_res' in mat:
            imagedata = np.atleast_3d(mat['MIP_res'])
            imagedata = imagedata.transpose((1,0,2,3))[::-1,::-1,:,:]
        if imagedata.ndim == 3:
            imagedata = imagedata.reshape(imagedata.shape + (1,))
        return imagedata

    def update_imagedata(self, imagedata, key=''):
        """
        Insert imagedata into self.data dictionary under the specified key.

        Update dataset metadata to be consistent with the imagedata.

        TODO: this assumes there is only a primary dataset. needs to rewritten for the
        case where there is primary and secondary data.

        Parameters
        ----------
        imagedata : array
            np array of voxel data
        key : str [default ''] (empty string is primary dataset)
            which key to use in data dict.

        """
        if imagedata.shape[0] != self.size[0] or imagedata.shape[1] != self.size[1]:
            log.warning('Image matrix discrepancy. Fixing the header, assuming imagedata is correct...')
            self.size = [imagedata.shape[0], imagedata.shape[1]]
            self.mm_per_vox_x = self.fov_x / self.size_y
            self.mm_per_vox_y = self.fov_y / self.size_y
        if imagedata.shape[2] != self.num_slices * self.num_bands:
            log.warning('Image slice count discrepancy. Fixing the header, assuming imagedata is correct...')
            self.num_slices = imagedata.shape[2]
        if imagedata.shape[3] != self.num_timepoints:
            log.warning('Image time frame discrepancy (header=%d, array=%d). Fixing the header, assuming imagedata is correct...'
                    % (self.num_timepoints, imagedata.shape[3]))
            self.num_timepoints = imagedata.shape[3]
        self.duration = self.num_timepoints * self.tr # FIXME: maybe need self.num_echos?
        self.data = {key: imagedata}

    def recon_hoshim(self, filepath, tempdir=None):
        log.debug('HOSHIM recon not implemented')
        self.is_non_image = True

    def recon_basic(self, filepath, tempdir=None):
        log.debug('BASIC recon not implemented')
        self.is_non_image = True

    def recon_spirec(self, filepath, tempdir=None):
        """
        Run spirec recon on the input filepath.

        Requires access to the spiral_recon submodule.

        Parameters
        ----------
        filepath : str
            path to P?????.7 to be reconstructed as spiral
        tempdir : str
            path to base of temporary directory

        """
        log.debug('SPIREC recon started')

        with tempfile.TemporaryDirectory(dir=tempdir) as temp_dirpath:
            log.debug('working in tempdir: %s' % temp_dirpath)
            if is_gzip(filepath):
                pfile_path = os.path.join(temp_dirpath, self.basename)
                with open(pfile_path, 'wb') as fd:
                    with gzip.open(self.filepath, 'rb') as gzfile:
                        fd.writelines(gzfile)
            else:
                pfile_path = filepath
            recon_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'spiral_recon'))
            basepath = os.path.join(temp_dirpath, 'recon')
            spirec_path = os.path.join(recon_path, 'spirec')
            cmd = '%s -l --rotate -90 --magfile --savefmap2 --b0navigator -r %s -t %s' % (spirec_path, pfile_path, basepath)
            log.debug(cmd)
            subprocess.call(shlex.split(cmd), cwd=temp_dirpath, stdout=open('/dev/null', 'w'))  # run spirec to generate .mag and fieldmap files
            log.debug(os.listdir(temp_dirpath))
            self.data = {'': np.fromfile(file=basepath+'.mag_float', dtype=np.float32).reshape([self.size[0],self.size[1],self.num_timepoints,self.num_echos,self.num_slices],order='F').transpose((0,1,4,2,3))}
            if os.path.exists(basepath+'.B0freq2') and os.path.getsize(basepath+'.B0freq2')>0:
                self.data['fieldmap'] = np.fromfile(file=basepath+'.B0freq2', dtype=np.float32).reshape([self.size[0],self.size[1],self.num_echos,self.num_slices],order='F').transpose((0,1,3,2))

            if self.num_echos == 2:
                # FIXME: Do a more robust test for spiralio!
                # FIXME: should do a quick motion correction here
                # Assume spiralio, so do a weighted average of the two echos.
                w_in = np.mean(self.data[''][:,:,:,:,0], 3)
                w_out = np.mean(self.data[''][:,:,:,:,1], 3)
                # self.data['w_in'] = w_in    # uncomment to save spiral_in
                # self.data['w_out'] = w_out  # uncomment to save spiral_out
                inout_sum = w_in + w_out
                w_in = w_in / inout_sum
                w_out = w_out / inout_sum
                avg = np.zeros(self.data[''].shape[0:4])
                for tp in range(self.data[''].shape[3]):
                    avg[:,:,:,tp] = w_in*self.data[''][:,:,:,tp,0] + w_out*self.data[''][:,:,:,tp,1]
                self.data[''] = avg

    def prep_convert(self):
        """
        Return criteria for identifying the additional aux file necessary to perform recon.

        This gives the processor a way of checking if additional files are necessary.  The processor has some knowledge about
        muxepi and how to identify a valid aux file.

        The main business logic of identifying candidate aux files has been moved into processor.py.
        """
        # FIXME: the following is a hack to get mux_epi2 SE-IR scans to recon properly. There *is* a more generic solution...
        # some mux scans don't have their own calibration, and require fetching calibration a scan with the same psd_name.
        if self.psd_type=='muxepi' and (self.num_mux_cal_cycle<2 or (self.psd_name=='mux_epi2' and self.ti>0)):
            aux_data = { 'psd': self.psd_name }
        else:
            aux_data = None
        return aux_data

    def recon_muxepi(self, filepath, tempdir=None, timepoints=[], octave_bin='octave'):
        """
        Do mux_epi image reconstruction and populate self.data.

        Always involves a tempdir.  If input is a pfile.tgz, the tempdir was created during
        unpacking will be re-used during recon. If input is a pfile.7, then the temp during
        will be created during this recon.

        Parameters
        ----------
        filepath : str
            input P?????.7 to be reconstructed
        tempdir : str
            path to base temporary directory
        timepoints : list
            if list is non-empty, restrict reconstruction to the listed timepoints
        octave_bin : str [default 'octave']
            path to octave executable, default assumes octave can be found in $PATH

        """
        start_sec = time.time()
        log.debug('MUXEPI recon of %s started at %d' % (filepath, start_sec))
        basename = os.path.basename(filepath)
        dirpath = os.path.dirname(filepath)

        fermi_filt = 1
        sense_recon = 0
        if self.recon_type == None:
            # set the recon type automatically, scans with mux>1, arc>1, caipi
            if self.is_dwi and self.num_bands>1 and self.phase_encode_undersample<1. and self.caipi:
                sense_recon = 1
        elif self.recon_type == 'sense':
            sense_recon = 1

        log.debug(self.aux_file)
        with tempfile.TemporaryDirectory(dir=tempdir) as temp_dirpath:
            # identification and extraction of the aux file has been moved into tempdir context manager
            # to allow extraction of the valid vrgf and ref within tempdir.
            # cal file must be either '', or dir/basename that is similar between ref.dat and vrgf.dat
            cal_file = ''  # XXX matlab/octave, base path determine names of ref.dat and vrgf.dat
            ref_file = os.path.join(dirpath, basename + '_ref.dat')
            vrgf_file = os.path.join(dirpath, basename + '_vrgf.dat')
            log.debug('checking num_mux_cal_cycle: %d' % self.num_mux_cal_cycle)
            # matlab code checks for num_mux_cal_cycles, so why don't we check the same thing
            # even scans without num_mux_cal_cycles will still have ref/vrgf files that exceed 64 bytes
            if self.num_mux_cal_cycle < 2:
                log.debug('num_mux_cal_cycle: %d. looking for calibration from a different acq.' % self.num_mux_cal_cycle)
                if self.aux_file:
                    with tarfile.open(self.aux_file) as aux_archive:
                        log.debug('inspecting aux file: %s' % self.aux_file)
                        aux_archive.extractall(path=temp_dirpath)
                    aux_subdir = os.listdir(temp_dirpath)[0]
                    aux_datadir = os.path.join(temp_dirpath, aux_subdir)
                    for f in os.listdir(aux_datadir):
                        if f.endswith('_ref.dat'):
                            cal_ref_file = os.path.join(aux_datadir, f)
                            log.debug('_ref.dat found, %s' % cal_ref_file)
                        elif f.endswith('_vrgf.dat'):
                            cal_vrgf_file = os.path.join(aux_datadir, f)
                            log.debug('_vrgf.dat found, %s' % cal_vrgf_file)
                        else:
                            pass
                if cal_ref_file and cal_vrgf_file:
                    if cal_ref_file.rsplit('_', 1)[0] == cal_vrgf_file.rsplit('_', 1)[0]:
                        cal_file = cal_ref_file.rsplit('_', 1)[0]
                    log.info('ref/vrgf.dat not found-- using calibration ref/vrgf from %s' % self.aux_file)
                else:
                    raise NIMSPFileError('ref.dat/vrgf.dat not found')

            # run the actual recon, spawning subprocess until all slices have been spawned.
            log.info('Running %d v-coil mux recon on %s in tempdir %s with %d jobs (sense=%d, fermi=%d, notch=%f).'
                    % (self.num_vcoils, filepath, temp_dirpath, self.num_jobs, sense_recon, fermi_filt, self.notch_thresh))
            if cal_file!='':
                log.info('Using calibration file: %s' % cal_file)
            recon_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mux_epi_recon'))
            outname = os.path.join(temp_dirpath, 'sl')

            mux_recon_jobs = []
            slice_num = 0
            while slice_num < self.num_slices:
                num_running_jobs = sum([job.poll()==None for job in mux_recon_jobs])
                if num_running_jobs < self.num_jobs:
                    # Recon each slice separately. Note the slice_num+1 to deal with matlab's 1-indexing.
                    # Use 'str' on timepoints so that an empty array will produce '[]'
                    cmd = ('%s --no-window-system -p %s --eval \'mux_epi_main("%s", "%s_%03d.mat", "%s", %d, %s, %d, 0, %s, %s, %s);\''
                        % (octave_bin, recon_path, filepath, outname, slice_num, cal_file, slice_num + 1, str(timepoints), self.num_vcoils, str(sense_recon), str(fermi_filt), str(self.notch_thresh)))
                    log.debug(cmd)
                    mux_recon_jobs.append(subprocess.Popen(args=shlex.split(cmd), stdout=open('/dev/null', 'w')))
                    slice_num += 1
                else:
                    time.sleep(1.)

            # Now wait for all the jobs to finish
            for job in mux_recon_jobs:
                job.wait()

            # Load the first slice to initialize the image array
            img = self.load_imagedata_from_file("%s_%03d.mat" % (outname, 0))
            for slice_num in range(1, self.num_slices):
                new_img = self.load_imagedata_from_file("%s_%03d.mat" % (outname, slice_num))
                # Allow for a partial last timepoint. This sometimes happens when the user aborts.
                t = min(img.shape[-1], new_img.shape[-1])
                img[...,0:t] += new_img[...,0:t]

            img = img.astype(np.float32)
            self.update_imagedata(img)
            elapsed = time.time() - start_sec
            log.info('Mux recon of %s with %d v-coils finished in %0.2f minutes using %d jobs.'
                        % (self.filepath, self.num_vcoils,  elapsed/60., min(self.num_jobs, self.num_slices)))

    def recon_mrs(self, filepath, tempdir=None):
        """
        Load raw spectro data.

        Currently just loads raw spectro data into self.data dictionary, to
        prepare for writing to nifti.

        Parameters
        ----------
        filepath : str
            path to input file, can be .7, .7.gz.  cannot be 7.tgz.

        Returns
        -------
        None : NoneType
            loads spectro data into self.data['']

        """
        log.debug('MRS recon started')
        self.data = {'': self.get_rawdata(filepath).transpose([0,5,3,1,2,4])}

    def get_rawdata(self, filepath, slices=None, passes=None, coils=None, echos=None, frames=None):
        """
        Read and return a chunck of data from the p-file.

        Specify the slices, timepoints, coils, and echos that you want.
        None means you get all of them. The default of all Nones will
        return all data.
        (based on https://github.com/cni/MRS/blob/master/MRS/files.py)

        """
        n_frames = self._hdr.rec.nframes + self._hdr.rec.hnover
        n_echos = self._hdr.rec.nechoes
        n_slices = self._hdr.rec.nslices / self._hdr.rec.npasses
        n_coils = self.num_receivers
        n_passes = self._hdr.rec.npasses
        frame_sz = self._hdr.rec.frame_size

        if passes == None: passes = range(n_passes)
        if coils == None: coils = range(n_coils)
        if slices == None: slices = range(n_slices)
        if echos == None: echos = range(n_echos)
        if frames == None: frames = range(n_frames)

        # Size (in bytes) of each sample:
        ptsize = self._hdr.rec.point_size
        data_type = [np.int16, np.int32][ptsize/2 - 1]

        # This is double the size as above, because the data is complex:
        frame_bytes = 2 * ptsize * frame_sz

        echosz = frame_bytes * (1 + n_frames)
        slicesz = echosz * n_echos
        coilsz = slicesz * n_slices
        passsz = coilsz * n_coils

        # Byte-offset to get to the data:
        offset = self._hdr.rec.off_data
        fp = gzip.open(filepath, 'rb') if is_gzip(filepath) else open(filepath, 'rb')
        data = np.zeros((frame_sz, len(frames), len(echos), len(slices), len(coils), len(passes)), dtype=np.complex)
        for pi,passidx in enumerate(passes):
            for ci,coilidx in enumerate(coils):
                for si,sliceidx in enumerate(slices):
                    for ei,echoidx in enumerate(echos):
                        for fi,frameidx in enumerate(frames):
                            fp.seek(passidx*passsz + coilidx*coilsz + sliceidx*slicesz + echoidx*echosz + (frameidx+1)*frame_bytes + offset)
                            # Unfortunately, numpy fromfile doesn't like gzip file objects. But we can
                            # safely load each chunk into RAM, since frame_sz is never very big.
                            #dr = np.fromfile(fp, data_type, frame_sz * 2)
                            dr = np.fromstring(fp.read(frame_bytes), data_type)
                            dr = np.reshape(dr, (-1, 2)).T
                            data[:, fi, ei, si, ci, pi] = dr[0] + dr[1]*1j
        fp.close()
        return data
