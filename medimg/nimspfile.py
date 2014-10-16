#!/usr/bin/env python
#
# @author:  Gunnar Schaefer
#           Kevin S Hahn

"""
nimsdata.nimspfile
==================

This module provides functions, classes and errors for reading for
minimally parsing a pfile.  Additional modules can be used to
full parse the pfiles.

NIMSPfile provides parsing and identification of GE pfiles.

NIMSPfile requires the pfile module to completely read and perform conversions.

"""

import os
import bson
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


def uncompress(filepath, tempdir):
    """
    Uncompress the contents of a tgz, or gz, into a specified tempdir.

    This function might be totally useless.

    """
    newpath = os.path.join(tempdir, os.path.basename(filepath)[:-3])
    # The following with pigz is ~4x faster than the python code above (with gzip, it's about 2.5x faster)
    if os.path.isfile('/usr/bin/pigz'):
        subprocess.call('pigz -d -c %s > %s' % (filepath, newpath), shell=True)
    elif os.path.isfile('/usr/bin/gzip') or os.path.isfile('/bin/gzip'):
        subprocess.call('gzip -d -c %s > %s' % (filepath, newpath), shell=True)
    else:
        with open(newpath, 'wb') as fd:
            with gzip.open(filepath, 'rb') as gzfile:
                fd.writelines(gzfile)
    return newpath


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
    Read pfile data and/or header.

    This class reads the data and/or header from a pfile, runs k-space reconstruction,
    and generates a NIfTI object, including header information.

    NIMSPFile object can handle several different input types
    .tgz of directory containing Pfile, auxillary files such as ref.dat, vrgf.dat and tensor.dat
    .gz of single pfile

    tgz cannot be "full parsed".  setting full_parse=True, with an input tgz, will raise an exception.

    the old way is no longer supported.
    Pfile.gz + _Pfile.7_ref.dat + _PFile.7._vrgf.dat _PFile.7_refscan.7, _P02048.7_param.dat

    NOTE: does nimspfile.NIMSPFile need to be backward compatible with nimsraw.NIMSPFile?

    Example:
        import pfile
        pf = pfile.PFile(filename='P56832.7')
        pf.to_nii(outbase='P56832.7')

    """
    domain = u'mr'
    filetype = u'pfile'
    parse_priority = 5
    state = ['orig']

    def __init__(self, filepath, load_data=False, full_parse=False, num_jobs=8, num_virtual_coils=16, notch_thresh=0, recon_type=None):
        super(NIMSPFile, self).__init__(filepath, load_data)
        self.full_parse = False  # self.full_parse = False indicates file is not fully parsed
        log.debug('parsing pfile %s, full_parse = %s, load_data = %s' % (self.filepath, full_parse, load_data))
        self.dirpath = os.path.dirname(self.filepath)
        self.filename = os.path.basename(self.filepath)
        self.basename, _ = os.path.splitext(self.filename)  # discard the file extension
        self.is_localizer = None
        self.num_vcoils = num_virtual_coils
        self.notch_thresh = notch_thresh
        self.recon_type = recon_type

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
                        self.epoch_id = _hdr.get('epoch')
                        self.timestamp = _hdr.get('timestamp')
                        self.group_name = _hdr.get('group')
                        self.experiment_name = _hdr.get('experiment')
                        break
                else:
                    raise NIMSPFileError('no json file with header section found. bailing', log_level=logging.WARNING)
        else:  # .7 or .7.gz, doing it old world style
            try:
                self.version = get_version(filepath)
                self._full_parse(filepath) if full_parse else self._min_parse(filepath)  # full_parse arg indicates run full_parse
            except:
                raise NIMSPFileError('not a PFile')

        self.metadata_status = 'pending'

        if load_data:  # load_data arg indicates run load_data
            self.load_data(num_jobs)  # load_data checks self.full_parse to see if full_parse is needed

    def _min_parse(self, filepath):
        """
        Parse a minimum required set of metadata from a Pfile fileobj.

        Can only parse Pxxxxx.7 and Pxxxxx.7.gz files, not tgz.  Does not load a self._hdr object.

        Parameters
        ----------
        filepath : str
            path of file to be parsed

        """
        fileobj = gzip.open(filepath, 'rb') if is_gzip(self.filepath) else open(filepath, 'rb')
        log.debug('min_parse')

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
            fileobj.seek(148972); self.psd_name = os.path.basename(struct.unpack("33s", fileobj.read(struct.calcsize("33s")))[0]).split('\0', 1)[0]
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
            fileobj.seek(65374); self.psd_name = os.path.basename(struct.unpack("33s", flleobj.read(struct.calcsize("33s")))[0]).split('\0', 1)[0]

        if self.im_datetime > 0:
            self.timestamp = datetime.datetime.utcfromtimestamp(self.im_datetime)
        else:
            month, day, year = map(int, self.scan_date.split('\0', 1)[0].split('/'))
            hour, minute = map(int, self.scan_time.split('\0', 1)[0].split(':'))
            self.timestamp = datetime.datetime(year + 1900, month, day, hour, minute)  # GE's epoch begins in 1900

        dcm.mr.ge.infer_psd_type(self)  # kinda awkward
        if self.psd_type == 'spiral':
            self.num_timepoints = int(self.rec_user0)
        elif self.psd_type == 'basic':
            self.num_timepoints = (self.num_timepoints * self.num_echos - 6) / 2
        elif self.psd_type == 'muxepi':
            self.num_timepoints = self.num_timepoints + int(self.rec_user6) * self.ileaves * (int(self.rec_user7) - 1)
        self.prescribed_duration = self.num_timepoints * self.tr
        self.subj_code, self.group_name, self.experiment_name = medimg.parse_patient_id(self.patient_id, 'ex' + self.exam_no)

    def _full_parse(self, filepath):
        """
        Work on .7.gz and .7.

        Requires the pfile submodule.
        """
        log.debug('full_parse')
        if tarfile.is_tarfile(filepath):
            raise NIMSPFileError('Cannot Full Parse a tgz, need to load all data first. bailing.', log_level=logging.WARNING)
        try:
            pfile = getattr(__import__('pfile.pfile%d' % self.version, globals()), 'pfile%d' % self.version)
        except ImportError as e:
            raise ImportError('%s\nNo Valid PFile parser for v%d\n' % (str(e), self.version))

        with gzip.open(filepath, 'rb') if is_gzip(filepath) else open(filepath, 'rb') as fileobj:
            self._hdr = pfile.POOL_HEADER(fileobj)
            if not self._hdr:
                raise NIMSPFileError('no pfile read', log_level=logging.WARNING)

            self.data = {}      # empty dict to start
            # self.data = None
            # self.fm_data = None

            self.exam_no = self._hdr.exam.ex_no
            self.exam_uid = unpack_uid(self._hdr.exam.study_uid)
            self.patient_id = self._hdr.exam.patidff.split('\0', 1)[0]
            self.series_no = self._hdr.series.se_no
            self.series_desc = self._hdr.series.se_desc.split('\0', 1)[0]
            self.series_uid = unpack_uid(self._hdr.series.series_uid)
            self.acq_no = self._hdr.image.scanactno
            self.subj_code, self.group_name, self.experiment_name = medimg.parse_patient_id(self.patient_id, 'ex' + str(self.exam_no))

            self.psd_name = os.path.basename(self._hdr.image.psdname.partition('\x00')[0])
            self.scan_type = self._hdr.image.psd_iname.split('\0', 1)[0]
            self.pfilename = 'P%05d' % self._hdr.rec.run_int
            self.subj_firstname, self.subj_lastname = medimg.parse_patient_name(self._hdr.exam.patnameff.split('\0', 1)[0])
            self.subj_dob = medimg.parse_patient_dob(self._hdr.exam.dateofbirth.split('\0', 1)[0])
            self.subj_sex = ('male', 'female')[self._hdr.exam.patsex-1] if self._hdr.exam.patsex in [1, 2] else None
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
            self.scan_type = self._hdr.image.psd_iname.split('\0', 1)[0]
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
            dcm.mr.ge.infer_psd_type(self)
            # self.psd_type = dcm.mr.ge.infer_psd_type(self.psd_name)
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
                # self.mm_per_vox[0:2] = [self.fov[0] / self.size[0]] * 2
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
            self.dwi_bvalue = self._hdr.rec.user22 if self.version == 24 else self._hdr.rec.user1
            self.is_dwi = True if self.dwi_numdirs >= 6 else False
            # if bit 4 of rhtype(int16) is set, then fractional NEX (i.e., partial ky acquisition) was used.
            self.partial_ky = self._hdr.rec.scan_type & np.uint16(16) > 0
            # was pepolar used to flip the phase encode direction?
            self.phase_encode_direction = -1 if np.bitwise_and(self._hdr.rec.dacq_ctrl,4)==4 else 1
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
            # self.scan_type = dcm.mr.generic_mr.infer_scan_type(self)
            dcm.mr.generic_mr.infer_scan_type(self)
            log.debug(self.psd_name)
            self.aux_files = None

            self.full_parse = True

    @property
    def canonical_filename(self):
        return self.pfilename

    @property
    def priority(self):
        return int(bool(self.recon_func)) * 2 - 1   # return 1 if we can recon, else -1

    def get_bvecs_bvals(self):
        """
        Retrieve tensor file from within tgz archive.

        Tensor file, like other supplementary files, will be located within a tgz

        """
        # TODO: have this identify tensor files from WITHIn a tgz
        # and then save results to self.bvecs and self.bvals
        log.debug('collecting tensor info')

        tensor_file = os.path.join(self.dirpath, '_' + self.basename + '_tensor.dat')
        with tarfile.open(self.filepath) as archive:
            try:
                fp = archive.extract_file(tensor_file)      # identify the tensor file by name
            except KeyError:
                log.warning('%s; tensor file not found')
            else:
                try:
                    uid = fp.readline().rstrip()
                    ndirs = int('0' + fp.readline().rstrip())
                except:
                    fp.seek(0, 0)
                    uid = None
                    ndirs = int('0' + fp.readline().rstrip())
                bvecs = np.fromfile(fp, sep=' ')

            if uid and uid != self._hdr.series.series_uid:  # if uid is provided, and does not match
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
        log.debug(self.psd_type)
        if self.psd_type == 'spiral':
            return self.recon_spirec
        elif self.psd_type == 'muxepi':
            return self.recon_mux_epi
        elif self.psd_type == 'mrs':
            return self.recon_mrs
        elif self.psd_type == 'hoshim':
            return self.recon_hoshim
        elif self.psd_type == 'basic':
            return self.recon_basic
        else:
            return None

    def prep_convert(self):
        # FIXME: the following is a hack to get mux_epi2 SE-IR scans to recon properly. There *is* a more generic solution...
        if self.psd_type == 'muxepi' and (self.num_mux_cal_cycle < 2 or (self.psd_name == 'mux_epi2' and self.ti > 0)):
            # Mux scan without internal calibration-- request other mux scans be handed to convert
            # to see if we can find a suitable calibration scan.
            aux_data = {'psd': self.psd_name}
        else:
            aux_data = None
        return aux_data

    def load_data(self, num_jobs=8, tempdir=None, num_virtual_coils=16):
        self.num_jobs = num_jobs
        self.num_vcoils = 16
        # first need to know what type of file it is...
        if tarfile.is_tarfile(self.filepath):
            log.debug('load_data, tgz')
            # this tempdir will be used both for untaring the tgz, and
            # for placing the temporary sl000.mat intermediate files
            with tempfile.TemporaryDirectory(dir=tempdir) as temp_dirpath:
                # extract everything into tarfile
                with tarfile.open(self.filepath) as archive:
                    archive.extractall(path=temp_dirpath)

                dirpath = os.listdir(temp_dirpath)[0]
                files = os.listdir(os.path.join(temp_dirpath, dirpath))
                for f in files:
                    fpath = os.path.join(temp_dirpath, dirpath, f)
                    try:
                        self.version = get_version(fpath)
                    except Exception as e:
                        log.debug(e)
                    else:
                        self._full_parse(fpath)
                        self.pfile_path = fpath
                        break
                else:
                    raise NIMSPFileError('no pfile found?', log_level=logging.WARNING)

                # find calibration files, if any
                self.cal_ref_file = os.path.join(temp_dirpath, dirpath, self.pfilename + '.7_ref.dat')
                self.cal_vrgf_file = os.path.join(temp_dirpath, dirpath, self.pfilename + '.7_vrgf.dat')
                if not os.path.exists(self.cal_ref_file):
                    self.cal_ref_file = ''           # set to empty, bc bob does it
                if not os.path.exists(self.cal_vrgf_file):
                    self.cal_vrgf_file = ''          # set to empty, bc bob does it

                # find dti tensor files, if any
                if self.is_dwi:
                    self.get_bvecs_bvals()      # needs to know where to look

                if not self.data:  # data starts as empty dict
                    if self.recon_func:
                        log.debug('running recon')
                        self.recon_func(fpath, temp_dirpath, num_jobs)
                    else:
                        raise NIMSPFileError('Recon not implemented for this type of data')

                log.debug('closing tempdir %s' % temp_dirpath)
        else:
            log.debug('load_data, .7 or .7.gz')
            # handle .7 or .7.gz
            if not self.full_parse:
                self._full_parse(self.filepath)
            # handle lone P file

            if not self.data:
                if self.recon_func:
                    log.debug('running recon')
                    self.recon_func(self.filepath, tempdir=tempdir, num_jobs=num_jobs)
                else:
                    raise NIMSPFileError('Recon not implemented for this type of data')

        # common ground stuff
        # in some cases, more than one file is supposed to be output.
        # how to deal wih datasets that have both self.data and self.fm_data?
        # nimsnifti can only write a single file at a time.
        # TODO: pfile needs to be updated to return data in dictionary
        if self.data.get('') is not None:   # data dictionary is not empty
            if self.reverse_slice_order:
                self.data[''] = self.data[''][:,:,::-1,]
                if self.data.get('fieldmap') is not None:
                    self.data['fieldmap'] = self.data['fieldmap'][:,:,::-1,]
            if self.flip_lr:
                self.data[''] = self.data[''][::-1,:,:,]
                if self.data.get('fieldmap') is not None:
                    self.data['fieldmap'] = self.data['fieldmap'][::-1,:,:,]
            if self.psd_type == 'spiral' and self.num_echos == 2:
                # Uncomment to save spiral in/out
                # nimsnifti.NIMSNifti.write(self, self.imagedata[:,:,:,:,0], outbase + '_in')
                # nimsnifti.NIMSNifti.write(self, self.imagedata[:,:,:,:,1], outbase + '_out')
                # FIXME: Do a more robust test for spiralio!
                # Assume spiralio, so do a weighted average of the two echos.
                # FIXME: should do a quick motion correction here
                w_in = np.mean(self.data[''][:,:,:,:,0], 3)
                w_out = np.mean(self.data[''][:,:,:,:,1], 3)
                inout_sum = w_in + w_out
                w_in = w_in / inout_sum
                w_out = w_out / inout_sum
                avg = np.zeros(self.data[''].shape[0:4])
                for tp in range(self.data[''].shape[3]):
                    avg[:,:,:,tp] = w_in*self.data[''][:,:,:,tp,0] + w_out*self.data[''][:,:,:,tp,1]
                self.data[''] = avg
                log.debug('writing spiral IO')
            else:
                log.debug('writing spiral')

            # field map data is automatically written if it exists in self.data dictionary
            # if self.data.get('fieldmap') is not None:
            #     log.debug('writing fieldmaps')   # logging message is performed by nimsnifti

    def load_imagedata_from_file(self, filepath):
        """Load raw image data from a file and do some sanity checking on num slices, matrix size, etc."""
        # TODO: confirm that the voxel reordering is necessary. Maybe lean on the recon folks to standardize their voxel order?
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

    def update_imagedata(self, imagedata):
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
        self.data = {'': imagedata}

    def recon_hoshim(self, filepath, tempdir, num_jobs):
        log.debug('Cannot recon HO SHIM data')

    def recon_basic(self, filepath, tempdir, num_jobs):
        log.debug('Cannot recon BASIC data')

    def recon_spirec(self, tempdir, num_jobs):
        """Do spiral image reconstruction and populate self.imagedata."""
        log.debug('spiral recon')
        with tempfile.TemporaryDirectory(dir=tempdir) as temp_dirpath:
            if self.compressed:
                pfile_path = os.path.join(temp_dirpath, self.basename)
                with open(pfile_path, 'wb') as fd:
                    with gzip.open(self.filepath, 'rb') as gzfile:
                        fd.writelines(gzfile)
            else:
                pfile_path = self.filepath
            basepath = os.path.join(temp_dirpath, 'recon')
            recon_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'spiral_recon/spirec'))
            cmd = recon_path + ' -l --rotate -90 --magfile --savefmap2 --b0navigator -r %s -t %s' % (pfile_path, 'recon')
            log.debug(cmd)
            subprocess.call(shlex.split(cmd), cwd=temp_dirpath, stdout=open('/dev/null', 'w'))  # run spirec to generate .mag and fieldmap files

            self.imagedata = np.fromfile(file=basepath+'.mag_float', dtype=np.float32).reshape([self.size[0],self.size[1],self.num_timepoints,self.num_echos,self.num_slices],order='F').transpose((0,1,4,2,3))
            if os.path.exists(basepath+'.B0freq2') and os.path.getsize(basepath+'.B0freq2')>0:
                self.fm_data = np.fromfile(file=basepath+'.B0freq2', dtype=np.float32).reshape([self.size[0],self.size[1],self.num_echos,self.num_slices],order='F').transpose((0,1,3,2))

    # TODO: this function needs to be re-assessed.  Is it even useful?  Pfiles that have associated aux files
    # will always have their aux files inside the tgz with the Pfile....
    def find_mux_cal_file(self):
        log.debug('find mux cal')
        cal_file = []
        if self.aux_files!=None and len(self.aux_files)>0 and self.aux_files[0]!=None:
            if self.num_mux_cal_cycle>=2:
                candidates = [pf for pf in [(NIMSPFile(f),f) for f in self.aux_files] if pf[0].num_bands==1]
            else:
                candidates = [pf for pf in [(NIMSPFile(f),f) for f in self.aux_files] if pf[0].num_mux_cal_cycle>=2]
            if len(candidates)==1:
                cal_file = candidates[0][1].encode()
            elif len(candidates)>1:
                series_num_diff = np.array([c[0].series_no for c in candidates]) - self.series_no
                closest = np.min(np.abs(series_num_diff))==np.abs(series_num_diff)
                # there may be more than one. We prefer the prior scan:
                closest = np.where(np.min(series_num_diff[closest])==series_num_diff)[0][0]
                cal_file = candidates[closest][1].encode()
        if len(cal_file)>0:
            cal_compressed = is_compressed(cal_file)
            cal_basename = cal_file[:-3] if cal_compressed else cal_file
            cal_ref_file  = os.path.join(os.path.dirname(cal_basename), '_'+os.path.basename(cal_basename)+'.7_ref.dat')
            cal_vrgf_file = os.path.join(os.path.dirname(cal_basename), '_'+os.path.basename(cal_basename)+'.7_vrgf.dat')
        else:
            cal_compressed = False
            cal_ref_file = ''
            cal_vrgf_file = ''
        # Make sure we return an empty string when none is found.
        if not cal_file:
            cal_file = ''
        return cal_file,cal_ref_file,cal_vrgf_file,cal_compressed

    def recon_mux_epi(self, filepath, tempdir, num_jobs, timepoints=[], octave_bin='octave'):
        # DOES NOT WORK FOR MICA scans.  Currently, the MICA option causes use of rand_lcg_gcc.m, which relies
        # on a matlab function to work.  The only way to side-step this.
        log.debug('mux epi recon')
        start_sec = time.time()
        log.debug(start_sec)
        """Do mux_epi image reconstruction and populate self.data."""
        ref_file  = os.path.join(self.dirpath, '_'+self.basename+'_ref.dat')
        vrgf_file = os.path.join(self.dirpath, '_'+self.basename+'_vrgf.dat')
        # See if external calibration data files are needed:
        cal_file,cal_ref_file,cal_vrgf_file,cal_compressed = self.find_mux_cal_file()
        # The dat files might be missing or empty if the vendor recon was disabled. If so, try to use the cal dat file.
        # FIXME: if the p-file is not compressed, the cal dat file will not be used! We should refactor the recon
        # code so that the dat files are always explicitly specified.
        # if not os.path.isfile(ref_file) or os.path.getsize(ref_file)<64:
        #     if cal_ref_file:
        #         ref_file = cal_ref_file
        #     else:
        #         raise NIMSPFileError('ref.dat file not found')
        # if not os.path.isfile(vrgf_file) or os.path.getsize(vrgf_file)<64:
        #     if cal_vrgf_file:
        #         vrgf_file = cal_vrgf_file
        #     else:
        #         raise NIMSPFileError('vrgf.dat file not found')
        ref_file = self.cal_ref_file
        vrgf_file = self.cal_vrgf_file

        if self.recon_type == None:
            # set the recon type automatically
            # scans with mux>1, arc>1, caipi
            if self.is_dwi and self.num_bands>1 and self.phase_encode_undersample<1. and self.caipi:
                sense_recon = 1
            else:
                sense_recon = 0
        elif self.recon_type == 'sense':
            sense_recon = 1
        else:
            sense_recon = 0

        fermi_filt = 1

        # already have a tempdir open...
        # with tempfile.TemporaryDirectory(dir=tempdir) as temp_dirpath:
        if True:
            pfile_path = filepath
            temp_dirpath = tempdir
            log.info('Running %d v-coil mux recon on %s in tempdir %s with %d jobs (sense=%d, fermi=%d, notch=%f).'
                    % (self.num_vcoils, filepath, tempdir, num_jobs, sense_recon, fermi_filt, self.notch_thresh))
            if cal_file!='':
                log.info('Using calibration file: %s.' % cal_file)
            recon_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'mux_epi_recon'))
            outname = os.path.join(temp_dirpath, 'sl')

            # Spawn the desired number of subprocesses until all slices have been spawned
            mux_recon_jobs = []
            slice_num = 0
            while slice_num < self.num_slices:
                num_running_jobs = sum([job.poll()==None for job in mux_recon_jobs])
                if num_running_jobs < num_jobs:
                    # Recon each slice separately. Note the slice_num+1 to deal with matlab's 1-indexing.
                    # Use 'str' on timepoints so that an empty array will produce '[]'
                    cmd = ('%s --no-window-system -p %s --eval \'mux_epi_main("%s", "%s_%03d.mat", "%s", %d, %s, %d, 0, %s, %s, %s);\''
                        % (octave_bin, recon_path, pfile_path, outname, slice_num, cal_file, slice_num + 1, str(timepoints), self.num_vcoils, str(sense_recon), str(fermi_filt), str(self.notch_thresh)))
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
                      % (self.filepath, self.num_vcoils,  elapsed/60., min(num_jobs, self.num_slices)))

    def recon_mrs(self, filepath, tempdir, num_jobs):
        """Currently just loads raw spectro data into self.imagedata so that we can save it in a nifti."""
        log.debug('mrs recon')
        # Reorder the data to be in [frame, num_frames, slices, passes (repeats), echos, coils]
        # This roughly complies with the nifti standard of x,y,z,time,[then whatever].
        # Note that the "frame" is the line of k-space and thus the FID timeseries.
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
