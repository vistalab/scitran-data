# @author:  Reno Bowen
#           Gunnar Schaefer
#           Bob Dougherty
#           Kevin S Hahn


import os
import abc
import dicom
import logging
import datetime
import dcmstack
import cStringIO
import dcmstack.extract

import nimspng
import nimsnifti
import nimsmrdata

log = logging.getLogger('nimsdicom')

dicom.config.enforce_valid_values = False


class NIMSStack(dcmstack.DicomStack):

    """
    Does not check for ImageOrietnationPatient consistency, affects localizers
    """

    def __init__(self):
        super(NIMSStack, self).__init__()

    def _chk_congruent(self, meta):
        """do not check for ImageOrientationPatient"""
        is_dummy = 'Rows' not in meta or 'Columns' not in meta
        if is_dummy and not self._allow_dummies:
            raise dcmstack.IncongruentImageError('Missing Rows/Columns')

        if self._ref_input is not None:
            self._chk_equal(('PixelSpacing', ), meta, self._ref_input)
            if not is_dummy:
                self._chk_equal(('Rows', 'Columns'), meta, self._ref_input)
        elif len(self._dummies) != 0:
            self._chk_equal(('PixelSpacing', ), meta, self._dummies[0][0])
        return is_dummy

    def get_shape(self):
        """get shape needs to artificially adjust slice spacing..."""
        # this is where a lot of edge case handling must go
        # localizer; build as one nifti
        # incomplete stack, insert dummy slices to make stack complete
        return super(NIMSStack, self).get_shape()

    def to_nifti_wrapper(self, voxel_order='LAS'):
        return super(NIMSStack, self).to_nifti_wrapper(voxel_order=voxel_order)


class NIMSDicomError(nimsmrdata.NIMSMRDataError):
    pass


class NIMSDicom(nimsmrdata.NIMSMRData):

    __metaclass__ = abc.ABCMeta

    filetype = u'dicom'
    priority = 0
    parse_priority = 9

    @abc.abstractmethod
    def __init__(self, archive):
        super(NIMSDicom, self).__init__()

        for ti in (ti for ti in archive if ti.isreg()):
            try:
                self.metadata._hdr = dicom.read_file(cStringIO.StringIO(archive.extractfile(ti).read()), stop_before_pixels=True)
                log.debug('dicom found')
            except Exception:
                pass
            else:
                break
        if not self.metadata._hdr:
            log.debug('no readable dicom?')
            raise NIMSDicomError('no readable data?')

        log.debug('parsing standard dicom tags')
        self.metadata.exam_no = self.getelem(self.metadata._hdr, 'StudyID', int)
        self.metadata.patient_id = self.getelem(self.metadata._hdr, 'PatientID')
        self.metadata.subj_code, self.metadata.group_name, self.metadata.experiment_name = self.parse_patient_id(self.metadata.patient_id, 'ex' + str(self.metadata.exam_no))

        def acq_date(header):
            if 'AcquisitionDate' in header:     return header.AcquisitionDate       #noqa
            elif 'StudyDate' in header:         return header.StudyDate             #noqa
            else:                               return '19000101'                   #noqa

        def acq_time(header):
            if 'AcquisitionTime' in header:     return header.AcquisitionTime       #noqa
            elif 'StudyTime' in header:         return header.StudyTime             #noqa
            else:                               return '000000'                     #noqa

        self.metadata.exam_uid = self.getelem(self.metadata._hdr, 'StudyInstanceUID')
        self.metadata.series_uid = self.getelem(self.metadata._hdr, 'SeriesInstanceUID')
        self.metadata.acq_no = self.getelem(self.metadata._hdr, 'AcquisitionNumber', int, 0)
        self.metadata.timestamp = datetime.datetime.strptime(acq_date(self.metadata._hdr) + acq_time(self.metadata._hdr)[:6], '%Y%m%d%H%M%S')
        self.metadata.series_no = self.getelem(self.metadata._hdr, 'SeriesNumber', int)
        self.metadata.series_desc = self.getelem(self.metadata._hdr, 'SeriesDescription')
        self.metadata.protocol_name = self.getelem(self.metadata._hdr, 'ProtocolName', None, 'unknown')
        self.metadata.operator = self.getelem(self.metadata._hdr, 'OperatorsName', None, '')
        self.metadata.subj_firstname, self.metadata.subj_lastname = self.parse_subject_name(self.getelem(self.metadata._hdr, 'PatientName', None, ''))
        self.metadata.subj_dob = self.parse_subject_dob(self.getelem(self.metadata._hdr, 'PatientBirthDate', None, ''))
        self.metadata.subj_sex = {'M': 'male', 'F': 'female'}.get(self.getelem(self.metadata._hdr, 'PatientSex'))
        self.metadata.manufacturer = self.getelem(self.metadata._hdr, 'Manufacturer', None, '')
        self.metadata.scanner_name = '%s %s'.strip() % (self.getelem(self.metadata._hdr, 'InstitutionName', None, ''), self.getelem(self.metadata._hdr, 'StationName', None, ''))
        self.metadata.scanner_type = '%s %s'.strip() % (self.metadata.manufacturer, self.getelem(self.metadata._hdr, 'ManufacturerModelName', None, ''))
        self.metadata.acquisition_type = self.getelem(self.metadata._hdr, 'MRAcquisitionType', None, 'unknown')
        self.metadata.image_type = self.getelem(self.metadata._hdr, 'ImageType', None, [])
        self.metadata.tr = self.getelem(self.metadata._hdr, 'RepetitionTime', float, 0.) / 1000.0
        self.metadata.ti = self.getelem(self.metadata._hdr, 'InversionTime', float, 0.) / 1000.0
        self.metadata.te = self.getelem(self.metadata._hdr, 'EchoTime', float, 0.) / 1000.0
        self.metadata.flip_angle = self.getelem(self.metadata._hdr, 'FlipAngle', float, 0.)
        self.metadata.pixel_bandwidth = self.getelem(self.metadata._hdr, 'PixelBandwidth', float, 0.)
        self.metadata.phase_encode = int(self.getelem(self.metadata._hdr, 'InPlanePhaseEncodingDirection', None, '') == 'COL')
        self.metadata.num_averages = self.getelem(self.metadata._hdr, 'NumberOfAverages', int, 1)
        self.metadata.num_echos = self.getelem(self.metadata._hdr, 'EchoNumbers', int, 1)

        # extract/calculate if possible
        self.metadata.size = [self.getelem(self.metadata._hdr, 'Columns', int, 0), self.getelem(self.metadata._hdr, 'Rows', int, 0)]
        self.metadata.mm_per_vox = self.getelem(self.metadata._hdr, 'PixelSpacing', float, [1., 1.]) + [self.getelem(self.metadata._hdr, 'SpacingBetweenSlices', float, self.getelem(self.metadata._hdr, 'SliceThickness', float, 1.))]

        # set up 'null' values
        self.metadata.bvecs = None
        self.metadata.bvals = None
        self.metadata.num_bands = 0
        self.metadata.notes = ''
        self.metadata.slice_order = nimsmrdata.SLICE_ORDER_UNKNOWN
        self.metadata.mt_offset_hz = 0              # GE
        self.metadata.effective_echo_spacing = 0    # GE
        self.metadata.slice_encode_undersample = 1  # GE
        self.metadata.phase_encode_undersample = 1  # GE

    def _postinit(self):
        """common methods that must be carried out at the end of init"""
        # TODO: higher level inferences will be done by "job" (not within nimsdata)
        self.adjust_fov()
        self.metadata.psd_type = nimsmrdata.infer_psd_type(self.metadata.psd_name)
        self.metadata.scan_type = self.infer_scan_type()

    def adjust_fov(self):
        """adjusts acquisition matrix and fov. returns None."""
        if self.metadata.phase_encode == 1:
            # The Acquisition matrix field includes four values: [freq rows, freq columns, phase rows, phase columns].
            # E.g., for a 64x64 image, it would be [64,0,0,64] if the image row axis was the frequency encoding axis or
            # [0,64,64,0] if the image row was the phase encoding axis.
            self.metadata.acquisition_matrix = self.getelem(self.metadata._hdr, 'AcquisitionMatrix', None, [0, 0, 0, 0])[0:4:3]
            self.metadata.fov[1] /= (self.getelem(self.metadata._hdr, 'PercentPhaseFieldOfView', float, 0.) / 100.) if 'PercentPhaseFieldOfView' in self.metadata._hdr else 1.
        else:
            # We want the acq matrix to always be ROWS,COLS, so we flip the order for the case where the phase encode is the first dim:
            self.metadata.acquisition_matrix = self.getelem(self.metadata._hdr, 'AcquisitionMatrix', None, [0, 0, 0, 0])[2:0:-1]
            self.metadata.fov[0] /= (self.getelem(self.metadata._hdr, 'PercentPhaseFieldOfView', float, 0.) / 100.) if 'PercentPhaseFieldOfView' in self.metadata._hdr else 1.

    @abc.abstractmethod
    def load_data(self, archive):
        """load data from archive into list of dicom objects"""
        super(NIMSDicom, self).load_data(archive)
        log.debug('loading data')
        self.metadata._dcm_list = []
        for ti in (ti for ti in archive if ti.isreg()):
            try:
                dcm = dicom.read_file(cStringIO.StringIO(archive.extractfile(ti).read()), stop_before_pixels=False)
                log.debug('read dicomfile %s' % ti)
            except Exception as e:
                log.debug(e)
                pass
            else:
                self.metadata._dcm_list.append(dcm)
        self.metadata._dcm_list.sort(key=lambda dcm: dcm.InstanceNumber)
        log.debug('data loaded')

    def convert(self, outbase):
        """write to file"""
        TYPE_SCREEN =   ['DERIVED', 'SECONDARY', 'SCREEN SAVE']     # FIXME; need a good catch all for screen saves, need to examine Siemens Screen Save

        if not self.metadata.image_type:
            log.warning('dicom conversion has failed for %s: ImageType not set in dicom header' % os.path.basename(outbase))
            return

        result = (None, None)
        if self.metadata.image_type == TYPE_SCREEN:
            for i, dcm in enumerate(self.metadata._dcm_list):
                result = ('bitmap', nimspng.NIMSPNG.write(self, dcm.pixel_array, outbase + '_%d' % (i+1)))
        elif 'PRIMARY' in self.metadata.image_type:
            imagedata = self.data
            metadata = self.metadata
            result = ('nifti', nimsnifti.NIMSNifti.write(metadata, imagedata, outbase, self.metadata.notes))

        if result[0] is None:
            log.warning('dicom conversion failed for %s: no applicable conversion defined' % os.path.basename(outbase))

        return result

    @staticmethod
    def getelem(hdr, tag, type_=None, default=None):
        """dicom specific element getter"""
        try:
            value = getattr(hdr, tag) if isinstance(tag, basestring) else hdr[tag].value
            if type_ is not None:
                value = [type_(x) for x in value] if isinstance(value, list) else type_(value)
        except (AttributeError, KeyError, ValueError):
            value = default
        return value
