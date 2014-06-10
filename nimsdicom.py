# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S Hahn
"""
nimsdata.nimsdicom
==================

nimsdicom stuff, extends nimsmrdata, and adds dicom specific parsing routines.

compatible with any subclass of NIMSMRWriter

"""

import dicom
import logging
import tarfile
import datetime
import dcmstack
import cStringIO
import collections
import dcmstack.extract
import nibabel.nicom.csareader

import numpy as np

import nimsmrdata

log = logging.getLogger(__name__)
dicom.config.auto_convert_VR_mismatch = True

MAX_LOC_DCMS = 150           # maximum number of dicoms allowed in a "localizer"


# Some of the DTI from Davis and UI have a badly formed header field that results in errors that cause an entire
# segment of the siemens header to not get parsed.  I, ksh, am not certain if this is an error with a Siemens Protocol,
# or an error in custom sequence(s).
# Currently, this problem is local to just one set of tags within the Csa Series Header, 'sWiPMemBlock.tFree'
def parse_phoenix_prot(prot_key, prot_val):
    """
    Parse a siemens MrProtocol or MrPhoenixProtocol string.

    modified from dcmstack, to not bail on an Exception.
    """
    if prot_key == 'MrPhoenixProtocol':         # syngo B
        str_delim = '""'
    elif prot_key == 'MrProtocol':              # syngo A
        str_delim = '"'
    else:
        raise ValueError('Unknown protocol key: %s' % prot_key)
    ascconv_start = prot_val.find('### ASCCONV BEGIN ###')
    ascconv_end = prot_val.find('### ASCCONV END ###')
    ascconv = prot_val[ascconv_start:ascconv_end].split('\n')[1:-1]

    result = collections.OrderedDict()
    for line in ascconv:
        try:
            parse_result = dcmstack.extract._parse_phoenix_line(line, str_delim)    # TODO: wrap in try/except, ALL fails result in continue
        except Exception:
            pass
        else:
            if parse_result:
                result[parse_result[0]] = parse_result[1]

    return result


# this function had to be exposed to call to the parse_phoenix_prot above
# however, it is basically EXACTLY the same as the dcmstack source
def csa_series_trans_func(elem):
    """Function for parsing the CSA series sub header element by element."""
    csa_dict = dcmstack.extract.simplify_csa_dict(nibabel.nicom.csareader.read(elem.value))

    # If there is a phoenix protocol, parse it and dump it into the csa_dict
    phx_src = None
    if 'MrPhoenixProtocol' in csa_dict:
        phx_src = 'MrPhoenixProtocol'
    elif 'MrProtocol' in csa_dict:
        phx_src = 'MrProtocol'

    if phx_src is not None:
        phoenix_dict = parse_phoenix_prot(phx_src, csa_dict[phx_src])           # THIS
        del csa_dict[phx_src]
        for key, val in phoenix_dict.iteritems():
            new_key = '%s.%s' % ('MrPhoenixProtocol', key)
            csa_dict[new_key] = val

    return csa_dict

csa_series_trans = dcmstack.extract.Translator(
        'CsaSeries',
        dicom.tag.Tag(0x29, 0x1020),
        'SIEMENS CSA HEADER',
        csa_series_trans_func)


# do not ignore private tags
# dcmstack.extract.MetaExtractor uses ignore_non_ascii_bytes, and ignore_private_tags
NIMSMetaExtractor = dcmstack.extract.MetaExtractor(
        ignore_rules=[dcmstack.extract.ignore_non_ascii_bytes],
        translators=[dcmstack.extract.csa_image_trans, csa_series_trans]
        )


class NIMSStack(dcmstack.DicomStack):

    """subclass of dcmstack.DicomStack with slightly altered behavior."""

    def __init__(self, time_order=None, vector_order=None, allow_dummies=False, meta_filter=None):
        super(NIMSStack, self).__init__(time_order, vector_order, allow_dummies, meta_filter)
        # add sorting item for Siemens multicoil
        self.sort_guesses = self.sort_guesses + ['CsaImage.ImaCoilString']


class GEDicom(object):

    """
    load all data from a set of GE Dicoms, almost all information is available via non-private tags.

    Cannot be instantiated.  This object is merely a container for GE specific processing functions.

    GE saves screenshot of the graphical prescription for each scan, the image contains some useful metadata.
    therefore, any missing values should remain MISSING, regardless if they are needed in calculations, to
    attempt to ensure correctness.

    """

    GEMS_TYPE_ORIG = ['ORIGINAL', 'PRIMARY', 'OTHER']
    GEMS_TYPE_SCREEN = ['DERIVED', 'SECONDARY', 'SCREEN SAVE']
    GEMS_DERIVED_RFMT = ['DERIVED', 'SECONDARY', 'REFORMATTED', 'AVERAGE']
    # GEMS derived reformatted average is some sort of sub selection of images from a previous scan.  The dicoms have
    # partial metadata, but may have variation in PixelSpacing that prevents dcmstack from creating a stack.  Identify
    # these types of files early on, and mark as non-image.

    TAG_BVALUE = (0x0043, 0x1039)                                       # CSA_BVALUE = 'Slop_int_6...Slop_int_9'
    TAG_BVEC = [(0x0019, 0x10bb), (0x0019, 0x10bc), (0x0019, 0x10bd)]   # CSA_BVEC = ['UserData20', 'UserData21', 'UserData22']

    @classmethod
    def parse_one(cls, ds):
        """
        Parse all metadata that can be parsed from a single dicom.

        Called by NIMSData init, if dicom manufacturer is GE Medical Sytems.

        """
        ds.psd_name = ds.getelem(ds._hdr, 'PulseSequenceName')
        ds.psd_iname = ds.getelem(ds._hdr, 'InternalPulseSequenceName')
        ds.fov_x, ds.fov_y = 2 * [ds.getelem(ds._hdr, 'ReconstructionDiameter', float)]
        ds.receive_coil_name = ds.getelem(ds._hdr, 'ReceiveCoilName')
        ds.mt_offset_hz = ds.getelem(ds._hdr, 'OffsetFrequency', float)
        effective_echo_spacing = ds.getelem(ds._hdr, 'EffectiveEchoSpacing', float)                     # default to None
        ds.effective_echo_spacing = effective_echo_spacing / 1e6 if effective_echo_spacing else None      # scale the value
        # some of the nordahl/app data is from GE Signa HDxt machine, which stores some information
        # slightly differently than the GE MR750. HDxt stores Asset_R is stored as a string '1\1',
        # MR750 stores Asset_R as list [1, 1].  The code below attempts to adjust the asset_r if it is a unicode string.
        asset_r = ds.getelem(ds._hdr, 'AssetRFactors')
        if isinstance(asset_r, unicode) and '\\' in asset_r:
            asset_r = map(int, asset_r.split('\\'))
        ds.phase_encode_undersample, ds.slice_encode_undersample = asset_r if asset_r else (None, None)
        # some very old Ge systems will output dicoms that don't define Locations in Acquition, or define it in a way
        # that is weird.  It may incorrectly label the value type as OB, but not be able to translate the value, resulting
        # in the NIMSMetaExtractor excluding it from the it's output metadata.
        ds.num_slices = ds.getelem(ds._hdr, 'LocationsInAcquisition', int)
        ds.total_num_slices = ds.getelem(ds._hdr, 'ImagesInAcquisition', int)
        ds.num_timepoints = ds.getelem(ds._hdr, 'NumberOfTemporalPositions', int)
        # slice check could end up wrong, if both total_num_slices and num_slices are None
        # could force num_slices and total_num_slices into different ORs, to prevent matching if both are None
        # thus only when they are both defined, AND not equal, can this test pass
        if (ds.total_num_slices or 1) == (ds.num_slices or 0):
            ds.total_num_slices = (ds.num_slices or 1) * (ds.num_timepoints or 1)
            log.debug('adjusted total_num_slices from %3d to %3d' % (ds.num_slices, ds.total_num_slices))  # num_slices == 'old' total_num
        # some localizer don't have header field to indicate the number of slices
        # per acquisition.  If the total_number of slices is set, and the num_timepoints is 1
        # then the number of slices should be equal to total number of slices
        if not ds.num_slices and (ds.num_timepoints or 1) == 1:
            ds.num_slices = ds.total_num_slices

        prescribed_duration = (ds.tr or 0) * (ds.num_timepoints or 0) * (ds.num_averages or 1)  # FIXME: only works for fMRI, not anatomical
        if prescribed_duration != 0:
            ds.prescribed_duration = prescribed_duration
            ds.duration = prescribed_duration
        else:
            ds.prescribed_duration = None
            ds.duration = None

        # handle some preliminary identification
        # is screenshot
        if ds.image_type == cls.GEMS_TYPE_SCREEN:
            ds.is_screenshot = True

        dwi_dirs = ds.getelem(ds._hdr, 'UserData24{#DTIDiffusionDir.,Release10.0&Above}', float)
        ds.dwi_dirs = int(dwi_dirs) if dwi_dirs else None  # convert to int if not None
        if ds.image_type == cls.GEMS_TYPE_ORIG and (ds.dwi_dirs or 0) >= 6:
            ds.is_dwi = True
            ds.num_timepoints = 1

    @classmethod
    def parse_all(cls, ds):
        """
        Parse all metadata that requires all dicoms.

        Called by NIMSDicom load_data, if dicom manufacturer is GE Medical System.

        """
        # is non-image
        if ds.image_type == cls.GEMS_DERIVED_RFMT:
            ds.is_non_image = True
        # is localizer
        if ds.total_num_slices and ds.total_num_slices < MAX_LOC_DCMS:
            ornts = set([tuple(d.get('ImageOrientationPatient', [0.] * 6)) for d in ds._dcm_list])
            num_ornts = len(ornts)
            log.debug('%d orientations' % num_ornts)
            ds.is_localizer = True if num_ornts > 1 else None
        # is dwi
        # depending on the version, dwi_dirs could come in as a float in a string ('30.000'),
        # as a regular integer (30), or integer in a string ('30').  casting to float works for all
        # cases, which can then be cast into an int.
        if ds.is_dwi:
            ds.bvals = np.array([float(ds.getelem(d, cls.TAG_BVALUE)[0]) for d in ds._dcm_list[0::ds.num_slices]])
            ds.bvecs = np.array([[ds.getelem(d, cls.TAG_BVEC[i], float) for i in range(3)] for d in ds._dcm_list[0::ds.num_slices]]).transpose()
        # multicoils have one volume worth of 'extra' slices, for the combined images
        # total_num_slices != Reps * num_slices * averages
        # localizers have num_averages = 0; and thus TRs * slices * averages * dti_dirs always will = 0
        # TODO: is there a way to identify a multicoil that does not rely on determining other special cases
        if not ds.is_dwi:
            if (ds.num_timepoints or 1) > 1 and (ds.num_timepoints or 1) * (ds.num_slices or 1) * (ds.num_averages or 1) != (ds.total_num_slices or 1):
                ds.is_multicoil = True
                log.debug('MULTICOIL')
                # GE multicoil is organized such that the combined coil is always after each individual coil
                # files are ordered as slice 1, 1,2,n,combined, slice 2 1,2,n,combined, slice 3 1,2,n,combined
                # this current method of grouping will probably not work for DWI or time series data.
                ds.num_receivers = ds.total_num_slices / ds.num_slices
                ds._dcm_groups = [ds._dcm_list[x::ds.num_receivers] for x in xrange(0, ds.num_receivers)]
                log.debug('groups: %3d; %3d coils + %3d combined' % (len(ds._dcm_groups), ds.num_receivers, 1))
        # attempt to calculate trigger times and slice duration, if the first dicom reports trigger time
        ds.slice_duration = None
        if ds.total_num_slices >= ds.num_slices and ds.getelem(ds._dcm_list[0], 'TriggerTime', float) is not None:
            log.debug('using trigger times to calculate slice order and slice duration')
            trigger_times = np.array([ds.getelem(d, 'TriggerTime', float) for d in ds._dcm_list[0:ds.num_slices]])
            if ds.reverse_slice_order:
                trigger_times = trigger_times[::-1]
            trigger_times_from_first_slice = trigger_times[0] - trigger_times
            if ds.num_slices > 2:
                ds.slice_duration = float(min(abs(trigger_times_from_first_slice[1:]))) / 1000    # msec to sec
                if trigger_times_from_first_slice[1] < 0:
                    ds.slice_order = nimsmrdata.SLICE_ORDER_SEQ_INC if trigger_times[2] > trigger_times[1] else nimsmrdata.SLICE_ORDER_ALT_INC
                else:
                    ds.slice_order = nimsmrdata.SLICE_ORDER_ALT_DEC if trigger_times[2] > trigger_times[1] else nimsmrdata.SLICE_ORDER_SEQ_DEC
            else:
                ds.slice_duration = trigger_times[0]
                ds.slice_order = nimsmrdata.SLICE_ORDER_SEQ_INC


class SiemensDicom(object):

    """
    Load all data from a set of siemens dicoms.

    Not all information required is available via non-private tags.

    Cannot be instantiated.  This is merely a container for Siemens specific functions.

    however, the CSA Series Header and CSA Image Header are stored in non-private tag groups (even group #).
    csa header tags are preferred over standard dicom tags. Although some information overlaps between
    the dicom standard tags and CSA headers, patient information is not included in the CSA.  Therefore,
    anonymization doesn't need to fiddle with the CSA headers.

    """

    # there are several varieties of CSAPARALLEL images, MPR, PROJECTION IMAGE.
    # the only common theme so far is that they all have 'CSAPARALLEL' in imagetype, and are invalid stacks
    # because the dcm metadata are inconsistent
    SIEMENS_TYPE_CSA = [u'ORIGINAL', u'PRIMARY', u'OTHER', u'CSA REPORT']
    SIEMENS_TYPE_TENSOR = [u'DERIVED', u'PRIMARY', u'DIFFUSION', u'TENSOR', u'ND']
    SIEMENS_TYPE_CSAPARALLEL = [u'DERIVED', u'SECONDARY', u'MPR', u'CSA MPR', u'', u'CSAPARALLEL', u'M', u'ND']
    SIEMENS_TYPE_DIFF_FA = [u'DERIVED', u'PRIMARY', u'DIFFUSION', u'FA', u'ND']
    SIEMENS_TYPE_DIFF_FA_NORM = [u'DERIVED', u'PRIMARY', u'DIFFUSION', u'FA', u'ND', u'NORM']
    SIEMENS_TYPE_ASL_TTEST = [u'DERIVED', u'PRIMARY', u'PERFUSION', u'ASL', u'ND', u'NORM', u'FILTERED', u'MOCO', u'SUB', u'TTEST', u'MOSAIC']
    SIEMENS_TYPE_RELCBF_TTEST = [u'DERIVED', u'PRIMARY', u'PERFUSION', u'ASL', u'RELCBF', u'ND', u'NORM', u'FILTERED', u'MOCO', u'SUB', u'TTEST', u'MOSAIC']
    SIEMENS_TYPE_CSA_3D = [u'DERIVED', u'SECONDARY', u'OTHER', u'CSA 3D EDITOR']
    SIEMENS_TYPE_POSDISP = [u'DERIVED', u'SECONDARY', u'POSDISP', u'M', u'ND', u'NORM', u'CSA RESAMPLED']
    # posdisp?  resample images of the previous anatomical scan, to view from 3 angles to estimate positioning?  not really sure
    # but the data was not acquired as such.

    # the Dicom SOPClass Secondaru Capture Storage (SC) are not valid MR Images.  They are missing
    # ImagePatientPosition information necessary for creating the affine.
    SIEMENS_SOPCLASS_SC = '1.2.840.10008.5.1.4.1.1.7'
    TAG_BVALUE = 'CsaSeries.MrPhoenixProtocol.sDiffusion.alBValue[1]'  # B_value
    TAG_BVEC = 'CsaImage.DiffusionGradientDirection'

    @classmethod
    def parse_one(cls, ds):
        """
        Parse all metadata that can be parsed from a single dicom.

        Called by NIMSDicom init, if dicom manufacturer is Siemens.

        """
        ds.psd_name = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.tSequenceFileName')
        ds.psd_iname = ds.getelem(ds._hdr, 'SeriesDescription')
        # FIXME: is fov [phase, readout], or [readout, phase]
        fov_phase = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.asSlice[0].dPhaseFOV', float)
        fov_readout = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.asSlice[0].dReadoutFOV', float)
        ds.fov = [fov_phase, fov_readout]
        ds.receive_coil_name = ds.getelem(ds._hdr, 'CsaImage.ImaCoilString')
        slice_duration = ds.getelem(ds._hdr, 'CsaImage.SliceMeasurementDuration', float, 0.)
        ds.slice_duration = slice_duration / 1e6 if slice_duration else None
        ds.prescribed_duration = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.lScanTimeSec')   # FIXME
        ds.duration = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.lTotalScanTimeSec')         # FIXME: not guaranteed
        ds.acq_no = None

        # preliminary identification
        if ds.image_type == cls.SIEMENS_TYPE_CSA:
            ds.is_non_image = True
        # is non-image
        if ds.getelem(ds._hdr, 'SOPClassUID') == cls.SIEMENS_SOPCLASS_SC:
            # Seconary Storage Class Images always are missing ImagePatientPositions
            # .: they never have enough information to be a valid MR dicom
            ds.is_non_image = True
        if ds.image_type == cls.SIEMENS_TYPE_POSDISP:
            ds.is_non_image = True
        if ds.image_type == cls.SIEMENS_TYPE_CSA_3D:
            # Siemens 3D Storage. what is this?
            ds.is_non_image = True
        if ds.image_type in [cls.SIEMENS_TYPE_ASL_TTEST or ds.image_type, cls.SIEMENS_TYPE_RELCBF_TTEST]:
            # CSA header contains WAY TOO MUCH info...ttest results?
            ds.is_non_image = True
        if ds.image_type in [cls.SIEMENS_TYPE_DIFF_FA_NORM, cls.SIEMENS_TYPE_DIFF_FA]:
            # Diffusion FA NORM dicoms count as is_non_image because they lack the information
            # necessary to create the affine xform for the nifti. Also, CSA data instead of pixel data.
            ds.is_non_image = True
        if u'CSAPARALLEL' in ds.image_type:
            # csaparallel images are secondary derived dataset
            # these will rarely be valid, because variation in the dicoms PixelSpacing
            # the first dicom will accurately hold PixelSpacing from the parent scan
            # all remaining dicoms will have PixelSpacing of [1, 1]
            ds.is_non_image = True
        if ds.getelem(ds._hdr, 'PrivateCreator_0x29_0x10') == 'SIEMENS CSA NON-IMAGE':
            ds.is_non_image = True
        if ds.getelem(ds._hdr, 'VariablePixelData') == 'SIEMENS CSA NON-IMAGE':
            ds.is_non_image = True
        if ds.image_type == cls.SIEMENS_TYPE_TENSOR:
            ds.is_non_image = True
        ds.dwi_dirs = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.sDiffusion.lDiffDirections', int, None)
        if 'ORIGINAL' in ds.image_type and (ds.dwi_dirs or 0) >= 6:
            ds.is_dwi = True
            ds.num_timepoints = 1

    @classmethod
    def parse_all(cls, ds):
        """
        Parse all metadata that requires all dicoms.

        Called by NIMSDicom load_data, if dicom manufacturer is Siemens.

        """
        if 'MOSAIC' not in ds.image_type:
            log.debug('SIEMENS SINGLE SLICE DICOM')
            ds.num_slices = len(set([tuple(ds.getelem(d, 'ImagePositionPatient', None, [0., 0., 0.])) for d in ds._dcm_list]))
            ds.total_num_slices = len(ds._dcm_list)
            ds.num_timepoints = ds.total_num_slices / ds.num_slices
        elif 'MOSAIC' in ds.image_type:           # explicit for readability
            log.debug('SIEMENS MOSAIC')
            ds.num_slices = ds.getelem(ds._hdr, 'CsaImage.NumberOfImagesInMosaic', int, ds.getelem(ds._hdr, 'NumberOfImagesInMosaic', int))
            ds.num_timepoints = len(ds._dcm_list)
            mosaic_dim = int(ds.num_slices ** 0.5)
            if mosaic_dim ** 2 < ds.num_slices:
                mosaic_dim += 1
            ds.size = [x / mosaic_dim for x in ds.size]
            ds.fov = [x / mosaic_dim for x in ds.fov]
            ds.total_num_slices = ds.num_slices * ds.num_timepoints

        ds.duration = ds.num_timepoints * (ds.num_averages or 1) * ds.tr if ds.num_timepoints and ds.tr else None

        # CsaSeries.MrPhoenixProtocol.sSliceArray.ucMode indicates siemens slice order:
        # Siemens; 1 Ascending, 2 Descending, and 4 Interleaved Ascending.
        # NIFTI;   1 Ascending, 2 Descending, and 4 Interleaving Descending
        # Siemens may output interleaved data one of two ways:
        # - nifti slice order 3: interleave_asc, odd first, odd num slices, interleave INC
        # - nifti slice order 5: interleave_asc, even first, even num slices, interleave INC 2
        ds.slice_order = ds.getelem(ds._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.ucMode', None, 0)
        if ds.slice_order == 4:   # don't try to guess if num_slices can't be determined
            if ds.num_slices % 2 != 0:
                ds.slice_order = 3  # interleaved ascending, odd first
            else:
                ds.slice_order = 5  # interleaved ascending, even first

        # is localizer
        if ds.total_num_slices < MAX_LOC_DCMS:
            ornts = [tuple(d.get('ImageOrientationPatient', [0.] * 6)) for d in ds._dcm_list]
            ds.is_localizer = bool(len(set(ornts)) > 1)
        if ds.is_dwi:
            ds.bvals = np.array([NIMSMetaExtractor(d).get(cls.TAG_BVALUE, 0.) for d in ds._dcm_list[0:ds.num_slices]])
            ds.bvecs = np.array([NIMSMetaExtractor(d).get(cls.TAG_BVEC, [0., 0., 0.]) for d in ds._dcm_list[0:ds.num_slices]]).transpose()
        # is multicoil?
        ds.num_receivers = len([ds._hdr[key] for key in ds._hdr if key.endswith('sCoilElementID.tCoilID')])

class NIMSDicomError(nimsmrdata.NIMSMRDataError):
    pass


class NIMSDicom(nimsmrdata.NIMSMRReader):

    """
    should be able to read either manufacturer dicoms.

    if load_data is set to True, the load_data method will be called at the end of init. scan identification
    is done during load_data. although some types can be determined with minimal data, identification is
    done at the same place across venders to provide consistency.

    parse gives a ball park idea of "what it is," without loading everything.

    parameters
    ----------
    path : str
        path to input file
    load_data : bool [default False]
        indicate if all data should be loaded by invoking load_data() at the end of init.

    """

    filetype = u'dicom'
    priority = 0
    parse_priority = 9

    def __init__(self, path, load_data=False):
        """
        Parse a single file from the input.

        Performs manufacturer specific single dicom parsing.

        """
        super(NIMSDicom, self).__init__(path, load_data)
        self.compressed = True      # v1 compat: indicates file is compressed, informing scheduler to not compress.
        with tarfile.open(path) as archive:
            for ti in archive:
                if not ti.isreg():
                    continue
                try:
                    self._hdr = NIMSMetaExtractor(dicom.read_file(cStringIO.StringIO(archive.extractfile(ti).read()), stop_before_pixels=True))
                except Exception:
                    pass
                else:
                    break
            else:
                raise NIMSDicomError('no readable dicom found during __init__')
        del ti
        del archive

        self.image_type = self.getelem(self._hdr, 'ImageType', None, [])
        if not self.image_type:
            log.warning('dicom parsing failed for %s: ImageType not set in dicom header' % self.filepath)
            return  # raise an error?

        log.debug('parsing standard manufacturer-agnostic tags')
        self.exam_no = self.getelem(self._hdr, 'StudyID')                       # siemens may have no studyID; find ex
        self.exam_uid = self.getelem(self._hdr, 'StudyInstanceUID')
        self.patient_id = self.getelem(self._hdr, 'PatientID')                  # siemens may not have patientID. find ex
        self.series_no = self.getelem(self._hdr, 'SeriesNumber', int)
        self.series_desc = self.getelem(self._hdr, 'SeriesDescription')
        self.series_uid = self.getelem(self._hdr, 'SeriesInstanceUID')
        self.subj_code, self.group_name, self.experiment_name = self.parse_patient_id(self.patient_id, 'ex' + self.exam_no)
        # Siemens and GE have different notions of 'acquisition'
        # GE: one "epoch series" may have multiple acquisitions (think 'add to series').
        # Siemens: each entire volume is a "acquisition". for mosaics, each mosaic is own acquisition
        self.acq_no = self.getelem(self._hdr, 'AcquisitionNumber', int, 1)         # wrong for siemens dicom, until load_data
        self.study_date = self.getelem(self._hdr, 'StudyDate')
        self.study_time = self.getelem(self._hdr, 'StudyTime')
        self.acq_date = self.getelem(self._hdr, 'AcquisitionDate')
        self.acq_time = self.getelem(self._hdr, 'AcquisitionTime')
        self.study_datetime = self.study_date and self.study_time and datetime.datetime.strptime(self.study_date + self.study_time[:6], '%Y%m%d%H%M%S')
        self.acq_datetime = self.acq_date and self.acq_time and datetime.datetime.strptime(self.acq_date + self.acq_time[:6], '%Y%m%d%H%M%S')
        self.timestamp = self.acq_datetime or self.study_datetime
        self.subj_firstname, self.subj_lastname = self.parse_subject_name(self.getelem(self._hdr, 'PatientName'))
        self.subj_dob = self.parse_subject_dob(self.getelem(self._hdr, 'PatientBirthDate', str))
        self.subj_sex = {'M': 'male', 'F': 'female'}.get(self.getelem(self._hdr, 'PatientSex'))
        tr = self.getelem(self._hdr, 'RepetitionTime', float)
        self.tr = tr / 1000.0 if tr else None
        ti = self.getelem(self._hdr, 'InversionTime', float)
        self.ti = ti / 1000.0 if ti else None
        te = self.getelem(self._hdr, 'EchoTime', float)
        self.te = te / 1000.0 if te else None
        self.flip_angle = self.getelem(self._hdr, 'FlipAngle', float)
        self.pixel_bandwidth = self.getelem(self._hdr, 'PixelBandwidth', float)
        phase_encode = self.getelem(self._hdr, 'InPlanePhaseEncodingDirection')
        self.phase_encode = int(phase_encode == 'COL') if phase_encode else None
        self.num_averages = self.getelem(self._hdr, 'NumberOfAverages', int)    # is this always an int? some GE scans show this as 0.5?
        self.num_echos = self.getelem(self._hdr, 'EchoNumbers', int)
        self.operator = self.getelem(self._hdr, 'OperatorsName')
        self.protocol_name = self.getelem(self._hdr, 'ProtocolName')
        self.scanner_name = '%s %s'.strip() % (
                self.getelem(self._hdr, 'InstitutionName', None, ''),
                self.getelem(self._hdr, 'StationName', None, '')
                )
        self.manufacturer = self.getelem(self._hdr, 'Manufacturer')
        self.manufacturer_model = self.getelem(self._hdr, 'ManufacturerModelName')
        self.scanner_type = '%s %s'.strip() % (
                self.manufacturer,
                self.manufacturer_model or ''
                )
        self.acquisition_type = self.getelem(self._hdr, 'MRAcquisitionType')
        self.mm_per_vox_x, self.mm_per_vox_y = self.getelem(self._hdr, 'PixelSpacing', float) or (None, None)
        self.mm_per_vox_z = self.getelem(self._hdr, 'SpacingBetweenSlices', float) or self.getelem(self._hdr, 'SliceThickness', float)
        self.size_x = self.getelem(self._hdr, 'Columns', int)
        self.size_y = self.getelem(self._hdr, 'Rows', int)
        # REMINDERS! make sure to do these in mfr specific
        self.reverse_slice_order = False
        self.slice_order = nimsmrdata.SLICE_ORDER_UNKNOWN
        # DEFAULT THINGS TO NONE!!!!!
        self.mt_offset_hz = None               # ge tag
        self.effective_echo_spacing = None     # ge tag, no siemens
        self.slice_encode_undersample = None   # ge tag, no siemens
        self.phase_encode_undersample = None   # ge tag, no siemens
        self.num_bands = None
        self.num_slices = None
        self.num_receivers = None
        self.total_num_slices = None
        self.num_timepoints = None
        # special cases
        self.is_localizer = None
        self.is_screenshot = None          # ge specific
        self.is_dwi = None
        self.is_multicoil = None
        self.is_non_image = None

        if not self.manufacturer:
            # check if any part of the image type starts with 'CSA'
            for sub in self.image_type:
                if sub.startswith('CSA'):
                    self.manufacturer = 'SIEMENS'

        if self.manufacturer == 'GE MEDICAL SYSTEMS':
            GEDicom.parse_one(self)
        elif self.manufacturer == 'SIEMENS':
            SiemensDicom.parse_one(self)
        elif self.manufacturer is not None:
            raise NIMSDicomError('dicom manufacturer %s is not supported' % self.manufacturer)
        else:
            log.warning('dicom manufacturer not defined. cannot parse any further')

        # manufacturer agnostic
        self.acquisition_matrix_x = self.acquisition_matrix_y = None
        acquisition_matrix = self.getelem(self._hdr, 'AcquisitionMatrix')
        if self.phase_encode == 1:
            # The Acquisition matrix field includes four values: [freq rows, freq columns, phase rows, phase columns].
            # E.g., for a 64x64 image, it would be [64,0,0,64] if the image row axis was the frequency encoding axis or
            # [0,64,64,0] if the image row was the phase encoding axis.
            if acquisition_matrix:
                self.acquisition_matrix_x, self.acquisition_matrix_y = acquisition_matrix[0:4:3]
            if self.fov != (None, None):
                fov_x, fov_y = self.fov
                fov_y /= (self.getelem(self._hdr, 'PercentPhaseFieldOfView', float, 0.) / 100.) if 'PercentPhaseFieldOfView' in self._hdr else 1.
                self.fov = fov_x, fov_y
        else:
            # We want the acq matrix to always be ROWS,COLS, so we flip the order for the case where the phase encode is the first dim:
            if acquisition_matrix:
                self.acquisition_matrix_x, self.acquisition_matrix_y = acquisition_matrix[2:0:-1]
            else:
                self.acquisition_matrix_x = self.acquisition_matrix_y = None
            if self.fov != (None, None):
                fov_x, fov_y = self.fov
                self.fov_x /= (self.getelem(self._hdr, 'PercentPhaseFieldOfView', float, 0.) / 100.) if 'PercentPhaseFieldOfView' in self._hdr else 1.
                self.fov = fov_x, fov_y

        self.psd_type = nimsmrdata.infer_psd_type(self.manufacturer, self.psd_name)
        self.scan_type = self.infer_scan_type()

        if load_data:
            self.load_data()
        self.metadata_status = 'pending'

    def load_data(self):
        """
        Load all dicoms and obtain more metadata.

        Performs manufacturer specific loads. Operates through side-effect, does not return anything.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        log.debug('%20s - %s' % ('start loading dicoms', str(datetime.datetime.now())))
        super(NIMSDicom, self).load_data()
        self._dcm_list = []
        with tarfile.open(self.filepath) as archive:
            for ti in archive:
                try:
                    dcm = dicom.read_file(cStringIO.StringIO(archive.extractfile(ti).read()), stop_before_pixels=False)
                    self._dcm_list.append(dcm)
                except Exception:
                    pass
        del ti
        del archive
        log.debug('%20s - %s' % ('done loading dicoms', str(datetime.datetime.now())))

        if self._dcm_list == []:
            raise NIMSDicomError('no dicoms loaded?')
        try:
            self._dcm_list.sort(key=lambda dcm: dcm.InstanceNumber)
        except AttributeError:
            log.debug('dicoms do not have InstanceNumber. cannot pre-sort')   # these are CSA non-image type
        log.debug(self.manufacturer)

        if self.manufacturer == 'GE MEDICAL SYSTEMS':
            GEDicom.parse_all(self)
        elif self.manufacturer == 'SIEMENS':
            SiemensDicom.parse_all(self)

        self.scan_type = self.infer_scan_type()

        log.debug('%20s - %s' % ('start recon', str(datetime.datetime.now())))
        if self.is_non_image:
            log.debug('file %s is a special non-image dicom' % self.filepath)
        elif self.is_screenshot:
            log.debug('screenshot recon')
            # there are some metadata that should never be set for screenshots
            self.slice_order = None  # slice order is not unknown, just Not Applicable.
            self.psd_type = None  # not a real PSD.
            self.data = {'': np.dstack([d.pixel_array for d in self._dcm_list])}
        elif self.is_localizer:
            log.debug('localizer recon')
            num_ornts = len(set([tuple(d.get('ImageOrientationPatient', [0.] * 6)) for d in self._dcm_list]))
            self.num_timepoints = num_ornts
            self.num_slices = self.total_num_slices / self.num_timepoints
            # determine cosines, slice norm and rotation
            cosines = self.getelem(self._hdr, 'ImageOrientationPatient', float, 6 * [np.nan])
            row_cosines = np.array(cosines[0:3])
            col_cosines = np.array(cosines[3:6])
            slice_norm = np.cross(row_cosines, col_cosines)
            rot = nimsmrdata.compute_rotation(row_cosines, col_cosines, slice_norm)
            # determine origin
            image_position = [tuple(self.getelem(d, 'ImagePositionPatient', float, [0., 0., 0.])) for d in self._dcm_list]
            origin = image_position[0] * np.array([-1, -1, 1])
            self.qto_xyz = nimsmrdata.build_affine(rot, self.mm_per_vox, origin)
            self.data = np.dstack([np.swapaxes(d.pixel_array, 0, 1) for d in self._dcm_list])
            dims = np.array((self.size_y, self.size_x, self.num_slices, self.num_timepoints))
            if np.prod(dims) == np.size(self.data):
                self.data = self.data.reshape(dims, order='F')
            self.data = {'': self.data}
            # FIXME: acq matrix reversed for some localizers?
        else:
            # check non-mosaic files for incomplete volumes or incomplete single slice dicom sets
            # AFAIK, mosaics only come in complete volumes, as such, are never missing single slices
            if 'MOSAIC' not in self.image_type:
                if self.total_num_slices < self.num_slices:                 # this might need 'or 0' syntax
                    # TODO: add comment to notes
                    raise NIMSDicomError('cannot reconstruct %s.\n%s' % (self.filepath, 'insufficient dicoms to recon 1 volume'))
                # check for and remove incomplete volumes.
                partial_vol_dcms = self.total_num_slices % self.num_slices  # this might need 'or 0' syntax
                if partial_vol_dcms:
                    log.debug('number of dicoms is not a integer multiple of number of unique slices positions. trimming.')
                    self._dcm_list = self._dcm_list[:-1 * partial_vol_dcms]
                    # TODO: add comment notes

            if self.is_multicoil:
                log.debug('multicoil recon')
                stacks = []
                group_id = 0
                for group in self._dcm_groups:
                    group_id += 1
                    log.debug('multicoil - %2s, %s dicom' % (str(group_id), str(len(group))))
                    num_positions = len(set([d.SliceLocation for d in group]))
                    if num_positions != self.num_slices:
                        raise NIMSDicomError('coil %s has %s unique positions; expected %s' % (group_id, num_positions, self.num_slices))
                    stack = NIMSStack()
                    for dcm in group:
                        meta = NIMSMetaExtractor(dcm)
                        stack.add_dcm(dcm, meta)
                    nii_wrp = stack.to_nifti_wrapper()
                    stacks.append(nii_wrp)
                try:
                    nii_wrp = dcmstack.dcmmeta.NiftiWrapper.from_sequence(stacks)
                except dcmstack.InvalidStackError as e:
                    raise NIMSDicomError('cannot reconstruct %s.\n%s' % (self.filepath, e), log_level=logging.ERROR)
                del self._dcm_groups
                del self._dcm_list
                del stacks
                del stack
                del dcm
            else:
                log.debug('standard recon')
                stack = NIMSStack()
                for dcm in self._dcm_list:
                    stack.add_dcm(dcm, NIMSMetaExtractor(dcm))
                try:
                    nii_wrp = stack.to_nifti_wrapper()
                except dcmstack.InvalidStackError as e:
                    raise NIMSDicomError('cannot reconstruct %s.\n%s' % (self.filepath, e), log_level=logging.ERROR)
                del self._dcm_list
                del stack
                del dcm

            nii = nii_wrp.nii_img
            self.data = {'': nii.get_data()}
            self.qto_xyz = nii.get_affine()
            del nii_wrp
            del nii

            # FIXME: is this a GE specific step?
            if self.is_dwi:
                self.bvecs, self.bvals = nimsmrdata.adjust_bvecs(self.bvecs, self.bvals, self.scanner_type, self.qto_xyz[0:3, 0:3])
            self.metadata_status = 'complete'

        log.debug('%20s - %s' % ('done recon', str(datetime.datetime.now())))

    @staticmethod
    def getelem(hdr, tag, type_=None, default=None):
        """
        Retrieve a dicom header element by its tag.

        Parameters
        ----------
        hdr : pydicom.dataset or csa_dict from dcmstack.extact.NIMSMetaExtractor
            Either a pydicom.dataset object, or a NestedMultiDict as returned by dcmstack.extract.NIMSMetaExtractor.
        tag : tuple or str
            A tuple or string that indicates which tag to retrieve.  The tuple tag can only be used with a pydicom.dataset.
        type_ : class
            Attempt to convert the tag's value to the specified type, such as str, int, float.  Default is None.
        default : any
            What value to return if the tag cannot be found, or if type conversion fails.

        Returns
        -------
        value : any
            the value held at the tag specificed, of type `type_` (if possible).  Or the default value if
            the tag does not exist, or could not be converted to `type_`.
        """
        value = default
        if isinstance(tag, tuple):
            # get by bytecode tuple
            elem = hdr.get(tag)
            if elem:
                value = elem.value
        else:
            # get by keyword
            value = hdr.get(tag, default)
        if type_ is not None:
            try:
                value = [type_(x) for x in value] if isinstance(value, list) else type_(value)
            except TypeError:   # TODO: check applicable Exceptions
                value = default
        return value
