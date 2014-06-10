# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S Hahn

"""
nimsdata.nimsmrdata
===================

Provides base class, and functions that apply to Magnetic Resonsnace (MR) data.

"""

import abc
import bson
import string
import logging
import datetime
import dcmstack

import numpy as np

import nimsdata

log = logging.getLogger(__name__)


scan_types = [
        'spectroscopy',
        'perfusion',
        'shim',
        'diffusion',
        'fieldmap',
        'functional',
        'calibration',
        'localizer',
        'anatomy_t1w',
        'anatomy_t2w',
        'anatomy',
        'unknown',
        ]
scan_types = type('Enum', (object,), dict(zip(scan_types, scan_types), all=scan_types))


# NIFITI1-style slice order codes:
SLICE_ORDER_UNKNOWN = 0
SLICE_ORDER_SEQ_INC = 1
SLICE_ORDER_SEQ_DEC = 2
SLICE_ORDER_ALT_INC = 3
SLICE_ORDER_ALT_DEC = 4
SLICE_ORDER_ALT_INC2 = 5        # interleaved, increased, starting at 2nd MRI slice
SLICE_ORDER_ALT_DEC2 = 6        # interleave, decreasing, starting at one before last MRI slice


def compute_rotation(row_cos, col_cos, slice_norm):
    """
    Compute rotation given row cosines, column cosines and slice norm.

    Parameters
    ----------
    row_cos : np.array
        numpy array of 3 floats
    col_cos : np.array
        numpy array of 3 floats
    slice_norm : np.array
        numpy array of 3 floats

    Returns
    -------
    rot : np.array

    """
    rot = np.zeros((3, 3))
    rot = np.matrix(((-row_cos[0], -col_cos[0], -slice_norm[0]),
                     (-row_cos[1], -col_cos[1], -slice_norm[1]),
                     (row_cos[2], col_cos[2], slice_norm[2])), dtype=float)
    return rot


def build_affine(rotation, scale, origin):
    """
    Compute affine matrix given rotation, scaling, and origin.

    Parameters
    ----------
    rotation : np.array
        rotation
    scale : np.array
        scale factor

    Returns
    -------
    aff : np.array [4x4]
        affine matrix

    """
    aff = np.zeros((4, 4))
    aff[0:3, 0:3] = rotation
    aff[:, 3] = np.append(origin, 1).T
    aff[0:3, 0:3] = np.dot(aff[0:3, 0:3], np.diag(scale))
    return aff


def adjust_bvecs(bvecs, bvals, vendor, rotation=None):
    """
    Compute bvec and bvals given bvecs, bvals and rotation.

    Parameters
    ----------
    bvecs : list of floats
        b vectors, extracted once per volume
    bvals : list of floats
        b value, extracted once per volume
    vendor : str
        manufacturer, exactly as it is in the dicom
    rotation: np.array
        rotation

    Returns
    -------
    (bvecs, bvals) : tuple(list of floats, list of floats)

    """
    bvecs, bvals = scale_bvals(bvecs, bvals)
    # TODO: Uncomment the following when we are ready to fix the bvec flip issue:
    # if vendor.lower().startswith('ge') and rotation != None:
    #   log.debug('rotating bvecs with image orientation matrix')
    #   bvecs,bvals = rotate_bvecs(bvecs, bvals, rotation)
    # else:
    #   bvecs,bvals = rotate_bvecs(bvecs, bvals, np.diag((-1.,-1.,1.)))
    return bvecs, bvals


def scale_bvals(bvecs, bvals):
    """
    Scale the b-values given non-unit-lengh bvecs.

    Scale the b-values in bvals given non-unit-length bvecs. E.g., if the magnitude a
    bvec is 0.5, the corresponding bvalue will be scaled by 0.5^2. The bvecs are also
    scaled to be unit-length. Returns the adjusted bvecs and bvals.

    Parameters
    ----------
    bvecs : list of floats
        b-vector list
    bvals : list of floats
        b-value list

    Returns
    -------
    (bvecs, bvals) : tuple(list of floats, list of floats)

    """
    if np.count_nonzero(bvecs) != 0 and np.count_nonzero(bvals) != 0:
        # if bvecs and bvals are all zeros, then there is no need to scale or rotate
        sqmag = np.array([bv.dot(bv) for bv in bvecs.T])
        # The bvecs are generally stored with 3 decimal values. So, we get significant fluctuations in the
        # sqmag due to rounding error. To avoid spurious adjustments to the bvals, we round the sqmag based
        # on the number of decimal values.
        # TODO: is there a more elegant way to determine the number of decimals used?
        try:
            num_decimals = np.nonzero([np.max(np.abs(bvecs - bvecs.round(decimals=d))) for d in range(9)])[0][-1] + 1
        except IndexError:
            # the bvecs don't have ANY decimals, and thus a rounding threshold cannot be set
            log.warning('there are no decimals, cannot intelligently scale bvecs/bvals')
            num_decimals = 1
        sqmag = np.around(sqmag, decimals=num_decimals - 1)
        bvals *= sqmag            # Scale each bval by the squared magnitude of the corresponding bvec
        sqmag[sqmag == 0] = np.inf  # Avoid divide-by-zero
        bvecs /= np.sqrt(sqmag)   # Normalize each bvec to unit length
    return bvecs, bvals


def rotate_bvecs(bvecs, bvals, rotation):
    """
    Rotate diffusion gradient directions (bvecs) based on the 3x3 rotation matrix.

    Returns the adjusted bvecs and bvals.

    Parameters
    ----------
    bvecs : list of floats
        b-vector list
    bvals : list of floats
        b-value list
    rotation : 3x3 np array
        rotation matrix

    Returns
    ------
    (bvecs, bvals) : tuple(list of floats, list of floats)

    """
    bvecs = np.array(np.matrix(rotation) * bvecs)
    # Normalize each bvec to unit length
    norm = np.sqrt(np.array([bv.dot(bv) for bv in bvecs.T]))
    norm[norm == 0] = np.inf  # Avoid divide-by-zero
    bvecs /= norm
    return bvecs, bvals


def infer_psd_type(manufacturer, psd_name):
    """
    Infer the psd type based on the manufacturer and psd name.

    This heuristic is entirely based on parsing the name.

    Parameters
    ----------
    manufacturer : str
        manufacturer, the same as it is contained within dicom
    psd_name : str
        psd_name, the sme as it is contained within the dicom

    Returns
    -------
    psd_type : str
        psd category

    """
    psd_type = 'unknown'
    if not psd_name:
        return psd_type

    psd_name = psd_name.lower()
    if manufacturer == 'GE MEDICAL SYSTEMS':
        if 'service' in psd_name:
            psd_type = 'service'
        elif psd_name == 'sprt':
            psd_type = 'spiral'
        elif psd_name == 'sprl_hos':
            psd_type = 'hoshim'
        elif psd_name == 'basic':
            psd_type = 'basic'
        elif 'mux' in psd_name:  # multi-band EPI!
            psd_type = 'muxepi'
        elif 'epi' in psd_name:
            psd_type = 'epi'
        elif psd_name in ['probe-mega', 'gaba_ss_cni']:
            psd_type = 'mrs'
        elif psd_name == 'asl':
            psd_type = 'asl'
        elif psd_name in ['bravo', '3dgrass']:
            psd_type = 'spgr'
        elif 'fgre' in psd_name:        # also want to catch efgre3d
            psd_type = 'gre'
        elif psd_name == 'ssfse':
            psd_type = 'fse'
        elif psd_name == 'cube':
            psd_type = 'cube'
        else:
            psd_type = 'unknown'
    elif manufacturer == 'SIEMENS':
        if psd_name == '%siemensseq%\\tse_vfl':
            psd_type = 'tse'
        elif psd_name == '%siemensseq%\\ep2d_diff':
            psd_type =  'epi'
        elif psd_name == '%siemensseq%\\ep2d_bold':
            psd_type =  'epi'
        elif psd_name == '%siemensseq%\\ep2d_asl':
            psd_type =  'asl'
        elif psd_name == '%siemensseq%\\gre':
            psd_type =  'gre'
        elif psd_name == '%siemensseq%\\tfl':
            psd_type =  'tfl'
        elif psd_name == '%siemensseq%\\gre_field_mapping':
            psd_type =  'gre'
        elif psd_name == '%serviceseq%\\rf_noise':
            psd_type =  'service'
        elif psd_name.startswith('%customerseq%\\ep2d_pasl'):
            psd_type = 'asl'
        elif psd_name.startswith('%customerseq%\\ep2d'):
            # %CustomerSeq%\ep2d_advdiff_511E
            psd_type = 'epi'
        elif psd_name == '%customerseq%\\wip711_moco\\tfl_multiecho_epinav_711':
            psd_type = 'tfl'
        else:
            psd_type = 'unknown'

    return psd_type


# TODO: infer_scan_type should be broken down by manufacturer for maintainability
# screensave, tensor and csa don't really matter here, mark them as 'nonstandard'?
def infer_scan_type(psd_type, num_timepoints, te, fov, mm_per_vox, is_dwi=None, is_localizer=None):
    """
    Infer scan type, based on metadata.

    Some scan types, such as localizer and DWI, can be positively identified if examining
    all the dicoms from a scan series.  Special params for positive identification can be used
    to "shortcut", these args default to None.

    Parameters
    ----------
    psd_type : str
        psd type, as output by infer_psd_type()
    num_timepoints : int
        repetitions, link to MR definition
    te : float
        echo time, link to MR definition
    fov : list, 2 item
    mm_per_vox : list, 3 item
        voxel size, link to MR definition
    is_dwi : boolean
        boolean indicate if this scan is diffusion weight image

    Returns
    -------
    scan_type : str
        category of scan type

    """
    scan_type = scan_types.unknown
    if is_dwi:
        scan_type = scan_types.diffusion
    elif is_localizer:
        scan_type = scan_types.localizer
    elif psd_type == 'mrs':
        scan_type = scan_types.spectroscopy
    elif psd_type == 'asl':
        scan_type = scan_types.perfusion
    elif psd_type == 'hoshim':
        scan_type = scan_types.shim
    elif psd_type == 'spiral' and num_timepoints == 2 and te < .05:
        scan_type = scan_types.fieldmap
    elif 'epi' in psd_type and te > 0.02 and te < 0.05 and num_timepoints > 2:
        scan_type = scan_types.functional
    elif (psd_type == 'gre' or psd_type == 'fse') and fov[0] >= 240. and fov[1] >= 240. and mm_per_vox[2] >= 4.5:
        # Could be a low-res calibration scan (e.g., ASSET cal).
        # this provides an alternative way to identify localizer from a single dcm.
        # is_localier requires all dicoms in a scan acquisition.
        if mm_per_vox[0] >= 2.:
            scan_type = scan_types.calibration
        else:
            scan_type = scan_types.localizer
    elif psd_type in ['spgr', 'tfl']:
        scan_type = scan_types.anatomy_t1w
    elif psd_type == 'cube':
        scan_type = scan_types.anatomy_t2w
    return scan_type


class NIMSMRDataError(nimsdata.NIMSDataError):
    pass


class NIMSMRReader(nimsdata.NIMSReader):

    """
    Base MR data reader class.

    A NIMSMRReader object will have a few special attributes that must be defined in order
    to proceed to writing.
    ds.data = {'imagedata': None'}      # must have at least this, which contains the primary imagedata

    in some cases, scans might have "supporting" or "auxillary" data that should also be writtten out.
    in these cases, the data is included in the dicionary under a different key name.

    for example, a scan that has both imagedata and fieldmap_data would define ds.data as
    `ds.data = {'imagedata': voxeldata, 'fieldmap_data': fm_voxeldata}`

    NIMSMRReaders assume that multiple data values have the same metadata; the metadata applies to
    all of the included voxel arrays in ds.data.

    This is especially useful in cases like multicoil data, images with included B0 (fieldmap).

    Parameters
    ----------
    path : str
        path to input file
    load_data : boolean
        indicate if a reader should attempt to immediately load all data. default False.

    """

    __metaclass__ = abc.ABCMeta

    _session_properties = {
        'exam_uid': {
            'field': 'exam_uid',
            'title': 'DICOM UID',
            'type': 'string',
        },
        'exam': {
            'field': 'exam_no',
            'title': 'Exam Number',
            'type': 'integer',
        },
        'patient_id': {
            'field': 'patient_id',
            'title': 'Patient ID',
            'type': 'string',
        },
        'subject': {
            'type': 'object',
            'properties': {
                'firstname': {
                    'field': 'subj_firstname',
                    'title': 'First Name',
                    'type': 'string',
                },
                'lastname': {
                    'field': 'subj_lastname',
                    'title': 'Last Name',
                    'type': 'string',
                },
                'dob': {
                    'field': 'subj_dob',
                    'title': 'Date of Birth',
                    'type': 'string',
                    'format': 'date',
                },
                'sex': {
                    'field': 'subj_sex',
                    'title': 'Sex',
                    'type': 'string',
                    'enum': ['male', 'female'],
                },
            },
        },
        'instrument': {
            'title': 'Instrument', # FIXME: should be surrounded by an array of multiple instruments
            'type': 'object',
            'properties': {
                'manufacturer': {
                    'field': 'manufacturer',
                    'title': 'Manufacturer',
                    'type': 'string',
                },
                'model': {
                    'field': 'manufacturer_model',
                    'title': 'Model',
                    'type': 'string',
                },
                'uid': {
                    'field': 'scanner_name',
                    'title': 'ID',
                    'type': 'string',
                },
            },
        },
    }
    session_properties = nimsdata.dict_merge(nimsdata.NIMSReader.session_properties, _session_properties)

    _epoch_properties = {
        'series_uid': {
            'field': 'series_uid',
            'title': 'DICOM UID',
            'type': 'string',
        },
        'series': {
            'field': 'series_no',
            'title': 'Series',
            'type': 'integer',
        },
        'acquisition': {
            'field': 'acq_no',
            'title': 'Acquisition',
            'type': 'integer',
        },
        'psd': {
            'field': 'psd_name',
            'title': 'PSD',
            'type': 'string',
            'maxLength': 64,
        },
        'tr': {
            'field': 'tr',
            'title': 'Tr',
            'type': 'number',
        },
        'te': {
            'field': 'te',
            'title': 'Te',
            'type': 'number',
        },
        'ti': {
            'field': 'ti',
            'title': 'Ti',
            'type': 'number',
        },
        'flip_angle': {
            'field': 'flip_angle',
            'title': 'Flip Angle',
            'type': 'integer',
        },
        'pixel_bandwidth': {
            'field': 'pixel_bandwidth',
            'title': 'Pixel Bandwidth',
            'type': 'number',
        },
        'num_averages': {
            'field': 'num_averages',
            'title': 'Averages',
            'type': 'integer',
        },
        'num_bands': {
            'field': 'num_bands',
            'title': 'Bands',
            'type': 'integer',
        },
        'num_echos': {
            'field': 'num_echos',
            'title': 'Echos',
            'type': 'integer',
        },
        'num_slices': {
            'field': 'num_slices',
            'title': 'Slices',
            'type': 'integer',
        },
        'num_timepoints': {
            'field': 'num_timepoints',
            'title': 'Time Points',
            'type': 'integer',
        },
        'rx_coil': {
            'field': 'receive_coil_name',
            'title': 'Coil',
            'type': 'string',
            'maxLength': 64,
        },
        'num_receivers': {
            'field': 'num_receivers',
            'title': 'Receivers',
            'type': 'integer',
        },
        'protocol': {
            'field': 'protocol_name',
            'title': 'Protocol',
            'type': 'string',
            'maxLength': 64,
        },
        'size': {
            'title': 'Size',
            'type': 'object',
            'properties': {
                'x': {
                    'field': 'size_x',
                    'title': 'X',
                    'type': 'integer',
                },
                'y': {
                    'field': 'size_y',
                    'title': 'Y',
                    'type': 'integer',
                },
            },
        },
        'fov': {
            'title': 'Field of View',
            'type': 'object',
            'properties': {
                'x': {
                    'field': 'fov_x',
                    'title': 'X',
                    'type': 'integer',
                },
                'y': {
                    'field': 'fov_y',
                    'title': 'Y',
                    'type': 'integer',
                },
            },
        },
        'mm_per_voxel': {
            'title': 'mm per Voxel',
            'type': 'object',
            'properties': {
                'x': {
                    'field': 'mm_per_vox_x',
                    'title': 'X',
                    'type': 'number',
                },
                'y': {
                    'field': 'mm_per_vox_x',
                    'title': 'Y',
                    'type': 'number',
                },
                'z': {
                    'field': 'mm_per_vox_z',
                    'title': 'Z',
                    'type': 'number',
                },
            },
        },
        'effective_echo_spacing': {
            'field': 'effective_echo_spacing',
            'title': 'Effective Echo Spacing',
            'type': 'number',
        },
        'duration': {
            'field': 'duration',
            'title': 'Duration',
            'type': 'number',
        },
        'prescribed_duration': {
            'field': 'prescribed_duration',
            'title': 'Prescribed Duration',
            'type': 'number',
        },
        'slice_encode_undersample': {
            'field': 'slice_encode_undersample',
            'title': 'Slice Encode Undersample',
            'type': 'integer',
        },
        'phase_encode_undersample': {
            'field': 'phase_encode_undersample',
            'title': 'Phase Encode Undersample',
            'type': 'integer',
        },
        'device': {
            'field': 'scanner_name',
            'title': 'Device',
            'type': 'string',
            'maxLength': 64,
        },
        'acquisition_matrix': {
            'title': 'Acquisition Matrix',
            'type': 'object',
            'properties': {
                'x': {
                    'field': 'acquisition_matrix_x',
                    'title': 'X',
                    'type': 'number',
                },
                'y': {
                    'field': 'acquisition_matrix_y',
                    'title': 'Y',
                    'type': 'number',
                },
            },
        },
    }
    epoch_properties = nimsdata.dict_merge(nimsdata.NIMSReader.epoch_properties, _epoch_properties)

    @abc.abstractmethod
    def __init__(self, path, load_data=False):
        super(NIMSMRReader, self).__init__(path, load_data)
        self.is_screenshot = False
        self.epoch_id = None
        self.study_datetime = None
        self.subj_code = None
        self.series_desc = None
        self.scan_type = None

    @abc.abstractmethod
    def load_data(self):
        super(NIMSMRReader, self).load_data()

    @property
    def nims_group_id(self):
        return self.group_name

    @property
    def nims_experiment(self):
        return self.experiment_name

    @property
    def nims_session_id(self):
        return self.exam_uid

    @property
    def nims_epoch_id(self):
        if self.epoch_id: # as in pfile json header
            return self.epoch_id
        series_uid = self.series_uid
        if self.is_screenshot:
            front, back = self.series_uid.rsplit('.', 1)
            series_uid = front + '.' + str(int(back) - 1)       # don't change real series_uid
        return series_uid + ('_' + str(self.acq_no) if self.acq_no is not None else '')

    @property
    def nims_type(self):
        return ('original', 'mri', self.filetype)

    @property
    def nims_filename(self):
        return self.nims_epoch_id + '_' + self.filetype

    @property
    def nims_file_ext(self):
        return '.tgz'

    @property
    def nims_timestamp(self):
        return self.timestamp.replace(tzinfo=bson.tz_util.FixedOffset(-7 * 60, 'pacific'))  # FIXME: use pytz

    @property
    def nims_timezone(self):
        pass  # FIXME

    @property
    def nims_session_label(self):
        return self.study_datetime and self.study_datetime.strftime('%Y-%m-%d %H:%M')

    @property
    def nims_session_subject(self):
        return self.subj_code

    @property
    def nims_session_type(self):
        pass  # FIXME

    @property
    def nims_epoch_label(self):
        return '%d.%d' % (self.series_no, self.acq_no) if self.acq_no is not None else str(self.series_no)

    @property
    def nims_epoch_description(self):
        return self.series_desc

    @property
    def nims_epoch_type(self):
        return self.scan_type

    def _get_mm_per_vox(self):
        return (self.mm_per_vox_x, self.mm_per_vox_y, self.mm_per_vox_z)
    def _set_mm_per_vox(self, value):
        self.mm_per_vox_x, self.mm_per_vox_y, self.mm_per_vox_z = value
    mm_per_vox = property(_get_mm_per_vox, _set_mm_per_vox)

    def _get_size(self):
        return (self.size_x, self.size_y)
    def _set_size(self, value):
        self.size_x, self.size_y = tuple(value)
    size = property(_get_size, _set_size)

    def _get_acquisition_matrix(self):
        return (self.acquisition_matrix_x, self.acquisition_matrix_y)
    def _set_acquisition_matrix(self, value):
        self.acquisition_matrix_x, self.acquisition_matrix_y = value
    acquisition_matrix = property(_get_acquisition_matrix, _set_acquisition_matrix)

    def _get_fov(self):
        return (self.fov_x, self.fov_y)
    def _set_fov(self, value):
        self.fov_x, self.fov_y = value
    fov = property(_get_fov, _set_fov)

    def parse_patient_id(self, patient_id, default_subj_code):
        """
        Parse a subject code, group name and experiment name from patient_id.

        If the patient id does not contain a subject code, rely on default (exam no).

        Expects formatting
        subjcode@group_name/experiment_name

        Parameters
        ----------
        patient_id : str
            patient_id string from dicom tag (0x10,0x20), 'PatientID'

        Returns
        -------
        subj_code : str
            string of subject identifer
        group_name : str
            string of group name
        experiment_name : str
            string of experiment name

        """
        subj_code = None
        group_name = None
        exp_name = None
        if patient_id is not None and default_subj_code is not None:
            subj_code, _, lab_info = patient_id.strip(string.punctuation + string.whitespace).lower().rpartition('@')
            group_name, _, exp_name = lab_info.partition('/')
        return subj_code or default_subj_code, group_name, exp_name

    def parse_subject_name(self, name):
        """
        Parse subject name.

        TODO: accepts a list of delimiters. if none are provided, default list is used

        expects "lastname" + "delimiter" + "firstname".

        Parameters
        ----------
        name : str
            string of subject first and last name, delimited by a '^' or ' '

        Returns
        -------
        firstname : str
            first name parsed from name
        lastname : str
            last name parsed from name

        """
        # TODO: make configurable or include more delimit chars
        firstname = None
        lastname = None
        if name:
            name = name.strip()
            if '^' in name:
                lastname, firstname = name.split('^', 1)        # GE
            elif ' ' in name:
                firstname, lastname = name.rsplit(None, 1)
            else:
                firstname, lastname = ('', name)
            firstname = firstname.title()
            lastname = lastname.title()
        return firstname, lastname

    def parse_subject_dob(self, dob):
        """
        Parse date string and sanity check.

        expects date string in YYYYMMDD format

        Parameters
        ----------
        dob : str
            dob as string YYYYMMDD

        Returns
        -------
        dob : datetime object

        """
        try:
            dob = datetime.datetime.strptime(dob, '%Y%m%d')
            if dob < datetime.datetime(1900, 1, 1):
                raise ValueError
        except (ValueError, TypeError):
            dob = None
        return dob

    def infer_scan_type(self):
        """
        Infer the scan type given common scan parameters.

        Invokes either infer_scan_type_one or infer_scan_type_all, if data is loaded.

        Parameters
        ----------
        None

        Returns
        -------
        scan_type : str

        """
        return infer_scan_type(self.psd_type, self.num_timepoints, self.te, self.fov, self.mm_per_vox, self.is_dwi, self.is_localizer)


class NIMSMRWriter(nimsdata.NIMSWriter):

    """
    Base MR data writer class.

    Cannot be instantiated.  This class serves as a container for static and class methods
    """

    @staticmethod
    def reorder_voxels(imagedata, qto_xyz, voxel_order):
        """
        Reorder voxel data to the specified voxel_order.

        Does not directly manipulate imagedata or qto_xyz.
        Returns manipulated copies of imagedata and qto_xyz.

        Parameters
        ----------
        imagedata : np.array
            single np.array of image data
        qto_xyz : np.matrix, 4x4
            patient space
        voxel_order : str, 3 char
            3 character voxel order string, such as 'LPS' or 'RAI'
        """
        new_data, new_qto_xyz, _, _ = dcmstack.reorder_voxels(imagedata, qto_xyz, voxel_order)
        return new_data, new_qto_xyz

    @nimsdata.abstractclassmethod
    def write(cls, metadata, imagedata, outbase, voxel_order=None):
        super(NIMSMRWriter, cls).write(metadata, imagedata, outbase)
