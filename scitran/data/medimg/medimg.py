# @author:  Kevin S. Hahn

"""
scitran.data.medimg.medimg
==========================

Contains generally useful functions to parse information from a medical image, such as dicom,
pfile, siemens raw, or nifti. Which allows all such medical images to use the same fxns to
parse data in a consistent way.

Contains base class MedImgReader, MedImgWriter, which provide schema information, properties
and functions for subclasses.

functions that have to do with parsing names, dob (which has nothing to do with MR data) are
located here.  MR fxns are contained with dcm.mr.generic_mr

"""

import abc
import string
import logging
import datetime
import dcmstack
import numpy as np

from .. import data
from .. import util

log = logging.getLogger(__name__)

# slice order codes. all medical imgaes will use nifti1 style slice order codes
# although these slice order codes "belong" to the nifti format, they can be
# used by any medical images
SLICE_ORDER_UNKNOWN = 0
SLICE_ORDER_SEQ_INC = 1
SLICE_ORDER_SEQ_DEC = 2
SLICE_ORDER_ALT_INC = 3
SLICE_ORDER_ALT_DEC = 4
SLICE_ORDER_ALT_INC2 = 5  # interleaved, increased, starting at 2nd MRI slice
SLICE_ORDER_ALT_DEC2 = 6  # interleave, decreasing, starting at one before last MRI slice

project_properties = data.project_properties

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
            'firstname_hash': {
                'field': 'subj_firstname',
                'type': 'string',
            },
            'lastname_hash': {
                'field': 'subj_lastname',
                'type': 'string',
            },
            'age': {
                'field': 'subj_age',
                'title': 'Age',
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
        'title': 'Instrument',  # FIXME: should be surrounded by an array of multiple instruments
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
session_properties = data.dict_merge(data.session_properties, _session_properties)

_acquisition_properties = {
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
                'field': 'mm_per_vox_y',
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
    'phase_encode_direction': {
        'field': 'phase_encode_direction',
        'title': 'Phase-encode Direction',
        'type': 'integer',
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
acquisition_properties = data.dict_merge(data.acquisition_properties, _acquisition_properties)


def get_slice_order(slice_order, num_slices):
    """
    Return the calculated slice order based on the nifti_slice_code and number of slices.

    See `brainder.org - nifti file format <http://brainder.org/2012/09/23/the-nifti-file-format/>`_ for
    more information.

    Parameters
    ----------
    slice_code : int
        nifti slice code, usually found as .slice_order.
    num_slices: int
        number of slices, usually found as .num_slices.

    Returns
    -------
    slice_order: list of int
        list of slice indexes, in the order acquired
    """
    if slice_order == SLICE_ORDER_UNKNOWN or slice_order is None:
        slice_order_array = None
    elif slice_order == SLICE_ORDER_SEQ_INC:
        slice_order_array = np.arange(0, num_slices)
    elif slice_order == SLICE_ORDER_SEQ_DEC:
        slice_order_array = np.arange(0, num_slices)[::-1]
    elif slice_order == SLICE_ORDER_ALT_INC:
        slice_order_array = np.hstack((np.arange(0, num_slices, 2), np.arange(1, num_slices, 2)))
    elif slice_order == SLICE_ORDER_ALT_DEC:
        slice_order_array = np.hstack((np.arange(0, num_slices, 2), np.arange(1, num_slices, 2)))[::-1]
    elif slice_order == SLICE_ORDER_ALT_INC2:
        slice_order_array = np.hstack((np.arange(1, num_slices, 2), np.arange(0, num_slices, 2)))
    elif slice_order == SLICE_ORDER_ALT_DEC2:
        slice_order_array = np.hstack((np.arange(1, num_slices, 2), np.arange(0, num_slices, 2)))[::-1]
    return slice_order_array


def parse_patient_id(patient_id, default_subj_code):
    """
    Parse a subject code, group name and project name from patient_id.

    If the patient id does not contain a subject code, rely on default (exam no).

    Expects formatting
    subjcode@group_name/project_name

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
    project_name : str
        string of project name

    """
    subj_code = group_name = exp_name = None
    if patient_id is not None and default_subj_code is not None:
        subj_code, _, lab_info = patient_id.strip(string.punctuation + string.whitespace).lower().rpartition('@')
        group_name, _, exp_name = lab_info.partition('/')

    log.debug([subj_code or default_subj_code, group_name, exp_name])
    return subj_code or default_subj_code, group_name, exp_name


def parse_patient_name(name):
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
    if '^' in name:
        lastname, _, firstname = name.partition('^')
    else:
        firstname, _, lastname = name.rpartition(' ')
    return firstname.strip().title(), lastname.strip().title()


def parse_patient_dob(dob):
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

    log.debug(dob)
    return dob

def parse_patient_age(age):
    """
    Parse patient age from string.
    convert from 70d, 10w, 2m, 1y to datetime.timedelta object.
    """
    conversion = {  # conversion to days
        'Y': 365,
        'M': 30,
        'W': 7,
        'D': 1,
    }
    scale = age[-1:]
    value = age[:-1]
    return datetime.timedelta(int(value) * conversion.get(scale))


class MedImgError(data.DataError):
    pass


class MedImgReader(data.Reader):

    """
    Base MR data reader class.

    A MedImgReader object will have a few special attributes that must be defined in order
    to proceed to writing.
    ds.data = {'imagedata': None'}      # must have at least this, which contains the primary imagedata

    in some cases, scans might have "supporting" or "auxillary" data that should also be writtten out.
    in these cases, the data is included in the dicionary under a different key name.

    for example, a scan that has both imagedata and fieldmap_data would define ds.data as
    `ds.data = {'imagedata': voxeldata, 'fieldmap_data': fm_voxeldata}`

    MedImgReaders assume that multiple data values have the same metadata; the metadata applies to
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
    project_properties = project_properties
    session_properties = session_properties
    acquisition_properties = acquisition_properties

    @abc.abstractmethod
    def __init__(self, path, load_data=False, timezone=None):
        super(MedImgReader, self).__init__(path, load_data, timezone)
        self.is_screenshot = False
        self.acquisition_id = None
        self.study_datetime = None
        self.subj_code = None
        self.series_desc = None
        self.scan_type = None
        self.is_localizer = None
        self.is_dwi = None
        self.is_fastcard = None
        self.is_non_image = None
        self.failure_reason = None

    @abc.abstractmethod
    def load_data(self):
        super(MedImgReader, self).load_data()

    @property
    def nims_metadata_status(self):
        return self.metadata_status

    @property
    def nims_group_id(self):
        return self.group_name

    @property
    def nims_project(self):
        return self.project_name

    @property
    def nims_session_id(self):
        return self.exam_uid

    @property
    def nims_session_label(self):
        return self.study_datetime and self.study_datetime.strftime('%Y-%m-%d %H:%M')

    @property
    def nims_session_subject(self):
        return self.subj_code

    @property
    def nims_acquisition_id(self):
        if self.acquisition_id:  # as in pfile json header
            return self.acquisition_id
        series_uid = self.series_uid
        if self.is_screenshot:
            front, back = self.series_uid.rsplit('.', 1)
            series_uid = front + '.' + str(int(back) - 1)  # don't change real series_uid
        return series_uid + ('_' + str(self.acq_no) if self.acq_no is not None else '')

    @property
    def nims_acquisition_label(self):
        return '%d.%d' % (self.series_no, self.acq_no) if self.acq_no is not None else str(self.series_no)

    @property
    def nims_acquisition_description(self):
        return self.series_desc

    @property
    def nims_file_name(self):
        if self.acquisition_id:  # as in pfile json header
            return self.acquisition_id + '_' + self.filetype
        return self.series_uid + ('_' + str(self.acq_no) if self.acq_no is not None else '') + '_' + self.filetype

    @property
    def nims_file_ext(self):
        return '.tgz'

    @property
    def nims_file_domain(self):
        return self.domain

    @property
    def nims_file_type(self):
        return self.filetype

    @property
    def nims_file_kinds(self):
        # this really SHOULDn't be scan_type, because not all medical images will set a scan type
        # pick a more general name. that is suitable for non-scan type
        # or define this property in EVERY class...
        return [self.scan_type]

    @property
    def nims_file_state(self):
        return self.state

    @property
    def nims_timestamp(self):
        return super(MedImgReader, self).nims_timestamp

    @property
    def nims_timezone(self):
        return super(MedImgReader, self).nims_timezone

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


class MedImgWriter(data.Writer):

    """
    Base MR data writer class.

    Cannot be instantiated.

    """

    __metaclass__ = abc.ABCMeta

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
        log.debug('reorienting to voxel order %s' % voxel_order)
        new_data, new_qto_xyz, _, _ = dcmstack.reorder_voxels(imagedata, qto_xyz, voxel_order)
        return new_data, new_qto_xyz

    @data.abstractclassmethod
    def write(cls, metadata, imagedata, outbase, voxel_order=None):
        super(MedImgWriter, cls).write(metadata, imagedata, outbase)
