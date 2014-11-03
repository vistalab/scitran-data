# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.mr.generic_mr
=================================

Generic MR Dicom functions to use for composing the NIMSDicom class.

These functions are meant to be imported by and used within
dcm.mr.ge or dcm.mr.siemens.

note: direct usage of these functions is an "advanced" topic.
"""

import logging
import dcmstack
import numpy as np

from ... import medimg
from ... import nimsdicom


log = logging.getLogger(__name__)


NIMSDicomError = nimsdicom.NIMSDicomError
MetaExtractor = nimsdicom.MetaExtractor
MAX_LOC_DCMS = nimsdicom.MAX_LOC_DCMS

# NIFITI1-style slice order codes:
SLICE_ORDER_UNKNOWN = 0
SLICE_ORDER_SEQ_INC = 1
SLICE_ORDER_SEQ_DEC = 2
SLICE_ORDER_ALT_INC = 3
SLICE_ORDER_ALT_DEC = 4
SLICE_ORDER_ALT_INC2 = 5  # interleaved, increased, starting at 2nd MRI slice
SLICE_ORDER_ALT_DEC2 = 6  # interleave, decreasing, starting at one before last MRI slice

# TODO scan types should be an external json file
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
    if vendor.lower().startswith('ge') and rotation != None:
        log.debug('rotating bvecs with image orientation matrix')
        bvecs,bvals = rotate_bvecs(bvecs, bvals, rotation)
    else:
        bvecs,bvals = rotate_bvecs(bvecs, bvals, np.diag((-1.,-1.,1.)))
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
            # log.debug('there are no decimals, cannot intelligently scale bvecs/bvals')
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


# TODO: infer_scan_type should be broken down by manufacturer for maintainability
def infer_scan_type(self):
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
    psd_type = self.psd_type
    num_timepoints = self.num_timepoints
    te = self.te
    fov = self.fov
    mm_per_vox = self.mm_per_vox
    is_dwi = self.is_dwi
    is_localizer = self.is_localizer

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
    # return scan_type
    log.debug(scan_type)
    self.scan_type = scan_type


def parse_standard_mr_tags(self):
    """
    Parse MR specific dicom tags.

    Given a dataset that contains parsed dicom header metadata
    in dataset._hdr, this will parse standard MR dicom header tags,
    and then updates the dataset metadata based on
    the parsed header metadata.

    Parameters
    ----------
    dataset : NIMSDicom

    Returns
    -------
    None : None

    """
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
    self.num_averages = self.getelem(self._hdr, 'NumberOfAverages', int)    # some GE scans show this as 0.5?
    self.num_echos = self.getelem(self._hdr, 'EchoNumbers', int)
    self.protocol_name = self.getelem(self._hdr, 'ProtocolName')
    self.acquisition_type = self.getelem(self._hdr, 'MRAcquisitionType')
    self.mm_per_vox_x, self.mm_per_vox_y = self.getelem(self._hdr, 'PixelSpacing', float) or (None, None)
    self.mm_per_vox_z = self.getelem(self._hdr, 'SpacingBetweenSlices', float) or self.getelem(self._hdr, 'SliceThickness', float)
    self.size_x = self.getelem(self._hdr, 'Columns', int)
    self.size_y = self.getelem(self._hdr, 'Rows', int)
    self.reverse_slice_order = False
    self.slice_order = medimg.SLICE_ORDER_UNKNOWN
    self.total_num_slices = None
    self.is_multiecho = None  # GE only, for now
    self.is_multicoil = None  # GE only
    self.is_non_image = None                            # does this ever still happen?
    self.mt_offset_hz = None                            # XXX: expected by nimsnifti from ALL mr dicoms. is GE specific


def adjust_fov_acqmat(self):
    """
    Parses acquisition matrix from dicom header.

    adjust FOV and acq_mat based on encoding direction
    and percent field of view.

    should be called at the end of parse_one.
    """
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
        if self.fov != (None, None):
            fov_x, fov_y = self.fov
            fov_x /= (self.getelem(self._hdr, 'PercentPhaseFieldOfView', float, 0.) / 100.) if 'PercentPhaseFieldOfView' in self._hdr else 1.
            self.fov = fov_x, fov_y


def non_image_handler(self):  # don't attempt recon.
    """
    Handle non-image mr dicoms.

    Currently, all non-image MR dicoms fail silently.
    This is sort of a placeholder until there
    is a better way to deal with varieties of non-image mr dicoms.

    Parameters
    ----------
    dataset : NIMSDicom

    Returns
    -------
    None : None

    """
    log.error('%s is non-image!?' % (self.filepath))
    # self.data = None      # implied.


def localizer_convert(self):
    """
    Reconstruct localizer from GE MR dicoms.

    Parameters
    ----------
    dataset : NIMSDicom

    Returns
    -------
    None : None

    """

    log.debug('localizer recon')
    num_ornts = len(set([tuple(d.get('ImageOrientationPatient', [0.] * 6)) for d in self._dcm_list]))
    self.num_timepoints = num_ornts
    self.num_slices = self.total_num_slices / self.num_timepoints
    # determine cosines, slice norm and rotation
    cosines = self.getelem(self._hdr, 'ImageOrientationPatient', float, 6 * [np.nan])
    row_cosines = np.array(cosines[0:3])
    col_cosines = np.array(cosines[3:6])
    slice_norm = np.cross(row_cosines, col_cosines)
    rot = compute_rotation(row_cosines, col_cosines, slice_norm)
    # determine origin
    image_position = [tuple(self.getelem(d, 'ImagePositionPatient', float, [0., 0., 0.])) for d in self._dcm_list]
    origin = image_position[0] * np.array([-1, -1, 1])
    self.qto_xyz = build_affine(rot, self.mm_per_vox, origin)
    # recon
    self.data = np.dstack([np.swapaxes(d.pixel_array, 0, 1) for d in self._dcm_list])
    dims = np.array((self.size_y, self.size_x, self.num_slices, self.num_timepoints))
    if np.prod(dims) == np.size(self.data):
        self.data = self.data.reshape(dims, order='F')
    self.data = {'': self.data}
    # XXX: acq matrix reversed for some localizers?


def partial_vol_check(self):        # AKA _pre_convert.
    """
    If not mosaic, AND know number of slices, then check if we have one complete volume.

    check that there are no dangling dicoms from an incomplete volume

    Parameters
    ----------
    dataset : NIMSDicom

    """
    if 'MOSAIC' not in self.image_type and self.num_slices:
        if len(self._dcm_list) < self.num_slices:                   # check if have at least 1 whole volume
            raise NIMSDicomError('cannot reconstruct. total dicoms < num slices')
            # TODO: add comment to notes
        partial_vol_dcms = len(self._dcm_list) % self.num_slices    # remove partial volumes
        if partial_vol_dcms:
            log.debug('number of dicoms is not a integer multiple of number of unique slices positions. trimming.')
            self._dcm_list = self._dcm_list[:-1 * partial_vol_dcms]
            # TODO: add comment notes


def post_convert(self):
    """Some things must happen after the standard reconstruction is complete."""
    # XXX: is this a GE specific step?
    if self.is_dwi:
        self.bvecs, self.bvals = adjust_bvecs(self.bvecs, self.bvals, self.scanner_type, self.qto_xyz[0:3, 0:3])


def standard_convert(self):
    """Standard reconstruction, works for most manufacturers and scans."""
    log.debug('standard recon')
    partial_vol_check(self)

    stack = dcmstack.DicomStack()
    for dcm in self._dcm_list:
        stack.add_dcm(dcm, MetaExtractor(dcm))
    try:
        nii_wrp = stack.to_nifti_wrapper()
    except dcmstack.InvalidStackError as e:
        raise NIMSDicomError('cannot reconstruct %s: %s' % (self.filepath, e))      # XXX FAIL! unexpect for recon to fail
    del self._dcm_list, stack, dcm

    nii = nii_wrp.nii_img
    self.data = {'': nii.get_data()}
    self.qto_xyz = nii.get_affine()
    del nii_wrp, nii

    post_convert(self)
