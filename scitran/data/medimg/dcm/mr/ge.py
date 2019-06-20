# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.mr.ge
=========================

load all data from a set of GE Dicoms, almost all information is available via non-private tags.

Cannot be instantiated.  This object is merely a container for GE specific processing functions.

GE saves screenshot of the graphical prescription for each scan, the image contains some useful metadata.
therefore, any missing values should remain MISSING, regardless if they are needed in calculations, to
attempt to ensure correctness.

"""

import logging
import dcmstack
import numpy as np

import mr
from .. import dcm

log = logging.getLogger(__name__)


GEMS_TYPE_ORIG = ['ORIGINAL', 'PRIMARY', 'OTHER']
GEMS_TYPE_DERIVED_RFMT = ['DERIVED', 'SECONDARY', 'REFORMATTED', 'AVERAGE']
TAG_RECON_FLAG = (0x0043, 0x107d)
TAG_BVALUE = (0x0043, 0x1039)                                       # CSA_BVALUE = 'Slop_int_6...Slop_int_9'
TAG_BVEC = [(0x0019, 0x10bb), (0x0019, 0x10bc), (0x0019, 0x10bd)]   # CSA_BVEC = ['UserData20', 'UserData21', 'UserData22']
MAX_LOC_DCMS = dcm.MAX_LOC_DCMS
MetaExtractor = dcm.MetaExtractor
DicomError = dcm.DicomError

def infer_psd_type(self):
    """
    Infer the psd type based on the manufacturer and psd name.

    This heuristic is entirely based on parsing the name.

    Parameters
    ----------
    self : NIMSDicom instance
        uses self.psd_type

    Returns
    -------
    None : NoneType
        sets self.psd_type

    """
    if not self.psd_name:
        self.psd_type = 'unknown'
    else:
        if 'service' in self.psd_name:
            self.psd_type = 'service'
        elif self.psd_name == 'sprt':
            self.psd_type = 'spiral'
        elif self.psd_name == 'sprl_hos':
            self.psd_type = 'hoshim'
        elif self.psd_name == 'basic':
            self.psd_type = 'basic'
        elif 'mux' in self.psd_name or 'mb_' in self.psd_name:
            self.psd_type = 'muxepi'
        elif 'epi' in self.psd_name:
            self.psd_type = 'epi'
        elif self.psd_name in ['probe-mega', 'gaba_ss_cni', 'gaba_special']:
            self.psd_type = 'mrs'
        elif self.psd_name.startswith('special_siam') or self.psd_name.startswith('mega_special'):
            self.psd_type = 'mrs'
        elif self.psd_name == 'asl':
            self.psd_type = 'asl'
        elif self.psd_name in ['bravo', '3dgrass']:
            self.psd_type = 'spgr'
        elif 'fgre' in self.psd_name:        # also want to catch efgre3d
            self.psd_type = 'gre'
        elif self.psd_name == 'ssfse':
            self.psd_type = 'fse'
        elif self.psd_name == 'cube':
            self.psd_type = 'cube'
        elif self.psd_name.endswith('b1map'):  # XXX general enough?
            self.psd_type = 'fieldmap'         # XXX is this correct "type"?
        else:
            self.psd_type = 'unknown'
    log.debug(self.psd_type)


def parse_one(self):
    """
    Composer function, parses all metadata that can be parsed from a single dicom.

    Called by NIMSData init, if dicom manufacturer is GE Medical Sytems.

    """

    mr.parse_standard_mr_tags(self)
    self.psd_name = self.getelem(self._hdr, 'PulseSequenceName', str, '').lower()
    self.psd_iname = self.getelem(self._hdr, 'InternalPulseSequenceName')
    self.fov_x, self.fov_y = 2 * [self.getelem(self._hdr, 'ReconstructionDiameter', float)]
    self.receive_coil_name = self.getelem(self._hdr, 'ReceiveCoilName')
    self.mt_offset_hz = self.getelem(self._hdr, 'OffsetFrequency', float)
    effective_echo_spacing = self.getelem(self._hdr, 'EffectiveEchoSpacing', float)
    self.effective_echo_spacing = effective_echo_spacing / 1e6 if effective_echo_spacing else None
    asset_r = self.getelem(self._hdr, 'AssetRFactors', None, [None, None])
    if isinstance(asset_r, unicode) and '\\' in asset_r:    # GE Signa HDxt stores asset as string '1\1'
        asset_r = map(int, asset_r.split('\\'))             # reformat to [1, 1] for consistency
    elif isinstance(asset_r, float):                        # asset_r can be single item float
        asset_r = [None, None]
    self.phase_encode_undersample, self.slice_encode_undersample = asset_r
    # some very old Ge systems will output dicoms that don't define Locations in Acquition, or define it in a way
    # that is weird.  It may incorrectly label the value type as OB, but not be able to translate the value, resulting
    # in the MetaExtractor excluding it from the it's output metadata.
    self.num_slices = self.getelem(self._hdr, 'LocationsInAcquisition', int)
    self.total_num_slices = self.getelem(self._hdr, 'ImagesInAcquisition', int)
    self.num_timepoints = self.getelem(self._hdr, 'NumberOfTemporalPositions', int)

    # slice check could end up wrong, if both total_num_slices and num_slices are None
    # could force num_slices and total_num_slices into different ORs, to prevent matching if both are None
    # thus only when they are both defined, AND not equal, can this test pass
    if (self.total_num_slices or 1) == (self.num_slices or 0):
        self.total_num_slices = (self.num_slices or 1) * (self.num_timepoints or 1)
        log.debug('adjusted total_num_slices from %3d to %3d' % (self.num_slices, self.total_num_slices))  # num_slices == 'old' total_num

    # some localizer don't have header field to indicate the number of slices
    # per acquisition.  If the total_number of slices is set, and the num_timepoints is 1
    # then the number of slices should be equal to total number of slices
    if not self.num_slices and (self.num_timepoints or 1) == 1:
        self.num_slices = self.total_num_slices

    prescribed_duration = (self.tr or 0) * (self.num_timepoints or 0) * (self.num_averages or 1)  # FIXME: only works for fMRI, not anatomical
    if prescribed_duration != 0:
        self.prescribed_duration = prescribed_duration
        self.duration = prescribed_duration
    else:
        self.prescribed_duration = None
        self.duration = None

    dwi_dirs = self.getelem(self._hdr, 'UserData24{#DTIDiffusionDir.,Release10.0&Above}', float)
    self.dwi_dirs = int(dwi_dirs) if dwi_dirs else None
    if self.image_type == GEMS_TYPE_ORIG and (self.dwi_dirs or 0) >= 6:
        self.is_dwi = True
        self.num_timepoints = 1

    if self.image_type == GEMS_TYPE_DERIVED_RFMT:
        self.is_non_image = True

    infer_psd_type(self)
    mr.adjust_fov_acqmat(self)
    mr.infer_scan_type(self)


def parse_all(self):
    """
    Parse all metadata that requires all dicoms.

    Called by NIMSDicom load_data, if dicom manufacturer is GE Medical System.

    """
    if self.total_num_slices < MAX_LOC_DCMS:
        slice_norms = [np.cross(np.matrix(d.get('ImageOrientationPatient')[0:3]), np.matrix(d.get('ImageOrientationPatient')[3:6]))[0] for d in self._dcm_list]
        norm_diff = [np.abs(np.dot(slice_norms[0], n)).round(2) for n in slice_norms]
        self.is_localizer = bool(len(set(norm_diff)) > 1)

    if self.is_dwi:
        # DTI scans could have 1+ non-DTI volume. num vols will be >= dwi_dirs + 1
        self.bvals = np.array([float(self.getelem(d, TAG_BVALUE)[0]) for d in self._dcm_list[0::self.num_slices]])
        self.bvecs = np.array([[self.getelem(d, TAG_BVEC[i], float) for i in range(3)] for d in self._dcm_list[0::self.num_slices]]).transpose()

    # TODO separatation of concern; identification vs dicom grouping
    recon_mode_flag = np.unique([self.getelem(d, TAG_RECON_FLAG, int, 0) for d in self._dcm_list])
    if self.psd_type not in ['fieldmap']:
        log.debug('recon_mode: 1, might be multicoil')
        vol_counter = 0
        ref_position = self.getelem(self._hdr, 'ImagePositionPatient')
        for d in self._dcm_list:
            if np.allclose(self.getelem(d, 'ImagePositionPatient'), ref_position):
                vol_counter += 1
        log.debug('found %d volumes; expected %d' % (vol_counter, (self.num_timepoints or 1) * (self.num_echos or 1)))
        if vol_counter > (self.num_timepoints or 1) * (self.num_echos or 1):
            self.is_multicoil = True
            self.num_receivers = (self.total_num_slices / self.num_slices) - 1   # actual #recv = -1 of num volumes
            self._dcm_groups = [self._dcm_list[x::self.num_receivers + 1] for x in xrange(0, self.num_receivers + 1)]
            log.debug('groups: %3d; %3d coils + 1 combined' % (len(self._dcm_groups), self.num_receivers))

    # attempt to calculate trigger times and slice duration, if the first dicom reports trigger time
    self.slice_duration = None
    if self.total_num_slices >= self.num_slices and self.getelem(self._dcm_list[0], 'TriggerTime', float) is not None:
        log.debug('using trigger times to calculate slice order and slice duration')
        trigger_times = np.array([self.getelem(d, 'TriggerTime', float) for d in self._dcm_list[0:self.num_slices]])
        if self.reverse_slice_order:
            trigger_times = trigger_times[::-1]
        # if trigger_times is mixed float and None, then dtype = 'O'
        if self.num_slices > 2 and len(trigger_times) > 2 and trigger_times.dtype == 'float64':
            trigger_times_from_first_slice = trigger_times[0] - trigger_times
            self.slice_duration = float(min(abs(trigger_times_from_first_slice[1:]))) / 1000    # msec to sec
            if trigger_times_from_first_slice[1] < 0:
                self.slice_order = mr.SLICE_ORDER_SEQ_INC if trigger_times[2] > trigger_times[1] else mr.SLICE_ORDER_ALT_INC
            else:
                self.slice_order = mr.SLICE_ORDER_ALT_DEC if trigger_times[2] > trigger_times[1] else mr.SLICE_ORDER_SEQ_DEC
        else:
            self.slice_duration = trigger_times[0]
            self.slice_order = mr.SLICE_ORDER_SEQ_INC

    # identify fastcard by looking at number of cardiac images as timepoints, with 5 dicoms per slice
    cardiac_images = self.getelem(self._hdr, 'CardiacNumberOfImages', int, 0)
    if self.total_num_slices == (cardiac_images * self.num_slices * 5):  # fastcard: 5 value per temporal X spatial location
        self.num_timepoints = cardiac_images
        self.is_fastcard = True
        self.velocity_encode_scale = self.getelem(self._hdr, 'VelocityEncodeScale', float)
        self.velocity_encoding = self.getelem(self._hdr, 'VelocityEncoding', int)

    log.debug(self.psd_name)
    log.debug(self.psd_type)

def fastcard_convert(self):
    """GE fast card conversion."""
    log.debug('fast card')
    stacks = []
    group_id = 0
    def _split_list(l, size):
        return [l[i:i+size] for i in range(0, len(l), size)]
    dcm_groups = _split_list(self._dcm_list, self.total_num_slices / 5)

    for group in dcm_groups:
        group_id += 1
        num_positions = len(set([d.SliceLocation for d in group]))
        if num_positions != self.num_slices:
            raise DicomError('volume %s has %s unique positions; expected %s' % (group_id, num_positions, self.num_slices))
        stack = dcmstack.DicomStack()
        for d in group:
            meta = MetaExtractor(d)
            stack.add_dcm(d, meta)
        nii_wrp = stack.to_nifti_wrapper()
        stacks.append(nii_wrp)
    try:
        nii_wrp = dcmstack.dcmmeta.NiftiWrapper.from_sequence(stacks)
    except dcmstack.InvalidStackError as e:
        raise DicomError('cannot reconstruct %s: %s' % (self.filepath, str(e)))
    nii = nii_wrp.nii_img
    self.data = {'': nii.get_data()}
    self.qto_xyz = nii.get_affine()
    mr.post_convert(self)

def multicoil_convert(self):
    """GE specific multicoil converesion. requires partial volume and missing slices check."""
    log.debug('multicoil recon')
    mr.partial_vol_check(self)

    stacks = []
    group_id = 0
    for group in self._dcm_groups:
        group_id += 1
        log.debug('multicoil - %2s, %s dicom' % (str(group_id), str(len(group))))
        num_positions = len(set([d.SliceLocation for d in group]))
        if num_positions != self.num_slices:
            raise DicomError('coil %s has %s unique positions; expected %s' % (group_id, num_positions, self.num_slices))
        stack = dcmstack.DicomStack()
        for dcm in group:
            meta = MetaExtractor(dcm)
            stack.add_dcm(dcm, meta)
        nii_wrp = stack.to_nifti_wrapper()
        stacks.append(nii_wrp)
    try:
        nii_wrp = dcmstack.dcmmeta.NiftiWrapper.from_sequence(stacks)
    except dcmstack.InvalidStackError as e:
        raise DicomError('cannot reconstruct %s: %s' % (self.filepath, e))      # XXX FAIL! unexpected for recon to fail
        # raise NIMSDicomError('cannot reconstruct %s: %s' % (self.filepath, e), log_level=logging.ERROR)
    del self._dcm_groups, self._dcm_list, stacks, stack, dcm
    nii = nii_wrp.nii_img
    self.data = {'': nii.get_data()}
    self.qto_xyz = nii.get_affine()
    del nii_wrp, nii

    mr.post_convert(self)

def convert(self):
    """
    Composer function, determines which convert function to use.

    Called by NIMSDicom load_data if dicom manufacturer is GE Medical Systems.

    """
    if self.is_non_image:
        mr.non_image_handler(self)
    elif self.is_localizer:
        mr.localizer_convert(self)
    elif self.is_multicoil:
        multicoil_convert(self)
    elif self.is_fastcard:
        fastcard_convert(self)
    else:
        mr.standard_convert(self)
