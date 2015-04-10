# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.mr.siemens
==============================

Load all data from a set of siemens dicoms.

Not all information required is available via non-private tags.

Cannot be instantiated.  This is merely a container for Siemens specific functions.

however, the CSA Series Header and CSA Image Header are stored in non-private tag groups (even group #).
csa header tags are preferred over standard dicom tags. Although some information overlaps between
the dicom standard tags and CSA headers, patient information is not included in the CSA.  Therefore,
anonymization doesn't need to fiddle with the CSA headers.

"""

import logging
import numpy as np

import mr
from .. import dcm

log = logging.getLogger(__name__)

# SIEMENS_TYPE_CSAPARALLEL = ['DERIVED', 'SECONDARY', 'MPR', 'CSA MPR', '', 'CSAPARALLEL', 'M', 'ND'] #
# SIEMENS_TYPE_POSDISP = ['DERIVED', 'SECONDARY', 'POSDISP', 'M', 'ND', 'NORM', 'CSA RESAMPLED']   # posdisp = position display
SIEMENS_TYPE_DIS2D = ['ORIGINAL', 'PRIMARY', 'M', 'RETRO', 'NORM', 'DIS2D', 'FM4_2', 'FIL']
# dis2d = distorted pixel and remapped, 2d distortion correct, retro = Retro image, or retrospective gating

TAG_BVALUE = 'CsaSeries.MrPhoenixProtocol.sDiffusion.alBValue[1]'  # B_value
TAG_BVEC = 'CsaImage.DiffusionGradientDirection'
MAX_LOC_DCMS = dcm.MAX_LOC_DCMS
MetaExtractor = dcm.MetaExtractor

def infer_psd_type(self):
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
    None : NoneType

    """
    if not self.psd_name:
        self.psd_type = 'unknown'
    else:
        if self.psd_name == 'siemensseq%\tse_vfl':
            self.psd_type = 'tse'
        elif self.psd_name == 'siemensseq%\ep2d_diff':
            self.psd_type =  'epi'
        elif self.psd_name == 'siemensseq%\ep2d_bold':
            self.psd_type =  'epi'
        elif self.psd_name == 'siemensseq%\ep2d_asl':
            self.psd_type =  'asl'
        elif self.psd_name == 'siemensseq%\gre':
            self.psd_type =  'gre'
        elif self.psd_name == 'siemensseq%\tfl':
            self.psd_type =  'tfl'
        elif self.psd_name == 'siemensseq%\gre_field_mapping':
            self.psd_type =  'gre'
        elif self.psd_name == 'serviceseq%\rf_noise':
            self.psd_type =  'service'
        elif self.psd_name.startswith('customerseq%\ep2d_pasl'):
            self.psd_type = 'asl'
        elif self.psd_name.startswith('customerseq%\ep2d_diff'):
            self.psd_type = 'epi'
        elif self.psd_name.startswith('customerseq%\ep2d'):
            self.psd_type = 'epi'
        elif self.psd_name == 'customerseq%\wip711_moco\tfl_multiecho_epinav_711':
            self.psd_type = 'tfl'
        else:
            self.psd_type = 'unknown'
    log.debug(self.psd_type)


def parse_one(self):
    """
    Composer function, parses all metadata that can be parsed from a single dicom.

    Called by NIMSDicom init, if dicom manufacturer is Siemens.

    """
    mr.parse_standard_mr_tags(self)

    self.psd_name = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.tSequenceFileName', str, '').lower().replace('%', '', 1)
    self.psd_iname = self.getelem(self._hdr, 'SeriesDescription')
    self.fov_x = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.asSlice[0].dPhaseFOV', float)    # XXX fov-x = PhaseFOV?
    self.fov_y = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.asSlice[0].dReadoutFOV', float)  # XXX fov-y = ReadoutFOV?
    self.receive_coil_name = self.getelem(self._hdr, 'CsaImage.ImaCoilString')
    slice_duration = self.getelem(self._hdr, 'CsaImage.SliceMeasurementDuration', float, 0.)
    self.slice_duration = slice_duration / 1e6 if slice_duration else None
    self.prescribed_duration = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.lScanTimeSec')   # FIXME
    self.duration = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.lTotalScanTimeSec')         # FIXME: not guaranteed
    self.acq_no = None      # siemens acq # indicates the brain volume instance. varies within one scan.

    self.dwi_dirs = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.sDiffusion.lDiffDirections', int, None)
    if (self.dwi_dirs or 0) > 1:
        self.is_dwi = True
        self.num_timepoints = 1

    # some siemens MR dicoms are not reconstructable
    if 'CSAPARALLEL' in self.image_type:
        self.is_non_image = True
    if 'POSDISP' in self.image_type:
        self.is_non_image = True
    if self.image_type == SIEMENS_TYPE_DIS2D:
        # had 2 image orientations, and other metdata that varied between dicoms, and is not a localizer. AFNI cannot reconstsruct.
        # this is a specific hack fix. A more general fix is to check each dicom to see if the orientation matches the first dicom
        # but checking EVERY dicom isn't very ideal solution (some datasets might have an ENORMOUS number of dicoms)
        self.is_non_image = True

    infer_psd_type(self)
    mr.adjust_fov_acqmat(self)
    mr.infer_scan_type(self)


def parse_all(self):
    """
    Composer function, parses all metadata that requires all dicoms.

    Called by NIMSDicom load_data, if dicom manufacturer is Siemens.

    """
    if 'MOSAIC' not in self.image_type:
        log.debug('SIEMENS SINGLE SLICE DICOM')
        self.total_num_slices = len(self._dcm_list)
        # num slices from CSA
        self.num_slices = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.lSize', int,
                                       self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.sGroupArray.asGroup[0].nSize', int, None)
                                       )
        self.num_timepoints = self.total_num_slices / self.num_slices if self.num_slices else None
    elif 'MOSAIC' in self.image_type:           # explicit for readability
        log.debug('SIEMENS MOSAIC')
        self.num_slices = self.getelem(self._hdr, 'CsaImage.NumberOfImagesInMosaic', int, self.getelem(self._hdr, 'NumberOfImagesInMosaic', int))
        self.num_timepoints = len(self._dcm_list)
        mosaic_dim = int(self.num_slices ** 0.5)
        if mosaic_dim ** 2 < self.num_slices:
            mosaic_dim += 1
        self.size = [x / mosaic_dim for x in self.size]
        self.fov = [x / mosaic_dim for x in self.fov]
        self.total_num_slices = self.num_slices * self.num_timepoints

    log.debug('num slices / vol: %s' % str(self.num_slices))  # stringify to be able to log NoneType

    self.duration = self.num_timepoints * (self.num_averages or 1) * self.tr if self.num_timepoints and self.tr else None

    # CsaSeries.MrPhoenixProtocol.sSliceArray.ucMode indicates siemens slice order:
    # Siemens; 1 Ascending, 2 Descending, and 4 Interleaved Ascending.
    # NIFTI;   1 Ascending, 2 Descending, and 4 Interleaving Descending
    # Siemens slice order 4 could be nifti slice order 3 or 5 depending on num_slices
    # - nifti slice order 3: interleave_asc, odd first, odd num slices, interleave INC
    # - nifti slice order 5: interleave_asc, even first, even num slices, interleave INC 2
    self.slice_order = self.getelem(self._hdr, 'CsaSeries.MrPhoenixProtocol.sSliceArray.ucMode', None, 0)
    if self.slice_order == 4:   # don't try to guess if num_slices can't be determined
        if self.num_slices % 2 != 0:
            self.slice_order = 3  # interleaved ascending, odd first
        else:
            self.slice_order = 5  # interleaved ascending, even first

    self.num_receivers = len([self._hdr[key] for key in self._hdr if key.endswith('sCoilElementID.tCoilID')])

    if self.total_num_slices < MAX_LOC_DCMS:
        slice_norms = [np.cross(np.matrix(d.get('ImageOrientationPatient')[0:3]), np.matrix(d.get('ImageOrientationPatient')[3:6]))[0] for d in self._dcm_list]
        norm_diff = [np.abs(np.dot(slice_norms[0], n)).round(2) for n in slice_norms]
        self.is_localizer = bool(len(set(norm_diff)) > 1)

    if self.is_dwi:
        self.bvals = np.array([MetaExtractor(d).get(TAG_BVALUE, 0.) for d in self._dcm_list[0:self.num_slices]])
        self.bvecs = np.array([MetaExtractor(d).get(TAG_BVEC, [0., 0., 0.]) for d in self._dcm_list[0:self.num_slices]]).transpose()


def convert(self):
    """
    Composer function, determines which convert function to use.

    Called by NIMSDicom load_data if dicom manufacturer is Siemens.

    """
    if self.is_non_image:
        mr.non_image_handler(self)
    elif self.is_localizer:
        mr.localizer_convert(self)
    else:
        mr.standard_convert(self)
