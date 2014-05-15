# @author: Kevin S Hahn
#

"""
Siemens Specific Implementation of NIMSDicom.

This module and it's classes focus on Siemens Specific Implementations of
private tag parsing, CSA header parsing, and computation of values that will be stored in the header
"""

import re
import logging
import dcmstack
import dcmstack.extract
import numpy as np

import nimspng
import nimsmrdata
import nimsdicom

log = logging.getLogger('nimsdicom.siemens')

# Siemens Private Tags
TAG_IMAGES_IN_MOSAIC =          (0x0019,0x100a)     #noqa
TAG_BVAL =                      (0x0019,0x100c)     #noqa
TAG_BVEC =                      (0x0019,0x100e)     #noqa
TAG_BMAT =                      (0x0019,0x1027)     #noqa
TAG_FOV =                       (0x0051,0x100c)     #noqa
TAG_TA =                        (0x0051,0x100a)     # not reliable across scan types #noqa
TAG_SLICE_DUR =                 (0x0019,0x100b)     # /1e6  #noqa
TAG_TIME_AFTER_START =          (0x0019,0x1016)     #noqa
TAG_COIL_ID =                   (0x0051,0x100f)     #noqa

# Siemens CSA labels
CSA_PROTOCOL_NAME =             'CsaSeries.MrPhoenixProtocol.tProtocolName'
CSA_DIFFUSION_DIRS =            'CsaSeries.MrPhoenixProtocol.sDiffusion.lDiffDirections'
CSA_SLICE_DUR =                 'CsaImage.SliceMeasurementDuration'     # /1e6
CSA_COIL_ID =                   'sCoilElementID.tCoilID'    # used in a 'string'.endswith(CSA_COIL_ID)
CSA_SCAN_TIME =                 'CsaSeries.MrPhoenixProtocol.lScanTimeSec'
CSA_TOTAL_SCAN_TIME =           'CsaSeries.MrPhoenixProtocol.lTotalScanTimeSec'
CSA_PARTITIONS =                'CsaSeries.MrPhoenixProtocol.sKSpace.lPartitions'
CSA_IMAGE_PER_SLAB =            'CsaSeries.MrPhoenixProtocol.sKSpace.lImagesPerSlab'
CSA_TAG_AFTER_START =           'CsaImage.TimeAfterStart'
CSA_IMAGES_IN_MOSAIC =          'CsaImage.NumberOfImagesInMosaic'
CSA_BVAL =                      'CsaSeries.MrPhoenixProtocol.sDiffusion.alBValue[1]'
CSA_BVEC =                      'CsaImage.DiffusionGradientDirection'

# maps siemens slice code to appropriate nifti slice code
SIEMENS_SLICE_CODE = {0: nimsmrdata.SLICE_ORDER_UNKNOWN,
                      1: nimsmrdata.SLICE_ORDER_SEQ_INC,
                      2: nimsmrdata.SLICE_ORDER_SEQ_DEC,
                      4: nimsmrdata.SLICE_ORDER_ALT_INC,
                      }


class NIMSDicomSiemensError(nimsdicom.NIMSDicomError):
    pass


class NIMSDicomSiemens(nimsdicom.NIMSDicom):

    """NIMSDicom Siemens specific parser"""

    def __init__(self, archive):
        super(NIMSDicomSiemens, self).__init__(archive)                         # self.metadata._hdr created
        if self.metadata.manufacturer != 'SIEMENS':
            raise NIMSDicomSiemensError('NIMSSiemensDicom received non-Siemens data. bailing out.')
        log.debug('SIEMENS')
        csa = dcmstack.extract.MetaExtractor()(self.metadata._hdr)

        self.metadata.psd_name = csa.get('CsaSeries.MrPhoenixProtocol.tSequenceFileName', 'unknown\\unknown').split('\\')[1]
        # self.metadata.psd_iname = csa.get(CSA_PROTOCOL_NAME, 'unknown')
        self.metadata.slice_duration = csa.get(CSA_SLICE_DUR, 0.) / 1e6
        self.metadata.fov = [int(x) for x in self.getelem(self.metadata._hdr, TAG_FOV, None, 'FoV 0*0').split()[1].split('*')]
        receive_channels = [csa[key] for key in csa if key.endswith('sCoilElementID.tCoilID')]
        self.metadata.receive_coil_name = csa.get('CsaImage.ImaCoilString')
        self.metadata.num_receivers = len(receive_channels)
        self.metadata.is_dwi = bool('ORIGINAL' in self.metadata.image_type and csa.get(CSA_DIFFUSION_DIRS, 0) >= 6)     # TODO: distinguish derived dti
        self.metadata.slice_order = SIEMENS_SLICE_CODE.get(csa.get('CsaSeries.MrPhoenixProtocol.sSliceArray.ucMode', 0), 0)

        if 'MOSAIC' in self.metadata.image_type:  # DTI or TimeSeries, where timepoints >1 or Value/pix >1
            log.debug('MOSAIC')
            self.metadata.num_slices = csa.get(CSA_IMAGES_IN_MOSAIC, 1)
            mosaic_dim = int(self.metadata.num_slices ** 0.5)
            if mosaic_dim ** 2 < self.metadata.num_slices:
                mosaic_dim += 1
            self.metadata.size = [x / mosaic_dim for x in self.metadata.size]
            self.metadata.fov = [x / mosaic_dim for x in self.metadata.fov]
            self.metadata.num_timepoints = 0                                    # TODO: fix num_timepoint estimation?
        else:
            # Anatomical scans, where num_timepoints=1 and value/pixel = 1. also multicoil
            log.debug('SINGLE IMAGE')
            self.metadata.num_timepoints = 1
            self.metadata.num_slices = csa.get(CSA_PARTITIONS, 0)               # TODO: is CSA_PARTITIONS 100% correct? is there a better alt?

        self.metadata.total_num_slices = self.metadata.num_slices * self.metadata.num_timepoints * self.metadata.num_averages
        self.metadata.prescribed_dur = csa.get(CSA_SCAN_TIME, 0)                # TODO: fix prescribed duration
        self.metadata.prescribed_duration = csa.get(CSA_TOTAL_SCAN_TIME, 0)     # TODO: fix prescribed duration
        self.metadata.duration = self.metadata.prescribed_duration              # TODO: fix observed duration

        self._postinit()

    def load_data(self, archive):
        """reconstructs dicom image data and parses metadata that requires information from all dicoms"""
        # TODO: refactor, do more "type identification" things up front, then specific blocks that deal with types
        super(NIMSDicomSiemens, self).load_data(archive)                        # all dicoms as self.metadata._dcm_list
        num_dcms = len(self.metadata._dcm_list)
        coils = [self.getelem(dcm, TAG_COIL_ID) for dcm in self.metadata._dcm_list]     # some dcms, like SRe or CSA Reports, do not have coil info
        is_multi_coil = bool(len(set(coils)) > 1)

        if 'MOSAIC' in self.metadata.image_type:
            log.debug('MOSAIC loading all data')
            if self.metadata.is_dwi:
                self.metadata.num_timepoints = 1
            else:
                self.metadata.num_timepoints = num_dcms
            # self.metadata.num_slices = self.getelem(self.metadata._hdr, TAG_IMAGES_IN_MOSAIC, int, 0)
            self.metadata.total_num_slices = self.metadata.num_slices * self.metadata.num_timepoints * self.metadata.num_averages
        else:
            log.debug('SINGLE IMAGES - ALL DATA LOADED')
            self.metadata.total_num_slices = num_dcms
            self.metadata.num_slices = len(set([tuple(dcm.ImagePositionPatient) for dcm in self.metadata._dcm_list]))  # number of positions
            self.num_timepoints = self.metadata.total_num_slices / self.metadata.num_slices

            self.metadata.duration = self.metadata.num_timepoints * self.metadata.num_averages * self.metadata.tr  # TODO: fix
            self.metadata.total_num_slices = self.metadata.num_slices * self.metadata.num_timepoints * self.metadata.num_averages  # TODO: fix

        if is_multi_coil and not self.metadata.is_dwi:
            log.debug('multicoil: ' + str(is_multi_coil))
            coil_names = list(set(coils))
            coil_names.sort(key=lambda coil: [int(c) if c.isdigit() else c for c in re.split('([0-9]+)', coil)])    # natural sort coil names
            if coil_names[0][0] != coil_names[1][0] and ':' in coil_names[0]:   # TODO: better check if combined coil is first
                coil_names.append(coil_names.pop(0))
            log.debug(coil_names)

            stack_list = []
            for coil in coil_names:
                log.debug('stacking %s' % coil)
                dcms = [dcm for dcm in self.metadata._dcm_list if self.getelem(dcm, TAG_COIL_ID) == coil]
                stack = nimsdicom.NIMSStack()
                for dcm in dcms:
                    stack.add_dcm(dcm)
                nii_wrp = stack.to_nifti_wrapper()
                stack_list.append(nii_wrp)

            log.debug('combining all stacks')
            merged = dcmstack.dcmmeta.NiftiWrapper.from_sequence(stack_list)
            nii = merged.nii_img
            self.data = nii.get_data()
            self.metadata.qto_xyz = nii.get_affine()  # for raw data, is sform == qform == affine?
            self.metadata.sform = nii.get_sform()
            self.metadata.qform = nii.get_qform()
        else:
            stack = nimsdicom.NIMSStack()
            for dcm in self.metadata._dcm_list:
                stack.add_dcm(dcm)
            nii = stack.to_nifti_wrapper().nii_img

            self.data = nii.get_data()
            self.metadata.qto_xyz = nii.get_affine()
            self.metadata.qform = nii.get_qform()
            self.metadata.sform = nii.get_sform()

        if self.metadata.is_dwi:            # TODO: consider restricting this to PRIMARY DIFUSSION MOSAIC types
            # TODO: some mosaics have no direction/matrix data, I think these are usually "derived" types
            rot = self.metadata.qto_xyz[0:3, 0:3]
            self.metadata.bvals = np.array([dcmstack.extract.MetaExtractor()(dcm).get(CSA_BVAL, 0.) for dcm in self.metadata._dcm_list[0:self.metadata.num_slices]])
            self.metadata.bvecs = np.array([dcmstack.extract.MetaExtractor()(dcm).get(CSA_BVEC, [0., 0., 0.]) for dcm in self.metadata._dcm_list[0:self.metadata.num_slices]]).transpose()
            self.metadata.bvecs, self.metadata.bvals = nimsmrdata.adjust_bvecs(self.metadata.bvecs, self.metadata.bvals, self.metadata.scanner_type, rot)
