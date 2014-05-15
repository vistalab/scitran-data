# @author:  Kevin S Hahn
#
#

import os
import dicom
import logging
import dcmstack
import numpy as np

import nimspng
import nimsdicom
import nimsmrdata

log = logging.getLogger('nimsdicom.ge')

dicom.config.enforce_valid_values = False

# GE SPECIFIC TYPES
TYPE_ORIGINAL = ['ORIGINAL', 'PRIMARY', 'OTHER']
TYPE_EPI =      ['ORIGINAL', 'PRIMARY', 'EPI', 'NONE']
TYPE_SCREEN =   ['DERIVED', 'SECONDARY', 'SCREEN SAVE']

# GE SPECIFIC TAGS
TAG_PSD_NAME =                      (0x0019, 0x109c)
TAG_PSD_INAME =                     (0x0019, 0x109e)
TAG_EPI_EFFECTIVE_ECHO_SPACING =    (0x0043, 0x102c)
TAG_PHASE_ENCODE_UNDERSAMPLE =      (0x0043, 0x1083)
TAG_SLICES_PER_VOLUME =             (0x0021, 0x104f)
TAG_DIFFUSION_DIRS =                (0x0019, 0x10e0)
TAG_BVALUE =                        (0x0043, 0x1039)
TAG_BVEC =                          [(0x0019, 0x10bb), (0x0019, 0x10bc), (0x0019, 0x10bd)]
TAG_MTOFF_HZ =                      (0x0043, 0x1034)
TAG_ACQ_DUR =                       (0x0019, 0x105a)


class NIMSDicomGEError(nimsdicom.NIMSDicomError):
    pass


class NIMSDicomGE(nimsdicom.NIMSDicom):

    """NIMSDicom GE Medical System specific parser"""

    def __init__(self, archive):
        super(NIMSDicomGE, self).__init__(archive)        # takes archive, reads dicom, creates self.metadata._hdr
        if self.metadata.manufacturer != 'GE MEDICAL SYSTEMS':
            raise NIMSDicomGEError('NIMSDicomGE recived non-GE data')
        log.debug('GE MEDICAL SYSTEMS')

        self.metadata.psd_name = os.path.basename(self.getelem(self.metadata._hdr, TAG_PSD_NAME, None, 'unknown'))
        # self.metadata.psd_iname = self.getelem(self.metadata._hdr, TAG_PSD_INAME, None, 'unknown')
        self.metadata.num_slices = self.getelem(self.metadata._hdr, TAG_SLICES_PER_VOLUME, int, 1)
        self.metadata.mt_offset_hz = self.getelem(self.metadata._hdr, TAG_MTOFF_HZ, float, 0.)
        self.metadata.slice_duration = 0.
        self.metadata.fov = 2 * [self.getelem(self.metadata._hdr, 'ReconstructionDiameter', float, 0.)]
        self.metadata.total_num_slices = self.getelem(self.metadata._hdr, 'ImagesInAcquisition', int, 0)
        self.metadata.num_timepoints = self.getelem(self.metadata._hdr, 'NumberOfTemporalPositions', int, self.metadata.total_num_slices / self.metadata.num_slices)
        if self.metadata.total_num_slices == self.metadata.num_slices:
            self.metadata.total_num_slices = self.metadata.num_slices * self.metadata.num_timepoints
        self.metadata.receive_coil_name = self.getelem(self.metadata._hdr, 'ReceiveCoilName', None, 'unknown')
        self.metadata.prescribed_duration = self.metadata.tr * self.metadata.num_timepoints * self.metadata.num_averages
        self.metadata.duration = self.metadata.prescribed_duration
        r = self.getelem(self.metadata._hdr, TAG_PHASE_ENCODE_UNDERSAMPLE, None, [1., 1.])
        self.metadata.phase_encode_undersample, self.metadata.slice_encode_undersample = [float(x) for x in (r.split('\\') if isinstance(r, basestring) else r)]
        self.metadata.num_bands = 1
        self.metadata.effective_echo_spacing = self.getelem(self.metadata._hdr, TAG_EPI_EFFECTIVE_ECHO_SPACING, float, 0.) / 1e6
        self.metadata.is_dwi = bool(self.metadata.image_type == TYPE_ORIGINAL and self.getelem(self.metadata._hdr, TAG_DIFFUSION_DIRS, int, 0) >= 6)
        self.metadata.is_multicoil = bool(self.getelem(self.metadata._hdr, 'NumberOfTemporalPositions') is None
                                          and self.metadata.total_num_slices / self.metadata.num_slices != 1
                                          and not self.metadata.is_dwi)
        self.metadata.num_receivers = 0                                 # FIXME: where to get this for GE? need to infer somehow?

        self._postinit()

    def load_data(self, archive):
        super(NIMSDicomGE, self).load_data(archive)
        log.debug('converting GEMS')
        num_dcms = len(self.metadata._dcm_list)
        self.metadata.total_num_slices = num_dcms

        if self.metadata.image_type == TYPE_SCREEN:     # FIXME; need a good catch all for screen saves
            # screen shots need to be dealt with differently, screensaves aren't the right sort of data for dcmstack
            self.data = np.dstack([np.swapaxes(dcm.pixel_array, 0, 1) for dcm in self.metadata._dcm_list])
        elif not self.metadata.is_multicoil:
            log.debug('not multicoil')
            stack = nimsdicom.NIMSStack()
            for dcm in self.metadata._dcm_list:
                stack.add_dcm(dcm)
            nii = stack.to_nifti_wrapper().nii_img

            self.data = nii.get_data()
            self.metadata.qto_xyz = nii.get_affine()
            self.metadata.sform = nii.get_sform()
            self.metadata.qform = nii.get_qform()
        else:
            log.debug('multicoil')
            # fix multicoil metadata fields
            self.metadata.num_receivers = self.metadata.num_timepoints
            log.debug('num recievers: ' + str(self.metadata.num_receivers))
            self.metadata.num_timepoints = 1    # TODO: think; is multi-coil always a single timeseries?
            # GE dicoms are sorted as slice1(coil 0..n, combined), slice2(coil 0..n, combined)...sliceX(coil 0..n, combined),
            groups = [self.metadata._dcm_list[x::self.metadata.num_receivers] for x in xrange(0, self.metadata.num_receivers)]
            log.debug('groups made: ' + str(len(groups)))
            stack_list = []
            group_id = 0
            for group in groups:
                group_id += 1
                if group_id > self.metadata.num_receivers:
                    group_id = 'combined'
                log.debug('stacking group %s, %s dicoms' % (str(group_id), str(len(group))))

                num_positions = len(set([dcm.SliceLocation for dcm in group]))
                if num_positions != self.metadata.num_slices:
                    raise NIMSDicomGEError('coil %s has %s unique positions, expected %s' % (group_id, num_positions, self.metadata.num_slices))

                stack = nimsdicom.NIMSStack()
                for dcm in group:
                    stack.add_dcm(dcm)
                nii_wrp = stack.to_nifti_wrapper()
                stack_list.append(nii_wrp)

            log.debug('combining stacks')
            merged = dcmstack.dcmmeta.NiftiWrapper.from_sequence(stack_list)
            nii = merged.nii_img
            self.data = nii.get_data()
            self.metadata.qto_xyz = nii.get_affine()
            self.metadata.sform = nii.get_sform()
            self.metadata.qform = nii.get_qform()

        # adjust dwi stuff, needs affine to proceed
        if self.metadata.is_dwi:
            rot = self.metadata.qto_xyz[0:3, 0:3]
            self.metadata.bvals = np.array([float(self.getelem(dcm, TAG_BVALUE)[0]) for dcm in self.metadata._dcm_list[0:self.metadata.num_slices]])
            self.metadata.bvecs = np.array([[self.getelem(dcm, TAG_BVEC[i], float) for i in range(3)] for dcm in self.metadata._dcm_list[0:self.metadata.num_slices]]).transpose()
            self.metadata.bvecs, self.metadata.bvals = nimsmrdata.adjust_bvecs(self.metadata.bvecs, self.metadata.bvals, self.metadata.scanner_type, rot)
