# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.sc.siemens
==============================

Not implemented.

"""

import logging

import sc


log = logging.getLogger(__name__)


SIEMENS_TYPE_DIFF_FA = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'FA', 'ND']
SIEMENS_TYPE_DIFF_TENSOR = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'TENSOR', 'ND']
SIEMENS_TYPE_DIFF_FA_NORM = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'FA', 'ND', 'NORM']
SIEMENS_TYPE_DIFF_TENSOR_NORM = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'TENSOR', 'ND', 'NORM']
SIEMENS_TYPE_CSA = ['ORIGINAL', 'PRIMARY', 'OTHER', 'CSA REPORT']
SIEMENS_TYPE_TENSOR = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'TENSOR', 'ND']
# there are several varieties of CSAPARALLEL images, MPR, PROJECTION IMAGE.
# the only common theme so far is that they all have 'CSAPARALLEL' in imagetype, and are invalid stacks
# because the dcm metadata are inconsistent
SIEMENS_TYPE_DIFF_FA = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'FA', 'ND']
SIEMENS_TYPE_DIFF_FA_NORM = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'FA', 'ND', 'NORM']
SIEMENS_TYPE_ASL_TTEST = ['DERIVED', 'PRIMARY', 'PERFUSION', 'ASL', 'ND', 'NORM', 'FILTERED', 'MOCO', 'SUB', 'TTEST', 'MOSAIC']
SIEMENS_TYPE_RELCBF_TTEST = ['DERIVED', 'PRIMARY', 'PERFUSION', 'ASL', 'RELCBF', 'ND', 'NORM', 'FILTERED', 'MOCO', 'SUB', 'TTEST', 'MOSAIC']
SIEMENS_TYPE_CSA_3D = ['DERIVED', 'SECONDARY', 'OTHER', 'CSA 3D EDITOR']
# SIEMENS_TYPE_POSDISP = ['DERIVED', 'SECONDARY', 'POSDISP', 'M', 'ND', 'NORM', 'CSA RESAMPLED']i  MR
# posdisp?  resample images of the previous anatomical scan, to view from 3 angles to estimate positioning?  not really sure
# but the data was not acquired as such.

SIEMENS_NON_IMAGE_TYPES = [SIEMENS_TYPE_DIFF_FA, SIEMENS_TYPE_DIFF_TENSOR, SIEMENS_TYPE_DIFF_FA_NORM, SIEMENS_TYPE_DIFF_TENSOR_NORM]


def parse_one(self):

    log.debug(self.image_type)
    self.is_non_image = True

    # preliminary identification
    # if self.image_type == SIEMENS_TYPE_CSA:
    #     self.is_non_image = True
    # # is non-image
    # if self.getelem(self._hdr, 'SOPClassUID') == SIEMENS_SOPCLASS_SC:
    #     # Seconary Storage Class Images always are missing ImagePatientPositions
    #     # .: they never have enough information to be a valid MR dicom
    #     self.is_non_image = True
    # if self.image_type == SIEMENS_TYPE_POSDISP:
    #     self.is_non_image = True
    # if self.image_type == SIEMENS_TYPE_CSA_3D:
    #     # Siemens 3D Storage. what is this?
    #     self.is_non_image = True
    # if self.image_type in [SIEMENS_TYPE_ASL_TTEST or self.image_type, SIEMENS_TYPE_RELCBF_TTEST]:
    #     # CSA header contains WAY TOO MUCH info...ttest results?
    #     self.is_non_image = True
    # if self.image_type in [SIEMENS_TYPE_DIFF_FA_NORM, SIEMENS_TYPE_DIFF_FA]:
    #     # Diffusion FA NORM dicoms count as is_non_image because they lack the information
    #     # necessary to create the affine xform for the nifti. Also, CSA data instead of pixel data.
    #     self.is_non_image = True
    # if u'CSAPARALLEL' in self.image_type:
    #     # csaparallel images are secondary derived dataset
    #     # these will rarely be valid, because variation in the dicoms PixelSpacing
    #     # the first dicom will accurately hold PixelSpacing from the parent scan
    #     # all remaining dicoms will have PixelSpacing of [1, 1]
    #     self.is_non_image = True
    # if self.getelem(self._hdr, 'PrivateCreator_0x29_0x10') == 'SIEMENS CSA NON-IMAGE':
    #     self.is_non_image = True
    # if self.getelem(self._hdr, 'VariablePixelData') == 'SIEMENS CSA NON-IMAGE':
    #     self.is_non_image = True
    # if self.image_type == SIEMENS_TYPE_TENSOR:


def parse_all(self):
    # set special metadata only available via having ALL the dicoms available
    pass


def convert(self):
    if self.is_screenshot:
        sc.convert_screenshot(self)
