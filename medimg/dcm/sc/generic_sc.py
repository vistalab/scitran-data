# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.sc.generic_sc
=================================

Generic functions that are shared by all Secondary Capture Dicoms.  These are manufacturer
agnostic.

"""

import logging
import numpy as np

log = logging.getLogger('dcm.sc.generic')


def convert_screenshot(self):
    """conver the voxel data from dcm list into np array."""
    log.debug('screenshot recon')
    self.slice_order = None  # slice order is not unknown, just Not Applicable.
    self.psd_type = None  # not a real PSD.
    self.qto_xyz = None  # screen capture has no affine
    self.data = {'': np.dstack([d.pixel_array for d in self._dcm_list])}
