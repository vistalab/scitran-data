# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.sc.ge
=========================

Composer functions for use with GE Secondary Capture Dicoms.

"""

import logging

import generic_sc


log = logging.getLogger(__name__)


GEMS_TYPE_SCREENSHOT = ['DERIVED', 'SECONDARY', 'SCREEN SAVE']
GEMS_TYPE_VXTL = ['DERIVED', 'SECONDARY', 'VXTL STATE']
GEMS_TYPE_RFMT = ['DERIVED', 'SECONDARY', 'REFORMATTED', 'AVERAGE']
GEMS_NON_IMAGE_TYPES = [GEMS_TYPE_VXTL, GEMS_TYPE_RFMT]


def parse_one(self):
    if self.image_type in [GEMS_TYPE_SCREENSHOT, GEMS_TYPE_VXTL]:
        self.is_screenshot = True
        self.scan_type = 'screenshot'
        self.nims_metadata_status = None

def parse_all(self):
    pass

def convert(self):
    if self.is_screenshot:
        generic_sc.convert_screenshot(self)
        self.nims_metadata_status = None
