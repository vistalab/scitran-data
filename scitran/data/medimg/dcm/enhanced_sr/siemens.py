# @author:  Kevin S. Hahn

"""
nimsdata.medimg.dcm.enhanced_sr.siemens
=======================================

Not implemented.

"""

import logging

from .. import dcm

log = logging.getLogger(__name__)


def parse_one(self):
    """
    Composer function.

    Parse one siemens ehanced sr dicom.

    """
    self.failure_reason = dcm.DicomError('enhanced sr/siemens has not been implemented')
    log.debug('enhanced SR not support yet, sorry!')

def parse_all(self):
    """
    Composer function.

    Parse a series of siemens enhanced sr dicoms.

    """
    self.failure_reason = dcm.DicomError('enhanced sr/siemens has not been implemented')
    log.debug('enhanced SR not support yet, sorry!')

def convert(self):
    """
    Composer function.

    Convert a series of siemens enhanced sr dicoms.

    """
    self.failure_reason = dcm.DicomError('enhanced sr/siemens has not been implemented')
    log.debug('enhanced SR not support yet, sorry!')
