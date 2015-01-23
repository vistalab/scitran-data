# @author:  Gunnar Schaefer
#           Bob Dougherty

"""
scitran.data.medimg.behavior
============================

Not implemented.

This will parse behavioral data that is related to medical imaging data.  This will potentially
parse outputs from stimulus presentation software that are commonly used in medical imaging. Such
as outputs from Matlab, Eprime2 and PsychoPy.

scitran.data.medimg.behavior will need to be paired with a writer that is capable of outputting to
text or csv.  Currently, there is no such writer.

"""

import logging
from .. import data
log = logging.getLogger(__name__)

class BehaviorError(data.DataError):
    pass


class Behavior(data.Reader):

    """Munges response data into an intermediate format."""

    domain = u'mr'
    filetype = u'behavior'
    state = ['orig']

    def __init__(self, filepath, load_data=False):
        super(Behavior, self).__init__(filepath, load_data)
        log.debug('behavior reader has not been implemented')
        self.failure_reason = BehaviorError('behavior reader has not been implemented')

    def load_data(self):
        super(Behavior, self).__init__()
        log.debug('behavior reader has not been implemented')
        self.failure_reason = BehaviorError('behavior reader has not been implemented')

    @property
    def nims_metadata_status(self):
        return self.metadata_status
