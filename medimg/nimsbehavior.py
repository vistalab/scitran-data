# @author:  Gunnar Schaefer
#           Bob Dougherty

"""
nimsdata.medimg.nimsbehavior
============================

Not implemented.

This will parse behavioral data that is related to medical imaging data.  This will potentially
parse outputs from stimulus presentation software that are commonly used in medical imaging. Such
as outputs from Matlab, Eprime2 and PsychoPy.

nimsdata.medimg.nimsbehavior will need to be paired with a writer that is capable of outputting to
text or csv.  Currently, there is no such writer.

"""

from .. import nimsdata


class NIMSBehaviorError(nimsdata.NIMSDataError):
    pass


class NIMSBehavior(nimsdata.NIMSReader):

    """Munges response data into an intermediate format."""

    domain = u'mr'
    filetype = u'behavior'
    state = ['orig']

    def __init__(self, filepath):
        raise NIMSBehaviorError('NIMSBehavior Reader class not yet implemented')
        super(NIMSBehavior, self).__init__()

    def load_data(self):
        raise NIMSBehaviorError('NIMSBehavior Reader class not yet implemented')
        super(NIMSBehavior, self).__init__()

    @property
    def nims_metadata_status(self):
        return self.metadata_status
