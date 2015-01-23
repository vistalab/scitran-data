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

from .. import data


class BehaviorError(data.DataError):
    pass


class Behavior(data.Reader):

    """Munges response data into an intermediate format."""

    domain = u'mr'
    filetype = u'behavior'
    state = ['orig']

    def __init__(self, filepath):
        raise BehaviorError('Behavior Reader class not yet implemented')
        super(Behavior, self).__init__()

    def load_data(self):
        raise BehaviorError('Behavior Reader class not yet implemented')
        super(Behavior, self).__init__()

    @property
    def nims_metadata_status(self):
        return self.metadata_status
