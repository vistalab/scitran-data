# @author:  Gunnar Schaefer
#           Bob Dougherty

"""
nimsdata.nimsbehavior
=====================

not implemented!

"""

from .. import nimsdata


class NIMSBehaviorError(nimsdata.NIMSDataError):
    pass


class NIMSBehavior(nimsdata.NIMSReader):

    domain = u'mr'
    filetype = u'behavior'
    state = ['orig']

    def __init__(self, filepath):
        raise NIMSBehaviorError('NIMSBehavior Reader class not yet implemented')
        super(NIMSBehavior, self).__init__()

    def load_data(self):
        raise NIMSBehaviorError('NIMSBehavior Reader class not yet implemented')
        super(NIMSBehavior, self).__init__()
