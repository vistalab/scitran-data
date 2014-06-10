# @author:  Gunnar Schaefer
#           Bob Dougherty

"""
nimsdata.nimsbehavior
=====================

not implemented!

"""

import nimsdata


class NIMSBehaviorError(nimsdata.NIMSDataError):
    pass


class NIMSBehavior(nimsdata.NIMSReader):

    def __init__(self, filepath):
        raise NIMSBehaviorError('NIMSBehavior Reader class not yet implemented')
        super(NIMSBehavior, self).__init__()

    def load_data(self):
        raise NIMSBehaviorError('NIMSBehavior Reader class not yet implemented')
        super(NIMSBehavior, self).__init__()
