import os
import numpy as np

from nose.plugins.attrib import attr
from numpy.testing.decorators import skipif
from nose.tools import ok_, eq_, raises, assert_raises

import nimsdata

DATADIR = os.path.join(os.path.dirname(__file__), 'nimsdata_testdata')
if not os.path.isdir(DATADIR):
    DATADIR = None


class Test_NIMSGEPhysio(object):

    @skipif(True, msg='not implemented')
    @skipif(not DATADIR, msg='no data')
    def test_class(self):
        pass

    @skipif(True, msg='not implemented')
    @skipif(not DATADIR, msg='no data')
    def test_parse(self):
        # MUST HAVE the following and their corresponding NIMS property
        # self.group
        # self.experiment
        # self.session
        # self.epoch
        # self.timestamp
        pass

    @skipif(True, msg='not implemented')
    @skipif(not DATADIR, msg='no data')
    def test_load_data(self):
        pass
