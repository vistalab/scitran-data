import os
import glob
import numpy as np

from nose.plugins.attrib import attr
from numpy.testing.decorators import skipif
from nose.tools import ok_, eq_, raises, assert_raises

import data
import data.tempdir as tempfile

# data is stored separately in nimsdata_testdata
# located at the top level of the testing directory
DATADIR = os.path.join(os.path.dirname(__file__), 'nimsdata_testdata')
if not os.path.isdir(DATADIR):
    DATADIR = None

class Test_PNG(object):

    @skipif(not DATADIR)
    def setUp(self):
        # TODO: remake ge_dcm_screenshot with json
        self.ds = data.parse(os.path.join(DATADIR, 'ge_dcm_sc_screenshot.tgz'), load_data=True, filetype='dicom', ignore_json=True)

    @skipif(not DATADIR)
    def test_write(self):
        with tempfile.TemporaryDirectory() as tempdir:
            outbase = os.path.join(tempdir, 'trashme')
            data.write(self.ds, self.ds.data, outbase=outbase, filetype='png')
            assert len(glob.glob(outbase + '*')) >= 1

