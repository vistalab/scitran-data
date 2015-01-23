import os
import numpy as np

from nose.plugins.attrib import attr
from numpy.testing.decorators import skipif
from nose.tools import ok_, eq_, raises, assert_raises

import data
import data.tempdir as tempfile

# data is stored separately in data_testdata
# located at the top level of the testing directory
DATADIR = os.path.join(os.path.dirname(__file__), 'data_testdata')
if not os.path.isdir(DATADIR):
    DATADIR = None


class Test_Montage(object):

    @skipif(not DATADIR)
    def setUp(self):
        self.ds = data.parse(os.path.join(DATADIR, 'ge_dcm_mr_localizer.tgz'), load_data=True)

    @skipif(not DATADIR)
    def test001_write_png(self):
        """flat png montage"""
        with tempfile.TemporaryDirectory() as tempdir:
            outbase = os.path.join(tempdir, 'trashme')
            ok_(data.write(self.ds, self.ds.data, outbase=outbase, filetype='montage', mtype='png'))

    @skipif(not DATADIR)
    def test002_write_dir(self):
        """directory jpeg montage"""
        with tempfile.TemporaryDirectory() as tempdir:
            outbase = os.path.join(tempdir, 'trashme')
            ok_(data.write(self.ds, self.ds.data, outbase=outbase, filetype='montage', mtype='dir'))

    @skipif(not DATADIR)
    def test003_write_pyrdb(self):
        with tempfile.TemporaryDirectory() as tempdir:
            outbase = os.path.join(tempdir, 'trashme')
            ok_(data.write(self.ds, self.ds.data, outbase=outbase, filetype='montage', mtype='sqlite'))

    @skipif(not DATADIR)
    def test004_get_size(self):
        with tempfile.TemporaryDirectory() as tempdir:
            outbase = os.path.join(tempdir, 'trashme')
            data.write(self.ds, self.ds.data, outbase=outbase, filetype='montage')
            outfile = os.path.join(tempdir, os.listdir(tempdir)[0])
            ok_(data.medimg.montage.get_info(outfile))              # gets made?

    @skipif(not DATADIR)
    def test005_get_tile(self):
        with tempfile.TemporaryDirectory() as tempdir:
            outbase = os.path.join(tempdir, 'trashme')
            data.write(self.ds, self.ds.data, outbase=outbase, filetype='montage')
            outfile = os.path.join(tempdir, os.listdir(tempdir)[0])
            ok_(data.medimg.montage.get_tile(outfile, 0, 0, 0))     # all montage have 0, 0, 0
