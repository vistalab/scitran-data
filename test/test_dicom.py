# @author:  Kevin S Hahn
"""Tests for nimsdata.nimsdata."""

import os
import numpy as np

from nose.plugins.attrib import attr
from numpy.testing.decorators import skipif
from nose.tools import ok_, eq_, raises, assert_raises

import nimsdata

# data is stored separately in nimsdata_testdata
# located at the top level of the testing directory
DATADIR = os.path.join(os.path.dirname(__file__), 'nimsdata_testdata')
if not os.path.isdir(DATADIR):
    DATADIR = None


class Test_NIMSDicom(object):

    """
    Test NIMSDicom primary interfaces.

    This suite of tests will large, this should try several combinations of files
    to confirm that certain data loads work correctly.

    """

    @skipif(not DATADIR)
    def test_parse(self):
        """
        Parse one dicom as header.

        NIMSDicom parses a single dicom, using nimsdicom.MetaExtractor, upon initialization.
        """
        # test 1 - keeps private tags
        # this test is rather indirect, nimsdicom uses the MetaExtractor
        # to set the dataset._hdr, testing it this way, allows re-using one of the test files
        ds = nimsdata.parse(os.path.join(DATADIR, 'ge_dcm_mr_localizer.tgz'))
        ok_(ds._hdr.get('PrivateCreator_0X9_0X0'))  # try getting the first PrivateCreator, tag (0x0009,0x0000)

    @skipif(not DATADIR)
    def test_composer(self):
        """
        nimsdicom composer.

        Check that nimsdicom composes based on different input dicoms.  Needs access
        to several different individual dicoms.  Should at least test dcm.mr.ge and dcm.mr.siemens.

        """
        pass

    @skipif(not DATADIR)
    def test_load_data(self):
        """Load pixels into self.data"""
        pass

    @skipif(not DATADIR, msg='no test data')
    def test_dcm_mr_ge(self):
        """dcm.mr.ge"""
        pass

    @skipif(not DATADIR, msg='no test data')
    def test_dcm_mr_siemens(self):
        """dcm.mr.siemens"""
        pass

    @skipif(not DATADIR, msg='no test data')
    def test_dcm_sc_ge(self):
        """dcm.sc.ge"""
        pass

    @skipif(not DATADIR, msg='no test data')
    def test_dcm_sc_siemens(self):
        """dcm.sc.siemens"""
        pass

    @skipif(not DATADIR, msg='no test data')
    def test_dcm_syngo_siemens(self):
        """dcm.syngo.siemens"""
        pass

    @skipif(not DATADIR, msg='no test data')
    def test_dcm_sr_ge(self):
        """dcm.sr.ge"""
        pass
