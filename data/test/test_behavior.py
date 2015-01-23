# @author:  Kevin S Hahn

import os
import numpy as np

from nose.tools import ok_, eq_, raises, assert_raises
from numpy.testing.decorators import skipif

import nibabel
import dicom

import data

DATADIR = os.path.join(os.path.dirname(__file__), 'nimsdata_testdata')
if not os.path.isdir(DATADIR):
    DATADIR = None


class Test_Behavior(object):

    @skipif(True, msg='not implemented')
    @skipif(not DATADIR, msg='no data')
    def test_parse(self):
        """behavior parse"""
        pass

    @skipif(True, msg='not implemented')
    @skipif(not DATADIR, msg='no data')
    def test_load(self):
        """behavior load."""
        pass
