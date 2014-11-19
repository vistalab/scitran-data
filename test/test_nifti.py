import os
import sys
import numpy as np

from nose.tools import ok_, eq_, raises, assert_raises
from nose.plugins.attrib import attr

from numpy.testing.decorators import skipif

import dicom
import nibabel

testdat_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

import nimsdata


class Test_NIMSNifti(object):

    @skipif(True, msg='nifti reader not implemented')
    def test_reading(self):
        """nimsnifti reader not implemented"""
        # raise error at init if input file is not a nifti

        # raise error at load_data if input file is not a nifti

        pass

    def test_writing(self):
        """
        Write pixeldata and metadata to nifti

        This should also test voxel_order param.
        """
        pass
