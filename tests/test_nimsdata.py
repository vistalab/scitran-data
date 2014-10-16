# @author:  Kevin S. Hahn
"""
tests for nimsdata package main API.

Relies on having valid MR dicoms to use for testing.

Download the test data from here! <linky plz>
"""

import os
import sys
import numpy as np

from nose.tools import ok_, eq_, raises, assert_raises
from nose.plugins.attrib import attr

from numpy.testing.decorators import skipif

import dicom
import nibabel
import inspect

testdat_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

import nimsdata


class API(object):

    """Tests the main nimsdata interfaces."""

    def parse(self):
        pass

    def write(self):
        pass

    def properties_by_type_list(self):
        pass

    def dict_merge(self):
        pass

    def NIMSDataError(self):
        pass
