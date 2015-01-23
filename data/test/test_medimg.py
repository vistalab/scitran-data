# @Author: Kevin S Hahn

"""Test nimsdata package."""

import os
import glob
import datetime
import numpy as np

from nose.plugins.attrib import attr
from numpy.testing.decorators import skipif
from nose.tools import ok_, eq_, raises, assert_raises

import data.medimg
import data.tempdir as tempfile

# data is stored separately in nimsdata_testdata
# located at the top level of the testing directory
DATADIR = os.path.join(os.path.dirname(__file__), 'nimsdata_testdata')
if not os.path.isdir(DATADIR):
    DATADIR = None

class test_get_slice_order(object):
    def test_seq_inc(self):
        pass
    def test_seq_dec(self):
        pass
    def test_alt_inc(self):
        pass
    def test_ald_dec(self):
        pass
    def test_alt_inc2(self):
        pass
    def test_alt_dec2(self):
        pass


class test_parse_patient_id(object):
    default_subj_code = 'default'

    def test_parse_valid_patient_id(self):
        patient_id = '1234@scitran/data'
        s, g, p = data.medimg.parse_patient_id(patient_id, self.default_subj_code)
        eq_(s, '1234')
        eq_(g, 'scitran')
        eq_(p, 'data')

    def test_parse_invalid_patient_id(self):
        patient_id = 'scitran'
        s, g, p = data.medimg.parse_patient_id(patient_id, self.default_subj_code)


class test_parse_patient_name(object):
    def test_patient_name_delimiters(self):
        f, l = data.medimg.parse_patient_name('Last^First')
        eq_((f, l), ('First', 'Last'))
        f, l = data.medimg.parse_patient_name('First Last')
        eq_((f, l), ('First', 'Last'))
        f, l = data.medimg.parse_patient_name('FirstLast')
        eq_((f, l), ('', 'Firstlast'))

class test_parse_patient_dob(object):
    def test_parse_valid_patient_dob(self):
        eq_(datetime.datetime(2011, 1, 1, 0, 0),  data.medimg.parse_patient_dob('20110101'))

    def test_parse_invalid_patient_dob(self):
        eq_(None,  data.medimg.parse_patient_dob('invalid'))
        eq_(None,  data.medimg.parse_patient_dob('18990101'))
