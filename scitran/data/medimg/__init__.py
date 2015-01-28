# @author:  Kevin S. Hahn

"""
scitran.data.medimg
===================

The sctran.data.medimg module provides reading and writing capabilities for medical image file types.
Reading from dicoms and niftis(in-progress), and writing to niftis.

Currently medimg module includes response data that is associated with medical images such as
physiological recordings or behavioral responses collected during an acquisition.

MedImgReaders should assign their data to self.data as a dictionary of numpy.darrays.  The primary
dataset should be assigned to the key "" (empty string), while any additional secondary data
should be assigned to some user-defined key.  For each dictionary item in self.data, the value
will be written to a file suffixed with the key.

.. code:: json

    self.data = {
        "": numpy.ndarray(),
        "fieldmap": numpy.ndarray()
    }

The same metadata will be written for each array of data in the self.data dictionary.

"""

from . import medimg

MedImageReader = medimg.MedImgReader
MedImageWriter = medimg.MedImgWriter
acquisition_properties = medimg.acquisition_properties
session_properties = medimg.session_properties
project_properties = medimg.project_properties
get_slice_order = medimg.get_slice_order
parse_patient_id = medimg.parse_patient_id
parse_patient_name = medimg.parse_patient_name
parse_patient_dob = medimg.parse_patient_dob
