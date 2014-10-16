# @author:  Kevin S. Hahn

"""
nimsdata.medimg
===============

The nimsdata.medimg module provides reading and writing capabilities for medical image file types.
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
epoch_properties = medimg.epoch_properties
session_properties = medimg.session_properties
experiment_properties = medimg.experiment_properties
