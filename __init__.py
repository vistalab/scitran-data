# @author:  Gunnar Schaefer
#           Kevin S Hahn

"""
nimsdata
========

The nimsdata package provides two main interfaces, nimsdata.parse() for reading an input file, and
nimsdata.write() for writing to an output file.  Nimsdata can read dicoms, and GE P-files.
Support for additional data domains and files types is being actively developed.

NIMSdata provides read and write capabilities for scientific data.  It is implemented in Python,
using pydicom, nibabel, PIL and dcmstack.  nimsdata is open source, released under the MIT License
(see LICENSE.txt for details).

Readers and Writers are associated with their labels in external json files, readers.json and
writers.json.  New readers and writers can be added to a by defining the class and adding an item
in the appropriate json file.

"""

import nimsdata

parse = nimsdata.parse
write = nimsdata.write
dict_merge = nimsdata.dict_merge
NIMSDataError = nimsdata.NIMSDataError
epoch_properties_by_type_list = nimsdata.epoch_properties_by_type_list
session_properties_by_type_list = nimsdata.session_properties_by_type_list
experiment_properties_by_type_list = nimsdata.experiment_properties_by_type_list
