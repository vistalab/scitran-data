# @author:  Gunnar Schaefer
#           Kevin S Hahn

"""
scitran.data
============

The scitran data package provides two main interfaces, scitran.data.parse() for reading an input file, and
scitran.data.write() for writing to an output file.  Scitran data can read dicoms, and GE P-files.
Support for additional data domains and files types is being actively developed.

Scitran data provides read and write capabilities for scientific data.  It is implemented in Python,
using pydicom, nibabel, PIL and dcmstack.  Scitran data is open source, released under the MIT License
(see LICENSE.txt for details).

Readers and Writers are associated with their labels in external json files, readers.json and
writers.json.  New readers and writers can be added to a by defining the class and adding an item
in the appropriate json file.

"""

import data

parse = data.parse
write = data.write
get_handler = data.get_handler
dict_merge = data.dict_merge
DataError = data.DataError
acquisition_properties_by_type_list = data.acquisition_properties_by_type_list
session_properties_by_type_list = data.session_properties_by_type_list
project_properties_by_type_list = data.project_properties_by_type_list
