# @author:  Gunnar Schaefer
#           Kevin S Hahn

"""
nimsdata
========

The nimsdata package provides two interfaces; one for reading from file, and the other for writing to file.

Readers and Writers are defined in external json files, readers.json and writers.json.  New readers and writers can be added to
a running system by defining it in the appropriate json file.

code is organized into a three tier structure - domain, filetype and file-type-variant.
- medical images            (domain)
- - dicoms                  (filetype)
- - - mr dicoms             (variant)
- - - - ge mr dicoms        (variant)
- - - - siemens mr dicoms   (variant)
- - - - genenric mr dicoms  (variant)
- - - sc dicoms             (variant)
- - - - ge mr dicoms        (variant)
- - - - siemens sc dicoms   (variant)

"""

import nimsdata

parse = nimsdata.parse
write = nimsdata.write
dict_merge = nimsdata.dict_merge
NIMSDataError = nimsdata.NIMSDataError
module_by_type = nimsdata.module_by_type
epoch_properties_by_type_list = nimsdata.epoch_properties_by_type_list
session_properties_by_type_list = nimsdata.session_properties_by_type_list
experiment_properties_by_type_list = nimsdata.experiment_properties_by_type_list
