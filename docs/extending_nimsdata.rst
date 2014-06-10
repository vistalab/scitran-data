Extending NIMSData
==================

Creating a custom NIMSReader
----------------------------
Readers generally have three specific roles

1. be capable of reading a specific filetype, e.g. dicom, nifti.
2. be capable of parsing metadata necessary for identification and sorting, defined in __init__
3. be capable of converting the data to a specified intermediate format, defined in load_data

Readers may manipulate the file they are working on, e.g. NIMSdicom anonymizing patient information of it's own dicoms.

Readers, being readers, should generally avoid:

1. writing a filetype other than their own, e.g. NIMSDicom attempting to write a nifti
2. other restriction...

NIMSReader should be subclasses to add support for additional scientific data domains/types.

For Example, NIMSMRReader and NIMSMRWriters are subclasses of NIMSReader and NIMSWriter, respectively.  The subclasses
add MR data specific properties and functions.  NIMSDicom and NIMSPFile are subclasses of NIMSMRReader, that add
more specific functions with respect to their datatypes.


Creating a custom NIMSWriter
----------------------------
Subclasses of NIMSWriter should stick to the following conventions:

1. NIMSWriter appends its own file extension.  The writer can optionally append a name suffix
to differentiate different types of outputs, or differentiate multiple outputs.

NIMS intermediate formats are always a dictionary of data.  each intermediate format is relatively
unique to each data domain.  see NIMSMRReader and NIMSMRWriter...

add instructions for subclassing NIMSWriter.


Creating a data domain
----------------------

Is 'data domain' okay name for categories of data that can share an intermediate data format, such as "MR", or  "Genetics"?

These data domains would define compatibility between readers and writers.

