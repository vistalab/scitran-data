"""
nimsdata.medimg.dcm
===================

TODO: docs? here? yes.

Currently supported dicoms include the following Manufacturers and SOPs

SUPPORTED_MFR
    - 'GE MEDICAL SYSTEMS'
    - 'SIEMENS'

SUPPORTED_SOP
    - '1.2.840.10008.5.1.4.1.1.4'       MR Image
    - '1.2.840.10008.5.1.4.1.1.7'       Secondary Capture
    - '1.3.12.2.1107.5.9.1'             Private Syngo CSA Non-Image; SIEMENS ONLY
    - '1.2.840.10008.5.1.4.1.1.88.22'   Enhanced SR
    - '1.2.840.10008.5.1.4.1.1.128'     PET
    - '1.2.840.10008.5.1.4.1.1.130'     Enhanced PET
    - '1.2.840.10008.5.1.4.1.1.128.1'   Legacy Enhanced PET

"""

import dcm

Dicom = dcm.Dicom
DicomError = dcm.DicomError
MAX_LOC_DCMS = dcm.MAX_LOC_DCMS
MetaExtractor = dcm.MetaExtractor
timestamp = dcm.timestamp
parse_patient_name = dcm.parse_patient_name
parse_patient_id = dcm.parse_patient_id
parse_patient_dob = dcm.parse_patient_dob
