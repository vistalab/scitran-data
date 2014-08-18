CLI Tutorial
============

nimsdata can be used as a command line utility.


`nimsdata.py` can take several other optional arguments
    - *outbase*
    - *--parser_kwarg*
    - *--writer_kwarg*
    - *ignore_json*, *-i*
    - *--verbose*, *-v*


Basic Usage
-----------

convert tgz of dicoms, input_dicoms.tgz, to a nifti file, output_nifti.nii.gz, without voxel reordering.

.. code-block:: sh

    nimsdata.py -p dicom /path/to/input_dicoms.tgz -w nifti output_nifti.nii.gz

convert tgz of dicoms, input_dicoms.tgz, to a nifti file, output_nifti_LPS.nii.gz, with voxel reordering to LPS.

.. code-block:: sh

    nimsdata.py -p dicom /path/to/input_dicoms.tgz -w nifti output_nifti_LPS.nii.gz --write_kwarg voxel_order=LPS


