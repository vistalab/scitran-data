CLI Tutorial
============

Nimsdata can be used as a command line utility to perform basic conversions.

`nimsdata.py` requires three arguments; what file to parse, which parser to use, what writer to use.
`nimsdata.py -p <parser> -w <writer> input_file.tgz`

`nimsdata.py` will also accept an outbase prefix, parser_kwarg and writer_kwarg.
`nimsdata.py -p <parser> -w <writer> input_file.tgz <outprefix> --writer_kwarg voxel_order=LPS`

TODO: it would be nice to be able to see all the parser and writers that are available, along with all of the
arguments they accept.  This would make the CLI much more useful. However, once nimsproc is ready, people should be
using nimsproc, which relies on nimsdata, rather than using nimsdata directly.

parser_kwarg and writer_kwarg can be called multiple times to set multi keyword arguments for the parser and/or writer.
Currently, parser_kwarg and writer_kwarg are specific to each parser and writer. For example, niftis have an 'voxel_order' kwarg, that determines the voxel order of the output nifti.  pfiles have a 'load_all' kwarg that determines if nimsdata should attempt to parse the entire file,


Basic Usage
-----------

convert tgz of dicoms, input_dicoms.tgz, to a nifti file, output_nifti.nii.gz, without voxel reordering.

.. code-block:: sh

    nimsdata.py -p dicom /path/to/input_dicoms.tgz -w nifti output_nifti.nii.gz

convert tgz of dicoms, input_dicoms.tgz, to a nifti file, output_nifti_LPS.nii.gz, with voxel reordering to LPS.

.. code-block:: sh

    nimsdata.py -p dicom /path/to/input_dicoms.tgz -w nifti output_nifti_LPS.nii.gz --write_kwarg voxel_order=LPS


Advanced Usage
--------------
Currently the CLI is limited to basic conversion.  Not much "advanced usage" to go over.
