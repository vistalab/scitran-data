Python Tutorial
===============

Basic Usage
-----------
Performing conversions from Python requires a few extra steps compared to the command line usage, but offers
great flexibility.

convert dcm to nifti, do not reorient voxels

.. code-block:: python

    ds = data.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    # ds, the dataset contains metadata, such as affine info
    # ds.data contains a dictionary of np arrays of pixel data
    data.write(ds, ds.data, outpath, filetype='nifti')


convert dcm to nifti, reorient voxels to 'LAS'

.. code-block:: python

    ds = data.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    data.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LAS')


convert dcm to nifti, backward compatible with data v1 (reorient to 'LPS')

.. code-block:: python

    ds = data.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    data.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LPS')


convert pfile to nifti, where there is primary data, and secondary field map data

.. code-block:: python

    ds = data.parse(pfile_tgz, load_data=True, ignore_json=False, filetype='pfile')
    # ds.data is {'': primary_data_array, 'fieldmap': secondary_data_array}
    data.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LPS')
    # creates outpath.nii.gz and outpath_fieldmap.nii.gz


Advanced Usage
--------------
Reading and writing are decoupled to create a user-accesible space between reading and writing.
This separation allows for conditional writing based upon the input file type
conditional writing

.. code-block:: python

    ds = data.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    if ds.is_screenshot:
        # PNG default voxel_order=LPS
        data.write(ds, ds.data, outpath, filetype='png')
    else:
        # Montage default voxel_order=LPS
        data.write(ds, ds.data, outpath, filetype='montage')
        # Nifti default voxel_order=None
        data.write(ds, ds.data, outpath, filetype='nifti')


.. code-block:: python

    import data

    ds = data.dcm.Dicom(input_tgz, load_data=True)
    if ds.is_screenshot:
        # PNG default voxel_order=LPS
        data.write(ds, ds.data, outpath, filetype='png')
    else:
        # Montage default voxel_order=LPS
        data.write(ds, ds.data, outpath, filetype='montage')
        # Nifti default voxel_order=None
        data.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LPS')
