Python Tutorial
===============

Basic Usage
-----------
Performing conversions from Python requires a few extra steps compared to the command line usage.

**dcm to nifti, no reorientation**

.. code-block:: python

    import nimsdata

    ds = nimsdata.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    # ds, the dataset contains metadata, such as affine info
    # ds.data contains np array of pixel data
    nimsdata.write(ds, ds.data, outpath, filetype='nifti')


**dcm to nifti, reorient to 'LAS'**

.. code-block:: python

    import nimsdata

    ds = nimsdata.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    nimsdata.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LAS')


**dcm to nifti, backward compatible with nimsdata v1 (reorient to 'LPS')**

.. code-block:: python

    import nimsdata

    ds = nimsdata.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    nimsdata.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LPS')


Advanced Usage
--------------
**conditional writing**

.. code-block:: python

    ds = nimsdata.parse(input_tgz, load_data=True, ignore_json=False, filetype='dicom')
    if ds.is_screenshot:
        # PNG default voxel_order=LPS
        nimsdata.write(ds, ds.data, outpath, filetype='png')
    else:
        # Montage default voxel_order=LPS
        nimsdata.write(ds, ds.data, outpath, filetype='montage')
        # Nifti default voxel_order=None
        nimsdata.write(ds, ds.data, outpath, filetype='nifti')


**directly access parser and writer**

.. code-block:: python

    import nimsdata.nimspng
    import nimsdata.nimsdicom
    import nimsdata.nimsnifti
    import nimsdata.nimsmontage

    ds = nimsdata.nimsdicom.NIMSDicom(input_tgz, load_data=True)
    if ds.is_screenshot:
        # PNG default voxel_order=LPS
        nimsdata.nimspng.NIMSPNG.write(ds, ds.data, outpath)
    else:
        # Montage default voxel_order=LPS
        nimsdata.nimsmontage.NIMSMontage.write(ds, ds.data, outpath)
        # Nifti default voxel_order=None
        nimsdata.nimsnifti.NIMSNifti.write(ds, ds.data, outpath, filetype='nifti', voxel_order='LPS')



