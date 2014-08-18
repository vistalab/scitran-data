Introduction
============

NIMSData is open-source software for converting scientific data from a raw format to some common output.

Provides out-of-the-box support for MRI data, capable of parsing dicoms and writing niftis, image montages, or png image files.

Provides abstract base classes for generic readers and writers that can be subclassed to expand the types of data supported.

**NIMSProcessor will provide a higher level interface, build on nimsdata. Recommend using NIMSProcessor.**


How to use the documentation
----------------------------
Documentation is available as docstrings provided within the code, and as a separate `web documentation <https://github.com>`_.


Dependencies
------------
Requires Python 2.7.

- nibabel       `pip install git+https:`
- numpy         `pip install numpy`
- dcmstack      `pip install git+https:`
- pydicom       `pip install hg+https:`
- PIL           `pip install PIL --allow-unverified PIL --allow-external PIL`


Installation
------------
- Download the latest from github. `git clone https://github.com/scitran/nimsdata.git`
- Create a vitrualenv (optional, but highly recommended). `cd /var/local; mkvirtualenv nims_env`
- install depenencies.
- run tests.


Basic Conversion
----------------
The software consists of the python package, *nimsdata*, with a single command line interface `nimsdata.py`.

NIMSData has a semi-standard input filetype, a tgz, that contains raw data and a json that indicates the raw data filetype, as well as metadata.

The `nimsdata.py` command line interface expects at least 3 options, *<input.tgz>*, *--parser <filetype>*, *--writer <filetype>*.
The following will take the *dicom* input *input.tgz* and convert it to nifti, *outprefix.nii.gz*.

.. code-block:: sh

    $ nimsdata.py -p dicom input.tgz -w nifti outprefix.nii.gz


And the equivelant command in python.

.. code-block:: python

    import nimsdata
    ds = nimsdata.parse('/path/to/input.gz', filetype='dicom')
    ds.load_data()
    nimsdata.write(ds, ds.data, 'outprefix', filetype='nifti')


For more information on using NIMSData in bash, see the :doc:`cli_tutorial`.
For more information on using NIMSData in python see the :doc:`python_tutorial`.
