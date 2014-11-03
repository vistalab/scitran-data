NIMSData |travis_badge|
=======================

.. |travis_badge| image:: https://travis-ci.org/scitran/nimsdata.svg?branch=master
    :target: https://travis-ci.org/scitran/nimsdata

Data handling utilities.

Provides class and functions for parsing, identifying and writing scientific datasets. Currently,
medical images (such as various flavors of MRI, dicoms, niftis, and some raw-files) are supported.

Also provides base classes that can be used to extend this package to cover additional data-domains.

Note that these utilites also have command-line interfaces, so you don't need a full NIMS
instance to use them. E.g., you can run our image reconstruction pipeline.


How to use the documentation
----------------------------
Documentation is available as docstrings provided within the code, and
`additional narrative documentation <https://scitran.github.io/nimsdata>`_.


Dependencies
------------

development version is currently using the following:

================ ====================
package          version
================ ====================
Pillow           2.6.1
numpy            1.9.0
nibabel          1.4.0dev
pydicom          0.9.9
dcmstack         0.7.0dev
================ ====================


installation on fresh ubuntu 14.04
----------------------------------
- install dependencies, **had to chmod o+w /var/local to allow regular user to write; not ideal**
    - pillow must be compiled with JPEG support, see `Python Image Library fails with message “decoder JPEG not available” - PIL
      <http://stackoverflow.com/questions/8915296/python-image-library-fails-with-message-decoder-jpeg-not-available-pil>`_.

.. code:: bash

    sudo apt-get update
    sudo apt-get upgrade
    sudo apt-get install python-dev python-virtualenv libjpeg-dev git
    sudo ln -s /usr/lib/x86_64-linux-gnu/libjpeg.so /usr/lib

- create and activate virtualenv

.. code:: bash

    virtualenv /var/local/nimsdata_env
    source /var/local/nimsdata_env/bin/bin/activate

- install dependencies, these are VERY specific version. once these versions are available
  through pypi, these installation commands will change.  (subject to change)

.. code:: bash

    pip install numpy==1.9.0
    pip install git+https://github.com/scitran/pydicom.gitmirror.git@value_mismatch
    pip install git+https://github.com/nipy/nibabel.git
    pip install git+https://github.com/moloney/dcmstack
    pip install pillow
    pip install pymongo

- clone nimsdata from github

.. code:: bash

    git clone https://github.com/scitran/nimsdata.git /var/local/nimsdata
    cd /var/local/nimsdata
    git submodule init  # optional
    git submodule update  # optional


Basic Conversion
----------------
The software consists of the python package, *nimsdata*, with a single command line interface
`nimsdata.py`.

NIMSData has a semi-standard input filetype, a tgz, that contains raw data and a json that
indicates the raw data filetype, header data, and metadata corrections.

`nimsdata.py` expects at least 3 options, *<input.tgz>*, *--parser <filetype>*, *--writer <filetype>*.

The following shell command will take the *dicom* input *input.tgz* and convert it to nifti, *outprefix.nii.gz*.

.. code-block:: sh

    nimsdata.py -p dicom input.tgz -w nifti outprefix.nii.gz


And the equivelant command in python.

.. code-block:: python

    import nimsdata
    ds = nimsdata.parse('/path/to/input.gz', filetype='dicom')
    ds.load_data()
    nimsdata.write(ds, ds.data, 'outprefix', filetype='nifti')


For more information on using NIMSData in bash, see `CLI tutorial <https://scitran.github.io/cli_tutorial.html>`_.

For more information on using NIMSData in python see `Python tutorial <https://scitran.github.io/nimsdata/python_tutorial.html>`_.


Developer Notes
---------------

To generate the docs locally, you will need sphinx, and numpydoc.

.. code:: bash

    pip install sphinx numpydoc


numpy 1.9 changes how numpy.unique() behaves when given an array of arrays.  Pre 1.9, np.unique
would return each unique array. Post 1.9, np.unique returns unique items from the arrays. dcmstack
is compatible with numpy 1.9, but numpy throws some FutureWarnings.  The current version of
dcmstack (0.7.0dev) may not be compatible with future version of numpy.

run the following git config commands to enable a git filter for the branch name.

.. code:: bash

    git config filter.brancher.smudge "./git_branch_filter.py smudge"
    git config filter.brancher.clean "./git_branch_filter.py clean"

Combined with .gitattributes, the smudge and clean filters will
replace 'branch=\_\_BRANCH\_\_' to indicate the current branch.
