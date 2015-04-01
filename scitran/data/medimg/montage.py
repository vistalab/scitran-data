# @author:  Bob Dougherty
#           Gunnar Schaefer
#           Kevin S Hahn

"""
scitran.data.medimg.montage
===========================

Montage provides montage writing capabilities for MR datasets read by any subclass of MedImgReader.

Provides a MedImgWriter subclass for creating image pyramids.

"""

import os
import math
import logging
import cStringIO
import subprocess
import numpy as np
from PIL import Image

import medimg

log = logging.getLogger(__name__)


def get_tile(tiff, z, x, y, size=256):
    """
    Return an image tile as string from PNG.

    If the requested tile does not exist, return an empty (blacked out) tile.
    """
    def null_tile(size=256):
        """Create a blacked out tile."""
        tile = cStringIO.StringIO()
        img = Image.fromarray(np.ones((size, size), dtype=np.int8), mode='L')
        img.save(tile, 'png')
        tile.seek(0)
        return tile.read()

    # figure out which tile to get
    tiff_info = get_info(tiff)
    z = (len(tiff_info) - (z + 2))  # er. what?
    rows, cols = tiff_info.get(z)
    if x > rows or y > cols:
        return null_tile()

    # fetch the tile
    with Image.open(tiff) as i:
        i.seek(z)
        index = (y * (rows+1)) + x
        fd = cStringIO.StringIO()
        crop_param = i.tile[index][1]
        i.tile = [i.tile[index]]
        cropped = i.crop(crop_param)
        resized = cropped.resize((size, size))
        resized.save(fd, 'png')
        fd.seek(0)
        return fd.read()

def get_info(tiff):
    """
    Return info about the pyramidal tiff.

    Returns a dictionary of zoom level keys, with their max tile coordinates as values.
    """
    # TODO: reverse zoom level order, 0 = fully zoomed out, n = fully zoomed ini
    # to be consistent with how d3tiles requests zoom levels
    tiff_info = {}
    with Image.open(tiff) as i:
        for page in range(i.ifd.named().get('PageNumber')[1]):
            log.info('parsing page %d' % page)
            try:
                i.seek(page)
            except Image.EOFError:
                break  # no more zoom levels
            tags = i.ifd.named()
            tsize = tags.get('TileWidth') + tags.get('TileLength')
            rows = i.size[0] / tsize[0]
            cols = i.size[1] / tsize[1]
            tiff_info[page] = (rows, cols)
    return tiff_info

def generate_montage(imagedata, timepoints=[], bits16=False):
    """Generate a montage."""
    # Figure out the image dimensions and make an appropriate montage.
    # NIfTI images can have up to 7 dimensions. The fourth dimension is
    # by convention always supposed to be time, so some images (RGB, vector, tensor)
    # will have 5 dimensions with a single 4th dimension. For our purposes, we
    # can usually just collapse all dimensions above the 3rd.
    # TODO: we should handle data_type = RGB as a special case.
    # TODO: should we use the scaled data (getScaledData())? (We do some auto-windowing below)

    # This transpose (usually) makes the resulting images come out in a more standard orientation.
    # TODO: we could look at the qto_xyz to infer the optimal transpose for any dataset.
    # TODO: every tile should be a square.  resize array as necessary.
    # possibly: data.resize((np.max(data.shape[:2]), np.max(data.shape[:2])), data.shape[2]) ??
    data = imagedata.transpose(np.concatenate(([1, 0], range(2, imagedata.ndim))))
    num_images = np.prod(data.shape[2:])

    if data.ndim < 2:
        raise MontageError('NIfTI file must have at least 2 dimensions')
    elif data.ndim == 2:
        # a single slice: no need to do anything
        num_cols = 1
        data = np.atleast_3d(data)
    elif data.ndim == 3:
        # a simple (x, y, z) volume- set num_cols to produce a square(ish) montage.
        rows_to_cols_ratio = float(data.shape[0])/float(data.shape[1])
        num_cols = int(math.ceil(math.sqrt(float(num_images)) * math.sqrt(rows_to_cols_ratio)))
    elif data.ndim >= 4:
        # timeseries (x, y, z, t) or more
        num_cols = data.shape[2]
        data = data.transpose(np.concatenate(([0, 1, 3, 2], range(4, data.ndim)))).reshape(data.shape[0], data.shape[1], num_images)
        if len(timepoints) > 0:
            data = data[..., timepoints]

    num_rows = int(np.ceil(float(data.shape[2])/float(num_cols)))
    montage = np.zeros((data.shape[0] * num_rows, data.shape[1] * num_cols), dtype=data.dtype)
    for im_num in range(data.shape[2]):
        slice_r, slice_c = im_num / num_cols * data.shape[0], im_num % num_cols * data.shape[1]
        montage[slice_r:slice_r + data.shape[0], slice_c:slice_c + data.shape[1]] = data[:, :, im_num]

    # montage = montage.copy()        # is this necessary? need a deep copy?
    if montage.dtype == np.uint8 and bits16:
        montage = np.cast['uint16'](data)
    elif montage.dtype != np.uint8 or (montage.dtype != np.uint16 and bits16):
        montage = montage.astype(np.float32)  # do scaling/clipping with floats
        clip_vals = np.percentile(montage, (20.0, 99.0))   # auto-window the data by clipping
        montage = montage.clip(clip_vals[0], clip_vals[1]) - clip_vals[0]
        if bits16:
            montage = np.cast['uint16'](np.round(montage/(clip_vals[1]-clip_vals[0])*65535))
        else:
            montage = np.cast['uint8'](np.round(montage/(clip_vals[1]-clip_vals[0])*255.0))
    return montage

def generate_flat(imagedata, filepath):
    """Generate a flat png montage."""
    montage = generate_montage(imagedata)
    Image.fromarray(montage).convert('L').save(filepath, optimize=True)

    if not os.path.exists(filepath):
        raise MontageError('montage (flat png) not generated')
    else:
        log.debug('generated %s' % os.path.basename(filepath))
        return filepath


class MontageError(medimg.MedImgError):
    pass


class Montage(medimg.MedImgReader, medimg.MedImgWriter):

    domain = u'mr'
    filetype = u'montage'
    state = ['orig']

    def __init__(self, filepath, load_data=False):
        super(Montage, self).__init__(load_data=load_data)
        self.data = None        # contains montage
        if load_data:
            self.load_data(preloaded=True)

    def load_data(self, preloaded=False):
        super(Montage, self).load_data(preloaded=preloaded)
        log.debug('loading %s')
        if not preloaded:
            # read the data
            pass

    def get_tile(self, x, y, z):
        get_tile(self.filepath, x, y, z)

    def get_info(self):
        get_info(self.filepath)

    @classmethod
    def write(cls, metadata, imagedata, outbase, voxel_order=None, multi=False):
        """
        Write the metadata and imagedata to image montage pyramid.

        Parameters
        ----------
        metadata : object
            fully loaded instance of a Reader.
        imagedata : dict
            dictionary of np.darrays. label suffix as keys, with np.darrays as values.
        outbase : str
            output name prefix.
        voxel_order : str [default None]
            three character string indicating the voxel order, ex. 'LPS'.
        multi : bool [default False]
            True indicates to write multiple files. False only writes primary data in imagedata['']

        Returns
        -------
        results : list
            list of files written.

        Raises
        ------
        DataError
            metadata or data is None.

        """
        super(Montage, cls).write(metadata, imagedata, outbase, voxel_order)

        results = []
        for data_label, data in imagedata.iteritems():
            if not multi and data_label is not '':
                continue

            if data is None:
                continue
            data = imagedata.get(data_label)
            outname = outbase + data_label

            if voxel_order:
                data, _ = cls.reorder_voxels(data, metadata.qto_xyz, voxel_order)

            log.debug('type: flat png')
            png_result = generate_flat(data, outname + '.png')
            tiff_result = outname + '.tiff'
            x, y = data.shape[:2]
            convert_cmd = 'convert %s -compress LZW -define tiff:tile-geometry=%dx%d ptif:%s' % (png_result, x, y, tiff_result)
            subprocess.check_call(convert_cmd.split())
            if os.path.exists(tiff_result):
                result = tiff_result
            os.remove(png_result)  # remove the intermediate png

            results.append(result)
        return results

write = Montage.write
