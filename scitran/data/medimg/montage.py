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
import json
import math
import logging
import sqlite3
import zipfile
import cStringIO
import numpy as np
from PIL import Image

import medimg

log = logging.getLogger(__name__)


def get_tile(montagezip, z, x, y):
    """Get a specific image tile from an sqlite db."""
    try:
        with zipfile.ZipFile(montagezip, 'r') as zf:
            for tile in zf.namelist():
                if tile.endswith('z%03d/x%03d_y%03d.jpg' % (z, x, y)):
                    return zf.open(tile).read()
            else:
                raise IndexError
    except zipfile.BadZipfile:
        log.error('bad zip file')
    except IndexError:
        log.error('tile does not exist')


def get_info(montagezip):
    """Return the tile_size, x_size, and y_size from the sqlite pyramid db."""
    try:
        with zipfile.ZipFile(montagezip, 'r') as zf:
            info = json.loads(zf.open(zf.namelist()[0]).read())
    except zipfile.BadZipfile:
        log.error('bad zip file')
    except ValueError:
        log.error('tileinfo.json could not be found')
    return info['tile_size'], info['zoom_levels']


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


def generate_pyramid(montage, tile_size):
    """
    Slice up a NIfTI file into a multi-res pyramid of tiles.

    We use the file name convention suitable for d3tiles
    The zoom level (z) is an integer between 0 and n, where 0 is fully zoomed out and n is zoomed in.
    E.g., z=0 is for 1 tile covering the whole world, z=1 is for 2x2=4 tiles, ... z=n is the original resolution.

    """
    montage_image = Image.fromarray(montage, 'L')
    montage_image = montage_image.crop(montage_image.getbbox())  # crop away edges that contain only zeros
    sx, sy = montage_image.size
    if sx * sy < 1:
        raise MontageError('degenerate image size (%d, %d): no tiles will be created' % (sx, sy))
    if sx < tile_size and sy < tile_size:  # Panojs chokes if the lowest res image is smaller than the tile size.
        tile_size = max(sx, sy)

    pyramid = {}
    pyramid_meta = {
        'tile_size': tile_size,
        'mimetype': 'image/jpeg',
        'zoom_levels': {},
    }
    divs = max(1, int(np.ceil(np.log2(float(max(sx, sy))/tile_size))) + 1)
    for z in range(divs):
        # flip the z label to be d3 friendly
        level = (divs - 1) - z
        ysize = int(round(float(sy)/pow(2, z)))
        xsize = int(round(float(ysize)/sy*sx))
        xpieces = int(math.ceil(float(xsize)/tile_size))
        ypieces = int(math.ceil(float(ysize)/tile_size))
        log.debug('level %s, size %dx%d, splits %d,%d' % (level, xsize, ysize, xpieces, ypieces))
        # TODO: we don't need to use 'thumbnail' here. This function always returns a square
        # image of the requested size, padding and scaling as needed. Instead, we should resize
        # and chop the image up, with no padding, ever. panojs can handle non-square images
        # at the edges, so the padding is unnecessary and, in fact, a little wrong.
        im = montage_image.copy()
        im.thumbnail([xsize, ysize], Image.ANTIALIAS)
        im = im.convert('L')    # convert to grayscale
        for x in range(xpieces):
            for y in range(ypieces):
                tile = im.copy().crop((x*tile_size, y*tile_size, min((x+1)*tile_size, xsize), min((y+1)*tile_size, ysize)))
                buf = cStringIO.StringIO()
                tile.save(buf, 'JPEG', quality=85)
                pyramid[(level, x, y)] = buf
        pyramid_meta['zoom_levels'][level] = (xpieces, ypieces)
    return pyramid, montage_image.size, pyramid_meta


def generate_dir_pyr(imagedata, outbase, tile_size=256):
    """Generate a panojs image pyramid directory."""
    montage = generate_montage(imagedata)
    pyramid, pyramid_size, pyramid_meta = generate_pyramid(montage, tile_size)

    # write directory pyramid
    image_path = os.path.join(outbase, 'images')
    if not os.path.exists(image_path):
        os.makedirs(image_path)
        for idx, tile_buf in pyramid.iteritems():
            with open(os.path.join(image_path, ('%03d_%03d_%03d.jpg' % idx)), 'wb') as fp:
                fp.write(tile_buf.getvalue())

    # check for one image, pyramid file
    if not os.path.exists(os.path.join(outbase, 'images', '000_000_000.jpg')):
        raise MontageError('montage (flat png) not generated')
    else:
        log.debug('generated %s' % outbase)
        return outbase

def generate_zip_pyr(imagedata, outbase, tile_size=256):
    montage = generate_montage(imagedata)
    pyramid, pyramid_size, pyramid_meta = generate_pyramid(montage, tile_size)
    zip_name = outbase + '.zip'
    with zipfile.ZipFile(zip_name, 'w', compression=zipfile.ZIP_STORED) as zf:
        metaname = os.path.join(os.path.basename(outbase), 'tileinfo.json')
        zf.writestr(metaname, json.dumps(pyramid_meta))
        for idx, tile_buf in pyramid.iteritems():
            tilename = 'z%03d/x%03d_y%03d.jpg' % idx
            arcname = os.path.join(os.path.basename(outbase), tilename)
            zf.writestr(arcname, tile_buf.getvalue())
    return zip_name

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
    def write(cls, metadata, imagedata, outbase, voxel_order=None, mtype='zip', tilesize=512, multi=False):
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
        mtype : str [default 'sqlite']
            type of montage to create. can be 'sqlite', 'dir', or 'png'.
        tilesize : int [default 512]
            tilesize for generated sqlite or directory pyramid. Has no affect on mtype 'png'.
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
            if mtype == 'png':
                log.debug('type: flat png')
                result = generate_flat(data, outname + '.png')
            elif mtype == 'dir':
                log.debug('type: directory')
                result = generate_dir_pyr(data, outname, tilesize)
            elif mtype == 'zip':
                log.debug('type: zip of tiles')
                result = generate_zip_pyr(data, outname, tilesize)
            else:
                raise MontageError('montage mtype must be sqlite, dir or png. not %s' % mtype)

            results.append(result)
        return results

write = Montage.write
