# @author:  Gunnar Schaefer
#           Kevin S Hahn
"""
NIMSMRPNG provides PNG image writing capabilities for MR datasets read by any subclass of NIMSMRReader.

NIMSPNG is a subclass of NIMSMRWriter, accepts data as a list of data components.
"""

import os
import Image
import logging
import numpy as np

import nimsmrdata

log = logging.getLogger(__name__)


class NIMSPNGError(nimsmrdata.NIMSMRDataError):
    pass


class NIMSPNG(nimsmrdata.NIMSMRWriter):

    datakind = u'derived'
    datatype = u'bitmap'
    filetype = u'png'

    @classmethod
    def write(cls, metadata, imagedata, outbase, voxel_order='LPS'):
        """Create png files for each image in a list of pixel data."""
        super(NIMSPNG, cls).write(metadata, imagedata, outbase, voxel_order)
        # nimsmrdata.NIMSMRWriter.write(metadata, imagedata, outbase, voxel_order)
        results = []
        for data_label, data in imagedata.iteritems():
            if data is None:
                continue
            # cannot reorder voxels if there is no affine
            # if voxel_order:
            #     data, qto_xyz = cls.reorder_voxels(data, metadata.qto_xyz, voxel_order)
            # else:
            #     qto_xyz = metadata.qto_xyz
            outname = outbase + data_label

            # cut the depth array
            data = np.dsplit(data, len(metadata._dcm_list))

            # squeeze arrays, remove axis that have 1 value
            data = [image.squeeze() for image in data]

            for i, data in enumerate(data):
                filepath = outname + '_%d' % (i + 1) + '.png'
                if data.ndim == 2:
                    data = data.astype(np.int32)
                    data = data.clip(0, (data * (data != (2**15 - 1))).max())   # -32768->0; 32767->brain.max
                    data = data * (2**8 - 1) / data.max()                             # scale to full 8-bit range
                    Image.fromarray(data.astype(np.uint8), 'L').save(filepath, optimize=True)
                elif data.ndim == 3:
                    data = data.reshape((data.shape[1], data.shape[2], data.shape[0]))
                    Image.fromarray(data, 'RGB').save(filepath, optimize=True)
                log.info('generated %s' % os.path.basename(filepath))

            log.debug('returning:  %s' % filepath)

            results.append(filepath)
        return results

write = NIMSPNG.write
