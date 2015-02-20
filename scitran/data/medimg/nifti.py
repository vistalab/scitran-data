# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S Hahn

"""
scitran.data.medimg.nifti
=========================

Nifti provide NIfti writing capabilities for MR datasets read by any subclass of MedImgReader.

Provides nifti specifics, inherits from MedImgReader, MedImgWriter.

"""

import os
import bson
import logging
import nibabel

import numpy as np

import medimg

log = logging.getLogger(__name__)

# NIFITI1-style slice order codes:
SLICE_ORDER_UNKNOWN = 0
SLICE_ORDER_SEQ_INC = 1
SLICE_ORDER_SEQ_DEC = 2
SLICE_ORDER_ALT_INC = 3
SLICE_ORDER_ALT_DEC = 4
SLICE_ORDER_ALT_INC2 = 5  # interleave, ascending, starting at 2nd slice
SLICE_ORDER_ALT_DEC2 = 6  # interleave, decreasing, starting at 2nd to last slice


class NiftiError(medimg.MedImgError):
    pass


class Nifti(medimg.MedImgReader, medimg.MedImgWriter):

    """
    Read elements of a MedImgReader subclass.

    Dataset must have a 'data' attribute that contains voxel data in a np.darray. Dataset must
    also have contain metadata attributes to define required attributes.

    Parameters
    ----------
    path : str
        filepath of input nifti, in .nii or .nii.gz format.
    load_data : bool [default False]
        attempt to load all data. has no affect.

    Returns
    -------
    None : NoneType

    Raises
    ------
    NiftiError
        TODO: explain plz.

    """

    domain = u'mr'
    filetype = u'nifti'
    state = ['orig']

    def __init__(self, path, load_data=False):
        super(Nifti, self).__init__(path, load_data)
        try:
            # TODO: load header only, unless load_data = True
            self.nifti = nibabel.load(path)
        except Exception as e:
            raise NiftiError(e)
        # parse some header info from the nifti
        # self.metadata._hdr = get header
        # first simple parse of _hdr

    def load_data(self, preloaded):
        super(Nifti, self).load_data(preloaded)
        if not preloaded:
            try:
                nifti = nibabel.load(self.filepath)
            except Exception as e:
                raise NiftiError(e)

        # TODO: nibabel nifti header reader
        # self.metadata.group
        # self.metadata.project
        # self.metadata.exam_uid
        self.data = nifti.imagedata.squeeze()
        self.qto_xyz = nifti.get_affine()
        self.sform = nifti.get_sform()
        self.qform = nifti.get_qform()
        self.image_type = ['derived', 'nifti', self.filetype]
        self.scan_type = 'unknown'   # FIXME

    @property
    def nims_group(self):
        return self.metadata.group

    @property
    def nims_project(self):
        return self.metadata.project

    @property
    def nims_session(self):
        return self.metadata.exam_uid.replace('.', '_')

    @property
    def nims_session_name(self):
        return self.metadata.timestamp.strftime('%Y-%m-%d %H:%M') if self.metadata.series_no == 1 and self.metadata.acq_no == 1 else None

    @property
    def nims_session_subject(self):
        return self.metadata.subj_code

    @property
    def nims_acquisition(self):
        return self.metadata.acquisition

    @property
    def nims_acquisition_name(self):
        pass

    @property
    def nims_acquisition_description(self):
        pass

    @property
    def nims_file_name(self):
        return self.nims_acquisition + '_' + self.filetype

    @property
    def nims_file_ext(self):
        return '.tgz'

    @property
    def nims_file_domain(self):
        return self.domain

    @property
    def nims_file_type(self):
        return self.filetype

    @property
    def nims_file_kinds(self):
        return self.scan_type

    @property
    def nims_file_state(self):
        return self.state

    @property
    def nims_timestamp(self):  # FIXME: should return UTC time and timezone
        return self.timestamp.replace(tzinfo=bson.tz_util.FixedOffset(-7 * 60, 'pacific'))  # FIXME: use pytz

    @property
    def nims_timezone(self):
        return None

    @classmethod
    def write(cls, metadata, imagedata, outbase, voxel_order=None):
        """
        Write the metadata and imagedata to niftis.

        Constructs a description from the metadata, and applies as much metadata into the nifti
        header as is applicable.  Creates .bvec and .bval files if bvecs and bvals exist.

        Parameters
        ----------
        metadata : object
            fully loaded instance of a MedImgReader.
        imagedata : dict
            dictionary of np.darrays. label suffix as keys, with np.darrays as values.
        outbase : str
            output name prefix.
        voxel_order : str [default None]
            three character string indicating the voxel order, ex. 'LPS'.

        Returns
        -------
        results : list
            list of files written.

        Raises
        ------
        DataError
            metadata or data is None.

        """
        super(Nifti, cls).write(metadata, imagedata, outbase, voxel_order)  # XXX FAIL! unexpected imagedata = None
        results = []
        for data_label, data in imagedata.iteritems():
            if data is None:
                continue
            if voxel_order:
                data, qto_xyz = cls.reorder_voxels(data, metadata.qto_xyz, voxel_order)
            else:
                qto_xyz = metadata.qto_xyz
            outname = outbase + data_label

            log.debug('creating nifti for %s' % data_label)

            # TODO: medimg.adjust_bvecs to use affine from after reorient
            if metadata.is_dwi and metadata.bvals is not None and metadata.bvecs is not None:
                filepath = outbase + '.bval'
                with open(filepath, 'w') as bvals_file:
                    bvals_file.write(' '.join(['%0.1f' % value for value in metadata.bvals]))
                log.debug('generated %s' % os.path.basename(filepath))
                filepath = outbase + '.bvec'
                with open(filepath, 'w') as bvecs_file:
                    bvecs_file.write(' '.join(['%0.4f' % value for value in metadata.bvecs[0, :]]) + '\n')
                    bvecs_file.write(' '.join(['%0.4f' % value for value in metadata.bvecs[1, :]]) + '\n')
                    bvecs_file.write(' '.join(['%0.4f' % value for value in metadata.bvecs[2, :]]) + '\n')
                log.debug('generated %s' % os.path.basename(filepath))

            # write actually nifti
            nifti = nibabel.Nifti1Image(data, None)
            nii_header = nifti.get_header()
            nifti.update_header()               # XXX are data and header ever "non-harmonious"
            num_slices = data.shape[2]          # Don't trust metatdata.num_slices; might not match the # acquired.
            nii_header.set_xyzt_units('mm', 'sec')
            nii_header.set_qform(qto_xyz, 'scanner')
            nii_header.set_sform(qto_xyz, 'scanner')
            nii_header.set_dim_info(*([1, 0, 2] if metadata.phase_encode == 0 else [0, 1, 2]))
            nii_header['slice_start'] = 0
            nii_header['slice_end'] = num_slices - 1

            nii_header.set_slice_duration(metadata.slice_duration)
            nii_header['slice_code'] = metadata.slice_order
            if np.iscomplexobj(data):
                clip_vals = np.percentile(np.abs(data), (10.0, 99.5))
            else:
                clip_vals = np.percentile(data, (10.0, 99.5))
            nii_header.structarr['cal_min'] = clip_vals[0]
            nii_header.structarr['cal_max'] = clip_vals[1]
            nii_header.set_data_dtype(data.dtype)

            # Stuff some extra data into the description field (max of 80 chars)
            # Other unused fields: nii_header['data_type'] (10 chars), nii_header['db_name'] (18 chars),
            te = 0 if not metadata.te else metadata.te
            ti = 0 if not metadata.ti else metadata.ti
            flip_angle = 0 if not metadata.flip_angle else metadata.flip_angle
            effective_echo_spacing = 0. if not metadata.effective_echo_spacing else metadata.effective_echo_spacing
            acquisition_matrix = [0, 0] if metadata.acquisition_matrix == (None, None) else metadata.acquisition_matrix
            mt_offset_hz = 0. if not metadata.mt_offset_hz else metadata.mt_offset_hz
            phase_encode_undersample = 1 if not metadata.phase_encode_undersample else metadata.phase_encode_undersample
            slice_encode_undersample = 1 if not metadata.slice_encode_undersample else metadata.slice_encode_undersample
            nii_header['descrip'] = 'te=%.2f;ti=%.0f;fa=%.0f;ec=%.4f;acq=[%s];mt=%.0f;rp=%.1f;' % (
                    te * 1000.,
                    ti * 1000.,
                    flip_angle,
                    effective_echo_spacing * 1000.,
                    ','.join(map(str, acquisition_matrix)),
                    mt_offset_hz,
                    1. / phase_encode_undersample,
                    )
            if '3D' in (metadata.acquisition_type or ''):
                nii_header['descrip'] = str(nii_header['descrip']) + 'rs=%.1f' % (1. / slice_encode_undersample)
            if metadata.phase_encode_direction != None:
                nii_header['descrip'] = str(nii_header['descrip']) + 'pe=%d' % (metadata.phase_encode_direction)
            if metadata.is_fastcard:
                nii_header['descrip'] = str(nii_header['descrip']) + 've=%f' % (metadata.velocity_encode_scale) or 0.

            nii_header['pixdim'][4] = metadata.tr   # XXX pixdim[4] = TR, even when non-timeseries. not nifti compliant

            filepath = outname + '.nii.gz'
            nibabel.save(nifti, filepath)
            log.debug('generated %s' % os.path.basename(filepath))
            results.append(filepath)

        return results

write = Nifti.write
