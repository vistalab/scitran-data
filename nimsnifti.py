# @author:  Bob Dougherty
#           Gunnar Schaefer
#           Kevin S Hahn

import os
import bson
import logging
import nibabel

import numpy as np

import nimsdata

log = logging.getLogger('nimsnifti')


class NIMSNiftiError(nimsdata.NIMSDataError):
    pass


class NIMSNifti(nimsdata.NIMSData):

    filetype = u'nifti'
    priority = 0
    parse_priority = 9

    # TODO: implement NIMSnifti as a parser
    def __init__(self, archive):
        super(NIMSNifti, self).__init__()
        try:
            nifti = nibabel.load(archive)
            self.imagedata = nifti.get_data()
            super(NIMSNifti, self).__init__()
        except Exception as e:
            raise NIMSNiftiError(str(e))

    @property
    def nims_group(self):
        return self.metadata.group

    @property
    def nims_experiment(self):
        return self.metadata.experiment

    @property
    def nims_session(self):
        return self.metadata.exam_uid.replace('.','_')

    @property
    def nims_epoch(self):
        return self.metadata.epoch

    @property
    def nims_type(self):
        return ('derived', 'nifti', self.filetype)

    @property
    def nims_filename(self):
        return self.nims_epoch + '_' + self.filetype

    @property
    def nims_timestamp(self):   # FIXME: should return UTC time and timezone
        return self.timestamp.replace(tzinfo=bson.tz_util.FixedOffset(-7*60, 'pacific'))  # FIXME: use pytz

    @property
    def nims_timezone(self):
        return None

    def load_data(self):
        pass

    def convert(self):
        pass

    @staticmethod
    def write(metadata, imagedata, outbase, notes=''):
        """Create a nifti file and possibly bval and bvec files from an ordered list of pixel data."""
        # write primary nifti header
        if notes != '':
            filepath = outbase + '_README.txt'
            with open(filepath, 'w') as fp:
                fp.write(notes)
            log.debug('generated %s' % os.path.basename(filepath))

        # if getattr(metadata, 'bvals', None) and getattr(metadata, 'bvecs', None):
        if metadata.bvals is not None and metadata.bvecs is not None:
            log.debug('dwi metadata found')
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

        # create image WITHOUT a header....
        nifti = nibabel.Nifti1Image(imagedata, None)
        nii_header = nifti.get_header()  # get the empty header and update the shizzles out of it
        num_slices = imagedata.shape[2]  # Don't trust metatdata.num_slices; # of resulting slices might not match the # acquired.
        nii_header.set_xyzt_units('mm', 'sec')
        nii_header.set_qform(metadata.qform, 'scanner')
        nii_header.set_sform(metadata.sform, 'scanner')
        nii_header.set_dim_info(*([1, 0, 2] if metadata.phase_encode == 0 else [0, 1, 2]))
        nii_header['pixdim'][4] = metadata.tr
        nii_header['slice_start'] = 0
        nii_header['slice_end'] = num_slices - 1
        nii_header.set_slice_duration(metadata.slice_duration)
        nii_header['slice_code'] = metadata.slice_order
        if np.iscomplexobj(imagedata):
            clip_vals = np.percentile(np.abs(imagedata), (10.0, 99.5))
        else:
            clip_vals = np.percentile(imagedata, (10.0, 99.5))
        nii_header.structarr['cal_min'] = clip_vals[0]
        nii_header.structarr['cal_max'] = clip_vals[1]
        nii_header.set_data_dtype(imagedata.dtype)
        # Stuff some extra data into the description field (max of 80 chars)
        # Other unused fields: nii_header['data_type'] (10 chars), nii_header['db_name'] (18 chars),
        nii_header['descrip'] = 'te=%.2f;ti=%.0f;fa=%.0f;ec=%.4f;acq=[%s];mt=%.0f;rp=%.1f;' % (
               metadata.te * 1000.,
               metadata.ti * 1000.,
               metadata.flip_angle,
               metadata.effective_echo_spacing * 1000.,
               ','.join(map(str, metadata.acquisition_matrix)),
               metadata.mt_offset_hz,
               1. / metadata.phase_encode_undersample,
               )
        if '3D' in metadata.acquisition_type:   # for 3D acquisitions, add the slice R-factor
            nii_header['descrip'] = str(nii_header['descrip']) + 'rs=%.1f' % (1. / metadata.slice_encode_undersample)

        nifti.update_header()

        filepath = outbase + '.nii.gz'
        nibabel.save(nifti, filepath)
        log.debug('generated %s' % os.path.basename(filepath))
        return filepath
