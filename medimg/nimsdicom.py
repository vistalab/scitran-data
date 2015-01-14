# @author:  Gunnar Schaefer
#           Bob Dougherty
#           Kevin S Hahn

"""
nimsdata.medimg.nimsdicom
=========================

nimsdicom stuff, and adds dicom specific parsing routines.
nimsdicom uses a composer design pattern to dynamically incoporate functions based on
the input data.

"""
import dicom
import types
import logging
import tarfile
import datetime
import dcmstack
import cStringIO
import collections
import dcmstack.extract
import nibabel.nicom.csareader

import medimg

log = logging.getLogger(__name__)


dicom.config.auto_convert_VR_mismatch = True
nibabel.nicom.csareader.MAX_CSA_ITEMS = 300
dcmstack.DicomStack.sort_guesses.append('CsaImage.ImaCoilString')
try:
    dcmstack.DicomStack.sort_guesses.remove('InversionTime')
    dcmstack.DicomStack.minimal_keys.remove('PixelSpacing')
except ValueError:
    pass

MAX_LOC_DCMS = 150  # maximum number of dicoms allowed in a "localizer"

SUPPORTED_MFR = {
    'GE MEDICAL SYSTEMS': 'ge',
    'SIEMENS': 'siemens',
    }
SUPPORTED_SOP = {
    '1.2.840.10008.5.1.4.1.1.4': 'mr',  # MR Image
    '1.2.840.10008.5.1.4.1.1.7': 'sc',  # Secondary Capture
    '1.3.12.2.1107.5.9.1': 'syngo_csa',  # Private Syngo CSA Non-Image; SIEMENS ONLY
    '1.2.840.10008.5.1.4.1.1.88.22': 'enhanced_sr',
    '1.2.840.10008.5.1.4.1.1.128': 'pet',
    # '1.2.840.10008.5.1.4.1.1.130': 'enhanced_pet',
    # '1.2.840.10008.5.1.4.1.1.128.1': 'legacy_enhanced_pet',
    }


# XXX; not compatible with dcmstack.extract.MetaExtractor warn_on_ex parameter.
# catching the error prevents error from getting to MetaExtractor class where
# it was being converted to a warning.
# this change is necessary because some of the DTI from Davis and UI have a badly
# formed header field.  The raised exception prevents the remaining
# sections of mr phoenix protocol to not get parsed.
# this problem is located in Csa Series Header, 'sWiPMemBlock.tFree'
class Extractor(dcmstack.extract.MetaExtractor):

    """
    Override the default dcmstack.extract.MetaExtractor.

    If an elements value cannot be translated using the supplied Value Representation (VR),
    then return an empty string.  dcmstack.extract.MetaExtractor would normally raise a ValueError
    upon value-VR mismatch.

    """
    def _get_elem_value(self, elem):
        """
        Get the value for any non-translated elements.

        If element is not translateable with its own VR, then return value=None.
        """
        try:
            value = super(Extractor, self)._get_elem_value(elem)
        except ValueError:
            value = ''
        return value


def _parse_phoenix_prot(prot_key, prot_val):
    """Parse a siemens MrProtocol or MrPhoenixProtocol string."""
    if prot_key == 'MrPhoenixProtocol':         # syngo B
        str_delim = '""'
    elif prot_key == 'MrProtocol':              # syngo A
        str_delim = '"'
    else:
        raise ValueError('Unknown protocol key: %s' % prot_key)
    ascconv_start = prot_val.find('### ASCCONV BEGIN ###')
    ascconv_end = prot_val.find('### ASCCONV END ###')
    ascconv = prot_val[ascconv_start:ascconv_end].split('\n')[1:-1]
    result = collections.OrderedDict()
    for line in ascconv:
        # added try except. does not stop parsing when exception is raised.
        try:
            parse_result = dcmstack.extract._parse_phoenix_line(line, str_delim)
        except Exception:       # TODO: specific exception?
            log.debug('error parsing line %s...' % str(line[:20]))
        else:
            if parse_result:
                result[parse_result[0]] = parse_result[1]
    return result


def _csa_series_trans_func(elem):
    """Function for parsing the CSA series sub header element by element."""
    csa_dict = dcmstack.extract.simplify_csa_dict(nibabel.nicom.csareader.read(elem.value))
    # If there is a phoenix protocol, parse it and dump it into the csa_dict
    phx_src = None
    if 'MrPhoenixProtocol' in csa_dict:
        phx_src = 'MrPhoenixProtocol'
    elif 'MrProtocol' in csa_dict:
        phx_src = 'MrProtocol'
    if phx_src is not None:
        phoenix_dict = _parse_phoenix_prot(phx_src, csa_dict[phx_src])  # calls custom _parse_phoenix_prot
        del csa_dict[phx_src]
        for key, val in phoenix_dict.iteritems():
            new_key = '%s.%s' % ('MrPhoenixProtocol', key)
            csa_dict[new_key] = val
    return csa_dict


def _simplify_csa_dict(csa_dict):
    """
    Simplify the result of nibabel.nicom.csareader and normalize list content types.

    dcmstack.extract.csa_dict naively inserts values into the result without attempting
    to normalize list item types.  For example, if a value is a list, and contains mixed
    string and unicode types, the value will contain mixed string and unicode types. If the first
    item in the list is a string, and the remaining are unicode, dcmstack.extract.MetaExtractor.__call__
    will attempt to cast all values as unicode(in_str, 'utf-8'), which will raise exception.

    By converting all list items that are strings to unicode, the first list item will be recognized
    as unicode, and the second round (MetaExtactor.__call__) of unicode conversion will be skipped.

    Alternatively, could convert all list items that are unicode to string.  In this case, the first
    list item would be recognized as a str, and the second round unicode conversion
    (MetaExtractor.__call__) should succeeed.

    """
    if csa_dict is None:
        return None

    result = collections.OrderedDict()
    for tag in csa_dict['tags']:
        items = csa_dict['tags'][tag]['items']
        if len(items) == 0:
            continue
        elif len(items) == 1:
            result[tag] = items[0]
        else:
            # convert any strings to unicode,
            result[tag] = [(unicode(i) if isinstance(i, str) else i) for i in items]
    return result


def _csa_image_trans_func(elem):
    return _simplify_csa_dict(nibabel.nicom.csareader.read(elem.value))

_csa_series_trans = dcmstack.extract.Translator(
        'CsaSeries',
        dicom.tag.Tag(0x29, 0x1020),
        'SIEMENS CSA HEADER',
        _csa_series_trans_func)

_csa_image_trans = dcmstack.extract.Translator(
        'CsaImage',
        dicom.tag.Tag(0x29, 0x1010),
        'SIEMENS CSA HEADER',
        _csa_image_trans_func)

MetaExtractor = Extractor(
        ignore_rules=[dcmstack.extract.ignore_non_ascii_bytes],
        translators=[_csa_image_trans, _csa_series_trans]
        )


class NIMSDicomError(medimg.MedImgError):
    pass


class NIMSDicom(medimg.MedImgReader):

    """
    Parse a series of Dicoms.

    if load_data is set to True, the load_data method will be called at the end of init. scan identification
    is done during load_data. although some types can be determined with minimal data, identification is
    done at the same place across venders to provide consistency.

    parse gives a ball park idea of "what it is," without loading everything.

    NIMSDicom is a special composer class.  Upon instantiation, NIMSDicom evaluates the input
    file to determine it's image type, manufactuere and SOP class.  Using that information, NIMSDicom
    imports the appropriate group of functions, and assigns each function as a member to itself.

    NIMSDicom will initally parse the SOP Class UID, image type, and manufacturer, to determine
    what parsing needs to be done.  NIMSDicom will compose itself by selecting functions that are
    specific to a mfr.sop such as ge.mr, or siemens.syngo_csa.

    parameters
    ----------
    path : str
        path to input file
    load_data : bool [default False]
        indicate if all data should be loaded by invoking load_data() at the end of init.

    """

    domain = u'mr'  # for now assuming all dicoms are MR
    filetype = u'dicom'
    priority = 0
    parse_priority = 9
    state = ['orig']

    def __init__(self, path, load_data=False):
        """
        Parse a single file from the input.

        Performs manufacturer specific single dicom parsing.

        calls mfr sop specific parse_one method

        """
        super(NIMSDicom, self).__init__(path, load_data)
        with tarfile.open(path) as archive:
            for ti in archive:
                try:
                    self._hdr = MetaExtractor(dicom.read_file(cStringIO.StringIO(archive.extractfile(ti).read()), stop_before_pixels=True))
                except (dicom.filereader.InvalidDicomError, AttributeError, ValueError):
                    # dicom.filereader.InvalidDicomError, not a dicom
                    # AttributeError,
                    # ValueError, dcmstack.extract, tag value not parseable with declared VR
                    pass
                else:
                    break
            else:
                raise NIMSDicomError('no header could be extracted')  # XXX FAIL! unexpected to not extract header
        del ti, archive

        self.image_type = self.getelem(self._hdr, 'ImageType', None, [])
        self.manufacturer = self.getelem(self._hdr, 'Manufacturer')
        if not self.manufacturer:
            for sub in self.image_type:
                if sub.startswith('CSA'):
                    self.manufacturer = 'SIEMENS'
                    break
        self.sop_class_uid = self.getelem(self._hdr, 'SOPClassUID', str)

        if not self.manufacturer:
            log.warning('could not determine manufacturer from dicom header')
        if not self.sop_class_uid:
            log.warning('SOP Class UID not set in dicom header')

        mfr = SUPPORTED_MFR.get(self.manufacturer)
        sop = SUPPORTED_SOP.get(self.sop_class_uid)
        composer = '%s.%s.%s' % ('dcm', sop, mfr)

        try:
            _temp = __import__(composer, globals(), fromlist=['parse_one', 'parse_all', 'convert'])
        except (ImportError, AttributeError):
            log.warning('no composer match. parsing basic info only.')
        else:
            self.parse_one = types.MethodType(_temp.parse_one, self)  # types.MethodType(fxn, i), turn fxn into method of instance
            self.parse_all = types.MethodType(_temp.parse_all, self)
            self.convert = types.MethodType(_temp.convert, self)
            log.debug('composing from %s' % composer)

        # standard dicoms stuff
        self.exam_no = self.getelem(self._hdr, 'StudyID')
        self.exam_uid = self.getelem(self._hdr, 'StudyInstanceUID')
        self.patient_id = self.getelem(self._hdr, 'PatientID')
        self.series_no = self.getelem(self._hdr, 'SeriesNumber', int)
        self.series_desc = self.getelem(self._hdr, 'SeriesDescription')
        self.series_uid = self.getelem(self._hdr, 'SeriesInstanceUID')
        self.subj_code, self.group_name, self.project_name = medimg.parse_patient_id(self.patient_id, 'ex' + self.exam_no)
        self.acq_no = self.getelem(self._hdr, 'AcquisitionNumber', int, 1)  # wrong for siemens dicom, until load_data
        self.study_date = self.getelem(self._hdr, 'StudyDate')
        self.study_time = self.getelem(self._hdr, 'StudyTime')
        self.acq_date = self.getelem(self._hdr, 'AcquisitionDate')
        self.acq_time = self.getelem(self._hdr, 'AcquisitionTime')
        self.study_datetime = self.study_date and self.study_time and datetime.datetime.strptime(self.study_date + self.study_time[:6], '%Y%m%d%H%M%S')
        self.acq_datetime = self.acq_date and self.acq_time and datetime.datetime.strptime(self.acq_date + self.acq_time[:6], '%Y%m%d%H%M%S')
        self.timestamp = self.acq_datetime or self.study_datetime
        self.subj_firstname, self.subj_lastname = medimg.parse_patient_name(self.getelem(self._hdr, 'PatientName'))
        self.subj_dob = medimg.parse_patient_dob(self.getelem(self._hdr, 'PatientBirthDate', str))
        self.subj_sex = {'M': 'male', 'F': 'female'}.get(self.getelem(self._hdr, 'PatientSex'))
        self.scanner_name = '%s %s'.strip() % (
                self.getelem(self._hdr, 'InstitutionName', None, ''),
                self.getelem(self._hdr, 'StationName', None, ''),
                )
        self.manufacturer = self.getelem(self._hdr, 'Manufacturer')
        self.manufacturer_model = self.getelem(self._hdr, 'ManufacturerModelName')
        self.scanner_type = '%s %s'.strip() % (
                self.manufacturer,
                self.manufacturer_model or '',
                )
        self.operator = self.getelem(self._hdr, 'OperatorsName')
        self.phase_encode_direction = None  # FIXME: how do diff mfr's store phase_encode_direction??

        self.metadata_status = 'pending'
        self.parse_one()  # COMPOSED; parses mfr sop specific data

        if load_data:
            self.load_data()

    def parse_one(self):
        """No-op."""
        pass

    def parse_all(self):
        """No-op."""
        pass

    def convert(self):
        """No-op."""
        self.data = {}

    def load_data(self):
        """
        Load all dicoms and obtain more metadata.

        Performs manufacturer specific loads. Operates through side-effect, does not return anything.

        calls mfr sop specific parse_all and convert functions.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        super(NIMSDicom, self).load_data()
        self._dcm_list = []
        with tarfile.open(self.filepath) as archive:
            for ti in archive:
                try:
                    dcm = dicom.read_file(cStringIO.StringIO(archive.extractfile(ti).read()), stop_before_pixels=False)
                    if self.getelem(dcm, 'SOPClassUID') != self.sop_class_uid:  # mismatch SOP, do not attempt recon
                        log.error('dicoms have inconsistent SOP Class UIDs')  # XXX expected error
                        self.is_non_image = True
                    self._dcm_list.append(dcm)
                except (dicom.filereader.InvalidDicomError, AttributeError):
                    pass
        del ti, archive

        if self._dcm_list == []:
            raise NIMSDicomError('no dicoms loaded?')  # XXX FAIL! unexpected for no dicoms to be loaded
        try:
            self._dcm_list.sort(key=lambda dcm: dcm.InstanceNumber)
        except AttributeError:
            log.debug('dicoms do not have InstanceNumber. cannot pre-sort')

        self.parse_all()  # COMPOSED; parses mfr sop specifics
        self.metadata_status = 'complete'  # if parse_all completes, metadata is assumed to be completed

        try:
            self.convert()  # COMPOSED; converts mfr sop specific, may also do last round of metadata touch ups
        except Exception as e:
            log.debug('%s pixel data could not be loaded: %s' % (self.filepath, str(e)))
            self.data = None
            self.failure_reason = e


    @staticmethod
    def getelem(hdr, tag, type_=None, default=None):
        """
        Retrieve a dicom header element by its tag.

        Parameters
        ----------
        hdr : pydicom.dataset or csa_dict from dcmstack.extact.MetaExtractor
            Either a pydicom.dataset object, or a NestedMultiDict as returned by dcmstack.extract.MetaExtractor.
        tag : tuple or str
            A tuple or string that indicates which tag to retrieve.  The tuple tag can only be used with a pydicom.dataset.
        type_ : class
            Attempt to convert the tag's value to the specified type, such as str, int, float.  Default is None.
        default : any
            What value to return if the tag cannot be found, or if type conversion fails.

        Returns
        -------
        value : any
            the value held at the tag specificed, of type `type_` (if possible).  Or the default value if
            the tag does not exist, or could not be converted to `type_`.

        """
        value = default
        if isinstance(tag, tuple):
            # get by bytecode tuple
            elem = hdr.get(tag)
            if elem:
                value = elem.value
        else:
            # get by keyword
            value = hdr.get(tag, default)
        if type_ is not None:
            try:
                value = [type_(x) for x in value] if isinstance(value, list) else type_(value)
            except (TypeError, ValueError):   # TypeError, for Nones, ValueError for casting
                value = default
        return value
