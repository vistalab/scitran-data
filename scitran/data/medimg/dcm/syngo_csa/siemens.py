# dcm.syngo_csa.siemens

"""
nimsdata.medimg.dcm.syngo_csa.siemens
=====================================

Not implemented.

"""

import logging

log = logging.getLogger(__name__)


SIEMENS_TYPE_UNKNOWN = ['ORIGINAL', 'PRIMARY']
SIEMENS_TYPE_DIFF_FA = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'FA', 'ND']
SIEMENS_TYPE_DIFF_TENSOR = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'TENSOR', 'ND']
SIEMENS_TYPE_DIFF_FA_NORM = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'FA', 'ND', 'NORM']
SIEMENS_TYPE_DIFF_TENSOR_NORM = ['DERIVED', 'PRIMARY', 'DIFFUSION', 'TENSOR', 'ND', 'NORM']


def parse_one(self):
    log.debug(self.image_type)
    self.is_non_image = True
    # do a bunch of generic parsing things
    pass


def parse_all(self):
    # TODO:  wtf does a syngo_csa even hold?
    pass


def convert(self):
    # TODO: figure out how to convert some of these SC dicoms
    pass

