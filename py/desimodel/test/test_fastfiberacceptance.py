# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.fastfiberacceptance
"""

import unittest
import numpy as np
from desimodel.fastfiberacceptance import FastFiberAcceptance

class TestFiberAcceptance(unittest.TestCase):
    """Test FastFiberAcceptance
    """

    def test_fastfiberacceptance(self):
        """Test interface options for tiles2pix"""

        fa = FastFiberAcceptance()
        platescale=70. #um/arsec
        fwhm_arcsec_to_sigma_um = platescale/2.35

        # only a test of type conversion
        val = fa.value("POINT",sigmas=1.1*fwhm_arcsec_to_sigma_um)
        assert(np.isscalar(val))
        val = fa.value("POINT",sigmas=[1.1*fwhm_arcsec_to_sigma_um,1.4*fwhm_arcsec_to_sigma_um])
        assert(val.size==2)
        val = fa.value("POINT",sigmas=np.linspace(1.0,1.5,3)*fwhm_arcsec_to_sigma_um)
        assert(val.size==3)
        val = fa.value("POINT",sigmas=np.linspace(1.0,1.5,3)*fwhm_arcsec_to_sigma_um,offsets=np.zeros(3))
        assert(val.size==3)

def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)

if __name__ == '__main__':
    unittest.main()
