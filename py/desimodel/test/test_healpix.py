# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.healpix
"""
from __future__ import print_function, division

import unittest
import numpy as np
from desimodel.healpix import tiles2pix

class TestHealpix(unittest.TestCase):
    """Test desimodel.healpix.
    """
    def test_tiles2pix(self):
        """Test that the correct healpix tiles are returned for a couple values of nside.
        """
        tiles = np.zeros(3, dtype=[('RA', float), ('DEC', float)])
        tiles['RA'] = [ 335.03,  333.22,  332.35]
        tiles['DEC'] = [ 19.88,  14.84,  12.32]

        pix = tiles2pix(nside=8, tiles=tiles[:3], radius=1.6)
        print('\n--- nside 8 pix', pix)
        print('--- compare', pix == np.array([196,197,208,302,303]))
        print('--- all', np.all(pix == np.array([196,197,208,302,303])))
        self.assertTrue( np.all(pix == np.array([196,197,208,302,303])) )

        pix = tiles2pix(nside=16, tiles=tiles[:3], radius=1.6)
        self.assertTrue( np.all(pix == np.array([785,788,789,791,832,1209,1211,1214,1215])) )

def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
