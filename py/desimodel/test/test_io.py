# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.io.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import unittest
from .. import io
#
# Try to import specter.
#
specter_available = True
try:
    import specter
except ImportError:
    specter_available = False
#
#
#
class TestIO(unittest.TestCase):
    """Test desimodel.io.
    """
    @classmethod
    def setUpClass(cls):
        global specter_available
        cls.specter_available = specter_available
        cls.data_dir = os.path.dirname(__file__)
        try:
            cls.old_desimodel = os.environ['DESIMODEL']
        except KeyError:
            cls.old_desimodel = None
        os.environ['DESIMODEL'] = cls.data_dir

    @classmethod
    def tearDownClass(cls):
        if cls.old_desimodel is None:
            del os.environ['DESIMODEL']
        else:
            os.environ['DESIMODEL'] = cls.old_desimodel

    @unittest.skipUnless(specter_available,"The specter package was not detected.")
    def test_throughput(self):
        """Test loading of throughput files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_throughput(channel)

    @unittest.skipUnless(specter_available,"The specter package was not detected.")
    def test_psf(self):
        """Test loading of PSF files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_psf(channel)

    def test_desiparams(self):
        """Test loading of basic DESI parameters.
        """
        p = io.load_desiparams()
        self.assertTrue(isinstance(p, dict))

    def test_fiberpos(self):
        """Test loading of fiber positioner data.
        """
        fiberpos = io.load_fiberpos()
        self.assertEqual(len(fiberpos), 5000)
        #- Check workaround for astropy 1.0.x bug for lower -> upper col names
        for key in ('FIBER', 'POSITIONER', 'SPECTROGRAPH', 'X', 'Y', 'Z'):
            self.assertIn(key, fiberpos.dtype.names)
            x = fiberpos[key]

    def test_tiles(self):
        """Test loading of tile files.
        """
        t1 = io.load_tiles(onlydesi=False)
        t2 = io.load_tiles(onlydesi=True)
        self.assertLess(len(t2), len(t1))
        self.assertEqual(len(set(t2.TILEID) - set(t1.TILEID)), 0)

#- This runs all test* functions in any TestCase class in this file
if __name__ == '__main__':
    unittest.main()
