# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.io.
"""
import os
import unittest
from .. import io
#
# Try to import specter.
#
specter_available = True
specter_message = "The specter package was not detected."
try:
    import specter
except ImportError:
    specter_available = False
#
# Try to load the DESIMODEL environment variable
#
desimodel_available = True
desimodel_message = "The desimodel data set was not detected."
try:
    spam = os.environ['DESIMODEL']
except KeyError:
    desimodel_available = False
    specter_available = False
    specter_message = desimodel_message


class TestIO(unittest.TestCase):
    """Test desimodel.io.
    """
    @classmethod
    def setUpClass(cls):
        global specter_available, desimodel_available
        cls.specter_available = specter_available
        cls.desimodel_available = desimodel_available
        # cls.data_dir = os.path.dirname(__file__)
        # try:
        #     cls.old_desimodel = os.environ['DESIMODEL']
        # except KeyError:
        #     cls.old_desimodel = None
        # os.environ['DESIMODEL'] = cls.data_dir

    @classmethod
    def tearDownClass(cls):
        pass
        # if cls.old_desimodel is None:
        #     del os.environ['DESIMODEL']
        # else:
        #     os.environ['DESIMODEL'] = cls.old_desimodel

    @unittest.skipUnless(specter_available, specter_message)
    def test_throughput(self):
        """Test loading of throughput files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_throughput(channel)

    @unittest.skipUnless(specter_available, specter_message)
    def test_psf(self):
        """Test loading of PSF files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_psf(channel)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_desiparams(self):
        """Test loading of basic DESI parameters.
        """
        p = io.load_desiparams()
        self.assertTrue(isinstance(p, dict))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_fiberpos(self):
        """Test loading of fiber positioner data.
        """
        fiberpos = io.load_fiberpos()
        self.assertEqual(len(fiberpos), 5000)
        # Check workaround for astropy 1.0.x bug for lower -> upper col names.
        for key in ('FIBER', 'POSITIONER', 'SPECTROGRAPH', 'X', 'Y', 'Z'):
            self.assertIn(key, fiberpos.dtype.names)
            x = fiberpos[key]

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_tiles(self):
        """Test loading of tile files.
        """
        t1 = io.load_tiles(onlydesi=False)
        t2 = io.load_tiles(onlydesi=True)
        self.assertLess(len(t2), len(t1))
        self.assertEqual(len(set(t2.TILEID) - set(t1.TILEID)), 0)


if __name__ == '__main__':
    unittest.main()
