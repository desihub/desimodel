# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.trim.
"""
import unittest
from os.path import abspath, dirname
import numpy as np
from ..trim import (inout, rebin_image, trim_focalplane, trim_inputs,
                    trim_sky, trim_targets)

skipMock = False
try:
    from unittest.mock import patch, MagicMock
except ImportError:
    # Python 2
    skipMock = True


class TestTrim(unittest.TestCase):
    """Test desimodel.trim.
    """

    def test_inout(self):
        """Test pathname helper function.
        """
        self.assertEqual(inout('/in', '/out', 'file.txt'),
                         ('/in/file.txt', '/out/file.txt'))

    def test_rebin_image(self):
        """Test image rebinning.
        """
        image = np.ones((10, 10), dtype='i2')
        i1 = rebin_image(image, 5)
        i2 = np.array([[25, 25], [25, 25]], dtype='i2')
        self.assertTrue((i1 == i2).all())

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_focalplane(self):
        """Test trim_focalplane().
        """
        with patch('shutil.copytree') as copytree:
            trim_focalplane('/in/focalplane', '/out/focalplane')
            copytree.assert_called_with('/in/focalplane', '/out/focalplane')

    def test_trim_inputs(self):
        """Test trim_inputs().
        """
        # This function does nothing, so we just check that we can call it.
        trim_inputs('/in/focalplane', '/out/focalplane')

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_targets(self):
        """Test trim_targets().
        """
        with patch('shutil.copytree') as copytree:
            trim_targets('/in/targets', '/out/targets')
            copytree.assert_called_with('/in/targets', '/out/targets')

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_targets(self):
        """Test trim_targets().
        """
        with patch('os.makedirs'):
            with patch('shutil.copy') as copy:
                trim_sky('/in/sky', '/out/sky')
                copy.assert_called_with('/in/sky/solarspec.txt',
                                        '/out/sky/solarspec.txt')


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
