# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.trim.
"""
from os.path import abspath, dirname
import unittest
from .. import __version__ as desimodel_version
from ..trim import inout, rebin_image


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


if __name__ == '__main__':
    unittest.main()
