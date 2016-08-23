# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.trim.
"""
from os.path import abspath, dirname
import unittest
from .. import __version__ as desimodel_version
from ..trim import inout


class TestTrim(unittest.TestCase):
    """Test desimodel.trim.
    """

    def test_inout(self):
        """Test pathname helper function.
        """
        self.assertEqual(inout('/in', '/out', 'file.txt'),
                         ('/in/file.txt', '/out/file.txt'))


if __name__ == '__main__':
    unittest.main()
