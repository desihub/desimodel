# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.inputs
"""
import unittest
import os
from requests.auth import HTTPDigestAuth
from ..inputs import docdb, fiberpos, gfa, throughput

desimodel_available = True
desimodel_message = "The desimodel data set was not detected."
try:
    spam = os.environ['DESIMODEL']
except KeyError:
    desimodel_available = False

skipMock = False
try:
    from unittest.mock import patch, MagicMock
except ImportError:
    # Python 2
    skipMock = True


class TestInputs(unittest.TestCase):
    """Test desimodel.inputs
    """
    def setUp(self):
        pass

    def test_docdb__xls_col2int(self):
        self.assertEqual(docdb._xls_col2int('A'), 0)
        self.assertEqual(docdb._xls_col2int('AA'), 26)

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_docdb__auth(self):
        with patch('netrc.netrc') as netrc:
            n = netrc.return_value
            n.authenticators.return_value = None
            with self.assertRaises(ValueError) as e:
                h = docdb._auth('www.example.com')
            self.assertEqual(str(e.exception), 'Unable to get user/pass from $HOME/.netrc for www.example.com')
            n.authenticators.return_value = ('username', 'foo', 'password')
            h = docdb._auth('www.example.com')
            try:
                self.assertIsInstance(h, HTTPDigestAuth)
            except AttributeError:
                # Python 2
                pass
            n.authenticators.assert_called_with('www.example.com')

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_build_gfa_table(self):
        gfa.build_gfa_table(testdir='.')
        self.assertTrue(os.path.exists('gfa.ecsv'), "Test Failed to create a valid file!")
        os.remove('gfa.ecsv')


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
