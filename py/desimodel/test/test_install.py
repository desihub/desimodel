# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.install.
"""
import os
from subprocess import CalledProcessError
import unittest
from ..install import default_install_dir, assert_svn_exists, svn_export, install
import desimodel

skipMock = False
try:
    from unittest.mock import patch, MagicMock
except ImportError:
    # Python 2
    skipMock = True


class TestInstall(unittest.TestCase):
    """Test desimodel.install.
    """

    def setUp(self):
        self.orig_version = desimodel.__version__

    def tearDown(self):
        desimodel.__version__ = self.orig_version

    def test_default_install_dir(self):
        """Test setting default install directory.
        """
        d = os.path.dirname
        d1 = default_install_dir()
        d2 = d(d(__file__))
        self.assertEqual(os.path.abspath(d1), os.path.abspath(d2))

    def test_svn_export(self):
        """Test svn export command.
        """
        base_url = "https://desi.lbl.gov/svn/code/desimodel/{0}/data"

        desimodel.__version__ = '1.2.3.dev456'
        cmd = svn_export()
        self.assertEqual(cmd[2], base_url.format('trunk'))

        desimodel.__version__ = '7.8.9'
        cmd = svn_export()
        self.assertEqual(cmd[2], base_url.format('tags/7.8.9'))

        desimodel.__version__ = '7.10'
        cmd = svn_export()
        self.assertEqual(cmd[2], base_url.format('tags/7.10'))

        desimodel.__version__ = 'blatfoo'
        cmd = svn_export()
        self.assertEqual(cmd[2], base_url.format('trunk'))

        # code version is tag-like, but we purposefully want to install a different version
        desimodel.__version__ = '0.19.3'
        cmd = svn_export('branches/test-0.19')
        self.assertEqual(cmd[2], base_url.format('branches/test-0.19'))

        desimodel.__version__ = self.orig_version
        cmd = svn_export('trunk')
        self.assertEqual(cmd[2], base_url.format('trunk'))

        cmd = svn_export('branches/v3')
        self.assertEqual(cmd[2], base_url.format('branches/v3'))

        cmd = svn_export('1.2.3')
        self.assertEqual(cmd[2], base_url.format('tags/1.2.3'))

        cmd = svn_export('1.2')
        self.assertEqual(cmd[2], base_url.format('tags/1.2'))

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_assert_svn_exists(self):
        """Test the check for svn presence.
        """
        mock = MagicMock(return_value=0)
        with patch('subprocess.check_output', mock):
            assert_svn_exists()
        mock.assert_called_with(['svn', '--version'])
        mock = MagicMock(side_effect=OSError(12345, "Mock error"))
        with patch('subprocess.check_output', mock):
            with self.assertRaises(AssertionError) as e:
                assert_svn_exists()
            self.assertEqual(str(e.exception), "svn command is not executable. Install svn to use the install script. Original Error is: 'Mock error'.")
        mock = MagicMock(side_effect=CalledProcessError(1, 'svn', 'Mock stdout', 'Mock stderr'))
        with patch('subprocess.check_output', mock):
            with self.assertRaises(AssertionError) as e:
                assert_svn_exists()
            self.assertEqual(str(e.exception), "The svn command (svn) on this system does not work. Output is: 'Mock stdout'.")

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_install(self):
        """Test the install function.
        """
        with patch.dict('os.environ'):
            try:
                del os.environ['DESIMODEL']
            except KeyError:
                pass
            with patch('os.path.exists') as exists:
                exists.return_value = True
                with self.assertRaises(ValueError) as e:
                    install(desimodel='/opt/desimodel')
                self.assertEqual(str(e.exception), "/opt/desimodel/data already exists!")
            exists.assert_called_with('/opt/desimodel/data')
            with patch('os.path.exists') as exists:
                exists.return_value = True
                with self.assertRaises(ValueError) as e:
                    install()
                self.assertEqual(str(e.exception), "{0}/data already exists!".format(default_install_dir()))
            exists.assert_called_with('{0}/data'.format(default_install_dir()))
            os.environ['DESIMODEL'] = '/opt/desimodel'
            with patch('os.path.exists') as exists:
                exists.return_value = True
                with self.assertRaises(ValueError) as e:
                    install()
                self.assertEqual(str(e.exception), "/opt/desimodel/data already exists!")
            exists.assert_called_with('/opt/desimodel/data')
        with patch.dict('os.environ'):
            try:
                del os.environ['DESIMODEL']
            except KeyError:
                pass
            with patch('desimodel.install.assert_svn_exists'):
                with patch('os.chdir'):
                    with patch('subprocess.Popen') as Popen:
                        proc = Popen.return_value
                        proc.communicate.return_value = ('Mock stdout', 'Mock stderr')
                        proc.returncode = 1
                        with self.assertRaises(RuntimeError) as e:
                            install(desimodel='/opt/desimodel')
                        self.assertEqual(str(e.exception), "Mock stderr")
                    Popen.assert_called_with(['svn', 'export', 'https://desi.lbl.gov/svn/code/desimodel/trunk/data'], stderr=-1, stdout=-1)
