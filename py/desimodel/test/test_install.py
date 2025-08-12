# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.install.
"""
import os
import shutil
from tempfile import mkdtemp
from subprocess import CalledProcessError
import unittest
from unittest.mock import patch, call
from ..install import default_install_dir, assert_svn_exists, svn_export, install, add_files_to_record
import desimodel


class TestInstall(unittest.TestCase):
    """Test desimodel.install.
    """

    def setUp(self):
        self.orig_version = desimodel.__version__
        self.tmp_dir = mkdtemp()

    def tearDown(self):
        desimodel.__version__ = self.orig_version
        shutil.rmtree(self.tmp_dir)

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
        cmd = svn_export(svn_checkout=True)
        self.assertEqual(cmd[1], 'checkout')
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

    @patch('desimodel.install.check_output')
    def test_assert_svn_exists(self, mock_output):
        """Test the check for svn presence.
        """
        mock_output.return_value = 0
        assert_svn_exists()
        mock_output.assert_called_with(['svn', '--version'])
        mock_output.side_effect = OSError(12345, "Mock error")
        with self.assertRaises(AssertionError) as e:
            assert_svn_exists()
        self.assertEqual(str(e.exception), "svn command is not executable. Install svn to use the install script. Original Error is: 'Mock error'.")
        mock_output.side_effect = CalledProcessError(1, 'svn', 'Mock stdout', 'Mock stderr')
        with self.assertRaises(AssertionError) as e:
            assert_svn_exists()
        self.assertEqual(str(e.exception), "The svn command (svn) on this system does not work. Output is: 'Mock stdout'.")

    @patch('builtins.print')
    @patch('desimodel.install.default_install_dir')
    def test_add_files_to_record(self, mock_dir, mock_print):
        """Test updating the RECORD metadata file.
        """
        version = '1.2.3'
        mock_dir.return_value = os.path.join(self.tmp_dir, 'desimodel')
        os.makedirs(os.path.join(self.tmp_dir, 'desimodel', 'data', 'weather'))
        added = add_files_to_record('desimodel', version)
        self.assertEqual(len(added), 0)
        desimodel_meta = os.path.join(self.tmp_dir, f'desimodel-{version}.dist-info')
        os.makedirs(desimodel_meta)
        record_file = os.path.join(desimodel_meta, 'RECORD')
        with open(record_file, 'w') as RECORD:
            RECORD.write('')
        with open(os.path.join(self.tmp_dir, 'desimodel', 'data', 'desi.yaml'), 'w') as YAML:
            YAML.write("foo: bar\n")
        with open(os.path.join(self.tmp_dir, 'desimodel', 'data', 'weather', 'weather.csv'), 'w') as CSV:
            CSV.write("day,temp\r\n20201212,15.0\r\n")
        self.assertTrue(os.path.exists(record_file))
        added = add_files_to_record('desimodel', version)
        self.assertEqual(len(added), 2)
        self.assertEqual(added[0], 'desimodel/data/desi.yaml,sha256=HavE48u9aggYvUYPOmyYVb_pXVBsdHJrwPLtsK7LH04,9')
        self.assertEqual(added[1], 'desimodel/data/weather/weather.csv,sha256=Jbi1Ugess7xrprvN2BkNvPQSRYcG2WRxylKXgaKxCA0,25')
        #
        # Dry-run mode
        #
        added = add_files_to_record('desimodel', version, data='data', dry_run=True)
        mock_print.assert_has_calls([call('desimodel/data/desi.yaml,sha256=HavE48u9aggYvUYPOmyYVb_pXVBsdHJrwPLtsK7LH04,9'),
                                     call('desimodel/data/weather/weather.csv,sha256=Jbi1Ugess7xrprvN2BkNvPQSRYcG2WRxylKXgaKxCA0,25')])

    @patch('os.path.exists')
    def test_install(self, mock_exists):
        """Test the install function.
        """
        with patch.dict('os.environ'):
            try:
                del os.environ['DESIMODEL']
            except KeyError:
                pass
            mock_exists.return_value = True
            with self.assertRaises(ValueError) as e:
                install(desimodel='/opt/desimodel')
            self.assertEqual(str(e.exception), "/opt/desimodel/data already exists!")
            mock_exists.assert_called_once_with('/opt/desimodel/data')
            mock_exists.reset_mock()
            with self.assertRaises(ValueError) as e:
                install()
            self.assertEqual(str(e.exception), "{0}/data already exists!".format(default_install_dir()))
            mock_exists.assert_called_once_with('{0}/data'.format(default_install_dir()))
            mock_exists.reset_mock()
            os.environ['DESIMODEL'] = '/opt/desimodel'
            with self.assertRaises(ValueError) as e:
                install()
            self.assertEqual(str(e.exception), "/opt/desimodel/data already exists!")
            mock_exists.assert_called_once_with('/opt/desimodel/data')
            mock_exists.reset_mock()

        mock_exists.return_value = False
        with patch.dict('os.environ'):
            try:
                del os.environ['DESIMODEL']
            except KeyError:
                pass
            with patch('desimodel.install.assert_svn_exists'):
                with patch('os.chdir'):
                    #
                    # Test with version explicitly set.
                    #
                    for v in ('0.19.3', '0.20.1.dev1234', 'branches/test-0.19', 'trunk'):
                        svn_command = svn_export(v)
                        with patch('desimodel.install.Popen') as Popen:
                            proc = Popen.return_value
                            proc.communicate.return_value = ('Mock stdout', 'Mock stderr')
                            proc.returncode = 1
                            with self.assertRaises(RuntimeError) as e:
                                install(desimodel='/opt/desimodel', version=v)
                            self.assertEqual(str(e.exception), "Mock stderr")
                        Popen.assert_called_once_with(svn_command, stderr=-1, stdout=-1)
                    #
                    # Test with version obtained from desimodel.__version__.
                    #
                    for v in ('0.19.3', '0.20.1.dev1234'):
                        desimodel.__version__ = v
                        svn_command = svn_export()
                        with patch('desimodel.install.Popen') as Popen:
                            proc = Popen.return_value
                            proc.communicate.return_value = ('Mock stdout', 'Mock stderr')
                            proc.returncode = 1
                            with self.assertRaises(RuntimeError) as e:
                                install(desimodel='/opt/desimodel')
                            self.assertEqual(str(e.exception), "Mock stderr")
                        Popen.assert_called_once_with(svn_command, stderr=-1, stdout=-1)
                    #
                    # Test dry_run mode.
                    #
                    with patch('builtins.print') as mock_print:
                        with patch('desimodel.install.default_install_dir') as mock_dir:
                            mock_dir.return_value = 'foo'
                            version='1.2.3'
                            desimodel_code = os.path.join(self.tmp_dir, 'desimodel')
                            os.makedirs(desimodel_code)
                            added = install(desimodel=desimodel_code, version=version, svn_checkout=False, dry_run=True)
                            mock_print.assert_has_calls([call(f'Installing desimodel data tags/{version} to {desimodel_code}'),
                                                         call('Dry run, would have run "svn export ' +
                                                              f'https://desi.lbl.gov/svn/code/desimodel/tags/{version}/data"')])
