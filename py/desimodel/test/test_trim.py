# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.trim.
"""
import unittest
from os.path import abspath, dirname
import numpy as np
from ..trim import (inout, rebin_image, trim_focalplane, trim_footprint, trim_inputs,
                    trim_sky, trim_targets, trim_data, trim_specpsf, trim_spectra,
                    trim_throughput)

skipMock = False
try:
    from unittest.mock import call, patch, mock_open, MagicMock, DEFAULT
except ImportError:
    # Python 2
    skipMock = True


class DummyData(object):
    """Simple object with a 'data' attribute.
    """
    def __init__(self, d):
        self.data = d


def iterable_mock_open(*args, **kargs):
    """unittest.mock.mock_open doesn't support line iteration.
    """
    f_open = mock_open(*args, **kargs)
    f_open.return_value.__iter__ = lambda self : iter(self.readline, '')
    return f_open


class TestTrim(unittest.TestCase):
    """Test desimodel.trim.
    """

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_data(self):
        """Test the trim wrapper function.
        """
        with patch.multiple('desimodel.trim', trim_focalplane=DEFAULT,
                            trim_footprint=DEFAULT, trim_inputs=DEFAULT,
                            trim_sky=DEFAULT, trim_specpsf=DEFAULT, trim_spectra=DEFAULT,
                            trim_targets=DEFAULT, trim_throughput=DEFAULT) as desimodel_trim:
            with patch('os.path.exists') as exists:
                exists.return_value = False
                with patch('os.makedirs') as makedirs:
                    with patch('shutil.copy') as copy:
                        trim_data('/in', '/out')
            exists.assert_called_with('/out')
            makedirs.assert_called_with('/out')
            copy.assert_called_with('/in/desi.yaml', '/out/desi.yaml')
            desimodel_trim['trim_focalplane'].assert_called_with('/in/focalplane', '/out/focalplane')
            desimodel_trim['trim_footprint'].assert_called_with('/in/footprint', '/out/footprint')
            desimodel_trim['trim_inputs'].assert_called_with('/in/inputs', '/out/inputs')
            desimodel_trim['trim_sky'].assert_called_with('/in/sky', '/out/sky')
            desimodel_trim['trim_specpsf'].assert_called_with('/in/specpsf', '/out/specpsf')
            desimodel_trim['trim_spectra'].assert_called_with('/in/spectra', '/out/spectra')
            desimodel_trim['trim_targets'].assert_called_with('/in/targets', '/out/targets')
            desimodel_trim['trim_throughput'].assert_called_with('/in/throughput', '/out/throughput')
            with patch('os.path.exists') as exists:
                exists.return_value = True
                with patch('os.makedirs') as makedirs:
                    with patch.multiple('shutil', copy=DEFAULT, rmtree=DEFAULT) as shutilmock:
                        trim_data('/in', '/out', overwrite=True)
            exists.assert_called_with('/out')
            makedirs.assert_called_with('/out')
            shutilmock['rmtree'].assert_called_with('/out')
            shutilmock['copy'].assert_called_with('/in/desi.yaml', '/out/desi.yaml')

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

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_footprint(self):
        """Test trim_footprint().
        """
        data = np.zeros((100,), dtype=[('RA', np.float64), ('DEC', np.float64)])
        data['RA'] = np.linspace(0, 360, 100, dtype=np.float64)
        data['DEC'] = np.linspace(-90, 90, 100, dtype=np.float64)
        with patch('os.path.exists') as exists:
            exists.return_value = False
            with patch('os.makedirs') as makedirs:
                with patch('astropy.io.fits.open') as hdulist:
                    hdulist.return_value.__enter__.return_value = [None, DummyData(data)]
                    with patch('astropy.table.Table.write') as w:
                        with patch('desimodel.trim.pixweight') as pixweight:
                            trim_footprint('/in/footprint', '/out/footprint')
        exists.assert_called_with('/out/footprint')
        makedirs.assert_called_with('/out/footprint')
        hdulist.assert_called_with('/in/footprint/desi-tiles.fits')
        w.assert_has_calls([call('/out/footprint/desi-tiles.fits', format='fits'),
                            call('/out/footprint/desi-tiles.ecsv', format='ascii.ecsv')])

    def test_trim_inputs(self):
        """Test trim_inputs().
        """
        # This function does nothing, so we just check that we can call it.
        trim_inputs('/in/focalplane', '/out/focalplane')

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_sky(self):
        """Test trim_sky().
        """
        with patch('os.makedirs'):
            with patch('shutil.copy') as copy:
                trim_sky('/in/sky', '/out/sky')
                copy.assert_called_with('/in/sky/solarspec.txt',
                                        '/out/sky/solarspec.txt')

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_specpsf(self):
        with patch.multiple('desimodel.trim', trim_quickpsf=DEFAULT, trim_psf=DEFAULT) as desimodel_trim:
            with patch('os.path.exists') as exists:
                exists.return_value = False
                with patch('os.makedirs') as makedirs:
                    trim_specpsf('/in/specpsf', '/out/specpsf')
            exists.assert_called_with('/out/specpsf')
            makedirs.assert_called_with('/out/specpsf')
            desimodel_trim['trim_quickpsf'].assert_called_with('/in/specpsf', '/out/specpsf', 'psf-quicksim.fits')
            desimodel_trim['trim_psf'].assert_has_calls([call('/in/specpsf', '/out/specpsf', 'psf-b.fits'),
                                                         call('/in/specpsf', '/out/specpsf', 'psf-r.fits'),
                                                         call('/in/specpsf', '/out/specpsf', 'psf-z.fits')])

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_spectra(self):
        """Test trim_spectra().
        """
        spectrum = '\n'.join(['#foo'] + [str(x) for x in range(42)])+'\n'
        with patch('os.path.exists') as exists:
            exists.return_value = False
            with patch('os.makedirs') as makedirs:
                with patch('builtins.open', new_callable=iterable_mock_open, read_data=spectrum) as mo:
                    trim_spectra('/in/spectra', '/out/spectra')
        exists.assert_called_with('/out/spectra')
        makedirs.assert_called_with('/out/spectra')
        # print(mo.mock_calls)
        mo.assert_has_calls([call('/in/spectra/spec-sky.dat'), call('/out/spectra/ZenithExtinction-KPNO.dat', 'w')], any_order=True)
        handle = mo()
        # print(handle.mock_calls)
        handle.write.assert_has_calls([call('#foo\n'), call('19\n'), call('39\n')], any_order=True)
        # print(handle.write.mock_calls)

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_targets(self):
        """Test trim_targets().
        """
        with patch('shutil.copytree') as copytree:
            trim_targets('/in/targets', '/out/targets')
            copytree.assert_called_with('/in/targets', '/out/targets')

    @unittest.skipIf(skipMock, "Skipping test that requires unittest.mock.")
    def test_trim_throughput(self):
        """Test trim_throughput().
        """
        with patch('os.path.exists') as exists:
            exists.return_value = False
            with patch('os.makedirs') as makedirs:
                with patch('shutil.copy') as copy:
                    with patch.multiple('astropy.io.fits', open=DEFAULT, HDUList=DEFAULT, BinTableHDU=DEFAULT) as fits:
                        trim_throughput('/in/throughput', '/out/throughput')
        exists.assert_called_with('/out/throughput')
        makedirs.assert_called_with('/out/throughput')
        fits['open'].assert_has_calls([call('/in/throughput/thru-b.fits'),
                                       call('/in/throughput/thru-r.fits'),
                                       call('/in/throughput/thru-z.fits')],
                                       any_order=True)
        fits['HDUList']().writeto.assert_has_calls([call('/out/throughput/thru-b.fits'),
                                                    call('/out/throughput/thru-r.fits'),
                                                    call('/out/throughput/thru-z.fits')])
        copy.assert_has_calls([call('/in/throughput/fiberloss-elg.dat', '/out/throughput/fiberloss-elg.dat'),
                               call('/in/throughput/fiberloss-lrg.dat', '/out/throughput/fiberloss-lrg.dat'),
                               call('/in/throughput/fiberloss-perfect.dat', '/out/throughput/fiberloss-perfect.dat'),
                               call('/in/throughput/fiberloss-qso.dat', '/out/throughput/fiberloss-qso.dat'),
                               call('/in/throughput/fiberloss-sky.dat', '/out/throughput/fiberloss-sky.dat'),
                               call('/in/throughput/fiberloss-star.dat', '/out/throughput/fiberloss-star.dat'),
                               call('/in/throughput/DESI-0347_blur.ecsv', '/out/throughput/DESI-0347_blur.ecsv'),
                               call('/in/throughput/DESI-0347_offset.ecsv', '/out/throughput/DESI-0347_offset.ecsv'),
                               call('/in/throughput/DESI-0347_random_offset_1.fits', '/out/throughput/DESI-0347_random_offset_1.fits'),
                               call('/in/throughput/galsim-fiber-acceptance.fits', '/out/throughput/galsim-fiber-acceptance.fits')])


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
