# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.io.
"""
import os
import numpy as np
from astropy.table import Table
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

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        """Ensure that any desimodel.io caches are clear before running
        any test.
        """
        io._thru = dict()
        io._psf = dict()
        io._params = None
        io._fiberpos = None
        io._tiles = None

    def tearDown(self):
        pass

    @unittest.skipUnless(specter_available, specter_message)
    def test_load_throughput(self):
        """Test loading of throughput files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_throughput(channel)

    @unittest.skipUnless(specter_available, specter_message)
    def test_load_psf(self):
        """Test loading of PSF files.
        """
        for channel in ('b', 'r', 'z'):
            t = io.load_psf(channel)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_desiparams(self):
        """Test loading of basic DESI parameters.
        """
        p = io.load_desiparams()
        self.assertTrue(isinstance(p, dict))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_fiberpos(self):
        """Test loading of fiber positioner data.
        """
        fiberpos = io.load_fiberpos()
        self.assertEqual(len(fiberpos), 5000)
        # Check workaround for astropy 1.0.x bug for lower -> upper col names.
        for key in ('FIBER', 'POSITIONER', 'SPECTROGRAPH', 'X', 'Y', 'Z'):
            self.assertIn(key, fiberpos.dtype.names)
            x = fiberpos[key]

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_tiles(self):
        """Test loading of tile files.
        """
        self.assertIsNone(io._tiles)
        t1 = io.load_tiles(onlydesi=False)
        tile_cache_id1 = id(io._tiles)
        self.assertIsNotNone(io._tiles)
        t2 = io.load_tiles(onlydesi=True)
        tile_cache_id2 = id(io._tiles)
        self.assertEqual(tile_cache_id1, tile_cache_id2)
        self.assertIs(t1['OBSCONDITIONS'].dtype, np.dtype(np.uint16))
        self.assertIs(t2['OBSCONDITIONS'].dtype, np.dtype(np.uint16))
        self.assertLess(len(t2), len(t1))
        # All tiles in DESI are also in full set.
        self.assertTrue(np.all(np.in1d(t2['TILEID'], t1['TILEID'])))
        # I think this is the exact same test as above, except using set theory.
        self.assertEqual(len(set(t2.TILEID) - set(t1.TILEID)), 0)
        t3 = io.load_tiles(onlydesi=False)
        tile_cache_id3 = id(io._tiles)
        self.assertEqual(tile_cache_id1, tile_cache_id3)
        self.assertIs(t3['OBSCONDITIONS'].dtype, np.dtype(np.uint16))
        # Check for extra tiles.
        a = io.load_tiles(extra=False)
        self.assertEqual(np.sum(np.char.startswith(a['PROGRAM'], 'EXTRA')), 0)
        b = io.load_tiles(extra=True)
        self.assertGreater(np.sum(np.char.startswith(b['PROGRAM'], 'EXTRA')), 0)
        self.assertLess(len(a), len(b))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_tiles_consistency(self):
        """Test consistency of tile files.

        - Validate different tile loading schemes.
        - Test numerical consistency between FITS & ECSV versions.
        - Assure that PROGRAM column has no trailing whitespace.
        """
        fitstiles = io.findfile('footprint/desi-tiles.fits')
        ecsvtiles = io.findfile('footprint/desi-tiles.ecsv')
        # tf=TilesFits  tt=TilesTable  te=TilesEcsv
        tf = io.load_tiles(onlydesi=False, extra=True)
        tt = Table.read(fitstiles)
        te = Table.read(ecsvtiles, format='ascii.ecsv')
        self.assertEqual(sorted(tf.dtype.names), sorted(tt.colnames))
        self.assertEqual(sorted(tf.dtype.names), sorted(te.colnames))

        for col in tt.colnames:
            self.assertTrue(np.all(tf[col]==tt[col]), 'fits[{col}] != table[{col}]'.format(col=col))
            if np.issubdtype(tf[col].dtype, float):
                self.assertTrue(np.allclose(tf[col], te[col], atol=1e-4, rtol=1e-4), 'fits[{col}] != ecsv[{col}]'.format(col=col))
            else:
                self.assertTrue(np.all(tf[col]==te[col]), 'fits[{col}] != ecsv[{col}]'.format(col=col))

        self.assertTrue(not np.any(np.char.endswith(tf['PROGRAM'], ' ')))
        self.assertTrue(not np.any(np.char.endswith(tt['PROGRAM'], ' ')))
        self.assertTrue(not np.any(np.char.endswith(te['PROGRAM'], ' ')))

    def test_get_tile_radec(self):
        """Test grabbing tile information by tileID.
        """
        io_tile_cache = io._tiles
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'
        io._tiles = tiles
        ra, dec = io.get_tile_radec(1)
        self.assertEqual((ra, dec), (0.0, 0.0))
        ra, dec, = io.get_tile_radec(2)
        self.assertEqual((ra, dec), (1.0, -1.0))
        io._tiles = io_tile_cache

    def test_is_point_in_desi_mock(self):
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'

        ret = io.is_point_in_desi(tiles, 0.0, -2.0, return_tile_index=True)
        self.assertEqual(ret, (True, 0))

        ret = io.is_point_in_desi(tiles, (0.0,), (-2.0,), return_tile_index=True)
        self.assertEqual(ret, ([True], [0]))

        ret = io.is_point_in_desi(tiles, 0.0, -3.7, return_tile_index=True)
        self.assertEqual(ret, (False, 0))

        ret = io.is_point_in_desi(tiles, -3.0, -2.0, return_tile_index=True)
        self.assertEqual(ret, (False, 0))

        ret = io.is_point_in_desi(tiles, tiles['RA'], tiles['DEC'])
        self.assertEqual(len(ret), len(tiles))

    def test_find_tiles_over_point(self):
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [0.0, 1.0, 2.0, 3.0]
        tiles['DEC'] = [-2.0, -1.0, 1.0, 2.0]
        tiles['IN_DESI'] = [0, 1, 1, 0]
        tiles['PROGRAM'] = 'DARK'

        ret = io.find_tiles_over_point(tiles, 0.0, -2.0)
        self.assertEqual(ret, [0, 1])

        ret = io.find_tiles_over_point(tiles, 1.0, -1.0)
        self.assertEqual(ret, [0, 1])

        # outside
        ret = io.find_tiles_over_point(tiles, 0.0, -3.7)
        self.assertEqual(ret, [])

        # array input
        ret = io.find_tiles_over_point(tiles, (0.0,), (-3.7,))
        self.assertEqual(len(ret), 1)
        self.assertEqual(ret[0], [])

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_spatial_real_tiles(self):
        tiles = io.load_tiles()
        rng = np.random.RandomState(1234)
        ra = rng.uniform(0, 360, 1000)
        dec = rng.uniform(-90, 90, 1000)

        # some points must be in the sky,
        # https://github.com/desihub/desimodel/pull/37#issuecomment-270435938
        indesi = io.is_point_in_desi(tiles, ra, dec)
        self.assertTrue(np.any(indesi))

        # now assert the consistency between find_tiles_over_point and is_point_in_desi
        ret = io.find_tiles_over_point(tiles, ra, dec)
        self.assertEqual(len(ret), len(ra))
        indesi2 = [len(i) > 0 for i in ret]

        # FIXME: is there a better assertion function for arrays?
        self.assertEqual(list(indesi), indesi2)

        # Just interesting to see how many tiles overlap a random point?
        print(np.bincount([len(i) for i in ret]))

def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
