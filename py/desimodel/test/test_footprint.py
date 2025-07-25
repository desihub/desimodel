# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.footprint.
"""
import unittest
import os
from collections import Counter
import numpy as np
from astropy.table import Table

from .. import io
from .. import footprint

desimodel_available = os.path.isdir(io.datadir())
desimodel_message = "The desimodel data set was not detected."

class TestFootprint(unittest.TestCase):
    """Test desimodel.footprint.
    """

    def setUp(self):
        io.reset_cache()

    def test_pass2program(self):
        '''Test footprint.pass2program().
        '''
        programs = footprint.pass2program(list(range(8)))
        count_layers = Counter(programs)
        self.assertEqual(count_layers['DARK'], 4)
        self.assertEqual(count_layers['GRAY'], 1)
        self.assertEqual(count_layers['BRIGHT'], 3)
        self.assertEqual(set(programs), set(['DARK', 'GRAY', 'BRIGHT']))

        #- confirm that it works with multiple kinds of input
        p0 = footprint.pass2program(0)
        p1 = footprint.pass2program(1)
        p01a = footprint.pass2program([0,1])
        p01b = footprint.pass2program(np.array([0,1]))
        self.assertEqual([p0,p1], p01a)
        self.assertEqual([p0,p1], p01b)

        p011 = footprint.pass2program([0,1,1])
        self.assertEqual([p0,p1,p1], p011)

        with self.assertRaises(KeyError):
            footprint.pass2program(999)


    def test_program2pass(self):
        '''Test footprint.program2pass().
        '''
        self.assertEqual(len(footprint.program2pass('DARK')), 7)
        # ADM in the real survey data there is no GRAY program...
#        self.assertEqual(len(footprint.program2pass('GRAY')), 1)
        # ADM ...but there is a BACKUP program.
        self.assertGreaterEqual(len(footprint.program2pass('BACKUP')), 1)
        self.assertEqual(len(footprint.program2pass('BRIGHT')), 5)

        passes = footprint.program2pass(['DARK', 'BACKUP', 'BRIGHT'])
        self.assertEqual(len(passes), 3)
        self.assertEqual(len(passes[0]), 7)
        self.assertGreaterEqual(len(passes[1]), 1)
        self.assertEqual(len(passes[2]), 5)

        with self.assertRaises(ValueError):
            footprint.program2pass('BLAT')

        #- confirm it works with column inputs too
        tiles = io.load_tiles()
        passes = footprint.program2pass(tiles['PROGRAM'])
        self.assertEqual(len(passes), len(tiles))
        for p in passes:
            self.assertNotEqual(p, None)
        passes = footprint.program2pass(Table(tiles)['PROGRAM'])
        self.assertEqual(len(passes), len(tiles))
        for p in passes:
            self.assertNotEqual(p, None)

    def test_ecsv(self):
        """Test consistency of ecsv vs. fits tiles files"""
        datadir = io.datadir()
        t1 = Table.read(f'{datadir}/footprint/desi-tiles.fits')
        t2 = Table.read(f'{datadir}/footprint/desi-tiles.ecsv',
            format='ascii.ecsv')
        for colname in t2.colnames:
            self.assertIn(colname, t1.colnames)
            try:
                assert(np.all(t1[colname] == t2[colname]))
            except AssertionError:
                self.assertTrue(np.allclose(t1[colname], t2[colname]), f'{colname} mismatch')

    def test_get_tile_radec(self):
        """Test grabbing tile information by tileID.
        """
        tx = io.load_tiles()
        tilefile = list(io._tiles.keys())[0]
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
        io._tiles[tilefile] = tiles

        #- TILEID 1 should be filtered out as not in DESI
        with self.assertRaises(ValueError):
            ra, dec = footprint.get_tile_radec(1)

        #- But TILEID 2 should be there with correct redshift
        ra, dec, = footprint.get_tile_radec(2)
        self.assertEqual((ra, dec), (1.0, -1.0))

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

        ret = footprint.is_point_in_desi(tiles, 0.0, -2.0, return_tile_index=True)
        self.assertEqual(ret, (True, 0))

        ret = footprint.is_point_in_desi(tiles, (0.0,), (-2.0,), radius=1.605, return_tile_index=True)
        self.assertEqual(ret, ([True], [0]))

        ret = footprint.is_point_in_desi(tiles, 0.0, -3.7, radius=1.605, return_tile_index=True)
        self.assertEqual(ret, (False, 0))

        ret = footprint.is_point_in_desi(tiles, -3.0, -2.0, radius=1.605, return_tile_index=True)
        self.assertEqual(ret, (False, 0))

        ret = footprint.is_point_in_desi(tiles, tiles['RA'], tiles['DEC'], radius=1.605)
        self.assertEqual(len(ret), len(tiles))

    def test_embed_sphere(self):
        rng = np.random.RandomState(1234)
        ra = rng.uniform(0, 360, size=1000)
        dec = rng.uniform(-90, 90, size=1000)

        ret = footprint._embed_sphere(ra, dec)
        self.assertEqual(ret.shape, (1000, 3))
        norm = np.einsum('ij, ij->i', ret, ret)
        self.assertTrue(all(abs(norm - 1.0) < 1e-5))

    def test_find_points_in_tiles(self):
        rng = np.random.RandomState(1234)

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

        ra = rng.uniform(-10, 10, 100000)
        dec = np.degrees(np.arcsin(rng.uniform(-0.1, 0.1, 100000)))
        lists = footprint.find_points_in_tiles(tiles, ra, dec, radius=1.6058)

        # assert we've found roughly same number of objects per tile
        counts = np.array([len(i) for i in lists])
        self.assertLess(counts.std() / counts.mean(), 1e-2)

        # assert each list are indeed in the cell.
        for i, ii in enumerate(lists):
            xyz = footprint._embed_sphere(ra[ii], dec[ii])
            xyzc = footprint._embed_sphere(tiles['RA'][i], tiles['DEC'][i])
            diff = xyz - xyzc
            dist = np.einsum('ij, ij->i', diff, diff) ** 0.5
            self.assertLess(dist.max(), 2 * np.sin(np.radians(1.6058) * 0.5))

        # tiles overlapped, so we must have duplicates
        full = np.concatenate(lists)
        self.assertLess(len(np.unique(full)), len(full))

        list1 = footprint.find_points_in_tiles(tiles[0], ra, dec, radius=1.6058)
        self.assertEqual(sorted(list1), sorted(lists[0]))

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

        ret = footprint.find_tiles_over_point(tiles, 0.0, -2.0)
        self.assertEqual(ret, [0, 1])

        ret = footprint.find_tiles_over_point(tiles, 1.0, -1.0, radius=1.605)
        self.assertEqual(ret, [0, 1])

        # outside
        ret = footprint.find_tiles_over_point(tiles, 0.0, -3.7, radius=1.605)
        self.assertEqual(ret, [])

        # array input
        ret = footprint.find_tiles_over_point(tiles, (0.0,), (-3.7,), radius=1.605)
        self.assertEqual(len(ret), 1)
        self.assertEqual(ret[0], [])

    def test_find_points_radec(self):
        """Checks if the function is successfully finding points within a certain radius of a telra and teldec.
        """
        answer = footprint.find_points_radec(0, 0,
                                             np.array([1.5, 0, 1.9, 0]),
                                             np.array([0, 1.5, 0, 1.3]))
        self.assertEqual(answer, [0, 1, 3])

    def test_tiles2pix(self):
        tiles = Table()
        tiles['RA'] = [1.0, 2.0]
        tiles['DEC'] = [10.0, 22.0]
        pix = footprint.tiles2pix(16, tiles=tiles, per_tile=True)
        self.assertEqual(len(pix), len(tiles))
        allpix = np.unique(np.concatenate(pix))
        pix = footprint.tiles2pix(16, tiles=tiles, per_tile=False)
        self.assertEqual(set(pix), set(allpix))

        #- also works with TILERA, TILEDEC
        tiles.rename_column('RA', 'TILERA')
        tiles.rename_column('DEC', 'TILEDEC')
        pix2 = footprint.tiles2pix(16, tiles=tiles, per_tile=False)
        self.assertTrue(np.all(pix == pix2))

        #- also works with dict and numpy structured array
        pix3 = footprint.tiles2pix(16, tiles=np.array(tiles), per_tile=False)
        self.assertTrue(np.all(pix == pix3))
        tdict = dict(RA=list(tiles['TILERA']), DEC=list(tiles['TILEDEC']))
        pix4 = footprint.tiles2pix(16, tiles=tdict, per_tile=False)
        self.assertTrue(np.all(pix == pix4))

    def test_tileids2pix(self):
        tiles = io.load_tiles()
        pix = footprint.tileids2pix(16, tiles['TILEID'][0:3])
        self.assertGreater(len(pix), 0)
        n = np.max(tiles['TILEID'])
        with self.assertRaises(ValueError):
            pix = footprint.tileids2pix(16, [n, n+1, n+2])

    def test_partial_pixels(self):
        """Check weights assigned to HEALPixels that partially overlap tiles.
        """
        tiles = np.zeros((4,), dtype=[('TILEID', 'i2'),
                                      ('RA', 'f8'),
                                      ('DEC', 'f8'),
                                      ('IN_DESI', 'i2'),
                                      ('PROGRAM', (str, 6)),
                                  ])

        #ADM I found a full (170) partial (406) and empty (1000) HEALPixel at nside=256
        #ADM that are also full, partial, empty at nside=64.
        #ADM You can find these with, e.g.:
        #ADM pixweight256 = footprint.pixweight(256,tiles=tiles,radius=radius,precision=0.01)
        #ADM np.where(pixweight256 == 1)
        fullpix256 = np.array([170])
        partpix256 = np.array([406])
        emptypix256 = np.array([1000])
        #ADM In the nested scheme you can traverse from 256 to 64 by integer division
        #ADM by (256/64)**2 = 16
        fullpix64 = fullpix256//16
        partpix64 = partpix256//16
        emptypix64 = emptypix256//16

        #ADM I then used:
        #ADM full = footprint.pix2tiles(256,[180])
        #ADM part = footprint.pix2tiles(256,[406])
        #ADM empty = footprint.pix2tiles(256,[20])
        #ADM to determine the tile coordinates for these pixels
        #ADM (the first tile in  the array is full and the last tile is empty)
        tiles['TILEID'] = np.arange(4) + 1
        tiles['RA'] = [43.05, 47.41, 47.08, 45.38]
        tiles['DEC'] = [1.54, 3.12, 3.17, 0.0]
        tiles['IN_DESI'] = [1, 1, 1, 1]
        tiles['PROGRAM'] = 'DARK'

        #ADM The approximate radius of DESI tiles
        radius = 1.605

        #ADM determine the weight array for pixels at nsides of 64 and 256
        pixweight64 = footprint.pixweight(64,tiles=tiles,radius=radius,precision=0.01)
        pixweight256 = footprint.pixweight(256,tiles=tiles,radius=radius,precision=0.01)

        #ADM check that the full pixel is assigned a weight of 1 at each nside
        self.assertTrue(np.all(pixweight64[fullpix64]==1))
        self.assertTrue(np.all(pixweight256[fullpix256]==1))

        #ADM check that the empty pixel is assigned a weight of 0 at each nside
        self.assertTrue(np.all(pixweight64[emptypix64]==0))
        self.assertTrue(np.all(pixweight256[emptypix256]==0))

        #ADM check weights of partial pixels agree reasonably at different nsides
        hirespixels = partpix64*16+np.arange(16)
        hiresweight = np.mean(pixweight256[hirespixels])
        loresweight = pixweight64[partpix64]
        #ADM really they should agree to much better than 2%. As "precision" is not set to be
        #ADM very high, this is just to check for catastrophic differences
        #ADM I checked that at precision = 0.01 this doesn't fail after 1000 attempts
        #ADM (the largest difference I encountered was 0.015)
        self.assertTrue(np.all(np.abs(hiresweight-loresweight) < 0.02))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_spatial_real_tiles(self):
        """Test code on actual DESI tiles.
        """
        tiles = io.load_tiles()
        rng = np.random.RandomState(1234)
        ra = rng.uniform(0, 360, 1000)
        dec = rng.uniform(-90, 90, 1000)

        # some points must be in the sky,
        # https://github.com/desihub/desimodel/pull/37#issuecomment-270435938
        indesi = footprint.is_point_in_desi(tiles, ra, dec, radius=1.605)
        self.assertTrue(np.any(indesi))

        # now assert the consistency between find_tiles_over_point and is_point_in_desi
        ret = footprint.find_tiles_over_point(tiles, ra, dec, radius=1.605)
        self.assertEqual(len(ret), len(ra))
        indesi2 = [len(i) > 0 for i in ret]

        # FIXME: is there a better assertion function for arrays?
        self.assertEqual(list(indesi), indesi2)

        # Just interesting to see how many tiles overlap a random point?
        ### print(np.bincount([len(i) for i in ret]))
