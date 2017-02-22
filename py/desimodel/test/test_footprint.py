import unittest
import os

import numpy as np

from .. import io
from .. import footprint

desimodel_available = True
desimodel_message = "The desimodel data set was not detected."
try:
    spam = os.environ['DESIMODEL']
except KeyError:
    desimodel_available = False

class TestFootprint(unittest.TestCase):
    
    def setUp(self):
        pass
            
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
        ra, dec = footprint.get_tile_radec(1)
        self.assertEqual((ra, dec), (0.0, 0.0))
        ra, dec, = footprint.get_tile_radec(2)
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
        lists = footprint.find_points_in_tiles(tiles, ra, dec)

        # assert we've found roughly same number of objects per tile
        counts = np.array([len(i) for i in lists])
        self.assertLess(counts.std() / counts.mean(), 1e-2)

        # assert each list are indeed in the cell.
        for i, ii in enumerate(lists):
            xyz = footprint._embed_sphere(ra[ii], dec[ii])
            xyzc = footprint._embed_sphere(tiles['RA'][i], tiles['DEC'][i])
            diff = xyz - xyzc
            dist = np.einsum('ij, ij->i', diff, diff) ** 0.5
            self.assertLess(dist.max(), 2 * np.sin(np.radians(1.605) * 0.5))

        # tiles overlapped, so we must have duplicates
        full = np.concatenate(lists)
        self.assertLess(len(np.unique(full)), len(full))

        list1 = footprint.find_points_in_tiles(tiles[0], ra, dec, radius=1.605)
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

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_spatial_real_tiles(self):
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
        print(np.bincount([len(i) for i in ret]))
                
def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
