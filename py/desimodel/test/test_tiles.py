import unittest
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table

import desimodel.io

class TestTiles(unittest.TestCase):
    '''
    Test desi-tiles.fits
    '''

    def setUp(self):
        self.fitstiles = desimodel.io.findfile('footprint/desi-tiles.fits')
        self.ecsvtiles = desimodel.io.findfile('footprint/desi-tiles.ecsv')

        #- tf=TilesFits  tt=TilesTable  te=TilesEcsv
        self.tf = fits.getdata(self.fitstiles)
        self.tt = Table.read(self.fitstiles)
        self.te = Table.read(self.ecsvtiles, format='ascii.ecsv')

    def test_options(self):
        a = desimodel.io.load_tiles(onlydesi=True)
        self.assertTrue(np.all(a['IN_DESI'] > 0))
        b = desimodel.io.load_tiles(onlydesi=False)
        self.assertTrue(np.any(b['IN_DESI'] == 0))
        self.assertLess(len(a), len(b))
        #- All tiles in DESI are also in full set
        self.assertTrue(np.all(np.in1d(a['TILEID'], b['TILEID'])))

        a = desimodel.io.load_tiles(extra=False)
        self.assertEqual(np.sum(np.char.startswith(a['PROGRAM'], 'EXTRA')), 0)
        b = desimodel.io.load_tiles(extra=True)
        self.assertGreater(np.sum(np.char.startswith(b['PROGRAM'], 'EXTRA')), 0)
        self.assertLess(len(a), len(b))        

    #- Test consistency of FITS vs. ECSV ASCII formats, allowing for rounding
    #- differences between decimal ASCII and IEEE float32
    def test_consistency(self):
        self.assertEqual(sorted(self.tf.dtype.names), sorted(self.tt.colnames))
        self.assertEqual(sorted(self.tf.dtype.names), sorted(self.te.colnames))

        for col in self.tt.colnames:
            self.assertTrue(np.all(self.tf[col]==self.tt[col]), 'fits[{col}] != table[{col}]'.format(col=col))
            if np.issubdtype(self.tf[col].dtype, float):
                self.assertTrue(np.allclose(self.tf[col], self.te[col], atol=1e-4, rtol=1e-4), 'fits[{col}] != ecsv[{col}]'.format(col=col))
            else:
                self.assertTrue(np.all(self.tf[col]==self.te[col]), 'fits[{col}] != ecsv[{col}]'.format(col=col))

    #- Test that PROGRAM does not have trailing white space in input files
    def test_program(self):
        self.assertTrue(not np.any(np.char.endswith(self.tf['PROGRAM'], ' ')))
        self.assertTrue(not np.any(np.char.endswith(self.tt['PROGRAM'], ' ')))
        self.assertTrue(not np.any(np.char.endswith(self.te['PROGRAM'], ' ')))

    def test_get_tile_radec(self):
        """Test grabbing tile information by tileID.
        """
        io_tile_cache = desimodel.io._tiles
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
        desimodel.io._tiles = tiles
        ra, dec = desimodel.io.get_tile_radec(1)
        self.assertEqual((ra, dec), (0.0, 0.0))
        ra, dec, = desimodel.io.get_tile_radec(2)
        self.assertEqual((ra, dec), (1.0, -1.0))
        desimodel.io._tiles = io_tile_cache

if __name__ == '__main__':
    unittest.main()
