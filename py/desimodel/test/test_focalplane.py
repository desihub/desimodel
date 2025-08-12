# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.focalplane.
"""
import unittest
import numpy as np
from ..focalplane import *
from .. import io
from astropy.table import Table

class TestFocalplane(unittest.TestCase):
    """Test desimodel.focalplane.
    """

    def test_random_offsets(self):
        """Test generating random positioner offsets.
        """
        dx, dy = generate_random_centroid_offsets(1.0)
        self.assertEqual(dx.shape, dy.shape)
        dr = np.sqrt(dx**2 + dy**2)
        self.assertAlmostEqual(np.sqrt(np.average(dr**2)), 1.0)
        self.assertLess(np.max(dr), 7)

    def test_check_radec(self):
        """Test the RA, Dec bounds checking.
        """
        F = FocalPlane()
        with self.assertRaises(ValueError):
            F._check_radec(365.0, 0.0)
        with self.assertRaises(ValueError):
            F._check_radec(0.0, 100.0)

    def test_set_tele_pointing(self):
        """Test setting the RA, Dec by hand.
        """
        F = FocalPlane()
        self.assertEqual(F.ra, 0.0)
        self.assertEqual(F.dec, 0.0)
        F.set_tele_pointing(180.0, 45.0)
        self.assertEqual(F.ra, 180.0)
        self.assertEqual(F.dec, 45.0)

    def test_get_radius(self):
        """Tests converting x, y coordinates on the focal plane to a radius in degrees and mm
        """
        ps = io.load_platescale()
        for i in np.linspace(0, len(ps)-1, 10).astype(int):
            degree = get_radius_deg(0, ps['radius'][i])
            self.assertAlmostEqual(degree, ps['theta'][i])
            radius = get_radius_mm(ps['theta'][i])
            self.assertAlmostEqual(radius, ps['radius'][i])

        #- Rotational invariance
        d1 = get_radius_deg(0, ps['radius'][10])
        d2 = get_radius_deg(ps['radius'][10], 0)
        d3 = get_radius_deg(0, -ps['radius'][10])
        d4 = get_radius_deg(-ps['radius'][10], 0)
        self.assertAlmostEqual(d1, d2)
        self.assertAlmostEqual(d1, d3)
        self.assertAlmostEqual(d1, d4)

        #- lists and arrays should also work
        degrees = get_radius_deg([0,0,0,0], ps['radius'][0:4])
        self.assertTrue(np.allclose(degrees, ps['theta'][0:4]))

        radius = get_radius_mm(ps['theta'][100:110])
        self.assertTrue(np.allclose(radius, ps['radius'][100:110]))

    def test_FocalPlane(self):
        n = 20
        x, y = np.random.uniform(-100, 100, size=(2,n))
        f = FocalPlane(ra=100, dec=20)
        ra, dec = f.xy2radec(x, y)
        xx, yy = f.radec2xy(ra, dec)
        self.assertLess(np.max(np.abs(xx-x)), 1e-4)
        self.assertLess(np.max(np.abs(yy-y)), 1e-4)

    def test_xy2qs(self):
        x, y = qs2xy(0.0, 100)
        self.assertAlmostEqual(y, 0.0)
        x, y = qs2xy(90.0, 100)
        self.assertAlmostEqual(x, 0.0)
        x, y = qs2xy(180.0, 100)
        self.assertAlmostEqual(y, 0.0)
        x, y = qs2xy(270.0, 100)
        self.assertAlmostEqual(x, 0.0)

        q, s = xy2qs(0.0, 100.0)
        self.assertAlmostEqual(q, 90.0)
        q, s = xy2qs(100.0, 0.0)
        self.assertAlmostEqual(q, 0.0)
        q, s = xy2qs(-100.0, 100.0)
        self.assertAlmostEqual(q, 135.0)

        n = 100
        s = 410*np.sqrt(np.random.uniform(0, 1, size=n))
        q = np.random.uniform(0, 360, size=n)

        x, y = qs2xy(q, s)
        qq, ss = xy2qs(x, y)
        xx, yy = qs2xy(qq, ss)
        self.assertLess(np.max(np.abs(xx-x)), 1e-4)
        self.assertLess(np.max(np.abs(yy-y)), 1e-4)
        self.assertLess(np.max(np.abs(ss-s)), 1e-4)
        dq = (qq-q) % 360
        dq[dq>180] -= 360
        self.assertLess(np.max(np.abs(dq)), 1e-4)

    def test_xy2radec_roundtrip(self):
        """Tests the consistency between the conversion functions
        radec2xy and xy2radec. Also tests the accuracy of the xy2radec
        on particular cases.
        """
        n = 100
        r = 410 * np.sqrt(np.random.uniform(0, 1, size=n))
        theta = np.random.uniform(0, 2*np.pi, size=n)
        x = r*np.cos(theta)
        y = r*np.sin(theta)

        test_telra = [0.0, 0.0, 90.0, 30.0]
        test_teldec = [0.0, 90.0, 0.0, -30.0]
        test_telra = np.concatenate([test_telra, np.random.uniform(0, 360, size=5)])
        test_teldec = np.concatenate([test_teldec, np.random.uniform(-90, 90, size=5)])

        for telra, teldec in zip(test_telra, test_teldec):
            ra, dec = xy2radec(telra, teldec, x, y)
            xx, yy = radec2xy(telra, teldec, ra, dec)
            dx = xx-x
            dy = yy-y
            self.assertLess(np.max(np.abs(dx)), 1e-4)   #- 0.1 um
            self.assertLess(np.max(np.abs(dy)), 1e-4)   #- 0.1 um

    def test_xy2radec_orientation(self):
        """Test that +x = -RA and +y = +dec"""
        n = 5
        x = np.linspace(-400, 400, n)
        y = np.zeros_like(x)
        ii = np.arange(n, dtype=int)
        for teldec in [-30, 0, 30]:
            ra, dec = xy2radec(0.0, teldec, x, y)
            self.assertTrue(np.all(np.argsort((ra+180) % 360) == ii[-1::-1]))
            ra, dec = xy2radec(359.9, teldec, x, y)
            self.assertTrue(np.all(np.argsort((ra+180) % 360) == ii[-1::-1]))
            ra, dec = xy2radec(30.0, teldec, x, y)
            self.assertTrue(np.all(np.argsort(ra) == ii[-1::-1]))

        y = np.linspace(-400, 400, n)
        x = np.zeros_like(y)
        for telra, teldec in [(0, -30), (0, 0), (0, 80), (30, 30)]:
            ra, dec = xy2radec(0.0, teldec, x, y)
            self.assertTrue(np.all(np.argsort(dec) == ii))

    def test_area_arcsec2(self):
        """ Test area calculation
        """
        x = [0., 50., 100., 200]
        y = [-0., -50., -100., -200]

        #- runs with scalar, list, or array
        a1 = fiber_area_arcsec2(10, 20)
        a2 = fiber_area_arcsec2([10, 30], [20, 40])
        a3 = fiber_area_arcsec2(np.array([10, 30]), np.array([20, 40]))
        self.assertTrue(np.allclose(a2, a3))
        self.assertTrue(np.allclose(a1, a2[0]))
        self.assertTrue(np.allclose(a1, a3[0]))

        #- approximately correct values (depends on platescale data)
        a0, a400 = fiber_area_arcsec2([0, 400], [0, 0])
        self.assertAlmostEqual(a0, 1.9744, 2)
        self.assertAlmostEqual(a400, 1.6754, 2)

        #- area is constant vs. radius
        a = fiber_area_arcsec2([30,40,0,-50], [40, 30, 50, 0])
        self.assertTrue(np.allclose(a[0], a))

        #- area decreases with radius
        x = np.arange(0,401,10)
        y = np.zeros(len(x))
        da = np.diff(fiber_area_arcsec2(x, y))
        self.assertTrue(np.all(da<0))

    def test_on_gfa(self):
        import numpy as np
        telra, teldec = 10, 20
        gfa = GFALocations()
        x = gfa.gfatable['X']
        y = gfa.gfatable['Y']
        ra, dec = xy2radec(telra, teldec, x, y)

        # GFA corners are on GFAs if scale>1.0
        gfa_loc = on_gfa(telra, teldec, ra, dec, scale=1.1)
        self.assertEqual(len(gfa_loc), len(ra))
        self.assertTrue(np.all(gfa_loc >= 0))

        # Test cases that aren't on a GFA
        ra = np.random.uniform(telra-1, telra+1, size=10)
        dec = np.random.uniform(teldec-1, teldec+1, size=10)
        gfa_loc = on_gfa(telra, teldec, ra, dec, scale=1.1)
        self.assertEqual(len(gfa_loc), len(ra))
        self.assertTrue(np.all(gfa_loc == -1))

    def test_on_tile_gfa(self):
        #- Generate some targets near a tile; some should land on the GFAs
        tile = io.load_tiles()[0]
        targets = Table()
        ntargets = 1000
        ra = np.random.uniform(tile['RA']-1.61, tile['RA']+1.61, size=ntargets)
        dec = np.random.uniform(tile['DEC']-1.61, tile['DEC']+1.61, size=ntargets)
        targets['RA'] = (ra+360) % 360
        targets['DEC'] = dec
        gfa_targets = on_tile_gfa(tile['TILEID'], targets)
        self.assertGreater(len(gfa_targets), 0)
        self.assertIn('GFA_LOC', gfa_targets.dtype.names)

        #- Test higher level get_gfa_targets
        targets['FLUX_R'] = 2000
        gfatargets = get_gfa_targets(targets, rfluxlim=1000)
        self.assertGreater(len(gfatargets), 0)
        self.assertLess(np.max(gfatargets['GFA_LOC']), 10)
        self.assertGreaterEqual(np.min(gfatargets['GFA_LOC']), 0)
        self.assertTrue(tile['TILEID'] in gfatargets['TILEID'])

        #- Shift +60 deg and verify that none are found for that tile
        targets['RA'] = (targets['RA'] + 60) % 360
        gfatargets = on_tile_gfa(tile['TILEID'], targets)
        self.assertEqual(len(gfatargets), 0)

        # Test an empty table
        gfatargets = get_gfa_targets(targets, rfluxlim=3000)
        self.assertEqual(len(gfatargets), 0)

    def test_qs_on_gfa(self):
        n = 1000
        s = np.random.uniform(391, 407, size=n)
        q = np.random.uniform(0, 360, size=n)
        gfa = GFALocations()
        for gfaloc in range(10):
            ok = gfa.qs_on_gfa(gfaloc, q, s)
            self.assertTrue(np.any(ok))

    def test_gfa_coverage(self):
        '''Test corner cases for testing coverage'''
        gfatable = io.load_gfa().copy()
        f = GFALocations(gfatable=gfatable)
        if 'PETAL' in gfatable.colnames:
            gfatable.rename_column('PETAL', 'GFA_LOC')
        else:
            gfatable.rename_column('GFA_LOC', 'PETAL')
        f = GFALocations(gfatable=gfatable)
        x, y = np.random.uniform(-100, 100, (2,10))
        ok = f.xy_on_gfa(0, x, y)
        self.assertEqual(len(ok), len(x))
        self.assertFalse(np.any(ok))
        ok = f.xy_on_gfa(1, x[0], y[0])
        self.assertFalse(ok)

        with self.assertRaises(ValueError):
            ok = f.xy_on_gfa(99, x[0], y[0])    #- no GFA_LOC=99

        if 'PETAL' in gfatable.colnames:
            del gfatable['PETAL']

        if 'GFA_LOC' in gfatable.colnames:
            del gfatable['GFA_LOC']

        with self.assertRaises(ValueError):
            fx = GFALocations(gfatable=gfatable)

        #- targets must have RA,DEC or TARGET_RA,TARGET_DEC
        telra, teldec = 10, 20
        targets = Table()
        targets['RA'] = np.random.uniform(5, 15, size=10)
        targets['DEC'] = np.random.uniform(5, 15, size=10)
        ok = f.targets_on_gfa(telra, teldec, targets)
        targets.rename_column('RA', 'TARGET_RA')
        targets.rename_column('DEC', 'TARGET_DEC')
        ok = f.targets_on_gfa(telra, teldec, targets)

        targets['BLAT'] = np.arange(len(targets))
        del targets['TARGET_RA']
        del targets['TARGET_DEC']
        with self.assertRaises(ValueError):
            ok = f.targets_on_gfa(telra, teldec, targets)  #- no ra/dec
