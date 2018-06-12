# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.focalplane.
"""
import unittest
import numpy as np
from ..focalplane import FocalPlane, generate_random_centroid_offsets, xy2radec, get_radius_mm, get_radius_deg, radec2xy
from ..focalplane import fiber_area_arcsec2, on_gfa, on_tile_gfa, get_gfa_targets
from ..focalplane import GFALocations
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
        degree = get_radius_deg(333.738, 217.766)
        truedegree = 1.5731300326614939
        radius = get_radius_mm(1.5731300326614939)
        trueradius = 398.5010456698936
        testradii = [1.5731300326614939, 1.5]
        radius1 = get_radius_mm(testradii)
        trueradius1 = [398.50104567, 378.53678987]
        self.assertAlmostEqual(truedegree, degree, 5)
        self.assertAlmostEqual(trueradius, radius, 5)
        self.assertAlmostEqual(all(trueradius1), all(radius1), 5)
        
        #- Arrays should also work
        x = np.array([333.738, 333.738])
        y = np.array([217.766, 217.766])
        truedegrees = np.array([truedegree, truedegree])
        degrees = get_radius_deg(x, y)
        self.assertTrue(np.allclose(degrees, truedegrees))

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
        area = fiber_area_arcsec2(x, y)
        self.assertFalse(np.any(np.fabs(area-np.array([1.97314482,  1.96207294,  1.92970506,  1.81091119])) >1E-6))

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


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
