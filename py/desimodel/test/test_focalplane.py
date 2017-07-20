# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.focalplane.
"""
import unittest
import numpy as np
from ..focalplane import FocalPlane, generate_random_centroid_offsets, xy2radec, get_radius_mm, get_radius_deg


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
        self.assertAlmostEqual(truedegree, degree)
        self.assertAlmostEqual(trueradius, radius)
    
    def new_test_xy2radec(self):
        """Should test the consistency between the conversion functions
        radec2xy and xy2radec. Right now tests the accuracy of the xy2radec 
        on particular cases.
        """
        truera = 8.927313423598427
        truedec = -9.324956250231294
        newra, newdec = xy2radec(8.37, -10.65, -138.345, -333.179)
        self.assertEqual(truera, newra)
        self.assertEqual(truedec, newdec)
    

    def test_xy2radec(self):
        """Test the consistency between the conversion functions
        radec2xy and xy2radec.
        """

        ra_list = np.array([0.0, 39.0, 300.0, 350.0, 359.9999, 20.0])
        dec_list = np.array([0.0, 45.0, 89.9999, -89.9999, 0.0, 89.9999])

        F = FocalPlane(ra_list, dec_list)

        # First test to check that it places x, y at the center.
        ra_obj = ra_list.copy()
        dec_obj = dec_list.copy()

        x_obj, y_obj = F.radec2xy(ra_obj, dec_obj)
        self.assertFalse(np.any(np.fabs(x_obj) > 1E-6) |
                         np.any(np.fabs(y_obj) > 1E-6),
                         ("Test Failed to recover position center with 1E-6 " +
                          "precision."))

        # Second test to check that it recovers input ra_obj,dec_obj.
        ra_out, dec_out = F.xy2radec(x_obj, y_obj)
        self.assertFalse(np.any(np.fabs(ra_out-ra_obj) > 1E-6) |
                         np.any(np.fabs(dec_out-dec_obj) >1E-6),
                         ("Test Failed to recover the input RA, Dec with " +
                          "1E-6 precision"))


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
