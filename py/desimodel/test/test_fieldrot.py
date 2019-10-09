# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.focalplane.
"""
import unittest
import numpy as np
from ..focalplane.fieldrot import *

class TestFieldRot(unittest.TestCase):
    
    """Test desimodel.focalplane.fieldrot
    """
    def test_xyz_radec_transfo(self) :
        ra = np.array([33.,33.])
        dec = np.array([12.,12])
        xyz = getXYZ(ra,dec)
        tra,tdec = getLONLAT(xyz)
        assert(np.allclose(ra,tra))
        assert(np.allclose(dec,tdec))

    def test_equatorial_ecliptic_transfo(self) :
        ra = np.array([33.,33.])
        dec = np.array([12.,12])
        lon,lat=radec2eclip(ra,dec)
        tra,tdec=eclip2radec(lon,lat)
        assert(np.allclose(ra,tra))

    def test_angle_with_or_without_astropy(self) :
        mjd = 58577.4797 # 2019/04/03
        maxdiff=0.
        for ra in np.linspace(0,360,10) :
            for dec in np.linspace(-30,80,5) :
                angle1 = field_rotation_angle(ra,dec,mjd)
                angle2 = field_rotation_angle(ra,dec,mjd,use_astropy=True)
                print("RA={} Dec={} angle (home) = {:3.2f} deg, (astropy) = {:3.2f} deg, diff = {:3.2f} arcsec".format(ra,dec,angle1,angle2,(angle2-angle1)*3600.))
                maxdiff=max(maxdiff,np.abs(angle2-angle1)*3600.)
        print("Max. difference = {:3.2f} arcsec".format(maxdiff))
        assert(maxdiff<10.)

def test_suite():
    """Allows testing of only this module with the command::
       
       python setup.py test -m <modulename>
    """

    # usage, in base directory of desimodel :
    # python setup.py test -m desimodel.test.test_fieldrot
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
