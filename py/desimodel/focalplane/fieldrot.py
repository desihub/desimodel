# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
desimodel.focalplane.fieldrot
=============================

Routines to estimate the field rotation in the DESI focal plane.
The rotation angle is defined in DESI-5190.

1) The CS5 system is defined in DESI-481. It is a plane tangent to the focal plane surface, attached to the petals. It has the axis X pointing to the physical west (increasing HA, decreasing RA), Y to the physical north (increasing Dec), and Z towards the mirror.
2) The ICRS is the sky coordinate system of the DESI target catalogs, the same as Gaia astrometry. RA Dec will always be in the ICRS system.
3) The Field rotation angle 'Theta' measures the rotation of the star images in the CS5 plane (not a rotation of the instrument, opposite sign!), tangent to the focal plane, assumed for assigning fibers.
4) Theta=0 corresponds to zero rotation after the coordinate transformation from the ICRS RA,Dec sky coordinates to the CS5, with Y pointing to the north (increasing Dec).
5) Theta is increasing from North to East (or from X(West) to Y(North)).


The DESI field rotation is due to a combination of precession, aberration, and polar misalignement and general mount imperfections.

Here we will just consider the most important terms.


A good fraction of the code below if from Mike Lampton, imported in desimodel by Julien Guy.

"""

import numpy as np
from desimodel.focalplane.geometry import xy2radec,radec2xy

###################################################################################
# from http://asa.usno.navy.mil/static/files/2018/Astronomical_Constants_2018.pdf
OBLIQ                      = 23.439279444444445 # 23+26/60.+21.406/3600. , obliquity of the ecliptic, Initial  Values  at  J2000·0
DAYS_PER_YEAR              = 365.256363004 
PRECESSION_PERIOD_IN_YEARS = 25771.5753382206 # 360./(5028.796195/100./3600.) , Rates of precession at J2000·0 (IAU 2006) , General precession in longitude
ICRS_MJD                   = 51544.5 # 2451545.0-2400000.5 J2000
###################################################################################

def sind(a):
    return np.sin(np.radians(a))

def cosd(a):
    return np.cos(np.radians(a))

def put360(degrees): # Puts an angle into range 0 to 360.
    return np.fmod(720.+degrees, 360)

def arctan2d(y, x):
    return put360(np.degrees(np.arctan2(y, x)))

def arcsind(x):
    return np.degrees(np.arcsin(x))

def getXYZ(lon,lat):  # Convert spherical angles (in degrees) into xyz triplet
    return np.array([cosd(lon)*cosd(lat), sind(lon)*cosd(lat), sind(lat)])

def getNorm(xyz):
    return np.sqrt(np.sum(xyz**2,axis=0))
    
def getNormalized(xyz):
    return xyz/getNorm(xyz)

def getLONLAT(xyz): # Convert xyz into its spherical angles
    xyz = getNormalized(xyz)  # usually unnecessary
    return  arctan2d(xyz[1],xyz[0]) , arcsind(xyz[2])  # degrees

def vecX(xdeg):  # For positive xdeg=cwRoll: +y=>+z; +z=>-y.
    c=np.cos(np.radians(xdeg)); s=np.sin(np.radians(xdeg))
    return np.array([[1,0,0], [0,c,-s], [0,+s,c]])

def vecY(ydeg):  # For positive ydeg=-elev: +z=>+x; +x=>-z.
    # do not overlook this minus sign: positive ydeg pitches downward.
    c=np.cos(np.radians(ydeg)); s=np.sin(np.radians(ydeg))
    return np.array([[c,0,+s], [0,1,0], [-s,0,c]])
    
def vecZ(zdeg):  # For positive zdeg=+az: +x=>+y; +y=>-x.
    c=np.cos(np.radians(zdeg)); s=np.sin(np.radians(zdeg))
    return np.array([[c,-s,0], [+s,c,0], [0,0,1]])

def refX(xdeg):  # Rolls reference frame clockwise about X
    return vecX(-xdeg)

def refY(elev):  # Elevates reference frame about Y
    return vecY(+elev)
    
def refZ(azim):  # Rotates reference frame to new azimuth
    return vecZ(-azim)

def radec2eclip(ra,dec):  # same epoch
    equatorial_xyz = getXYZ(ra,dec)
    ecliptic_xyz = np.dot(refX(OBLIQ), equatorial_xyz)  # roll frame clockwise
    return getLONLAT(ecliptic_xyz)

def eclip2radec(lon,lat):  # same epoch
    ecliptic_xyz = getXYZ(lon,lat)
    equatorial_xyz = np.dot(refX(-OBLIQ), ecliptic_xyz)  # roll frame counterclockwise by obliq
    return getLONLAT(equatorial_xyz)

def precess(ra, dec, years) :
    # see DESI-4957
    # Equator and zero longitude point glide westward at 0.0139 deg/year, so..
    # Star ecliptic longitudes increase +0.0139 deg/year from combined lunar and solar torques on the Earth.
    # To precess any sta'’s {RA,DEC}:
    # 1. Convert to ecliptic coordinates {elon, elat}
    # 2. Add 0.0139 deg * number of years to elon
    # 3. Convert back to {RA,DEC}
    deltaELON = 360.* (years / PRECESSION_PERIOD_IN_YEARS) # degrees
    lon,lat=radec2eclip(ra,dec)
    xyz_ecliptic  = getXYZ(lon,lat)
    xyz_precessed = np.dot(vecZ(deltaELON), xyz_ecliptic) 
    lon,lat = getLONLAT(xyz_precessed)
    return eclip2radec(lon,lat)

def rotation_angle(x1,y1,x2,y2) :
    """
    returns the angle from (x1,y1) to (x2,y2), in deg
    """
    return np.mean( arcsind((x1*y2-x2*y1)/np.sqrt((x1**2+y1**2)*(x2**2+y2**2))) )

def field_rotation_angle(ra,dec,mjd,use_astropy=False) :
    """
    precessiom, etc: https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=4957

    Parameters
    ---------- 
    ra : Right ascension of center of focal plane in ICRS, in degrees
    dec : Declinaison  of center of focal plane in ICRS, in degrees
    mjd : Modified Julian Date, decimal number, of the observation
    
    Returns
    -------
    Field rotation angle, in degrees
    """

    # a cross 
    x1 = np.array([0,1,0,-1,0])
    y1 = np.array([0,0,1,0,-1])
    
    # ra dec of cross given telescope pointing
    ra1 , dec1 = xy2radec(ra, dec, x1, y1)
    
    # transformed coordinates
    if use_astropy : #  using astropy
        from astropy.coordinates import SkyCoord,FK5,GCRS
        import astropy.time
        fk5  = FK5(equinox=astropy.time.Time(mjd,format="mjd")) # precession
        c1   = SkyCoord(ra1,dec1, frame='icrs', unit='deg')
        c2   = c1.transform_to(fk5)
        ra2  = c2.ra.value
        dec2 = c2.dec.value
    else :
        # precession    
        years = (mjd-ICRS_MJD)/DAYS_PER_YEAR
        ra2, dec2 = precess(ra1, dec1, years)
    
    # back to x y
    x2 , y2 = radec2xy(ra2[0],dec2[0],ra2[1:],dec2[1:]) # first set of coords is pointing = field center

    # remove central point
    x1=x1[1:]
    y1=y1[1:]
    
    # compute angle
    return rotation_angle(x1,y1,x2,y2)
