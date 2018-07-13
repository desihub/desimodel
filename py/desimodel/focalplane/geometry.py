# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
desimodel.focalplane.geometry
=============================

Dimensions and coordinate system transforms for the DESI focal plane.
"""

import os
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
from astropy.table import Table

from ..footprint import find_points_in_tiles, find_points_radec, get_tile_radec
from ..io import (load_desiparams, load_fiberpos, load_platescale,
    load_tiles, load_deviceloc)

_tile_radius_deg = None
_tile_radius_mm = None

def get_tile_radius_mm():
    '''Returns radius in mm to the middle of the outermost positioner.
    '''
    global _tile_radius_mm
    if _tile_radius_mm is None:
        fp = load_fiberpos()
        p = load_desiparams()
        _tile_radius_mm = np.sqrt(np.max(fp['X']**2 + fp['Y']**2))  # + p['positioners']['radius_max']
    return _tile_radius_mm


def get_tile_radius_deg():
    '''Returns radius in degrees to the middle of the outermost positioner.
    '''
    global _tile_radius_deg
    if _tile_radius_deg is None:
        rmax = get_tile_radius_mm()
        platescale = load_platescale()
        fn = interp1d(platescale['radius'], platescale['theta'], kind='quadratic')
        _tile_radius_deg = float(fn(rmax))
    return _tile_radius_deg

def get_radius_mm(theta):
    """Returns an array of radii in mm given an array of radii in degrees using the platescale data
    relative to the center of the focal plane as (0,0). Supports scalar and vector inputs.

    Parameters
    ----------
    theta : :class:`float` or array-like
        An array that represents the angle from the center of the focal plane.

    Returns
    -------
    :class:`float` or array-like
        Radii in mm.
    """
    platescale = load_platescale()
    # Uses a quadratic one-dimensional interpolation to approximate the radius in degrees versus radius in mm
    fn = interp1d(platescale['theta'], platescale['radius'], kind = 'quadratic')
    radius = fn(theta)
    if(np.isscalar(theta)):
        return float(radius)
    else:
        return radius

def get_radius_deg(x, y):
    """Returns the radius in degrees given `x`, `y` coordinates using the
    platescale data.

    Parameters
    ----------
    x : :class:`float`
        The x coordinate in mm of a location on the focal plane
    y : :class:`float`
        The y coordinate in mm of a location on the focal plane

    Returns
    -------
    :class:`float`
        Radius corresponding to `x`, `y`.
    """
    radius = np.sqrt(x**2 + y**2)
    platescale = load_platescale()
    fn = interp1d(platescale['radius'], platescale['theta'],
                                    kind='quadratic')
    degree = fn(radius).astype(float)
    return degree

def _extrapolate_r_s(r, s):
    '''
    Utility function for xy2qs and qs2xy; returns new r, s with extrapolations
    to 0 and max(r)+10 mm.
    '''
    # quadratic extrapolation beyond rmax
    ii = (r>410)
    c = np.polyfit(r[ii], s[ii], 2)
    xr = np.linspace(np.max(r)+0.1, np.max(r)+10, 10)
    xs = np.polyval(c, xr)

    full_r = np.concatenate([[0.0,], r, xr])
    full_s = np.concatenate([[0.0,], s, xs])
    return full_r, full_s


def xy2qs(x, y):
    '''Focal tangent plane x,y -> angular q,s on curved focal surface

    Args:
        x, y: cartesian location on focal tangent plane in mm

    Returns (q, s) where q=angle in degrees; s=focal surface radial dist [mm]

    Notes: (x,y) are in the "CS5" DESI coordinate system tangent plane to
    the curved focal surface.  q is the radial angle measured counter-clockwise
    from the x-axis; s is the radial distance along the curved focal surface;
    it is *not* sqrt(x**2 + y**2).  (q,s) are the preferred coordinates for
    the DESI focal plane hardware engineering team.
    '''
    #- Load device locations, use just petal 0 for interpolation
    d = load_deviceloc()
    d = d[d['PETAL'] == 0]

    q = np.degrees(np.arctan2(y, x))
    q[q<0] += 360.0
    r = np.sqrt(x**2 + y**2)

    dr, ds = _extrapolate_r_s(np.sqrt(d['X']**2 + d['Y']**2), d['S'])
    ii = np.argsort(dr)

    fn = interp1d(dr[ii], ds[ii], kind='quadratic', assume_sorted=True)
    s = fn(r)

    return q, s

def qs2xy(q, s):
    '''angular q,s on curved focal surface -> focal tangent plane x,y

    Args:
        q: angle in degrees
        s: focal surface radial distance in mm

    Returns (x, y) cartesian location on focal tangent plane in mm

    Notes: (x,y) are in the "CS5" DESI coordinate system tangent plane to
    the curved focal surface.  q is the radial angle measured counter-clockwise
    from the x-axis; s is the radial distance along the curved focal surface;
    it is *not* sqrt(x**2 + y**2).  (q,s) are the preferred coordinates for
    the DESI focal plane hardware engineering team.
    '''
    #- Load device locations, use just petal 0 for interpolation
    d = load_deviceloc()
    d = d[d['PETAL'] == 0]

    dr, ds = _extrapolate_r_s(np.sqrt(d['X']**2 + d['Y']**2), d['S'])
    ii = np.argsort(ds)

    fn = interp1d(ds[ii], dr[ii], kind='quadratic')
    r = fn(s)
    x = r*np.cos(np.radians(q))
    y = r*np.sin(np.radians(q))

    return x, y

def xy2radec(telra, teldec, x, y):
    """Returns the new RA and Dec of an `x`, `y` position on the focal plane
    in the sky given an arbitrary telescope pointing in RA and Dec.

    Parameters
    ----------
    telra : :class:`float`
        The telescope's RA pointing in degrees.
    teldec : :class:`float`
        The telescope's Dec pointing in degrees.
    x : :class:`float`
        The x coordinate in mm of a location on the focal plane
    y : :class:`float`
        The y coordinate in mm of a location on the focal plane

    Returns
    -------
    tuple
        The RA, Dec corresponding to `x`, `y`.
    """
    # radial distance on the focal plane in radians
    r_rad = np.radians(get_radius_deg(x, y))

    # q signifies the angle the position makes with the +x-axis of focal plane
    q = np.degrees(np.arctan2(y, x))
    q_rad = np.radians(q)

    # Consider a unit sphere (x,y,z)
    # starting at (RA,dec) = (0,0) -> v0 = (1,0,0)
    # v0 = np.array([1.0, 0.0, 0.0])

    # The focal plane is oriented with +yfocal = +dec but +xfocal = -RA
    # Rotate clockwise around z by r_rad
    # zrotate = np.zeros(shape=(3,3))
    # zrotate[0] = [np.cos(r_rad), np.sin(r_rad), 0]
    # zrotate[1] = [-np.sin(r_rad), np.cos(r_rad), 0]
    # zrotate[2] = [0, 0, 1]
    # v1 = zrotate.dot(v0)

    x1 = np.cos(r_rad)      # y0=0 so drop sin(r_rad) term
    y1 = -np.sin(r_rad)     # y0=0 so drop cos(r_rad) term
    z1 = np.zeros_like(x1)

    # clockwise rotation around the x-axis
    # xrotate = np.zeros(shape=(3,3))
    # q_rad = np.radians(q)
    # xrotate[0] = [1, 0, 0]
    # xrotate[1] = [0, np.cos(q_rad), np.sin(q_rad)]
    # xrotate[2] = [0, -np.sin(q_rad), np.cos(q_rad)]

    x2 = x1
    y2 = y1*np.cos(q_rad)           # z1=0 so drop sin(q_rad) term
    z2 = -y1*np.sin(q_rad)          # z1=0 so drop cos(q_rad) term
    v2 = np.stack([x2, y2, z2])

    # Clockwise rotation around y axis by declination of the tile center
    decrotate = np.zeros(shape=(3,3))
    teldec_rad = np.radians(teldec)
    decrotate[0] = [np.cos(teldec_rad), 0, -np.sin(teldec_rad)]
    decrotate[1] = [0, 1, 0]
    decrotate[2] = [np.sin(teldec_rad), 0, np.cos(teldec_rad)]

    # Counter-clockwise rotation around the z-axis by the right ascension of the tile center
    rarotate = np.zeros(shape=(3,3))
    telra_rad = np.radians(telra)
    rarotate[0] = [np.cos(telra_rad), -np.sin(telra_rad), 0]
    rarotate[1] = [np.sin(telra_rad), np.cos(telra_rad), 0]
    rarotate[2] = [0, 0, 1]

    x3, y3, z3 = v3 = rarotate.dot(decrotate.dot(v2))

    ra_deg = np.degrees(np.arctan2(y3, x3)) % 360
    dec_deg = np.degrees((np.pi/2) - np.arccos(z3))

    return ra_deg, dec_deg


def radec2xy(telra, teldec, ra, dec):
    """Returns arrays of the x, y positions of given celestial objects
    on the focal plane given an arbitrary telescope pointing in RA and Dec and
    arrays of the `ra` and `dec` of celestial objects in the sky.

    Parameters
    ----------
    telra : :class:`float`
        The telescope's RA pointing in degrees.
    teldec : :class:`float`
        The telescope's Dec pointing in degrees.
    ra : array-like
        An array of RA values for locations in the sky.
    dec : array-like
        An array of Dec values for locations in the sky.

    Returns
    -------
    tuple
        The x, y positions corrsponding to `ra`, `dec`.

    Notes
    -----
    Implements the Haversine formula.
    """
    # Inclination is 90 degrees minus the declination in degrees
    dec = np.asarray(dec)
    inc = 90 - dec
    ra = np.asarray(ra)
    x0 = np.sin(np.radians(inc)) * np.cos(np.radians(ra))
    y0 = np.sin(np.radians(inc)) * np.sin(np.radians(ra))
    z0 = np.cos(np.radians(inc))
    coord = [x0, y0, z0]

    # Clockwise rotation around the z-axis by the right ascension of the tile center
    rarotate = np.zeros(shape=(3,3))
    telra_rad = np.radians(telra)
    rarotate[0] = [np.cos(telra_rad), np.sin(telra_rad), 0]
    rarotate[1] = [-np.sin(telra_rad), np.cos(telra_rad), 0]
    rarotate[2] = [0, 0, 1]

    # Counter-Clockwise rotation around y axis by declination of the tile center
    decrotate = np.zeros(shape=(3,3))
    teldec_rad = np.radians(teldec)
    decrotate[0] = [np.cos(teldec_rad), 0, np.sin(teldec_rad)]
    decrotate[1] = [0, 1, 0]
    decrotate[2] = [-np.sin(teldec_rad), 0, np.cos(teldec_rad)]

    coord1 = np.matmul(rarotate, coord)
    coord2 = np.matmul(decrotate, coord1)
    x = coord2[0]
    y = coord2[1]
    z = coord2[2]

    newteldec = 0
    newtelra = 0
    ra_rad = np.arctan2(y, x)
    dec_rad = (np.pi / 2) - np.arccos(z / np.sqrt((x**2) + (y**2) + (z**2)))
    radius_rad = 2 * np.arcsin(np.sqrt((np.sin((dec_rad - newteldec) / 2)**2) + ((np.cos(newteldec)) * np.cos(dec_rad) * (np.sin((ra_rad - newtelra) / 2)**2))))
    radius_deg = np.degrees(radius_rad)

    q_rad = np.arctan2(z, -y)

    radius_mm = get_radius_mm(radius_deg)

    x_focalplane = radius_mm * np.cos(q_rad)
    y_focalplane = radius_mm * np.sin(q_rad)

    return x_focalplane, y_focalplane


class FocalPlane(object):
    """A class for modeling the DESI focal plane and converting between
    focal plane coordinates (in mm) and RA, Dec on the sky (in degrees).
    Provides utility functions for mapping which positioners cover
    which (RA, Dec) or (x, y) locations and vice versa.

    Parameters
    ----------
    ra, dec : :class:`float`
        Initialize DESI focal plane model with the telescope pointing
        at (`ra`, `dec`) in degrees.

    NOTE: this class is deprecated (or should be further expanded), but
    I'm not removing it yet in order to not arbitrarily break code that
    might be using it.
    """

    def __init__(self, ra=0.0, dec=0.0):
        """
        """
        # Read $DESIMODEL/data/focalplane/fiberpos.fits and platescale.txt
        # to construct focal plane model.  May also need data/desi.yaml .
        self._check_radec(ra, dec)
        self.ra = ra
        self.dec = dec
        self._fiberpos_file = os.path.join(os.environ['DESIMODEL'],
                                           'data', 'focalplane',
                                           'fiberpos.fits')
        with fits.open(self._fiberpos_file) as hdulist:
            self.fiberpos = hdulist[1].data

    def _check_radec(self, ra, dec):
        """Raise ValueError if RA or dec are out of bounds.
        """
        if np.any( (ra < 0) | (ra >= 360) ):
            raise ValueError("RA must be 0 <= RA < 360")
        if np.any( (dec < -90) | (dec > +90) ):
            raise ValueError("Dec must be -90 <= dec <= 90")

    def set_tele_pointing(self, ra, dec):
        """Set telescope pointing to (RA, Dec) in degrees.

        Parameters
        ----------
        ra
        dec : :class:`float`
            Telescope pointing in degrees.
        """
        self._check_radec(ra, dec)
        self.ra = ra
        self.dec = dec

    def radec2xy(self, ra, dec):
        """Convert (RA, Dec) in degrees to (x, y) in mm on the focal plane
        given the current telescope pointing.

        If RA and Dec are floats, returns a tuple (x, y) of floats
        If RA and Dec are numpy arrays, returns a tuple (x, y) of numpy arrays

        Parameters
        ----------
        ra, dec : :class:`float` or :class:`numpy.ndarray`
            Sky position.

        Returns
        -------
        :func:`tuple`
            A tuple containing the (x, y) coordinates in mm.
        """
        self._check_radec(ra, dec)
        return radec2xy(self.ra, self.dec, ra, dec)

    def xy2radec(self, x, y):
        """Convert (x, y) in mm on the focal plane to (ra_object, dec_object)
        in degrees on the sky given the current telescope pointing towards
        (RA, Dec).

        `x`, `y` must be floats. This function is vectorized in xy2radec(),
        which doesn't appear to exist.

        Parameters
        ----------
        x, y : :class:`float`
            Position on the focal plane in mm.

        Returns
        -------
        :func:`tuple`
            Coordinates of object.
        """
        return xy2radec(self.ra, self.dec, x, y)

    def radec2pos(self, ra, dec):
        """Identify which positioners cover (`ra`, `dec`).

        If `ra`, `dec` are floats, return an array of positioner IDs that
        cover it. The array could be empty if no positioner covers that
        location.

        If `ra`, `dec` are numpy arrays, return a list of arrays.
        The ith element is an array of positioner IDs that cover
        (ra[i], dec[i]).

        .. warning:: This method is not implemented!
        """
        raise NotImplementedError

    def xy2pos(self, x, y):
        """Identify which positioners cover (`x`, `y`).

        .. warning:: This method is not implemented!
        """
        raise NotImplementedError


def fiber_area_arcsec2(x, y):
    '''Returns area of fibers at (`x`, `y`) in arcsec^2.
    '''
    params = load_desiparams()
    fiber_dia = params['fibers']['diameter_um']
    x = np.asarray(x)
    y = np.asarray(y)
    r = np.sqrt(x**2 + y**2)

    #- Platescales in um/arcsec
    ps = load_platescale()
    radial_scale = np.interp(r, ps['radius'], ps['radial_platescale'])
    az_scale = np.interp(r, ps['radius'], ps['az_platescale'])

    #- radial and azimuthal fiber radii in arcsec
    rr  = 0.5 * fiber_dia / radial_scale
    raz = 0.5 * fiber_dia / az_scale
    fiber_area = (np.pi * rr * raz)
    return fiber_area

