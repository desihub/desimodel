# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
desimodel.focalplane.gfa
========================

Guide-Focus-Alignment location routines
"""

import numpy as np
from matplotlib.path import Path
import astropy.table

from ..io import load_deviceloc, load_gfa, load_platescale, load_tiles
from ..footprint import find_points_in_tiles, find_points_radec, get_tile_radec
from .geometry import get_radius_deg, radec2xy, qs2xy

class GFALocations(object):
    def __init__(self, gfatable=None, scale=1.0):
        '''
        Wrapper class for GFA locations
        
        Options:
            gfatable: table with GFA corners in X, Y and either GFA_LOC or PETAL
            scale: scale factor for GFA size

        Notes:

          * GFA coordinates are identified by a positive integer GFA_LOC or PETAL
            and four rows of corners specified in the X,Y tangent plane to the focal surface.
          * Default loads the GFA locations from desimodel.io.load_gfa().
            scale>1 enables selecting targets that are just off of the GFA to assist
            with initial acquisition.
        '''
        from matplotlib.path import Path

        if gfatable is None:
            gfatable = load_gfa()
        
        if 'GFA_LOC' in gfatable.dtype.names:
            gfa_locations = gfatable['GFA_LOC']
        elif 'PETAL' in gfatable.dtype.names:
            gfa_locations = gfatable['PETAL']
        else:
            raise ValueError('gfatable must have GFA_LOC or PETAL column')
        
        self.gfatable = gfatable
        self.gfa_polygons = list()
        self.gfa_locations = np.array(sorted(set(gfa_locations)))
        self.max_radius_mm = 0.0
        for loc in self.gfa_locations:
            ii = (gfa_locations == loc)
            xcorners = gfatable['X'][ii]
            ycorners = gfatable['Y'][ii]
            
            if scale != 1.0:
                x0, y0 = np.mean(xcorners), np.mean(ycorners)
                xcorners = x0 + scale*(xcorners-x0)
                ycorners = y0 + scale*(ycorners-y0)
            
            r = np.sqrt(xcorners**2 + ycorners**2)
            self.max_radius_mm = max(self.max_radius_mm, np.max(r))
            
            self.gfa_polygons.append(Path(list(zip(xcorners, ycorners))))

        self.max_radius_deg = float(get_radius_deg(self.max_radius_mm, 0.0))

    def xy_on_gfa(self, gfa_loc, x, y):
        '''
        Return boolean array of whether points (x,y) are on GFA at GFA_LOC
        
        Args:
            gfa_loc: integer GFA location [0-9] = PETAL_LOC for DESI
            x, y: array of focal tangent plane locations in mm
        
        Returns boolean array of whether (x,y) is on GFA at location GFA_LOC
        '''
        if gfa_loc not in self.gfa_locations:
            raise ValueError('GFA_LOC {} not in {}'.format(gfa_loc, self.gfa_locations))

        scalar_inputs = np.isscalar(x)
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        xy = np.vstack([x, y]).T
        i = np.where(self.gfa_locations == gfa_loc)[0][0]
        inside = self.gfa_polygons[i].contains_points(xy)
        if scalar_inputs:
            return inside[0]
        else:
            return inside

    def qs_on_gfa(self, gfa_loc, q, s):
        '''
        Return boolean array of whether points (q,s) are on GFA at GFA_LOC
        
        Args:
            gfa_loc: integer GFA location ID
            q: angular coordinate in degrees
            s: radial coordinate along focal surface in mm
        
        Returns boolean array of whether (q,s) is on GFA at location GFA_LOC
        '''
        x, y = qs2xy(q, s)
        return self.xy_on_gfa(gfa_loc, x, y)
        
    def targets_on_gfa(self, telra, teldec, targets):
        '''
        Returns subset of targets table with new GFA_LOC column
        
        Args:
            telra, teldec: Telescope pointing (RA,dec) in degrees
            targets: table with columns RA, DEC
        
        Returns gfa_targets table: subset of targets with new GFA_LOC column
        '''
        #- Trim to targets that might be on this pointing
        ps = load_platescale()
        rmax = ps['theta'][-1]
        
        columns = targets.dtype.names
        if ('RA' in columns) and ('DEC' in columns):
            ra = targets['RA']
            dec = targets['DEC']
        elif ('TARGET_RA' in columns) and ('TARGET_DEC' in columns):
            ra = targets['TARGET_RA']
            dec = targets['TARGET_DEC']
        else:
            raise ValueError('input targets must have columns RA,DEC or TARGET_RA,TARGET_DEC; has columns {}'.format(targets.dtype.names))
        
        ii = find_points_radec(telra, teldec, ra, dec, radius=rmax)
        targets = targets[ii]
        ra = ra[ii]
        dec = dec[ii]

        #- Identify targets on GFAs
        x, y = radec2xy(telra, teldec, ra, dec)
        gfa_loc = np.zeros(len(targets), dtype=int) - 1
        keep = np.zeros(len(targets), dtype=bool)
        for loc in self.gfa_locations:
            ii = self.xy_on_gfa(loc, x, y)
            keep |= ii
            gfa_loc[ii] = loc
        
        results = astropy.table.Table(targets[keep], copy=False)
        results['GFA_LOC'] = gfa_loc[keep].astype(np.int16)
        
        assert np.all(results['GFA_LOC'] >= 0)
        assert np.all(np.in1d(results['GFA_LOC'], self.gfa_locations))
        
        return results
        
#-------------------------------------------------------------------------

def on_gfa(telra, teldec, ra, dec, scale=1.0):
    """Checks if a target is on any of the 10 GFAs given `telra`, `teldec`
    and an array of `ra` and `dec` pointings, as well as a parameter for
    degrees of tolerance one would like to allow.

    Parameters
    ----------
    telra : :class:`float`
        The telescope's arbitrary RA pointing.
    teldec : :class:`float`
        The telescope's arbitrary Dec pointing.
    ra : array-like
        An array of RA values for locations in the sky.
    dec : array-like
        An array of Dec values for locations in the sky.
    scale : :class:`float`, optional
        Scale factor for GFA size to allow slightly off the edge targets

    Returns
    -------
    array
        An array of the same length as input `ra`, giving the GFA_LOC that
        each (ra, dec) is on.  -1 means not on a GFA.
    """
    gfa = GFALocations(scale=scale)
    x, y = radec2xy(telra, teldec, ra, dec)
    gfaloc = np.full(len(x), -1, dtype=int)
    for loc in gfa.gfa_locations:
        ii = gfa.xy_on_gfa(loc, x, y)
        gfaloc[ii] = loc

    return gfaloc
    

def on_tile_gfa(tileid, targets, scale=1.0):
    """This function takes a tileid, a table of targets, and an optional
    buffer_arcsec parameter to return the indices of targets lying on the GFA
    as well as the GFA locations from 0-9.

    Parameters
    ----------
    tileid : :class:`int`
        DESI tile ID, used to lookup telescope (RA, dec).
    targets : Table
        Table with columns RA, DEC.
    scale : :class:`float`, optional
        Scale factor for GFA size to allow slightly off the edge targets

    Returns
    -------
    Table
        Subset of input targets with new GFA_LOC column
    """
    gfa = GFALocations(scale=scale)
    telra, teldec = get_tile_radec(tileid)
    return gfa.targets_on_gfa(telra, teldec, targets)


def get_gfa_targets(targets, rfluxlim=1000, tiles=None, scale=1.0):
    """This function takes a table of targets, as well as optional parameters
    including a minimum flux in the r-band, a list of tiles, and a buffer in arcseconds
    and returns a table of targets on the GFA satisfying a minimum flux_r

    Parameters
    ----------
    targets : Table
        Table columns RA, DEC, FLUX_R.
    rfluxlim : :class:`float`, optional
        r-band flux limit; default 1000 = rmag 15.
    tiles : Table, optional
        Table of tiles, default to :func:`desimodel.io.load_tiles`.
    scale : :class:`float`, optional
        Scale factor for GFA size to allow slightly off the edge targets

    Returns
    -------
    Table
        A subset of input `targets` with additional columns:
        TILEID: (integer) DESI tile ID; GFA_LOC: (integer) GFA location [0-9].

    Notes
    -----
    * The same target could be repeated with different TILEID, GFA_LOC.
    * The function returns an empty Table if no targets are on any GFAs or
      of sufficient brightness.
    """
    # Checks if the flux_r meets a minimum threshold
    targets = astropy.table.Table(targets[targets['FLUX_R'] > rfluxlim])
    if len(targets) == 0:
        targets['TILEID'] = np.zeros(0, dtype=np.int32)
        targets['GFA_LOC'] = np.zeros(0, dtype=np.int8)
        return targets

    if(tiles is None):
        tiles = load_tiles()
    
    gfa = GFALocations(scale=scale)

    points = find_points_in_tiles(tiles, targets['RA'], targets['DEC'],
        radius=gfa.max_radius_deg)

    target_tables = list()

    for i, pointlist in enumerate(points):
        if pointlist:
            tileid = tiles[i]['TILEID']
            tiletargets = on_tile_gfa(tileid, targets, scale=scale)
            if len(tiletargets) > 0:
                tiletargets['TILEID'] = tileid
                target_tables.append(tiletargets)
            
    return astropy.table.vstack(target_tables)
