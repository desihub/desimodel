# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.healpix
=================

Healpix utility functions for files in desimodel.

"""
import numpy as np
import healpy as hp

def tiles2pix(nside, tiles=None, radius=None):
    '''
    Returns sorted array of healpix pixels that overlap DESI tiles. 

    Args:
        nside: integer healpix nside, 2**k where 0 < k < 30

    Optional:
        tiles: table with RA,DEC columns;
            if None use all DESI tiles from desimodel.io.load_tiles()
        radius: tile radius in degrees;
            if None use desimodel.focalplane.get_tile_radius_deg()

    Returns:
        array of integer pixel numbers using nested numbering scheme
    '''
    if tiles is None:
        import desimodel.io
        tiles = desimodel.io.load_tiles()
    
    if radius is None:
        import desimodel.focalplane
        radius = desimodel.focalplane.get_tile_radius_deg()
    
    theta, phi = np.radians(90-tiles['DEC']), np.radians(tiles['RA'])
    vec = hp.ang2vec(theta, phi)
    ipix = [hp.query_disc(nside, vec[i], radius=np.radians(radius), inclusive=True, nest=True) for i in range(len(tiles))]
    return np.sort(np.unique(np.concatenate(ipix)))
