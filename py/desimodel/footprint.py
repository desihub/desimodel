#- Utility functions for working with the DESI footprint

import numpy as np
from . import focalplane
from . import io

def _embed_sphere(ra, dec):
    """ embed RA DEC to a uniform sphere in three-d """
    phi = np.radians(np.asarray(ra))
    theta = np.radians(90.0 - np.asarray(dec))
    r = np.sin(theta)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = np.cos(theta)
    return np.array((x, y, z)).T

def is_point_in_desi(tiles, ra, dec, radius=None, return_tile_index=False):
    """Return if points given by ra, dec lie in the set of _tiles.

    This function is optimized to query a lot of points.
    radius is in units of degrees.

    `tiles` is the result of load_tiles.

    If a point is within `radius` distance from center of any tile,
    it is in desi.

    The shape of ra, dec must match. The current implementation
    works only if they are both 1d vectors or scalars.

    If return_tile_index is True, return the index of the nearest tile in tiles array.

    default radius is from desimodel.focalplane.get_tile_radius_deg()
    """
    from scipy.spatial import cKDTree as KDTree
    
    if radius is None:
        radius = focalplane.get_tile_radius_deg()

    tilecenters = _embed_sphere(tiles['RA'], tiles['DEC'])
    tree = KDTree(tilecenters)
    # radius to 3d distance
    threshold = 2.0 * np.sin(np.radians(radius) * 0.5)
    xyz = _embed_sphere(ra, dec)
    d, i = tree.query(xyz, k=1)

    indesi = d < threshold
    if return_tile_index:
        return indesi, i
    else:
        return indesi

def find_tiles_over_point(tiles, ra, dec, radius=None):
    """Return a list of indices of tiles that covers the points.

    This function is optimized to query a lot of points.
    radius is in units of degrees. The return value is an array
    of list objects that are the indices of tiles that cover each point.

    The indices are not sorted in any particular order.

    if ra, dec are scalars, a single list is returned.

    default radius is from desimodel.focalplane.get_tile_radius_deg()
    """
    from scipy.spatial import cKDTree as KDTree

    if radius is None:
        radius = focalplane.get_tile_radius_deg()

    tilecenters = _embed_sphere(tiles['RA'], tiles['DEC'])
    tree = KDTree(tilecenters)

    # radius to 3d distance
    threshold = 2.0 * np.sin(np.radians(radius) * 0.5)
    xyz = _embed_sphere(ra, dec)
    indices = tree.query_ball_point(xyz, threshold)
    return indices

def find_points_in_tiles(tiles, ra, dec, radius=None):
    """Return a list of indices of points that are within each provided tile(s).

    This function is optimized to query a lot of points with relatively few tiles.

    radius is in units of degrees. The return value is an array
    of lists that contains the index of points that are in each tile.
    The indices are not sorted in any particular order.

    if tiles is a scalar, a single list is returned.

    default radius is from desimodel.focalplane.get_tile_radius_deg()
    """
    from scipy.spatial import cKDTree as KDTree

    if radius is None:
        radius = focalplane.get_tile_radius_deg()

    # check for malformed input shapes. Sorry we currently only
    # deal with vector inputs. (for a sensible definition of indices)

    assert ra.ndim == 1
    assert dec.ndim == 1

    points = _embed_sphere(ra, dec)
    tree = KDTree(points)

    # radius to 3d distance
    threshold = 2.0 * np.sin(np.radians(radius) * 0.5)
    xyz = _embed_sphere(tiles['RA'], tiles['DEC'])
    indices = tree.query_ball_point(xyz, threshold)
    return indices

#
#
#
def get_tile_radec(tileid):
    """Return (ra, dec) in degrees for the requested tileid.

    If tileid is not in DESI, return (0.0, 0.0)
    TODO: should it raise and exception instead?
    """
    tiles = io.load_tiles()
    if tileid in tiles['TILEID']:
        i = np.where(tiles['TILEID'] == tileid)[0][0]
        return tiles[i]['RA'], tiles[i]['DEC']
    else:
        return (0.0, 0.0)

