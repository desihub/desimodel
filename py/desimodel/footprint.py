#- Utility functions for working with the DESI footprint

import numpy as np
import os
from . import focalplane
from . import io

def radec2pix(nside, ra, dec):
    '''Convert ra,dec to nested pixel number

    Args:
        ra: float or array, Right Accention in degrees
        dec: float or array, Declination in degrees

    Returns:
        array of integer pixel numbers using nested numbering scheme

    Note: this is syntactic sugar around
    `hp.ang2pix(nside, ra, dec, lonlat=True, nest=True)`, but also works
    with older versions of healpy that didn't have `lonlat` yet.
    '''
    import healpy as hp
    theta, phi = np.radians(90-dec), np.radians(ra)
    return hp.ang2pix(nside, theta, phi, nest=True)

def tiles2pix(nside, tiles=None, radius=None, per_tile=False):
    '''
    Returns sorted array of pixels that overlap the tiles

    Args:
        nside: integer healpix nside, 2**k where 0 < k < 30

    Optional:
        tiles:
            array-like integer tile IDs; or
            integer tile ID; or
            Table-like with RA,DEC columns; or
            None to use all DESI tiles from desimodel.io.load_tiles()
        radius: tile radius in degrees;
            if None use desimodel.focalplane.get_tile_radius_deg()
        per_tile: if True, return a list of arrays of pixels per tile

    Returns pixels:
        integer array of pixel numbers that cover these tiles; or
        if per_tile is True, returns list of arrays such that pixels[i]
            is an array of pixel numbers covering tiles[i]
    '''
    import healpy as hp
    if tiles is None:
        import desimodel.io
        tiles = desimodel.io.load_tiles()

    if radius is None:
        import desimodel.focalplane
        radius = desimodel.focalplane.get_tile_radius_deg()

    theta, phi = np.radians(90-tiles['DEC']), np.radians(tiles['RA'])
    vec = hp.ang2vec(theta, phi)
    ipix = [hp.query_disc(nside, vec[i], radius=np.radians(radius),
                inclusive=True, nest=True) for i in range(len(tiles))]
    if per_tile:
        return ipix
    else:
        return np.sort(np.unique(np.concatenate(ipix)))

def tileids2pix(nside, tileids, radius=None, per_tile=False):
    '''
    Like tiles2pix, but accept integer tileid or list of tileids instead
    of table of tiles
    '''
    import desimodel.io
    tiles = desimodel.io.load_tiles()
    ii = np.in1d(tiles['TILEID'], tileids)
    if np.count_nonzero(ii) > 0:
        return tiles2pix(nside, tiles[ii], radius=radius, per_tile=per_tile)
    else:
        raise ValueError('TILEID(s) {} not in DESI footprint'.format(tileids))

def tiles2fracpix(nside=64, step=1, tiles=None, radius=None):
    '''
    Returns a sorted array of just the *fractional* pixels that overlap the tiles

    Optional Args:
        nside: integer healpix nside, 2**k where 0 < k < 30
        step: The number of integration steps around the edges of 
              a HEALPix pixel. step=1 means just the pixel vertices (e.g., see
              http://healpy.readthedocs.io/en/latest/generated/healpy.boundaries.html)
              step=2 means the vertices and the corners and the points halfway 
              between the vertices.
        tiles:
            Table-like with RA,DEC columns; or
            None to use all DESI tiles from desimodel.io.load_tiles()
        radius: tile radius in degrees;
            if None use desimodel.focalplane.get_tile_radius_deg()

    Returns fracpix:
        integer array of pixel numbers that cover these tiles, *excluding
        pixels that fully overlap the tiles* (i.e., just pixels that
        *partially* overlap the tiles). The integers are sorted.

    Notes:
        there are potentially malicious cases where a pixel just brushes
        a tile, such that there is a very small area where the pixel overlaps
        the tile. To guard against these case, call this function with
        progressively larger step values until it converges.
    '''
    #ADM set up healpy and set default tiles and radius
    import healpy as hp
    import desimodel
    if tiles is None:
        tiles = desimodel.io.load_tiles()

    if radius is None:
        radius = desimodel.focalplane.get_tile_radius_deg()

    #ADM obtain ALL pixels that overlap the tiles
    pix = desimodel.footprint.tiles2pix(nside, tiles=tiles, radius=radius) 

    #ADM the recovered number of pixels, and the total number of points
    #ADM that will be integrated around the boundary of the pixel
    npix = len(pix)
    nvertsperpix = 4*step

    #ADM find points around the boundary of all pixels in Cartesian coordinates
    xyzverts = hp.boundaries(nside,pix,step=step,nest=True)
    #ADM convert to RA/Dec
    theta, phi = hp.vec2ang(np.hstack(xyzverts).T) 
    ra, dec = np.degrees(phi), 90-np.degrees(theta) 

    #ADM calculate which boundary points are in the tiles
    verts_in = desimodel.footprint.is_point_in_desi(tiles,ra,dec,radius=radius) 

    #ADM reshape this into an array with nvertsperpix columns
    pix_verts_in = np.reshape(verts_in,(npix,nvertsperpix))
    #ADM any row with a column not in the tiles must be a fractional pixel
    isfracpix = ~np.all(pix_verts_in,axis=1)

    #ADM the pixel integers where pixels are fractional
    return pix[np.where(isfracpix)]


def pixweight(nside=256, tiles=None, radius=None, write=True):
    '''
    Create a rec array of the fraction of each pixel that overlaps the passed tiles

    Optional Args:
        nside: integer healpix nside, 2**k where 0 < k < 30
        tiles:
            Table-like with RA,DEC columns; or
            None to use all DESI tiles from desimodel.io.load_tiles()
        radius: tile radius in degrees;
            if None use desimodel.focalplane.get_tile_radius_deg()
        write: if True, then write the pixel->weight array to file in the

    Returns pixweight:
        a `~numpy.ndarray` where the first column is all possible pixels at the
        passed nside (called `HPXPIXEL`) and the second column (called `WEIGHT`)
        is the fracion of the pixel that overlaps the passed tiles:
        `WEIGHT=1` for the pixel is entirely contained in the tiles
        `WEIGHT=0` for the pixel is entirely outside of the tiles
        `0 < WEIGHT < 1` for a pixel that overlaps the tiles

    Notes:
        it's sufficient to create the weights at a suitably high nside, say
        nside=256 (0.052456 sq. deg. per pixel) as pixel numbers at
        lower nsides can be obtained by integer division by powers of 4, e.g.
        pix_@_nside_128 = pix@nside_256//4 and fractional weights at lower
        nsides are the mean of the 4 pixels at the higher nside
    '''

    
    if write:
        #ADM get path to DESIMODEL footprint directory, create output file name
        import desimodel.io
        outfile = os.path.join(desimodel.io.datadir(),
                               'footprint','desi-healpix-weights.fits')

        #ADM create rec array of pixels and weights and populate it
        outdata = np.empty(npix, dtype=[('HPXPIXEL', '>i8'), ('WEIGHT', '>f4')])
        outdata["HPXPIXEL"] = pix
        outdata["WEIGHT"] = weight

        #ADM write information indicating HEALPix setup to file header
        from desiutil import depend
        import fitsio
        hdr = fitsio.FITSHDR()
        depend.setdep(hdr, 'HPXNSIDE', nside)
        depend.setdep(hdr, 'HPXNEST', True)

        fitsio.write(outfile, outdata, extname='PIXWEIGHTS', header=hdr, clobber=True)


def pix2tiles(nside, pixels, tiles=None, radius=None):
    '''
    Returns subset of tiles that overlap the list of pixels

    Args:
        nside: integer healpix nside, 2**k with 1 <= k <= 30
        pixels: array of integer pixels using nested numbering scheme

    Optional:
        tiles:
            Table-like with RA,DEC columns; or
            None to use all DESI tiles from desimodel.io.load_tiles()
        radius: tile radius in degrees;
            if None use desimodel.focalplane.get_tile_radius_deg()

    Returns:
        table of tiles that cover these pixels

    TODO: add support for tiles as integers or list/array of integer TILEIDs
    '''
    import healpy as hp
    import desimodel.footprint

    if tiles is None:
        import desimodel.io
        tiles = desimodel.io.load_tiles()

    if radius is None:
        import desimodel.focalplane
        radius = desimodel.focalplane.get_tile_radius_deg()

    #- Trim tiles to ones that *might* overlap these pixels
    theta, phi = hp.pix2ang(nside, pixels, nest=True)
    ra, dec = np.degrees(phi), 90 - np.degrees(theta)
    pixsize = np.degrees(hp.nside2resol(nside))
    ii = desimodel.footprint.find_tiles_over_point(tiles, ra, dec, radius=radius+pixsize)
    if np.isscalar(pixels):
        tiles = tiles[ii]
    else:
        ii = np.unique(np.concatenate(ii))
        tiles = tiles[ii]

    #- Now check in detail
    theta, phi = np.radians(90-tiles['DEC']), np.radians(tiles['RA'])
    vec = hp.ang2vec(theta, phi)
    ii = list()
    for i in range(len(tiles)):
        tilepix = hp.query_disc(nside, vec[i], radius=np.radians(radius), inclusive=True, nest=True)
        if np.any(np.in1d(pixels, tilepix)):
            ii.append(i)
    return tiles[ii]

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

