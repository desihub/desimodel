# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
desimodel.footprint
===================

Utility functions for working with the DESI footprint.
"""
import os
from time import time
import numpy as np

from desiutil.log import get_logger

from .io import load_tiles
from . import __version__ as desimodel_version


log = get_logger()


_pass2program = None


def pass2program(tilepass, surveyops=False):
    '''Converts integer tile pass number to string program name.

    Args:
        tilepass (int or int array): tiling pass number.
    surveyops (bool): ``True`` to look for tiles in $DESI_SURVEYPOPS.

    Returns:
        Program name for each pass (str or list of str).
    '''
    global _pass2program
    # ADM this function isn't useful if looking in the DESI_SURVEYOPS
    # ADM directory as the real data is not ordered by pass.
    if surveyops == True:
        msg = "Function is not meaningful when looking in the DESI_SURVEYOPS "
        msg += "directory as the real, Main Survey, data is not ordered by pass"
        log.critical(msg)
        raise ValueError(msg)
    if _pass2program is None:
        tiles = load_tiles(surveyops=False)
        _pass2program = dict(set(zip(tiles['PASS'], tiles['PROGRAM'])))
    if np.isscalar(tilepass):
        return _pass2program[tilepass]
    else:
        return [_pass2program[p] for p in tilepass]


def program2pass(program):
    '''Convert string program name to tile passes for that program.

    Args:
        program (str for str array): program name, *e.g.* DARK, BRIGHT, or GRAY.

    Returns:
        List of integer passes that cover that program, or list of lists
        if input was array-like.
    '''
    tiles = load_tiles()

    if np.isscalar(program):
        passes = sorted(list(set(tiles['PASS'][tiles['PROGRAM'] == program])))
        if len(passes) > 0:
            return passes
        else:
            known_programs = set(tiles['PROGRAM'])
            msg = 'Unknown program {}; known programs are {}'.format(
                program, known_programs)
            raise ValueError(msg)
    else:
        program = np.asarray(program)
        passes = [None,] * len(program)
        for thisprogram in np.unique(program):
            thesepasses = program2pass(thisprogram)
            for i in np.where(program == thisprogram)[0]:
                passes[i] = thesepasses

        return passes


def radec2pix(nside, ra, dec):
    '''Convert `ra`, `dec` to nested pixel number.

    Args:
        nside (int): HEALPix `nside`, ``2**k`` where 0 < k < 30.
        ra (float or array): Right Accention in degrees.
        dec (float or array): Declination in degrees.

    Returns:
        Array of integer pixel numbers using nested numbering scheme.

    Notes:
        This is syntactic sugar around::

            hp.ang2pix(nside, ra, dec, lonlat=True, nest=True)

        but also works with older versions of healpy that didn't have
        `lonlat` yet.
    '''
    import healpy as hp
    theta, phi = np.radians(90-dec), np.radians(ra)
    if np.isnan(np.sum(theta)) :
        raise ValueError("some NaN theta values")

    if np.sum((theta < 0)|(theta > np.pi))>0 :
        raise ValueError("some theta values are outside [0,pi]: {}".format(theta[(theta < 0)|(theta > np.pi)]))

    return hp.ang2pix(nside, theta, phi, nest=True)


def tiles2pix(nside, tiles=None, radius=None, per_tile=False, fact=2**7):
    '''Returns sorted array of pixels that overlap the tiles.

    Args:
        nside (int): HEALPix `nside`, ``2**k`` where 0 < k < 30.
        tiles (Table-like, optional): tiles with columns RA,DEC or
            TILERA,TILEDEC, or None to load
            :func:`desimodel.io.load_tiles` (deprecated)
        radius (float, optional): tile radius in degrees;
            if ``None`` use :func:`desimodel.focalplane.get_tile_radius_deg`.
        per_tile (bool, optional): If ``True``, return a list of arrays of
            pixels per tile.
        fact (int, optional): Factor healpy uses to resolve pixel overlaps.
            When this is large there are fewer false positives at the expense
            of run time (although ``fact=2**8`` seems fast). Must be a
            power of 2.

    Returns:
        Integer array of pixel numbers that cover these tiles; or
        if per_tile is `True`, returns list of arrays such that ``pixels[i]``
        is an array of pixel numbers covering ``tiles[i]``.
    '''
    import healpy as hp
    from .focalplane import get_tile_radius_deg
    if tiles is None:
        log.warning('DESIMODEL tiles table is deprecated')
        tiles = load_tiles()

    if radius is None:
        radius = get_tile_radius_deg()

    if isinstance(tiles, dict):
        from astropy.table import Table
        tiles = Table(tiles)

    for col in ['RA', 'TILERA']:
        if col in tiles.dtype.names:
            ra = tiles[col]
            break
    else:
        raise ValueError('tiles table needs RA or TILERA')

    for col in ['DEC', 'TILEDEC']:
        if col in tiles.dtype.names:
            dec = tiles[col]
            break
    else:
        raise ValueError('tiles table needs DEC or TILEDEC')


    theta, phi = np.radians(90-dec), np.radians(ra)
    vec = hp.ang2vec(theta, phi)
    ipix = [hp.query_disc(nside, vec[i], radius=np.radians(radius),
                inclusive=True, nest=True, fact=fact) for i in range(len(ra))]
    if per_tile:
        return ipix
    else:
        return np.sort(np.unique(np.concatenate(ipix)))


def tileids2pix(nside, tileids, radius=None, per_tile=False):
    '''Like :func:`~desimodel.footprint.tiles2pix`, but accept integer
    tileid or list of tileids instead of table of tiles.
    '''
    tiles = load_tiles()
    ii = np.isin(tiles['TILEID'], tileids)
    if np.count_nonzero(ii) == np.asarray(tileids).size:
        return tiles2pix(nside, tiles[ii], radius=radius, per_tile=per_tile)
    else:
        extra = set(tileids) - set(tiles['TILEID'])
        raise ValueError('{}/{} TILEID(s) not in DESI footprint: {}'.format(
            len(extra), len(tileids), extra))

def tiles2fracpix(nside, step=1, tiles=None, radius=None, fact=2**7):
    '''Returns a sorted array of just the *fractional* pixels that overlap the
    tiles.

    Args:
        nside (int): HEALPix `nside`, ``2**k`` where 0 < k < 30.
        step (int, optional): The number of integration steps around the edges
            of a HEALPix pixel. ``step=1`` means just the pixel vertices.
            ``step=2`` means the vertices and the corners and the points halfway
            between the vertices.  See also the
            `HEALPix boundary document <http://healpy.readthedocs.io/en/latest/generated/healpy.boundaries.html>`_ .
        tiles (Table-like, optional): Table-like with RA,DEC columns; or
            ``None`` to use all DESI tiles from :func:`desimodel.io.load_tiles`.
        radius (float, optional): Tile radius in degrees;
            if ``None`` use :func:`desimodel.focalplane.get_tile_radius_deg`.
        fact (int, optional): Factor healpy uses to resolve pixel overlaps.
            When this is large there are fewer false positives at the expense
            of run time (although ``fact=2**8`` seems fast). Must be a
            power of 2.

    Returns:
        Integer array of pixel numbers that cover these tiles, *excluding
        pixels that fully overlap the tiles* (*i.e.*, just pixels that
        *partially* overlap the tiles). The integers are sorted.

    Notes:
        There are potentially malicious cases where a pixel just brushes
        a tile, such that there is a very small area where the pixel overlaps
        the tile. To guard against these case, call this function with
        progressively larger step values until it converges.
    '''
    #ADM set up healpy and set default tiles and radius
    import healpy as hp
    from .focalplane import get_tile_radius_deg

    if tiles is None:
        tiles = load_tiles()

    if radius is None:
        radius = get_tile_radius_deg()

    #ADM obtain ALL pixels that overlap the tiles (and perhaps a
    #ADM few more if fact is a small number
    pix = tiles2pix(nside, tiles=tiles, radius=radius, fact=fact)

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
    verts_in = is_point_in_desi(tiles, ra, dec, radius=radius)

    #ADM reshape this into an array with nvertsperpix columns
    pix_verts_in = np.reshape(verts_in,(npix,nvertsperpix))
    #ADM any row with a column not in the tiles must be a fractional pixel
    isfracpix = ~np.all(pix_verts_in,axis=1)

    #ADM the pixel integers where pixels are fractional
    return pix[np.where(isfracpix)]


def pixweight(nside, tiles=None, radius=None, precision=0.01, outfile=None, outplot=None):
    '''Create an array of the fraction of each pixel that overlaps the passed tiles.

    Args:
        nside (int): HEALPix `nside`, ``2**k`` where 0 < k < 30.
        tiles (Table-like, optional): Table-like with RA,DEC columns; or
            ``None`` to use all DESI tiles from :func:`desimodel.io.load_tiles`.
        radius (float, optional): Tile radius in degrees;
            if `None` use :func:`desimodel.focalplane.get_tile_radius_deg`.
        precision (float, optional): Approximate precision at which to
            calculate the area of pixels that partially overlap the footprint
            in SQUARE DEGREES (*e.g.* 0.01 means precise to
            0.01 sq. deg., or 36 sq. arcmin.). Lower numbers mean better precision.
        outfile (str, optional): Write the pixel->weight array to the file
            passed as `outfile` (could be full directory path + file).
        outplot (str, optional): Create a plot named `outplot`
           (pass a *name* for a plot in the current directory, a *full path*
           for a plot in a different directory). This is passed to
           matplotlib.pyplot's savefig routine.

    Returns pixweight:
        An array of the weight for each pixel at the passed nside. The
        weight is the fracion of the pixel that overlaps the passed tiles:
        `WEIGHT=1` for the pixel is entirely contained in the tiles;
        `WEIGHT=0` for the pixel is entirely outside of the tiles;
        `0 < WEIGHT < 1` for a pixel that overlaps the tiles.
        The index of the array is the HEALPixel integer.

    Notes:
        It is sufficient to create the weights at a suitably high nside, say
        nside=256 (0.052456 sq. deg. per pixel) as pixel numbers at
        lower nsides can be obtained by integer division by powers of 4, *e.g.*
        pix_@_nside_128 = pix@nside_256//4 and fractional weights at lower
        nsides are the mean of the 4 pixels at the higher nside
        :func:`desimodel.io.load_pixweight` can downsample the array to lower nsides.
    '''
    t0 = time()

    # ADM if tiles or radius is None, load the DESI model defaults.
    from .focalplane import get_tile_radius_deg
    if tiles is None:
        tiles = load_tiles()

    if radius is None:
        radius = get_tile_radius_deg()

    #ADM create an array that is zero for each integer pixel at this nside
    import healpy as hp
    npix = hp.nside2npix(nside)
    weight = np.zeros(npix,float)

    #ADM recover pixels that are likely to be in the DESI footprint and
    #ADM set their weight to one (it's the case, then, that anything that
    #ADM is *definitely outside of* the footprint has a weight of zero)

    pix = tiles2pix(nside, tiles=tiles, radius=radius, fact=2**8)
    weight[pix] = 1.

    #ADM loop through to find the "edge" (fractional) pixels, until convergence
    log.info('Start integration around partial pixels...')
    setfracpix = set([-1])
    #ADM only have a limited range, to prevent this running forever
    for i in range(20):
        log.info('Trying {} pixel boundary points (step={})...t={:.1f}s'
                 .format(4*2**i,2**i,time()-t0))
        #ADM find the fractional pixels at this step
        fracpix = tiles2fracpix(nside, step=2**i, tiles=tiles, radius=radius,
                                fact=2**8)
        log.info('...found {} fractional pixels...t={:.1f}s'
                 .format(len(fracpix),time()-t0))
        if set(fracpix) == setfracpix:
            break
        #ADM if we didn't converge, loop through again with the new
        #ADM set of fractional pixels
        setfracpix = set(fracpix)

    #ADM warn the user if the integration didn't converge at 4*2**20 boundary points
    if i == 20:
        log.warning('Integration around pixel boundaries did NOT converge!')

    #ADM create a mask that is True for fractional pixels, false for all other pixels
    mask = np.zeros(npix,bool)
    mask[fracpix] = True

    #ADM find the minimum and maximum dec of interest (there's no need to Monte Carlo
    #ADM integrate over declinations that lie beyond the fractional pixels)
    xyzverts = hp.boundaries(nside,fracpix,nest=True)
    theta, phi = hp.vec2ang(np.hstack(xyzverts).T)
    ra, dec = np.degrees(phi), 90-np.degrees(theta)
    decmin, decmax = np.min(dec), np.max(dec)
    sindecmin, sindecmax = np.sin(np.radians(decmin)), np.sin(np.radians(decmax))
    area = 360.*np.degrees(sindecmax-sindecmin)
    log.info('Populating randoms between {:.2f} and {:.2f} degrees, an area of {:.1f} sq. deg....t={:.1f}s'
             .format(decmin,decmax,area,time()-t0))

    #ADM determine the required precision for the area of interest
    nptpersqdeg = int((1./precision)**2)
    npt = int(nptpersqdeg * area)
    log.info('Generating {} random points...t={:.1f}s'.format(npt,time()-t0))

    #ADM loop over chunks (if npt > 1e7) to reach npt points while avoiding memory issues
    nchunk = int(1e7)
    pixinmask = []
    rainmask = []
    decinmask = []
    cnt = 0
    while cnt < npt:
        #ADM if a chunk would pass too many points (> npt), revert to the remaining number
        #ADM of points instead of creating a full chunk
        if nchunk + cnt > npt:
            nchunk = npt - cnt
        #ADM populate the portion of the sphere of interest with random points
        ra = np.random.uniform(0.,360.,nchunk)
        dec = np.degrees(np.arcsin(1.-np.random.uniform(1-sindecmax,1-sindecmin,nchunk)))
        #ADM convert the random points to pixel number
        pix = radec2pix(nside,ra,dec)
        #ADM retain random points for which the mask is True (i.e. just the fractional pixels)
        inmask = np.where(mask[pix])[0]
        decinmask.append(dec[inmask])
        rainmask.append(ra[inmask])
        pixinmask.append(pix[inmask])
        cnt += nchunk
        log.info('...generated {} random points...t={:.1f}s'
                 .format(cnt,time()-t0))

    #ADM collapse the 2-D chunks into a 1-D array
    from itertools import chain
    rainmask = np.array(list(chain.from_iterable(rainmask)))
    decinmask = np.array(list(chain.from_iterable(decinmask)))
    pixinmask = np.array(list(chain.from_iterable(pixinmask)))

    log.info('{} of the random points are in fractional pixels...t={:.1f}s'
             .format(len(pixinmask),time()-t0))

    #ADM find which random points in the fractional pixels are in the DESI footprint
    log.info('Start integration over fractional pixels at edges of DESI footprint...')
    indesi = is_point_in_desi(tiles,rainmask,decinmask)
    log.info('...{} of the random points in fractional pixels are in DESI...t={:.1f}s'
             .format(np.sum(indesi),time()-t0))

    #ADM assign the weights of the fractional pixels as the fraction of random points
    #ADM in the fractional pixels that are in the DESI footprint
    allinfracpix = np.histogram(pixinmask,bins=np.arange(npix))[0][fracpix]
    desiinfracpix = np.histogram(pixinmask[np.where(indesi)],bins=np.arange(npix))[0][fracpix]
    #ADM guard against integer division (for backwards-compatability with Python2)
    #ADM and create the final array of weights
    weight[fracpix] = desiinfracpix.astype('float64')/allinfracpix

    if outfile is not None:
        #ADM write information indicating HEALPix setup to file header
        #ADM include desimodel version as a check in case footprint changes
        import fitsio
        from desiutil import depend

        hdr = fitsio.FITSHDR()
        depend.setdep(hdr, 'desimodel', desimodel_version)
        hdr['PRECISE'] = precision
        hdr['HPXNSIDE'] = nside
        hdr['HPXNEST'] = True

        fitsio.write(outfile, weight, extname='PIXWEIGHTS', header=hdr, clobber=True)

    #ADM if outplot was passed, make a plot of the final mask in Mollweide projection
    if outplot is not None:
        import matplotlib.pyplot as plt
        hp.mollview(weight, nest=True)
        plt.savefig(outplot)

    log.info('Done...t={:.1f}s'.format(time()-t0))

    return weight


def pix2tiles(nside, pixels, tiles=None, radius=None):
    '''Returns subset of tiles that overlap the list of pixels.

    Args:
        nside (int): HEALPix `nside`, ``2**k`` where 0 < k < 30.
        pixels (array-like): Array of integer pixels using nested numbering scheme.
        tiles (Table-like, optional): Table-like with RA,DEC columns; or
            ``None`` to use all DESI tiles from :func:`desimodel.io.load_tiles`.
        radius (float, optional): Tile radius in degrees;
            if `None` use :func:`desimodel.focalplane.get_tile_radius_deg`.

    Returns:
        Table of tiles that cover these pixels.

    TODO: add support for tiles as integers or list/array of integer TILEIDs.
    '''
    import healpy as hp
    from .focalplane import get_tile_radius_deg

    if tiles is None:
        tiles = load_tiles()

    if radius is None:
        radius = get_tile_radius_deg()

    #- Trim tiles to ones that *might* overlap these pixels
    theta, phi = hp.pix2ang(nside, pixels, nest=True)
    ra, dec = np.degrees(phi), 90 - np.degrees(theta)
    pixsize = np.degrees(hp.nside2resol(nside))
    ii = find_tiles_over_point(tiles, ra, dec, radius=radius+pixsize)
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
        if np.any(np.isin(pixels, tilepix)):
            ii.append(i)
    return tiles[ii]


def _embed_sphere(ra, dec):
    """Embed `ra`, `dec` to a uniform sphere in three dimensions.
    """
    phi = np.radians(np.asarray(ra))
    theta = np.radians(90.0 - np.asarray(dec))
    r = np.sin(theta)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    z = np.cos(theta)
    return np.array((x, y, z)).T


def is_point_in_desi(tiles, ra, dec, radius=None, return_tile_index=False, workers = -1):
    """If a point (`ra`, `dec`) is within `radius` distance from center of any
    tile, it is in DESI.

    Args:
        tiles (Table-like): The output of :func:`desimodel.io.load_tiles`, or
            a similar Table.
        ra (scalar or array-like): Right Ascension in degrees.
        dec (scalar or array-like): Declination in degrees.  The size of `dec`
            must match the size of `ra`.
        radius (float, optional): Tile radius in degrees;
            if `None` use :func:`desimodel.focalplane.get_tile_radius_deg`.
        return_tile_index (bool, optional): If ``True``, return the index of
            the nearest tile in tiles array.
        workers (int, optional): Number of workers to pass to KDTree search. 
            Defaults to -1 i.e. all available threads.

    Returns:
        Return ``True`` if points given by `ra`, `dec` lie in the set of `tiles`.

    Notes:
        This function is optimized to query a lot of points.
    """
    from scipy.spatial import cKDTree as KDTree
    from .focalplane import get_tile_radius_deg

    if radius is None:
        radius = get_tile_radius_deg()

    tilecenters = _embed_sphere(tiles['RA'], tiles['DEC'])
    tree = KDTree(tilecenters)
    # radius to 3d distance
    threshold = 2.0 * np.sin(np.radians(radius) * 0.5)
    xyz = _embed_sphere(ra, dec)
    if not xyz.flags['C_CONTIGUOUS']:
        xyz = xyz.copy()

    d, i = tree.query(xyz, k=1, distance_upper_bound = threshold, workers = workers)

    #indesi = d < threshold
    indesi = np.isfinite(d)
    if return_tile_index:
        return indesi, i
    else:
        return indesi


def find_tiles_over_point(tiles, ra, dec, radius=None, workers = -1):
    """Return a list of indices of tiles that covers the points.

    This function is optimized to query a lot of points.
    radius is in units of degrees. The return value is an array
    of list objects that are the indices of tiles that cover each point.

    The indices are not sorted in any particular order.

    if ra, dec are scalars, a single list is returned.

    default radius is from desimodel.focalplane.get_tile_radius_deg()
    """
    from scipy.spatial import cKDTree as KDTree
    from .focalplane import get_tile_radius_deg

    if radius is None:
        radius = get_tile_radius_deg()

    tilecenters = _embed_sphere(tiles['RA'], tiles['DEC'])
    tree = KDTree(tilecenters)

    # radius to 3d distance
    threshold = 2.0 * np.sin(np.radians(radius) * 0.5)
    xyz = _embed_sphere(ra, dec)
    if not xyz.flags['C_CONTIGUOUS']:
        xyz = xyz.copy()

    indices = tree.query_ball_point(xyz, threshold, workers = workers)
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
    return find_points_radec(tiles['RA'], tiles['DEC'], ra, dec, radius)


def find_points_radec(telra, teldec, ra, dec, radius = None, workers = -1):
    """Return a list of indices of points that are within a radius of an arbitrary telra, teldec.

    This function is optimized to query a lot of points with a single telra and teldec.

    radius is in units of degrees. The return value is a list
    that contains the index of points that are in each tile.
    The indices are not sorted in any particular order.

    if tiles is a scalar, a single list is returned.

    default radius is from desimodel.focalplane.get_tile_radius_deg()

    Note: This is simply a modified version of find_points_in_tiles, but this function does not know about tiles.
    """
    from scipy.spatial import cKDTree as KDTree
    from .focalplane import get_tile_radius_deg

    if radius is None:
        radius = get_tile_radius_deg()

    # check for malformed input shapes. Sorry we currently only
    # deal with vector inputs. (for a sensible definition of indices)

    assert ra.ndim == 1
    assert dec.ndim == 1

    points = _embed_sphere(ra, dec)
    tree = KDTree(points)

    # radius to 3d distance
    threshold = 2.0 * np.sin(np.radians(radius) * 0.5)
    xyz = _embed_sphere(telra, teldec)
    if not xyz.flags['C_CONTIGUOUS']:
        xyz = xyz.copy()

    indices = tree.query_ball_point(xyz, threshold, workers = workers)
    return indices


def get_tile_radec(tileid):
    """Get the coordinates of a tile.

    Args:
        tileid (int): ID of a tile.

    Returns:
        tuple: (ra, dec) in degrees for the requested `tileid`.

    Raises:
        ValueError: If tileid is not in list of known tiles.
    """
    tiles = load_tiles()
    if tileid in tiles['TILEID']:
        i = np.where(tiles['TILEID'] == tileid)[0][0]
        return tiles[i]['RA'], tiles[i]['DEC']
    else:
        raise ValueError('Unknown tileid {}'.format(tileid))
