# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.io
============

I/O utility functions for files in desimodel.
"""
import os
import sys
import re
import warnings
from datetime import datetime

import yaml
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column

from desiutil.log import get_logger
log = get_logger()


_thru = dict()
def load_throughput(channel):
    """Returns specter Throughput object for the given channel 'b', 'r', or 'z'.

    Parameters
    ----------
    channel : {'b', 'r', 'z'}
        Spectrograph channel.
    """
    import specter.throughput
    channel = channel.lower()
    global _thru
    if channel not in _thru:
        thrufile = os.path.join(os.environ['DESIMODEL'],'data','throughput','thru-{0}.fits'.format(channel))
        _thru[channel] = specter.throughput.load_throughput(thrufile)
    return _thru[channel]
#
#
#
_psf = dict()
def load_psf(channel):
    """Returns specter PSF object for the given channel 'b', 'r', or 'z'.

    Parameters
    ----------
    channel : {'b', 'r', 'z'}
        Spectrograph channel.
    """
    import specter.psf
    channel = channel.lower()
    global _psf
    if channel not in _psf:
        psffile = os.path.join(os.environ['DESIMODEL'],'data','specpsf','psf-{0}.fits'.format(channel))
        _psf[channel] = specter.psf.load_psf(psffile)
    return _psf[channel]
#
#
#
_params = None
def load_desiparams():
    """Returns DESI parameter dictionary loaded from desimodel/data/desi.yaml.
    """
    global _params
    if _params is None:
        desiparamsfile = os.path.join(os.environ['DESIMODEL'],'data','desi.yaml')
        with open(desiparamsfile) as par:
            _params = yaml.safe_load(par)

        #- for temporary backwards compability after 'exptime' -> 'exptime_dark'
        if ('exptime' not in _params) and ('exptime_dark' in _params):
            _params['exptime'] = _params['exptime_dark']

        #- Augment params with wavelength coverage from specpsf files
        #- wavemin/max = min/max wavelength covered by *any* fiber on the CCD
        #- wavemin/max_all = min/max wavelength covered by *all* fibers
        for channel in ['b', 'r', 'z']:
            hdr = fits.getheader(findfile('specpsf/psf-{}.fits'.format(channel)), 0)
            _params['ccd'][channel]['wavemin'] = hdr['WAVEMIN']
            _params['ccd'][channel]['wavemax'] = hdr['WAVEMAX']
            _params['ccd'][channel]['wavemin_all'] = hdr['WMIN_ALL']
            _params['ccd'][channel]['wavemax_all'] = hdr['WMAX_ALL']

    return _params
#
#
#
# Added and still needs to be committed and pushed to desihub
_gfa = None
def load_gfa():
    """Returns GFA table from desimodel/data/focalplane/gfa.ecsv"""
    global _gfa
    from astropy.table import Table
    # os is imported already in the desimodel io.py
    import os
    if _gfa is None:
        gfaFile = os.path.join(os.environ['DESIMODEL'], 'data', 'focalplane', 'gfa.ecsv')
        _gfa = Table.read(gfaFile, format = 'ascii.ecsv')
    return _gfa
#
#
#
_deviceloc = None
def load_deviceloc():
    global _deviceloc
    from astropy.table import Table
    if _deviceloc is None:
        fiberposfile = os.path.join(os.environ['DESIMODEL'],'data','focalplane','fiberpos-all.fits')
        _deviceloc = Table.read(fiberposfile)

    #- Convert to upper case if needed
    #- Make copy of colnames b/c they are updated during iteration
    for col in list(_deviceloc.colnames):
        if col.islower():
            _deviceloc.rename_column(col, col.upper())

    return _deviceloc

_fiberpos = None
def load_fiberpos():
    """Returns fiberpos table from desimodel/data/focalplane/fiberpos.fits.
    """
    global _fiberpos
    from astropy.table import Table
    if _fiberpos is None:
        fiberposfile = os.path.join(os.environ['DESIMODEL'],'data','focalplane','fiberpos.fits')
        _fiberpos = Table.read(fiberposfile)
        _fiberpos.sort('FIBER')
        #- Convert to upper case if needed
        #- Make copy of colnames b/c they are updated during iteration
        for col in list(_fiberpos.colnames):
            if col.islower():
                _fiberpos.rename_column(col, col.upper())

        #- Temporary backwards compatibility for renamed columns
        if 'POSITIONER' in _fiberpos.colnames:
            warnings.warn('old fiberpos.fits with POSITIONER column instead of LOCATION; please update your $DESIMODEL checkout', DeprecationWarning)
            _fiberpos['LOCATION'] = _fiberpos['POSITIONER']
        else:
            _fiberpos['POSITIONER'] = _fiberpos['LOCATION']


        if 'SPECTROGRAPH' in _fiberpos.colnames:
            warnings.warn('old fiberpos.fits with SPECTROGRAPH column instead of SPECTRO; please update your $DESIMODEL checkout', DeprecationWarning)
            _fiberpos['SPECTRO'] = _fiberpos['SPECTROGRAPH']
        else:
            _fiberpos['SPECTROGRAPH'] = _fiberpos['SPECTRO']

    return _fiberpos
#
#
#
_tiles = dict()
def load_tiles(onlydesi=True, extra=False, tilesfile=None, cache=True):
    """Return DESI tiles structure from desimodel/data/footprint/desi-tiles.fits.

    Parameters
    ----------
    onlydesi : :class:`bool` (default True)
        If ``True``, trim to just the tiles in the DESI footprint.
    extra : :class:`bool`, (default False)
        If ``True``, include extra layers with PROGRAM='EXTRA'.
    tilesfile : (str)
        Name of tiles file to load; or None for default.
        Without path, look in $DESIMODEL/data/footprint, otherwise load file.
    cache : :class:`bool`, (default True)
        Use cache of tiles data.
    """
    global _tiles

    if tilesfile is None:
        # Use the default
        tilesfile = os.path.join(
            os.environ['DESIMODEL'], 'data', 'footprint', 'desi-tiles.fits')
    else:
        # If full path isn't included, check local vs $DESIMODEL/data/footprint
        tilepath, filename = os.path.split(tilesfile)
        if tilepath == '':
            have_local = os.path.isfile(tilesfile)
            checkfile = os.path.join(os.environ['DESIMODEL'],
                                     'data', 'footprint', tilesfile)
            have_dmdata = os.path.isfile(checkfile)
            if have_dmdata:
                if have_local:
                    msg = '$DESIMODEL/data/footprint/{} is shadowed by a local'\
                          ' file. Choosing $DESIMODEL file.'\
                          ' Use tilesfile="./{}" if you want the local copy'\
                          ' instead'.format(tilesfile, tilesfile)
                    warnings.warn(msg)
                tilesfile = checkfile

            if not (have_local or have_dmdata):
                msg = 'File "{}" does not exist locally or in '\
                      '$DESIMODEL/data/footprint/'.format(tilesfile)
                if sys.version_info.major == 2:
                    raise IOError(msg)
                else:
                    raise FileNotFoundError(msg)

    #- standarize path location
    tilesfile = os.path.abspath(tilesfile.format(**os.environ))
    log.debug('Loading tiles from %s', tilesfile)

    if cache and tilesfile in _tiles:
        tiledata = _tiles[tilesfile]
    else:
        with fits.open(tilesfile, memmap=False) as hdulist:
            tiledata = hdulist[1].data
        #
        # Temporary workaround for problem identified in
        # https://github.com/desihub/desimodel/issues/30
        #
        if any([c.bzero is not None for c in tiledata.columns]):
            foo = [_tiles[k].dtype for k in tiledata.dtype.names]

        #- Check for out-of-date tiles file
        if np.issubdtype(tiledata['OBSCONDITIONS'].dtype, np.unsignedinteger):
            warnings.warn('old desi-tiles.fits with uint16 OBSCONDITIONS; please update your $DESIMODEL checkout', DeprecationWarning)

        #- load cache for next time
        if cache:
            _tiles[tilesfile] = tiledata

    #- Filter to only the DESI footprint if requested
    subset = np.ones(len(tiledata), dtype=bool)
    if onlydesi:
        subset &= tiledata['IN_DESI'] > 0

    #- Filter out PROGRAM=EXTRA tiles if requested
    if not extra:
        subset &= ~np.char.startswith(tiledata['PROGRAM'], 'EXTRA')

    if np.all(subset):
        return tiledata
    else:
        return tiledata[subset]

_platescale = None
def load_platescale():
    '''
    Loads platescale.txt, returning structured array with columns

        radius: radius from center of focal plane [mm]
        theta: radial angle that has a centroid at this radius [deg]
        radial_platescale: Meridional (radial) plate scale [um/arcsec]
        az_platescale: Sagittal (azimuthal) plate scale [um/arcsec]
    '''
    global _platescale
    if _platescale is not None:
        return _platescale

    infile = findfile('focalplane/platescale.txt')
    columns = [
        ('radius', 'f8'),
        ('theta', 'f8'),
        ('radial_platescale', 'f8'),
        ('az_platescale', 'f8'),
        ('arclength', 'f8'),
    ]
    _platescale = np.loadtxt(infile, usecols=[0,1,6,7,8], dtype=columns)
    return _platescale


_focalplane = None
def load_focalplane(time=None):
    """Load the focalplane state that is valid for the given time.

    Options:
        time (datetime):  The time to query. default to current time.

    Returns:
        (tuple):  The (FP layout, exclusion polygons, state, time string).
            The FP layout is a Table.  The exclusion polygons are a dictionary
            indexed by names that are referenced in the state.  The state
            is a Table.  The time string is the resulting UTC ISO format
            time string used by the lookup.

    """
    if time is None:
        time = datetime.now()

    global _focalplane
    if _focalplane is None:
        # First call, load all data files.
        fpdir = os.path.join(datadir(), "focalplane")
        fppat = re.compile(r"desi-focalplane_(.*)\.ecsv")
        stpat = re.compile(r"desi-state_(.*)\.ecsv")
        expat = re.compile(r"desi-exclusion_(.*)\.yaml")
        fpraw = dict()
        msg = "Loading focalplanes from {}".format(fpdir)
        log.debug(msg)
        for root, dirs, files in os.walk(fpdir):
            for f in files:
                fpmat = fppat.match(f)
                if fpmat is not None:
                    dt = fpmat.group(1)
                    if dt not in fpraw:
                        fpraw[dt] = dict()
                    fpraw[dt]["fp"] = os.path.join(root, f)
                    continue
                stmat = stpat.match(f)
                if stmat is not None:
                    dt = stmat.group(1)
                    if dt not in fpraw:
                        fpraw[dt] = dict()
                    fpraw[dt]["st"] = os.path.join(root, f)
                    continue
                exmat = expat.match(f)
                if exmat is not None:
                    dt = exmat.group(1)
                    if dt not in fpraw:
                        fpraw[dt] = dict()
                    fpraw[dt]["ex"] = os.path.join(root, f)
            break
        # Check that we have all 3 files needed for each timestamp
        for ts, files in fpraw.items():
            for key in ["fp", "st", "ex"]:
                if key not in files:
                    msg = "Focalplane state for time {} is missing one of \
                          the 3 required files (focalplane, state, exclusion)"\
                          .format(ts)
                    raise RuntimeError(msg)
        # Now load the files for each time into our cached global variable.
        _focalplane = list()
        for ts in sorted(fpraw.keys()):
            dt = datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S")
            fp = Table.read(fpraw[ts]["fp"], format="ascii.ecsv")
            st = Table.read(fpraw[ts]["st"], format="ascii.ecsv")
            ex = None
            with open(fpraw[ts]["ex"], "r") as f:
                ex = yaml.safe_load(f)
            _focalplane.append((dt, fp, ex, st))

    # Search the list of states for the most recent time that is before our
    # requested time.  There should not be too many different states, or else
    # we are using the wrong format for storing these.  Therefore, a linear
    # search should be fast enough.
    fp_data = None
    excl_data = None
    fullstate = None
    tmstr = None
    for dt, fp, ex, st in _focalplane:
        if time > dt:
            fp_data = fp
            excl_data = ex
            fullstate = st
            try:
                tmstr = dt.isoformat(timespec="seconds")
            except TypeError:
                # This must be python < 3.6, with no timespec option.
                # Since the focalplane time is read from the file name without
                # microseconds, the microseconds should be zero and so the
                # default return string will be correct.
                tmstr = dt.isoformat()
        else:
            break

    if fullstate is None:
        msg = "Cannot find focalplane for time {}".format(time)
        raise RuntimeError(msg)

    # Now "replay" the state up to our requested time.
    locstate = dict()
    for row in range(len(fullstate)):
        tm = datetime.strptime(fullstate[row]["TIME"], "%Y-%m-%dT%H:%M:%S")
        if tm <= time:
            loc = fullstate[row]["LOCATION"]
            pet = fullstate[row]["PETAL"]
            dev = fullstate[row]["DEVICE"]
            st = fullstate[row]["STATE"]
            excl = fullstate[row]["EXCLUSION"]
            if loc not in locstate:
                locstate[loc] = dict()
            locstate[loc]["PETAL"] = pet
            locstate[loc]["DEVICE"] = dev
            locstate[loc]["STATE"] = st
            locstate[loc]["EXCLUSION"] = excl

    nloc = len(locstate)
    state_cols = [
        Column(name="PETAL", length=nloc, dtype=np.int32,
               description="Petal location [0-9]"),
        Column(name="DEVICE", length=nloc, dtype=np.int32,
               description="Device location on the petal"),
        Column(name="LOCATION", length=nloc, dtype=np.int32,
               description="Global device location (PETAL * 1000 + DEVICE)"),
        Column(name="STATE", length=nloc, dtype=np.uint32,
               description="State bit field (good == 0)"),
        Column(name="EXCLUSION", length=nloc, dtype=np.dtype("a9"),
               description="The exclusion polygon for this device"),
    ]
    state_data = Table()
    state_data.add_columns(state_cols)
    row = 0
    for loc in sorted(locstate.keys()):
        state_data[row]["PETAL"] = locstate[loc]["PETAL"]
        state_data[row]["DEVICE"] = locstate[loc]["DEVICE"]
        state_data[row]["LOCATION"] = loc
        state_data[row]["STATE"] = locstate[loc]["STATE"]
        state_data[row]["EXCLUSION"] = locstate[loc]["EXCLUSION"]
        row += 1

    return (fp_data, excl_data, state_data, tmstr)


def reset_cache():
    '''Reset I/O cache'''
    global _thru, _psf, _params, _gfa, _fiberpos, _tiles, _platescale,\
        _focalplane
    _thru = dict()
    _psf = dict()
    _params = None
    _gfa = None
    _fiberpos = None
    _tiles = dict()
    _platescale = None
    _focalplane = None

def load_target_info():
    '''
    Loads data/targets/targets.yaml and returns the nested dictionary

    This is primarily syntactic sugar to avoid end users constructing
    paths and filenames by hand (which e.g. broke when targets.dat was
    renamed to targets.yaml)
    '''
    targetsfile = os.path.join(datadir(),'targets','targets.yaml')
    if not os.path.exists(targetsfile):
        targetsfile = os.path.join(datadir(),'targets','targets.dat')

    with open(targetsfile) as fx:
        data = yaml.safe_load(fx)

    return data

def load_pixweight(nside, pixmap=None):
    '''
    Loads desimodel/data/footprint/desi-healpix-weights.fits

    Args:
        nside: after loading, the array will be resampled to the
            passed HEALPix nside

    Options:
        pixmap: input pixel weight map (optional, defaults to None)

    Returns healpix weight map for the DESI footprint at the requested nside
    '''
    import healpy as hp

    if pixmap is not None:
        log.debug('Using input pixel weight map of length {}.'.format(len(pixmap)))
    else:
        #ADM read in the standard pixel weights file
        pixfile = os.path.join(os.environ['DESIMODEL'],'data','footprint','desi-healpix-weights.fits')
        with fits.open(pixfile) as hdulist:
            pixmap = hdulist[0].data

    #ADM determine the file's nside, and flag a warning if the passed nside exceeds it
    npix = len(pixmap)
    truenside = hp.npix2nside(len(pixmap))
    if truenside < nside:
        log.warning("downsampling is fuzzy...Passed nside={}, but file {} is stored at nside={}"
                  .format(nside,pixfile,truenside))

    #ADM resample the map
    return hp.pixelfunc.ud_grade(pixmap,nside,order_in='NESTED',order_out='NESTED')

def findfile(filename):
    '''
    Return full path to data file $DESIMODEL/data/filename

    Note: this is a precursor for a potential future refactor where
    desimodel data would be installed with the package and $DESIMODEL
    would become an optional override.
    '''
    return os.path.join(datadir(), filename)

def datadir():
    '''
    Returns location to desimodel data

    if set, $DESIMODEL overrides data installed with the package
    '''
    if 'DESIMODEL' in os.environ:
        return os.path.abspath(os.path.join(os.environ['DESIMODEL'], 'data'))
    else:
        import pkg_resources
        return pkg_resources.resource_filename('desimodel', 'data')
