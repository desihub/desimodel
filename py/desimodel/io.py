# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.io
============

I/O utility functions for files in desimodel.
"""
import os
from astropy.io import fits
import yaml
import numpy as np
import warnings
#
#- PSF and throughput, which require specter
#
try:
    from specter.throughput import load_throughput as specter_load_throughput
except ImportError as e:
    warnings.warn(str(e))
    warnings.warn("Unable to import specter.throughput.load_throughput(); desimodel.io.load_throughput() won't work.")
try:
    from specter.psf import load_psf as specter_load_psf
except ImportError as e:
    warnings.warn(str(e))
    warnings.warn("Unable to import specter.psf.load_psf(); desimodel.io.load_psf() won't work.")
#
#
#
_thru = dict()
def load_throughput(channel):
    """Returns specter Throughput object for the given channel 'b', 'r', or 'z'.

    Parameters
    ----------
    channel : {'b', 'r', 'z'}
        Spectrograph channel.
    """
    channel = channel.lower()
    global _thru
    if channel not in _thru:
        thrufile = os.path.join(os.environ['DESIMODEL'],'data','throughput','thru-{0}.fits'.format(channel))
        _thru[channel] = specter_load_throughput(thrufile)
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
    channel = channel.lower()
    global _psf
    if channel not in _psf:
        psffile = os.path.join(os.environ['DESIMODEL'],'data','specpsf','psf-{0}.fits'.format(channel))
        _psf[channel] = specter_load_psf(psffile)
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
            _params = yaml.load(par)

    #- for temporary backwards compability after 'exptime' -> 'exptime_dark'
    if ('exptime' not in _params) and ('exptime_dark' in _params):
        _params['exptime'] = _params['exptime_dark']

    return _params
#
#
#
_fiberpos = None
def load_fiberpos():
    """Returns fiberpos table from desimodel/data/focalplane/fiberpos.fits.
    """
    global _fiberpos
    if _fiberpos is None:
        fiberposfile = os.path.join(os.environ['DESIMODEL'],'data','focalplane','fiberpos.fits')
        with fits.open(fiberposfile) as hdulist:
            _fiberpos = hdulist[1].data
        if 'FIBER' not in _fiberpos.dtype.names:
            #
            # File contains lower-case column names, but we want upper-case.
            #
            for i, key in enumerate(_fiberpos.dtype.names):
                _fiberpos.columns[i].name = key.upper()
    return _fiberpos
#
#
#
_tiles = None
def load_tiles(onlydesi=True, extra=False):
    """Return DESI tiles structure from desimodel/data/footprint/desi-tiles.fits

    Options
    -------
    onlydesi : :class:`bool` (default False)
        If ``True``, trim to just the tiles in the DESI footprint.
    extra : :class:`bool`, (default True)
        If ``True``, include extra layers with PROGRAM='EXTRA'.
    """
    global _tiles
    if _tiles is None:
        footprint = os.path.join(os.environ['DESIMODEL'],'data','footprint','desi-tiles.fits')
        with fits.open(footprint) as hdulist:
            _tiles = hdulist[1].data
        #
        # Temporary workaround for problem identified in
        # https://github.com/desihub/desimodel/issues/30
        #
        if any([c.bzero is not None for c in _tiles.columns]):
            foo = [_tiles[k].dtype for k in _tiles.dtype.names]

    #- Filter to only the DESI footprint if requested
    subset = np.ones(len(_tiles), dtype=bool)
    if onlydesi:
        subset &= _tiles['IN_DESI'] > 0

    #- Filter out PROGRAM=EXTRA tiles if requested
    if not extra:
        subset &= ~np.char.startswith(_tiles['PROGRAM'], 'EXTRA')

    if np.all(subset):
        return _tiles
    else:
        return _tiles[subset]
#
#
#
def get_tile_radec(tileid):
    """Return (ra, dec) in degrees for the requested tileid.

    If tileid is not in DESI, return (0.0, 0.0)
    TODO: should it raise and exception instead?
    """
    tiles = load_tiles()
    if tileid in tiles['TILEID']:
        i = np.where(tiles['TILEID'] == tileid)[0][0]
        return tiles[i]['RA'], tiles[i]['DEC']
    else:
        return (0.0, 0.0)

def findfile(filename):
    '''
    Return full path to data file $DESIMODEL/data/filename

    Note: this is a precursor for a potential future refactor where
    desimodel data would be installed with the package and $DESIMODEL
    would become an optional override.
    '''
    return os.path.join(os.getenv('DESIMODEL'), 'data', filename)
