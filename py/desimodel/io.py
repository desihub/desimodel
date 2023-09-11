# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.io
============

I/O utility functions for files in desimodel.
"""
import os
import re
import warnings
from datetime import datetime, timezone

import yaml
import gzip
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column

from desiutil.log import get_logger

log = get_logger()


_thru = dict()

# ADM raise a custom exception when an environment variable is missing.
class MissingEnvVar(Exception):
    pass


def load_throughput(channel):
    """Returns specter Throughput object for the given channel 'b', 'r', or 'z'.

    Parameters
    ----------
    channel : {'b', 'r', 'z'}
        Spectrograph channel.

    Returns
    -------
    Throughput
        A specter throughput object.
    """
    import specter.throughput

    channel = channel.lower()
    global _thru
    if channel not in _thru:
        thrufile = findfile("throughput/thru-{0}.fits".format(channel))
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

    Returns
    -------
    PSF
        A specter PSF object.
    """
    import specter.psf

    channel = channel.lower()
    global _psf
    if channel not in _psf:
        psffile = findfile("specpsf/psf-{0}.fits".format(channel))
        _psf[channel] = specter.psf.load_psf(psffile)
    return _psf[channel]


#
#
#
_params = None


def load_desiparams():
    """Returns DESI parameter dictionary loaded from ``$DESIMODEL/data/desi.yaml``.

    Returns
    -------
    :class:`dict`
        The parameters read from the YAML file.
    """
    global _params
    if _params is None:
        desiparamsfile = findfile("desi.yaml")
        with open(desiparamsfile) as par:
            _params = yaml.safe_load(par)

        # - for temporary backwards compability after 'exptime' -> 'exptime_dark'
        if ("exptime" not in _params) and ("exptime_dark" in _params):
            _params["exptime"] = _params["exptime_dark"]

        # - Augment params with wavelength coverage from specpsf files
        # - wavemin/max = min/max wavelength covered by *any* fiber on the CCD
        # - wavemin/max_all = min/max wavelength covered by *all* fibers
        for channel in ["b", "r", "z"]:
            hdr = fits.getheader(findfile("specpsf/psf-{}.fits".format(channel)), 0)
            _params["ccd"][channel]["wavemin"] = hdr["WAVEMIN"]
            _params["ccd"][channel]["wavemax"] = hdr["WAVEMAX"]
            _params["ccd"][channel]["wavemin_all"] = hdr["WMIN_ALL"]
            _params["ccd"][channel]["wavemax_all"] = hdr["WMAX_ALL"]

    return _params


#
#
#
# Added and still needs to be committed and pushed to desihub
_gfa = None


def load_gfa():
    """Returns GFA table from ``$DESIMODEL/data/focalplane/gfa.ecsv``.

    Returns
    -------
    :class:`~astropy.table.Table`
        The data from the ECSV file.
    """
    global _gfa
    if _gfa is None:
        gfaFile = findfile("focalplane/gfa.ecsv")
        _gfa = Table.read(gfaFile, format="ascii.ecsv")
    return _gfa


#
#
#
_deviceloc = None


def load_deviceloc():
    """Returns a table from ``$DESIMODEL/data/focalplane/fiberpos-all.fits``.

    Returns
    -------
    :class:`~astropy.table.Table`
        The data from the FITS file, with columns converted to uppercase.
    """
    global _deviceloc
    if _deviceloc is None:
        fiberposfile = findfile("focalplane/fiberpos-all.fits")
        _deviceloc = Table.read(fiberposfile)

    # - Convert to upper case if needed
    # - Make copy of colnames b/c they are updated during iteration
    for col in list(_deviceloc.colnames):
        if col.islower():
            _deviceloc.rename_column(col, col.upper())

    return _deviceloc


_fiberpos = None


def load_fiberpos():
    """Returns fiberpos table from ``$DESIMODEL/data/focalplane/fiberpos.fits``.

    Returns
    -------
    :class:`~astropy.table.Table`
        The data from the FITS file, sorted by ``FIBER``.
    """
    global _fiberpos
    if _fiberpos is None:
        fiberposfile = findfile("focalplane/fiberpos.fits")
        _fiberpos = Table.read(fiberposfile)
        _fiberpos.sort("FIBER")
        # - Convert to upper case if needed
        # - Make copy of colnames b/c they are updated during iteration
        for col in list(_fiberpos.colnames):
            if col.islower():
                _fiberpos.rename_column(col, col.upper())

        # - Temporary backwards compatibility for renamed columns
        if "POSITIONER" in _fiberpos.colnames:
            warnings.warn(
                "old fiberpos.fits with POSITIONER column instead of LOCATION; please update your $DESIMODEL checkout",
                DeprecationWarning,
            )
            _fiberpos["LOCATION"] = _fiberpos["POSITIONER"]
        else:
            _fiberpos["POSITIONER"] = _fiberpos["LOCATION"]

        if "SPECTROGRAPH" in _fiberpos.colnames:
            warnings.warn(
                "old fiberpos.fits with SPECTROGRAPH column instead of SPECTRO; please update your $DESIMODEL checkout",
                DeprecationWarning,
            )
            _fiberpos["SPECTRO"] = _fiberpos["SPECTROGRAPH"]
        else:
            _fiberpos["SPECTROGRAPH"] = _fiberpos["SPECTRO"]

    return _fiberpos


#
#
#
_tiles = dict()


def load_tiles(onlydesi=True, extra=False, tilesfile=None, cache=True, programs=None):
    """Return DESI tiles structure from ``$DESI_SURVEYOPS/trunk/ops/tiles-main.ecsv``.

    Parameters
    ----------
    onlydesi : :class:`bool`, optional
        If ``True``, trim to just the tiles in the DESI footprint.
    extra : :class:`bool`, optional
        If ``True``, include extra layers with ``PROGRAM='EXTRA'``.
    tilesfile : :class:`str`, optional
        Name of tiles file to load; or None for default.  See Notes for
        details.
    cache : :class:`bool`, optional
        If ``False``, force reload of data from tiles file, instead of
        using cached values.
    programs : :class:`list` or `str`, optional
        Pass a list of program names to restrict to only those programs,
        e.g. ["DARK", "BACKUP"].

    Returns
    -------
    :class:`~astropy.io.fits.FITS_rec`
        The data table portion of the FITS file.

    Raises
    ------
    :exc:`FileNotFoundError`
        If the value of `tilesfile` does not exist.

    Notes
    -----
    Keyword-based environment variable expansion is performed on the `tilesfile`
    value, so *e.g.*::

        tiles = load_tiles(tilesfile='{HOME}/my-tiles.fits')

    will be expanded with the value of :envvar:`HOME`.

    If the parameter `tilesfile` is set, this function uses the following
    search method:

    0. Paths corresponding to both $DESI_SURVEYOPS/trunk/ops and
       $DESI_SURVEYOPS/ops are always both checked, to cover different
       svn checkout approaches.
    1. If the value includes an explicit path, even ``./``, use that file.
    2. If the value does *not* include an explicit path, *and* the file 
       name is identical to a file in ``$DESI_SURVEYOPS/trunk/ops/``, use 
       the file in ``$DESI_SURVEYOPS/trunk/ops/`` and issue a warning.
    3. If no matching file can be found at all, raise an exception.
    """
    global _tiles

    if tilesfile is None:
        # Use the default
        tilesfile = findfile("tiles-main.ecsv", surveyops=True)
    else:
        # If full path isn't included, check local vs $DESI_SURVEYOPS/ops
        tilepath, filename = os.path.split(tilesfile)
        if tilepath == "":
            have_local = os.path.isfile(tilesfile)
            checkfile = findfile(tilesfile, surveyops=True)
            have_dmdata = os.path.isfile(checkfile)
            if have_dmdata:
                if have_local:
                    msg = (
                        "$DESI_SURVEYOPS/(trunk)/ops/{0} is shadowed by a local"
                        + " file. Choosing $DESI_SURVEYOPS file."
                        + ' Use tilesfile="./{0}" if you want the local copy'
                        + " instead."
                    ).format(tilesfile)
                    warnings.warn(msg)
                tilesfile = checkfile

            if not (have_local or have_dmdata):
                msg = (
                    'File "{}" does not exist locally or in '
                    + "$DESI_SURVEYOPS/(trunk)/ops/!"
                ).format(tilesfile)
                raise FileNotFoundError(msg)

    # - standarize path location
    tilesfile = os.path.abspath(tilesfile.format(**os.environ))
    log.debug("Loading tiles from %s", tilesfile)

    # ADM allow reading from either .fits or .ecsv files.
    # ADM guard against the possibility that the file is zipped.
    isfits = ".fits" in os.path.basename(tilesfile)

    if cache and tilesfile in _tiles:
        tiledata = _tiles[tilesfile]
    else:
        if isfits:
            with fits.open(tilesfile, memmap=False) as hdulist:
                tiledata = hdulist[1].data
            #
            # Temporary workaround for problem identified in
            # https://github.com/desihub/desimodel/issues/30
            #
            if any([c.bzero is not None for c in tiledata.columns]):
                foo = [_tiles[k].dtype for k in tiledata.dtype.names]

            # - Check for out-of-date tiles file
            if np.issubdtype(tiledata["OBSCONDITIONS"].dtype, np.unsignedinteger):
                warnings.warn(
                    "Old desi-tiles.fits with uint16 OBSCONDITIONS; please update your $DESIMODEL checkout.",
                    DeprecationWarning,
                )
        else:
            tiledata = Table.read(tilesfile)
        # - load cache for next time
        if cache:
            _tiles[tilesfile] = tiledata

    # - Filter to only the DESI footprint if requested
    subset = np.ones(len(tiledata), dtype=bool)
    if onlydesi:
        subset &= tiledata["IN_DESI"] > 0

    # - Filter out PROGRAM=EXTRA tiles if requested
    if not extra:
        subset &= ~np.char.startswith(tiledata["PROGRAM"], "EXTRA")

    # ADM filter to program names if requested.
    if programs is not None:
        # ADM guard against a single string being passed.
        programs = np.atleast_1d(programs)
        isprog = np.zeros(len(tiledata), dtype=bool)
        for program in programs:
            isprog |= tiledata["PROGRAM"] == program
        subset &= isprog

    if np.all(subset):
        return tiledata
    else:
        return tiledata[subset]


_platescale = None


def load_platescale():
    """Loads platescale.txt.

    Returns
    -------
    :class:`~numpy.recarray`
        The data table read from the file.

    Notes
    -----
    The returned object has these columns:

    radius
        Radius from center of focal plane [mm].

    theta
        Radial angle that has a centroid at this radius [deg].

    radial_platescale
        Meridional (radial) plate scale [um/arcsec].

    az_platescale:
        Sagittal (azimuthal) plate scale [um/arcsec].

    arclength:
        Unknown description.
    """
    global _platescale
    if _platescale is not None:
        return _platescale

    infile = findfile("focalplane/platescale.txt")
    columns = [
        ("radius", "f8"),
        ("theta", "f8"),
        ("radial_platescale", "f8"),
        ("az_platescale", "f8"),
        ("arclength", "f8"),
    ]
    try:
        _platescale = np.loadtxt(infile, usecols=[0, 1, 6, 7, 8], dtype=columns)
    except IndexError:
        # - no "arclength" column in this version of desimodel/data
        # - Get info from separate rzs file instead

        _platescale = np.loadtxt(infile, usecols=[0, 1, 6, 7, 7], dtype=columns)
        rzs = Table.read(findfile("focalplane/rzsn.txt"), format="ascii")

        from scipy.interpolate import interp1d
        from numpy.lib.recfunctions import append_fields

        arclength = interp1d(rzs["R"], rzs["S"], kind="quadratic")
        _platescale["arclength"] = arclength(_platescale["radius"])

    return _platescale


_focalplane = None


def load_focalplane(time=None):
    """Load the focalplane state that is valid for the given time.

    Parameters
    ----------
    time : :class:`~datetime.datetime`
        The time to query with explicit timezone. Default to current time
        (:meth:`~datetime.datetime.now`) with timezone UTC.

    Returns
    -------
    :class:`tuple`
        A tuple of (FP layout, exclusion polygons, state, time string).
        The FP layout is a Table.  The exclusion polygons are a dictionary
        indexed by names that are referenced in the state.  The state
        is a Table.  The time string is the resulting UTC ISO format
        time string for the creation date of the FP model.
    """
    if time is None:
        time = datetime.now(tz=timezone.utc)
    elif time.tzinfo is None:
        msg = "Requested focalplane time '{}' is not timezone-aware.  Assuming UTC.".format(time)
        log.warning(msg)
        time = time.replace(tzinfo=timezone.utc)

    # Convert requested time to UTC
    time = time.astimezone(tz=timezone.utc)

    global _focalplane
    if _focalplane is None:
        # First call, parse all data files.
        fpdir = os.path.join(datadir(), "focalplane")
        fppat = re.compile(r"^desi-focalplane_(.*)\.ecsv$")
        stpat = re.compile(r"^desi-state_(.*)\.ecsv$")
        expat = re.compile(r"^desi-exclusion_(.*)\.(?:yaml|json).*$")
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
                          the 3 required files (focalplane, state, exclusion)".format(
                        ts
                    )
                    raise RuntimeError(msg)
        # Now load the files for each time into our cached global variable.
        _focalplane = list()
        for ts in sorted(fpraw.keys()):
            dt = None
            file_dt = None
            try:
                file_dt = datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S%z")
                dt = file_dt
            except ValueError:
                # This is an old file with implicit UTC times (no offset)
                # Load it as-is and then set the time zone.
                file_dt = datetime.strptime(ts, "%Y-%m-%dT%H:%M:%S")
                dt = file_dt.replace(tzinfo=timezone.utc)
            _focalplane.append((
                dt, file_dt, {
                    "fp_file": fpraw[ts]["fp"],
                    "st_file": fpraw[ts]["st"],
                    "ex_file": fpraw[ts]["ex"],
                    "fp_data": None,
                    "st_data": None,
                    "ex_data": None,
                }
            ))

    # Search the list of states for the most recent time that is before our
    # requested time.  There should not be too many different states, or else
    # we are using the wrong format for storing these.  Therefore, a linear
    # search should be fast enough.
    focalplane_props = None
    file_tmstr = None
    for dt, file_dt, props in _focalplane:
        if time >= dt:
            focalplane_props = props
            try:
                file_tmstr = file_dt.isoformat(timespec="seconds")
            except TypeError:
                # This must be python < 3.6, with no timespec option.
                # Since the focalplane time is read from the file name without
                # microseconds, the microseconds should be zero and so the
                # default return string will be correct.
                file_tmstr = file_dt.isoformat()
        else:
            break

    if focalplane_props is None:
        msg = "Cannot find focalplane for time {}".format(time)
        raise RuntimeError(msg)

    # Load the data if this is the first time working with this focalplane
    if focalplane_props["fp_data"] is None:
        focalplane_props["fp_data"] = Table.read(
            focalplane_props["fp_file"],
            format="ascii.ecsv"
        )
        focalplane_props["st_data"] = Table.read(
            focalplane_props["st_file"],
            format="ascii.ecsv"
        )
        if (focalplane_props['ex_file'].endswith('.json') or
            focalplane_props['ex_file'].endswith('.json.gz')):
            import json
            loadroutine = json.load
        else:
            loadroutine = yaml.safe_load
        try:
            # First try to load uncompressed
            with open(focalplane_props["ex_file"], "r") as f:
                focalplane_props["ex_data"] = loadroutine(f)
        except:
            # Must be gzipped
            with gzip.open(focalplane_props["ex_file"], "rb") as f:
                focalplane_props["ex_data"] = loadroutine(f)

    # Now "replay" the state up to our requested time.
    st_data = focalplane_props["st_data"]
    locstate = dict()
    for row in range(len(st_data)):
        try:
            tm = datetime.strptime(st_data[row]["TIME"], "%Y-%m-%dT%H:%M:%S%z")
        except ValueError:
            # Old format with implicit UTC timezone
            tm = datetime.strptime(st_data[row]["TIME"], "%Y-%m-%dT%H:%M:%S")
            tm = tm.replace(tzinfo=timezone.utc)
        if tm <= time:
            loc = st_data[row]["LOCATION"]
            locstate[loc] = st_data[row]

    rows = list()
    for loc in sorted(locstate.keys()):
        rows.append(locstate[loc])
    state_data = Table(rows=rows, names=st_data.colnames)
    state_data.remove_column("TIME")

    return (
        focalplane_props["fp_data"],
        focalplane_props["ex_data"],
        state_data,
        file_tmstr
    )


def reset_cache():
    """Reset I/O cache."""
    global _thru, _psf, _params, _gfa, _fiberpos, _tiles, _platescale, _focalplane
    _thru = dict()
    _psf = dict()
    _params = None
    _gfa = None
    _fiberpos = None
    _tiles = dict()
    _platescale = None
    _focalplane = None


def load_target_info():
    """Loads data/targets/targets.yaml and returns the nested dictionary.

    This is primarily syntactic sugar to avoid end users constructing
    paths and filenames by hand (which *e.g.* broke when targets.dat was
    renamed to targets.yaml).

    Returns
    -------
    :class:`dict`
        The dictionary read from the YAML file.
    """
    targetsfile = findfile("targets/targets.yaml")
    if not os.path.exists(targetsfile):
        targetsfile = findfile("targets/targets.dat")

    with open(targetsfile) as fx:
        data = yaml.safe_load(fx)

    return data


def load_pixweight(nside, pixmap=None):
    """Loads ``$DESIMODEL/data/footprint/desi-healpix-weights.fits``.

    Parameters
    ----------
    nside : :class:`int`
        After loading, the array will be resampled to the passed HEALPix `nside`.
    pixmap : :class:`~astropy.io.fits.FITS_rec`, optional
        Input pixel weight map, already read from a weights file.

    Returns
    -------
    Weight
        HEALPix weight map for the DESI footprint at the requested `nside`.
    """
    import healpy as hp

    if pixmap is not None:
        log.debug("Using input pixel weight map of length {}.".format(len(pixmap)))
    else:
        # ADM read in the standard pixel weights file
        pixfile = findfile("footprint/desi-healpix-weights.fits")
        with fits.open(pixfile) as hdulist:
            pixmap = hdulist[0].data

    # ADM determine the file's nside, and flag a warning if the passed nside exceeds it
    npix = len(pixmap)
    truenside = hp.npix2nside(len(pixmap))
    if truenside < nside:
        log.warning(
            "downsampling is fuzzy...Passed nside={}, but file {} is stored at nside={}".format(
                nside, pixfile, truenside
            )
        )

    # ADM resample the map
    return hp.pixelfunc.ud_grade(pixmap, nside, order_in="NESTED", order_out="NESTED")


def findfile(filename, surveyops=False):
    """Return full path to data file ``$DESIMODEL/data/filename``.

    Parameters
    ----------
    filename : :class:`str`
        Name of the file, relative to the desimodel data directory.

    surveyops : :class:`bool`
        If ``True`` then find the relevant path for the $DESI_SURVEYOPS
        directory rather than the $DESIMODEL directory.

    Returns
    -------
    :class:`str`
        The full path.

    Notes
    -----
    This is a precursor for a potential future refactor where
    desimodel data would be installed with the package and :envvar:`DESIMODEL`
    would become an optional override.
    """
    return os.path.join(datadir(surveyops), filename)


def datadir(surveyops=False):
    """Returns location to desimodel data.

    Parameters
    ----------
    surveyops : :class:`bool`
        If ``True`` then find the relevant path for the $DESI_SURVEYOPS
        directory rather than the $DESIMODEL directory.

    If set, :envvar:`DESIMODEL` overrides data installed with the package.
    """
    if surveyops:
        if "DESI_SURVEYOPS" in os.environ:
            surveyops = os.environ["DESI_SURVEYOPS"]
            # ADM test whether surveyops directory was checked out to trunk.
            if os.path.isdir(os.path.join(surveyops, "trunk", "ops")):
                surveyops = os.path.join(surveyops, "trunk")
            return os.path.abspath(os.path.join(surveyops, "ops"))
        # ADM raise a custom exception if $DESI_SURVEYOPS is not set.
        else:
            raise MissingEnvVar(f"$DESI_SURVEYOPS is not set")
    else:
        if "DESIMODEL" in os.environ:
            return os.path.abspath(os.path.join(os.environ["DESIMODEL"], "data"))
        else:
            import pkg_resources

            return pkg_resources.resource_filename("desimodel", "data")
