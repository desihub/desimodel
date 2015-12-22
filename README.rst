=========
desimodel
=========

Introduction
------------

This product contains information about DESI hardware designs in machine
readable formats to be used by simulations.  It is intended to be correct
but non-authoritative.  If in question, the DocDB version of a design
parameter is correct.  This is intended to reflect the same information
while being more conveniently organized and formatted for simulations.

Data Files
----------

Files in data/inputs/ are copied from DocDB and are not intended to be used
directly.  These are parsed and reformatted to produce files in other data/
directories which are for use.

data/desi.yaml
    basic scalar parameters, organized in a nested tree

data/focalplane/
    informaton about positioner locations and platescale

data/specpsf/
    spectrograph point-spread-function (PSF) for specter
    CCD pixel-level simulations

data/throughput/
    throughput vs wavelength (also contained in specpsf)

data/footprint/
    DESI footprint with RA, dec of tiles

data/spectra/
    Example benchmark spectra

| Stephen Bailey
| Early 2014

License
-------

desimodel is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
