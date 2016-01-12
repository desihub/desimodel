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

Desimodel Data
--------------

Finding the Data
~~~~~~~~~~~~~~~~

The data that accompanies the desimodel code is not stored with the code.
Due to its size, it is kept in the DESI svn repository.  What follows is
a description of these data.

Data Files
~~~~~~~~~~

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

Branches
~~~~~~~~

There are a couple of permanent branches that were used for testing
alternative designs.  These will never be merged into trunk but we
will keep them around for the record:

altccd
    500 um vs. 250 um thick CCDs

bb
    recreating assumptions used during early BigBOSS projections

In addition to these historical branches, there is a permanent 'testing' branch
that contains smaller versions of the desimodel files.  This branch is
intended for use in desimodel unit tests.

Full Documentation
------------------

Please visit `desimodel on Read the Docs`_

.. image:: https://readthedocs.org/projects/desimodel/badge/?version=latest
    :target: http://desimodel.readthedocs.org/en/latest/
    :alt: Documentation Status

.. _`desimodel on Read the Docs`: http://desimodel.readthedocs.org/en/latest/

Travis Build Status
-------------------

.. image:: https://img.shields.io/travis/desihub/desimodel.svg
    :target: https://travis-ci.org/desihub/desimodel
    :alt: Travis Build Status


Test Coverage Status
--------------------

.. image:: https://coveralls.io/repos/desihub/desimodel/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/desihub/desimodel?branch=master
    :alt: Test Coverage Status

License
-------

desimodel is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
