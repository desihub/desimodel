=========
desimodel
=========


|Actions Status| |Coveralls Status| |Documentation Status|

.. |Actions Status| image:: https://github.com/desihub/desimodel/workflows/CI/badge.svg
    :target: https://github.com/desihub/desimodel/actions
    :alt: GitHub Actions CI Status

.. |Coveralls Status| image:: https://coveralls.io/repos/desihub/desimodel/badge.svg
    :target: https://coveralls.io/github/desihub/desimodel
    :alt: Test Coverage Status

.. |Documentation Status| image:: https://readthedocs.org/projects/desimodel/badge/?version=latest
    :target: https://desimodel.readthedocs.io/en/latest/
    :alt: Documentation Status


Introduction
------------

This product contains information about DESI hardware designs in machine
readable formats to be used by simulations.  It is intended to be correct
but non-authoritative.  If in question, the DocDB version of a design
parameter is correct.  This is intended to reflect the same information
while being more conveniently organized and formatted for simulations.

**PLEASE KEEP THIS FILE IN SYNC WITH THE EQUIVALENT FILE IN SVN.**

Desimodel Code
--------------

The svn product described below contains only the *data* associated with
desimodel. The code is in Github: https://github.com/desihub/desimodel.

Desimodel Data
--------------

Getting the Data
~~~~~~~~~~~~~~~~

The data that accompanies the desimodel code is not stored with the code.
Due to its size, it is kept in the DESI svn repository.  The public, read-only
URL for svn access is https://desi.lbl.gov/svn/code/desimodel, with the usual
trunk/, tags/ and branches/ directories.

Once you have installed this package, using either ``pip`` or ``desiInstall``,
the associated data must also be installed.  **For most every case, you
should install the tag in svn that corresponds to the same tag in git.**

Installing this package will create the command-line script
``install_desimodel_data``.  It should appear in your ``PATH`` assuming
a successful install.  ``install_desimodel_data --help`` will show you
how to use this script.

You can also call the function ``desimodel.install.install()`` from
inside other Python code.

When desimodel is used with other downstream code, it might become
necessary to set the environment variable ``DESIMODEL``.
You should set the ``DESIMODEL`` environment variable to point to the directory
containing the data/ directory.

``install_desimodel_data`` with default options, may install the data in
a relatively obscure location. To find the value to set for ``DESIMODEL``,
you can do the following::

    python -c "from importlib.resources import files; print(str(files('desimodel')))"


Data Files
~~~~~~~~~~

Files in data/inputs/ are copied from DocDB and are not intended to be used
directly.  These are parsed and reformatted to produce files in other data/
directories which are for use.

data/desi.yaml
    Basic scalar parameters, organized in a nested tree.

data/focalplane/
    Information about positioner locations and platescale.

data/footprint/
    The areas of the sky that will be observed by DESI, with RA, Dec of tiles.

data/sky/
    Sample sky spectra.

data/specpsf/
    Spectrograph point-spread-function (PSF) for specter_
    CCD pixel-level simulations.

data/spectra/
    Example benchmark spectra.

data/targets/
    Expected n(z) information per target class for cosmology projections.

data/throughput/
    Throughput *versus* wavelength (also contained in specpsf).

data/weather/
    Historical weather data.

.. _specter: https://github.com/desihub/specter

Branches
~~~~~~~~

There are several permanent branches that were used for testing
alternative designs.  Some of these will never be merged into trunk but we
will keep them around for the record:

altccd
    500 micron *versus* 250 micron thick CCDs.

bb
    Recreating assumptions used during early BigBOSS projections.

newtiles
    An improved tiling dither pattern from Eddie Schlafly, intended
    to be merged prior to the start of the DESI survey.

update_inputs
    This branch *might* be present during updates to the inputs to
    the desimodel data files.  See the `updating desimodel inputs`_ document
    for further details.

In addition to these historical branches, there are a set of ``test-*`` branches
that contain smaller versions of the desimodel files.  These branches are
intended for use in desimodel unit tests.  See the `desimodel testing`_
document for further details.

.. _`desimodel testing`: https://desimodel.readthedocs.io/en/latest/testing.html
.. _`updating desimodel inputs`: https://desimodel.readthedocs.io/en/latest/update_inputs.html

Tagging
-------

If *either* the data *or* the code changes, a new tag should be created in
both git and svn.

Full Documentation
------------------

Please visit `desimodel on Read the Docs`_

.. _`desimodel on Read the Docs`: https://desimodel.readthedocs.io/en/latest/

License
-------

desimodel is free software licensed under a 3-clause BSD-style license. For details see
the ``LICENSE.rst`` file.
