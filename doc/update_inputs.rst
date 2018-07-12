===============================
How to update inputs from DocDB
===============================

Context: Inputs from DocDB that require reformatting before using are kept
in data/inputs/ and then reformatted using scripts in bin/ .  Inputs from
DocDB that can be used as-is are directly added to their final location
under data/.

Basic Setup
===========

Make branches of both desimodel GitHub code and svn data::

    git checkout -b update_inputs
    base=https://desi.lbl.gov/svn/code/desimodel
    svn copy $base/trunk $base/branches/update_inputs
    svn checkout $base/branches/update_inputs/data
    export DESIMODEL=`pwd`

Add entries to $HOME/.netrc to enable downloading DocDB files without
having to enter a password every time::

    machine desi.lbl.gov
    login StephenBailey
    password NotMyRealPassword

The code in :mod:`desimodel.inputs.docdb` requires `requests`_
(for communicating with DocDB) and `xlrd`_ (for reading Microsoft Excel spreadsheets).
Both of these are available via Anaconda.

.. _`requests`: http://docs.python-requests.org/en/master/
.. _`xlrd`: http://www.python-excel.org/

Inputs to update
================

DESI-0530-v13 Excel spreadsheet to .ecsv file for GFA locations
—--------------------------------------------------------------

This writes the .ecsv file containing the GFA data, which
is pulled from the "GFALocation" tab on the DESI-0530-v13 Excel spreadsheet
and from rows 16-23 and columns A-I. The function
:func:`~desimodel.inputs.gfa.build_gfa_table` writes the file in the current directory.

::

    from desimodel.inputs.gfa import build_gfa_table
    build_gfa_table()

Positioner to Fiber Mapping
---------------------------

This updates the mapping of device locations on the focal plane to
spectrograph fiber numbers using DESI-0530 and DESI-2721.

::

    import desimodel.inputs.fiberpos
    desimodel.inputs.fiberpos.update()

Throughput
----------

This updates the throughput model from DESI-0347 and DESI-0344.

::

    import desimodel.inputs.throughput
    desimodel.inputs.throughput.update()

To update the version of DESI-347 that is used, change the default value of
``desi347_version`` in the ``update()`` function.  Similiarly for DESI-344.

Only three rows of the throughput spreadsheet from DESI-347 are used, with
hard-coded row numbers.  There are some simple checks that these are correct,
using the ``specthru_row`` and ``thru_row`` arguments to ``load_throughput()``,
but check the outputs carefully if you think the spreadsheet structure might
have changed.

To Do
=====

Update methodology and document how to update the following:

  * PSF model from DESI-0334
  * PSF spots -> PSF for quicksim
  * Fiber input loss calculations
  * desi.yaml
  * desimodel/data/focalplane/platescale.txt
  * trimming the full files into a small test set (:mod:`desimodel.trim`)
  * Optical distortions

      * data/inputs/throughput/raytracing.txt
      * data/throughput/DESI-0347_blur.ecsv
      * data/throughput/DESI-0347_offset.ecsv
