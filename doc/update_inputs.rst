===============================
How to update inputs from DocDB
===============================

Context: Inputs from DocDB that require reformatting before using are kept
in data/inputs/ and then reformatted using scripts in bin/ .  Inputs from
DocDB that can be used as-is are directly added to their final location
under data/.

Basic Setup
===========

Make branches of both desimodel github code and svn data::

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

Inputs to update
================

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

To Do
=====

Update methodology and document how to update the following:

  * PSF model from DESI-0334
  * PSF spots -> PSF for quicksim
  * Fiber input loss calculations
  * desi.yaml
  * desimodel/data/focalplane/platescale.txt
  * trimming the full files into a small test set (desimodel.trim)
  * Optical distortions

      * data/inputs/throughput/raytracing.txt
      * data/throughput/DESI-0347_blur.ecsv
      * data/throughput/DESI-0347_offset.ecsv

