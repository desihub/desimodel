========================================
Components to the throughput calculation
========================================

Introduction
============

THIS DESCRIPTION IS OUT OF DATE.

The various contributions to the throughput are stored in a binary table
following the throughput format used by Specter:

https://github.com/sbailey/specter/blob/dev/doc/datamodel/throughput.md

It is basically a binary table with columns:

* wavelength : in Angstroms
* extinction : atmospheric extinction in magnitudes per airmass
* fiberinput : geometrical loss at fiber input
* throughput : all other throughput terms, e.g. mirrors, spectrograph, CCDs

Different types of sources are affected by different combinations of
throughput terms.

========== ====== === =====
Term       OBJECT SKY CALIB
========== ====== === =====
EXTINCTION yes    yes no
FIBERINPUT yes    no  no
THROUGHPUT yes    yes yes
========== ====== === =====

Source types:

* OBJECT: astronomical objects, affected by all sources of throughput loss
* SKY: sky spectra do not have a geometrical loss term for the fiber
  input.  Positioner misalignments and changes to the atmospheric
  PSF still get the same amount of sky light down the fiber.
* CALIB: calibration lamps internal to the dome do not see atmospheric
  extinction or fiber geometric loss terms.

The sections below detail the input data used to generate this table.

Atmospheric Extinction
======================

Affects astronomical objects and sky spectra, but not calibration exposures.

Depends upon airmass; extinction curve from ZenithExtinction-KPNO.fits
included in this product.

Fiber Input
===========

Affects astronomical objects, but not sky or calibration spectra.

Includes both PSF/aperture losses and pointing/guiding mis-alignment.

From:

* DESI-0347v2
  * row 16 "PSF and Aperture efficiency"
  * row 27 "Lateral errors"
  * row 52 "Fiber Defocus overall budget"
  * Multiplied then interpolated with a cubic spline

All other telescope, fiber, and instrument throughputs
======================================================

Affects all object types.

From:

* DESI-0347v2
  * row 6 "Telescope to fiber input"
  * row 87 "Fiber System other"
  * row 95 Fiber slit block
  * row 100 VPHG: Substrate row 100
  * but *not* the other Spectrograph entries in DESI-0347v2 that are
    represented at higher resolution in DESI-0334v1 and DESI-0336v3
* Spectrograph: DESI-0334v1
  * DESI-0334-blue-thru.txt "total" column
  * DESI-0334-red-thru.txt "total" column
  * DESI-0334-NIR-thru.txt "total" column
* CCD: DESI-0336v3
  * Blue: e2vqe.txt
  * Red: lbnl250qe.txt
  * NIR: lbnl500qe.txt

Example
=======

::

    python $DESIMODEL/bin/combine_throughputs -o $DESIMODEL/data/throughput


