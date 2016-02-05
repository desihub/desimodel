=======================
desimodel Release Notes
=======================

1.0.0 (unreleased)
------------------

* "First" post-separation tag.
* Add code to create trimmed versions of the data.
* No significant changes to *data* from 0.4.2.

0.4.2 (2016-02-04)
------------------

* Improved svn download instructions in the README file.
* Changes to data on svn side

  * updated desi.yaml with dark vs. bright exptime
  * updated targets.dat to include MWS placeholders

* io.load_params() adds 'exptime' -> 'exptime_dark' key for temporary
  backwards compatibility
* Removed deprecated fibers module
* Use `ci-helpers`_ to handle most of the dirty work of Travis build scripts.
* Make `specter`_ import errors more verbose.

.. _`ci-helpers`: https://github.com/astropy/ci-helpers
.. _`specter`: https://github.com/desihub/specter

0.4.1 (2016-01-25)
------------------

* Last tag prior to separating desimodel into code (GitHub) and data (svn)
  repositories.
* pip install support (BAW).
* Replace fitsio dependency with astropy.io.fits.

0.4 (2015-12-14)
----------------

* Added tile file for the bright time survey.

0.3.8 (2015-10-30)
------------------

* Adds python io library (SJB).

0.3.7 (2015-04-16)
------------------

* Tag to support dogwood production (SJB).

0.3.6 (2015-01-30)
------------------

* Adds ``focalplane/FocalPlane.xy2radec()`` from Jaime (SJB).

0.3.5 (2014-12-28)
------------------

data/targets/targets.dat
    added fractions for sky and stdstar fibers (SJB).

py/desimodel/focalplane.py
    bug fixes for transformations (SJB).

0.3.4 (2014-09-23)
------------------

* Fix a simple import error (BAW).

0.3.3 (2014-09-23)
------------------

* Fix a simple version error (BAW).

0.3.2 (2014-09-23)
------------------

* Change how version is set (BAW).
* Updated target numbers.

0.3.1 (2014-07-23)
------------------

* Also updated quicksim sn-spec* file output, using IDL version which is slightly
  more optimistic than the python version (diff is dark current?) (SJB).

0.3 (2014-07-23)
----------------

* Updated throughput files for real.
* Added initial "compare_versions.py" script to make it easier to visualize
  differences in versions.  This script should grow as various parameters
  change; right now it only makes a thoughput difference plot (SJB).
* Updated throughput files from 0334v3 (spectro) and 0347v5 (system throughput)
  Correction: thoughput files didn't make it into that change (SJB, 2014-07-08).
* Updated psf-b.fits and psf-quicksim.fits to match new npix_y for blue
  STA/ITL CCDs (SJB, 2014-07-08).

0.2 (2014-07-08)
----------------

2014-07-07 SJB
~~~~~~~~~~~~~~

* Added ELG spectrum with continuum and multiple emission lines

2014-07-07 David Kirkby
~~~~~~~~~~~~~~~~~~~~~~~

Python quicksim

* add readnoise contributions in quadrature during the downsampling
* Refactor for speed, results now named ndarray, updated plots
* Allow different base directories

2014-07-02 DJS
~~~~~~~~~~~~~~

* Put sky back to dimmer UVES sky model

0.1 (2014-07-01)
----------------

2014-06-29 SJB
~~~~~~~~~~~~~~

* Extended fiberloss range from 3500-10000 instead of 3600-10000
* Added data/throughput/fiberloss-qso.dat (same as fiberloss-star.dat)

2014-06-27 SJB
~~~~~~~~~~~~~~

* Updated data/focalplane/platescale.txt with latest from DESI-0329v14.
  This includes a new "theta" column.
* Updated desi.yaml from DESI-0347v4.  This removes the FWHM and wavemin/max
  params which are not derived quantities associated with the PSFs.
* Updated throughput files with new numbers from DESI-0347v4.
* Updated spectrograph throughput files with new numbers from DESI-0334v2.
* Updated py/fiberloss.py -> bin/fiberloss.py .  Biggest change is ELG
  half light radius 0.35" -> 0.45" which drops us below 7-sigma.
* Updated data/throughput/fiberloss-\*.dat files with calculation based
  upon fiberloss.py
* bin/psf2quicksim.py extracted PSF parameters needed for quicksim.
    - pro/desi_quicksim.pro updated, but it still treats FWHM as constant
      rather than wavelength dependent.
    - python quicksim will be broken until it is updated to use new inputs.
* Reorganized data/inputs/throughput/
* spots2psf.py: leftover spot mirroring bug removed, PSFs updated

2014-06-12 SJB
~~~~~~~~~~~~~~

* Updated throughputs to not double count central obscuration.
* Updated PSF files to remove throughputs to avoid possible inconsistency.
* Added wavemin_all, wavemax_all to desi.yaml with min/max wavelength
  seen by all spectra

2014-06-06 SJB
~~~~~~~~~~~~~~

* Updated CCD pixel dimensions and regenerated PSFs to match.
