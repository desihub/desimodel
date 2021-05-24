=======================
desimodel Release Notes
=======================

0.15.1 (unreleased)
-------------------

* Update ``etc/desimodel_sync_kpno_cron.sh`` for automatically syncing the
  focalplane model to the latest DB dump (PR `#148`_).

.. _`#148`: https://github.com/desihub/desimodel/pull/148

0.15.0 (2021-04-19)
-------------------

Code in GitHub:

* Use UTC time everywhere in the focalplane model (PR `#147`_).
  This is backwards compatible with old files, but new FP models will not
  be readable by previous code tags.

Data in svn:

* DB sync 2021-04-03T23:53:23 appended to ``desi-state-2021-03-17T23:20:01.ecsv``.
* DB sync 2021-04-10T20:00:39 appended to ``desi-state_2021-03-17T23:20:01.ecsv``.
* DB sync 2021-04-13T20:00:30 appended to ``desi-state_2021-03-17T23:20:01.ecsv``.

.. _`#147`: https://github.com/desihub/desimodel/pull/147

0.14.2 (2020-03-31)
-------------------

Data changes to svn, no code changes:

* Added LYA TSNR2 templates.
* Focalplane model updated 2021-03-17.

0.14.1 (2021-03-18)
-------------------

* Add fastfiberacceptance code originally in specsim (PR `#145`_).

.. _`#145`: https://github.com/desihub/desimodel/pull/145

0.14.0 (2021-02-10)
-------------------

* Code (in GitHub):

  * travis test fixes for old astropy (PR `#141`_).
  * move command line scripts from svn to git (PR `#142`_).
  * add option to exclude petals from restricted reach (PR `#144`_).

* Data (in svn):

  * added Template Signal-to-Noise (TSNR) ensembles
  * added pre-calculated Noise Equivalent Area (NEA) from PSF model
  * Corrected restricted reach focalplane model (@135002)
  * Focal plan model with full reach for petal locs 0,2,4,5 (@135236)

.. _`#141`: https://github.com/desihub/desimodel/pull/141
.. _`#142`: https://github.com/desihub/desimodel/pull/142
.. _`#144`: https://github.com/desihub/desimodel/pull/144

0.13.1 (2020-08-03)
-------------------

* New tag of data+code since data had been erroneously pre-tagged 0.13.0

0.13.0 (2020-08-03)
-------------------

* Fix py3.8 syntax warnings (PR `#140`_).
* Fix corner cases in generating and using focalplane models (PR `#139`_).
* Use DESI-5501 (as built) instead of DESI-334 (design) for spectrograph throughput (PR `#137`_).

.. _`#140`: https://github.com/desihub/desimodel/pull/140
.. _`#139`: https://github.com/desihub/desimodel/pull/139
.. _`#137`: https://github.com/desihub/desimodel/pull/137

0.12.0 (2020-03-13)
-------------------

* update platescale to as-built DESI-4037v5 (PR `#136`_).
* update desi-focalplane model for limited phi range 20200306 (svn data).
* fix bug in generating focalplane model from old fiberpos files (PR `#139`_).
* use >= not > when comparing runtime to focalplane model `#139`_).

.. _`#136`: https://github.com/desihub/desimodel/pull/136
.. _`#139`: https://github.com/desihub/desimodel/pull/139

0.11.0 (2020-03-13)
-------------------

* Updated ``data/footprint/desi-tiles.fits`` and
  ``desi-healpix-weights.fits`` with new dither pattern; see DESI-0717.
  Layers 0=GRAY, 1-4=DARK instead of 0-3=DARK, 4=GRAY. (PR `#135`_).
* Update documentation for :mod:`desimodel.io`; use
  :func:`desimodel.io.findfile` consistently throughout the module (PR `#133`_).
* Update README file and Travis tests (PR `#132`_).
* Include S (curved focal surface arc length) vs. R (CS5 xy radius)
  table from DESI-0530 (PR `#130`_ and `#135`_).

.. _`#130`: https://github.com/desihub/desimodel/pull/130
.. _`#132`: https://github.com/desihub/desimodel/pull/132
.. _`#133`: https://github.com/desihub/desimodel/pull/133
.. _`#135`: https://github.com/desihub/desimodel/pull/135

0.10.3 (2019-12-20)
-------------------

* Pass multiple sets of exclusion polygons (PR `#128`_).
* Propagate existing focalplane state to new focalplanes (PR `#129`_).

.. _`#128`: https://github.com/desihub/desimodel/pull/128
.. _`#129`: https://github.com/desihub/desimodel/pull/129

0.10.2 (2019-10-31)
-------------------

* Improve focalplane creation code (PR `#127`_).

.. _`#127`: https://github.com/desihub/desimodel/pull/127

0.10.1 (2019-10-17)
-------------------

* Workaround upstream bugs in positioner locations (PR `#118`_).
* Added `desimodel.focalplate.fieldrot.field_rotation_angle` with
  field rotation CS5 vs. ICRS due to precession (PR `#119`_).
* Add focalplane model documentation (PR `#125`_).

.. _`#118`: https://github.com/desihub/desimodel/pull/118
.. _`#119`: https://github.com/desihub/desimodel/pull/119
.. _`#125`: https://github.com/desihub/desimodel/pull/125

0.10.0 (2019-09-25)
-------------------

* Store petal and gfa keepouts in the focalplane model (PR `#112`_).
* When generating a focalplane, check for device locations assigned to the
  same slitblock and fiber (PR `#113`_).
* Fix support for ``datetime.isoformat()`` in Python 3.5 (PR `#114`_).
* Update tests and documentation to be consistent with latest desiutil versions (PR `#115`_).

.. _`#112`: https://github.com/desihub/desimodel/pull/112
.. _`#113`: https://github.com/desihub/desimodel/pull/113
.. _`#114`: https://github.com/desihub/desimodel/pull/114
.. _`#115`: https://github.com/desihub/desimodel/pull/115


0.9.12 (2019-08-09)
-------------------

* Support for time-varying focal plane state (*e.g.* broken fibers) (PR `#105`_).
* Documentation about CI weather *versus* model (PR `#107`_).
* Fix :func:`~desimodel.footprint.find_points_radec` for scipy 1.3 (PR `#109`_).
* Replace deprecated ``yaml.load`` with ``yaml.safe_load`` (PR `#110`_).

.. _`#105`: https://github.com/desihub/desimodel/pull/105
.. _`#107`: https://github.com/desihub/desimodel/pull/107
.. _`#109`: https://github.com/desihub/desimodel/pull/109
.. _`#110`: https://github.com/desihub/desimodel/pull/110

0.9.11 (2019-05-30)
-------------------

* Added data/footprint/ci-tiles-v7.fits, data/focalplane/ci-corners.ecsv
  to svn and docs to GitHub (PR `#103`_).

.. _`#103`: https://github.com/desihub/desimodel/pull/103

0.9.10 (2019-02-28)
-------------------

* ``io.load_tiles(tilesfile)`` warns if local copy exists, but :envvar:`DESIMODEL`
  version wins (PR `#98`_ and `#101`_).
* Update default tile radius (max radius, not typical outer pos radius)
  (PR `#102`_).

.. _`#98`: https://github.com/desihub/desimodel/pull/98
.. _`#101`: https://github.com/desihub/desimodel/pull/101
.. _`#102`: https://github.com/desihub/desimodel/pull/102

0.9.9 (2018-09-27)
------------------

* Change default healpy pixel overlap factor from 4 to 128 (PR `#93`_).

.. _`#93`: https://github.com/desihub/desimodel/pull/93

0.9.8 (2018-09-05)
------------------

* Implement :func:`~desimodel.weather.dome_close_fractions` to replay daily Mayall weather history (PR `#92`_).
* Run tests using new svn branch test-0.9.8.
* Bug fix for GFA target selection when no targets overlap a GFA (PR `#91`_).

.. _`#91`: https://github.com/desihub/desimodel/pull/91
.. _`#92`: https://github.com/desihub/desimodel/pull/92

0.9.7 (2018-07-30)
------------------

* Create DESI-3977 in doc/tex/desi3977/ to track ELG SNR with changes to the DESI model.
* Add accompanying notebook doc/nb/ELG_SNR.ipynb.

0.9.6 (2018-07-18)
------------------

* Update data and associated code to reflect changes in DESI-347-v13 (PR `#89`_):
  * ``data/throughput/thru-[brz].fits``: new corrector coatings.
  * ``data/throughput/DESI-0347_blur.ecsv``: new achromatic blurs.
  * ``data/desi.yaml``: new read noise and dark currents.
  * ``data/focalplane/gfa.ecsv``: replace ``RADIUS_MM`` with ``S``.
  * ``data/throughput/DESI-0347_static_[123].fits``: replace random offset files (RMS=10.886um) with static offset files (RMS=8.0um).
* Use a new svn branch test-0.9.6 for travis tests (was test-0.9.3).

.. _`#89`: https://github.com/desihub/desimodel/pull/89

0.9.5 (2018-06-27)
------------------

* Increase test coverage, especially for :mod:`desimodel.trim` (PR `#82`_).
* Reorganize :mod:`desimodel.focalplane` and add more GFA selection code (PR `#85`_).
* Allow an environment variable in the tilesfile filename (PR `#87`_).

.. _`#82`: https://github.com/desihub/desimodel/pull/82
.. _`#85`: https://github.com/desihub/desimodel/pull/85
.. _`#87`: https://github.com/desihub/desimodel/pull/87

0.9.4 (2018-03-29)
------------------

* Download script will create :envvar:`INSTALL_DIR` if it doesn't exist (PR `#80`_).

.. _`#80`: https://github.com/desihub/desimodel/pull/80

0.9.3 (2018-03-14)
------------------

* Fix some installation bugs, and update to latest versions on various
  dependencies (PR `#77`_).
* Ensure that desimodel tests are compatible with Astropy 2 and 3, and with
  other DESI packages (PR `#78`_).
* Add ``footprint/desi-healpix-weights.fits`` and
  ``throughput/galsim-fiber-acceptance.fits`` to the trimmed test data set
  (PR `#79`_).

.. _`#77`: https://github.com/desihub/desimodel/pull/77
.. _`#78`: https://github.com/desihub/desimodel/pull/78
.. _`#79`: https://github.com/desihub/desimodel/pull/79


0.9.2 (2018-02-27)
------------------

* Update LyA S/N calculation (PR `#73`_).
* Optionally use an input pixel weight map in :func:`~desimodel.io.load_pixweight`
  (PR `#74`_).

.. _`#73`: https://github.com/desihub/desimodel/pull/73
.. _`#74`: https://github.com/desihub/desimodel/pull/74

0.9.1 (2017-11-10)
------------------

* Extracts wavelength coverage from specpsf files into params dictionary
  (PR `#68`_).
* Added :func:`~desimodel.footprint.program2pass` and
  :func:`~desimodel.footprint.pass2program` to convert between
  tiling integer pass number and string program name (PR `#67`_).

.. _`#67`: https://github.com/desihub/desimodel/pull/67
.. _`#68`: https://github.com/desihub/desimodel/pull/68

0.9.0 (2017-09-19)
------------------

* Added desimodel.focalplane.radec2xy, which converts RA, Dec coordinates to x, y coordinates on the focal plane, which accepts vector inputs.
* Added desimodel.focalplane.on_gfa() and its respective helper functions to check if a target is on a GFA of arbitrary telescope pointing
* Added desimodel.focalplane.on_tile_gfa() to check return a list of indices of targets on a specific tile
* Added desimodel.focalplane.get_gfa_targets() to return a table with added columns GFA_LOC and TILEID that consists of all targets on any GFA on any tile satisfying a minimum flux in the r-band.
* Unittests for the desimodel.focalplane functions were updated accordingly.
* Added desimodel.footprint.find_points_in_tel_range() to return a list of indices withnin a radius of an arbitray telescope pointing, unaware of tiles (Added respective unittest)
* Adds desimodel.focalplane.fiber_area_arcsec2()
* Updates tests to work with trimmed data subset

0.8.0 (2017-08-07)
------------------

* Add new weather module to specify assumed atmospheric seeing and transparency
  distributions at KPNO, with accompanying DESI-doc and jupyter notebook.
* Remove seeing module, which is superseded by new weather module.
* Added `desimodel.footprint.pixweight()` in :mod:`desimodel.footprint` to create an array of what fraction
  of every HEALPixel at a given nside overlaps the DESI footprint
* Also added `desimodel.footprint.tiles2fracpix()` to estimate which HEALPixels overlap the footprint edges
* Added `desimodel.io.load_pixweight()` in :mod:`desimodel.io` to load the array created by
  `desimodel.footprint.pixweight()` and resample it to any HEALPix nside
* Modified path to Lya SNR spectra files used in desi_quicklya.py, used in Lya Fisher forecast.
* Added desimodel.inputs.build_gfa_table and its helper functions to write a .ecsv file for GFA data
* Added desimodel.io.load_gfa to return the GFA data table
* Added desimodel.focalplane.xy2radec, which converts x,y coordinates on the focal plane to RA, Dec coordinates
* don't print warnings in desimodel.io if specter isn't installed

0.7.0 (2017-06-15)
------------------

* Added desimodel.footprint.tiles2pix and .pix2tiles for mapping healpix
  to DESI tiles.
* fixed psf-quicksim.fits units to be astropy-friendly
* added `desimodel.io.load_target_info()`

0.6.0 (2017-03-27)
------------------

* Add desimodel.seeing module with functions that model the expected DESI
  zenith seeing at 6355A, with an accompanying jupyter notebook.
* Altered xy offset RMS calculation in focalplane.py to scale the distribution
  RMS rather than the sample standard deviation.
* Update focal plane to positioner mapping
* z-channel 250 um CCD instead of 500 um CCD
* Update DocDB -> desimodel update method for fiberpos and throughput

0.5.1 (2016-12-01)
------------------

* By default, desimodel.io.load_tiles now excludes PROGRAM=EXTRA layers
* Adds desi-tiles.* tests

0.5.0 (2016-11-21)
------------------

* Moved test of focalplane code into the actual test suite.
* Preparing for Python 3.
* Changed default svn version to trunk and added error handling to
  :command:`install_desimodel_data`.
* Update template module file to reflect DESI+Anaconda infrastructure.
* Add code to generate random centroid offsets in :mod:`desimodel.focalplane`.
* Add jupyter notebook documenting new throughput files of `PR#29`_.
* Use Astropy-recommended method of reading FITS data tables.
* Remove reference to Travis scripts in MANIFEST.in.

.. _`PR#29`: https://github.com/desihub/desimodel/pull/29

0.4.5 (2016-07-15)
------------------

* Fixed a minor bug that made the help message for :command:`install_desimodel_data`
  garbled.
* Add additional files to lightweight test data to work with quickgen

0.4.4 (2016-03-15)
------------------

* Allow :command:`desiInstall` to download and install the data from svn.
* No changes to data in svn.

0.4.3 (2016-03-10)
------------------

* "First" post-separation tag.
* Added :func:`desimodel.trim.trim_data` for trimming a data directory into a
  lightweight version for testing.
* svn data includes targets.dat: preliminary numbers for MWS and BGS densities
  (Still waiting upon supporting technote).

0.4.2 (2016-02-04)
------------------

* Improved svn download instructions in the README file.
* Changes to data on svn side

  * updated desi.yaml with dark vs. bright exptime
  * updated targets.dat to include MWS placeholders

* :func:`desimodel.io.load_desiparams` adds 'exptime' -> 'exptime_dark' key for temporary
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

* Adds :meth:`desimodel.focalplane.FocalPlane.xy2radec` from Jaime (SJB).

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
