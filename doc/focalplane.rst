========================================
Focalplane Hardware Model
========================================

Overview
============

The DESI focalplane consists of petals with individual positioners and each
positioner is a unique device connected to a fiber.  The fibers travel to the
spectrographs and are connected to a position on a slitblock.  Each positioner
has a range of angular motion along two axes (theta and phi).

The canonical focalplane calibration is stored in the ICS database,
but we want to copy it into desimodel for offline analysis and to be
able to tag exactly what version was used for fiber assignment and
data processing runs.

The focalplane model currently uses the CS5 coordinate system for the X/Y
locations of devices on the focalplane.

Tracking Changes
====================

Although the generated focalplane models are checked into svn, we need to be
able to get the state of the hardware at any time, across the full history of
the DESI project.  We enable this feature through the use of our file format
(see section below on the details of the format). Since most changes are small
(a positioner breaks or gets stuck, etc), we keep a running log of these small
changes.  A completely new model is only generated for large events (e.g. a
petal is swapped out).  When a focalplane is loaded, the most recent model for
a given time is found and the events in this log are replayed up to the
requested time.  It is worth emphasizing the previous text again.  When a
focalplane is loaded for a particular timestamp:

1.  Each focalplane model has a starting time, and remains valid until it is superceded by a newer model.  The most recent set of 3 focalplane files which comes before the requested timestamp is read from disk.  The static focalplane properties are kept as a Table and the exclusion polygons are read into a dictionary.

2.  The state log (one of the 3 files) is parsed line by line.  If the timestamp for that line comes before the requested time, then the event in that line is applied.  All positioners have an initial state specified in the state log with a timestamp that matches the starting time of the focalplane.  Subsequent events are appended to the state log with a timestamp.

The file formats used are text-based (ECSV and YAML).  **However**, these are
intended to be modified by the included scripts, which can ensure that the
formatting is correct.  The risk of typos and subtle errors if hand-editing
these files is large.  If you find that you frequently need to edit these
files, then please open an issue to document your use case.


Inputs
============

Focalplane models are generated and updated from ICS database dumps at KPNO in
``/data/focalplane/calibration/*.ecsv``.

Additionally, the following DocDB files are downloaded and parsed to
get the mapping of fiber focalplane
location to slitblock location (i.e. fiber number on the spectrographs).

==============   ===========
DocDB Number     Purpose
==============   ===========
4042             Petal verification for petal ID 02
4043             Petal verification for petal ID 03
4807             Petal verification for petal ID 04
4808             Petal verification for petal ID 05
4809             Petal verification for petal ID 06
4190             Petal verification for petal ID 07
4806             Petal verification for petal ID 08
4810             Petal verification for petal ID 09
4042             Petal verification for petal ID 10
4042             Petal verification for petal ID 11
==============   ===========


Generating a New Model
==========================

Under normal circumstances focalplane updates are done daily by a KPNO cronjob
running ``etc/desimodel_sync_kpno_cron.sh`` (update existing model)
or ``etc/desimodel_sync_kpno_force.sh`` (use ``--reset`` to make a new model).
The rest of this section documents what those are doing "under the hood",
but it should not be necessary to run by hand.

A new focalplane model is generated with the ``desi_sync_focalplane`` script:

.. code-block:: console

    usage: desi_sync_focalplane [-h] --calib_file CALIB_FILE [--test] [--reset]
                                [--simulate_good] [--debug_dir DEBUG_DIR]
                                [--commit]

    optional arguments:
      -h, --help            show this help message and exit
      --calib_file CALIB_FILE
                            The ECSV database dump file
      --test                Go through the process of updating the focalplane, but
                            do not actually write new files.
      --reset               Create a new focalplane model from the calib file,
                            ignoring all previous state information
      --simulate_good       Create a focalplane model for simulations. Non-broken
                            fibers set to good
      --debug_dir DEBUG_DIR
                            Override the output directory for debugging.
      --commit              Commit updated focalplane model to svn.


Note that the ``--reset`` option generates a new focalplane model,
while without that option it updates the state ledger of the current
focalplane model for only the positioners that changed.  In addition to
the ``--calib_file`` input with the latest ICS database focalplane dump,
this script automatically downloads the necessary DocDB entries listed above,
which requires you to have DESI DocDB credentials stored in your
``$HOME/.netrc`` file.

See desimodel tags 0.16.0 and prior for documentation of an older script
``desi_generate_focalplane`` which uses lab-measured focalplane metrology
from DESI SVN ``code/focalplane/fp_settings/pos_settings`` to create
a new focalplane model.  This has been removed from
newer versions of desimodel in favor of loading the metrology as measured
in-situ at KPNO.


Updating the State of a Model
================================

After a focalplane model is created, the state can be updated by rerunning
``desi_sync_focalplane`` with a new ICS database dump without using the
``--reset`` option.  If needed, one can override the database to
update the state of an individual
positioner with the following command line tool:

.. include:: update_focalplane_state.inc


Updating the Exclusion Polygons in a Model
================================================

After a focalplane model is created, one can update the available exclusion
polygons with the following command line tool:

.. include:: update_focalplane_excl.inc

Existing exclusion polygons with the same name as any input files will be
replaced.  New polygons will be appended.


File Format and Loading
============================

A single focalplane model (with a starting valid datetime) consist of 3 files
on disk.  These files contain matching date / time stamps that correspond the
the first valid time for that focalplane model.  For example::

    desi-focalplane_2019-09-16T00:00:00.ecsv
    desi-exclusion_2019-09-16T00:00:00.yaml
    desi-state_2019-09-16T00:00:00.ecsv

The first is an ECSV text file containing the static information such as the
positioner locations and angle range, the mapping of device locations to
fibers, etc.

============   ===========   =============
Column         Data Type     Description
============   ===========   =============
PETAL          int32         Petal location [0-9]
DEVICE         int32         Device location on the petal
LOCATION       int32         PETAL * 1000 + DEVICE
PETAL_ID       int32         The physical petal ID
DEVICE_ID      string        The physical device ID string
DEVICE_TYPE    string        The device type (POS, ETC, FIF)
SLITBLOCK      int32         The slit block where this fiber goes
BLOCKFIBER     int32         The fiber index within the slit block
CABLE          int32         The cable ID
CONDUIT        string        The conduit
FIBER          int32         PETAL * 500 + SLITBLOCK * 25 + BLOCKFIBER
FWHM           float64       FWHM at f/3.9
FRD            float64       FRD Throughput
ABS            float64       ABS Throughput
OFFSET_X       float64       X location of positioner center
OFFSET_Y       float64       Y location of positioner center
OFFSET_T       float64       THETA zero point angle
OFFSET_P       float64       PHI zero point angle
LENGTH_R1      float64       Length of THETA arm
LENGTH_R2      float64       Length of PHI arm
MAX_T          float64       Maximum THETA angle relative to OFFSET_T
MIN_T          float64       Minimum THETA angle relative to OFFSET_T
MAX_P          float64       Maximum PHI angle relative to OFFSET_P
MIN_P          float64       Minimum PHI angle relative to OFFSET_P
============   ===========   =============

The second file is a YAML format file which contains one or more exclusion
polygons for the positioners.  Each named exclusion entry actually has multiple
polygons:  for the GFA, petal boundary, theta arm and phi arm.  These define
the shape of the polygon at the origin, which is then translated and rotated
differently for every positioner based on the arm length, etc.  Exclusion
polygons are specified in terms of lists of circles and line segments.

The third file is another ECSV format file that contains the *state log* for
the focalplane model.  This is the running log of *events* that happen which
modify the instantaneous state of the focalplane.

============   ===========   =============
Column         Data Type     Description
============   ===========   =============
TIME           string        The timestamp of the event (UTC, ISO format)
PETAL          int32         Petal location [0-9]
DEVICE         int32         Device location on the petal
LOCATION       int32         PETAL * 1000 + DEVICE
STATE          uint32        State bit field (good == 0)
EXCLUSION      string        The exclusion polygon for this device
============   ===========   =============

The file formats above are documented for completeness, but you should not
generally read these manually.  Instead, one calls the `load_focalplane()`
function:

.. autofunction:: desimodel.io.load_focalplane
    :noindex:

The state table returned by this function contains the instantaneous state of
the focalplane at the input time (i.e. the log of events has already been
replayed and the state at the requested time is returned)


TO-DO
==========

There are several small features needed:

- When marking fibers as broken or stuck, their current X/Y or theta/phi
  location should be marked.  See
  https://github.com/desihub/desimodel/issues/122

- We should build this focalplane model from the online instrument DB.  See
  https://github.com/desihub/desimodel/issues/124
