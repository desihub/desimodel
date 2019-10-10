========================================
Focalplane Hardware Model
========================================

Overview
============

The DESI focalplane consists of petals with individual positioners and each
positioner is a unique device connected to a fiber.  The fibers travel to the
spectrographs and are connected to a position on a slitblock.  Each positioner
has a range of angular motion along two axes (theta and phi).

As of October, 2019, much of this information is scattered across many
documents.  We want to gather it all into one place for use in offline
analysis.  Eventually this information will be consolidated into a database and
our focalplane model will be based on that database instead of the many
individual pieces of upstream information.

Although the generated focalplane models are checked into svn, we need to be
able to get the state of the hardware at any time, across the full history of
the DESI project.  We enable this feature through the use of our file format.
Since most changes are small (a positioner breaks or gets stuck, etc), we keep
a running log of these small changes.  A completely new model is only generated
for large events (e.g. a petal is swapped out).  When a focalplane is loaded,
the most recent model for a given time is found and the events in this log are
replayed up to the requested time.

Inputs
============

The current inputs require a local svn checkout of the fp_settings repo.  This
includes a "pos_settings" directory of python ConfigObj files, one per device.
It also includes a collision_settings directory with some ConfigObj files
representing the exclusion polygons as line segments.

Additionally, the following DocDB files are downloaded and parsed:

==============   ===========
DocDB Number     Purpose
==============   ===========
0530             X / Y locations of positioner centers in petal coordinates
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

A new focalplane model is generated with the desi_generate_focalplane script:

.. include:: generate_focalplane.inc

This script will download and cache some files from DocDB.  It also requires
that you have a local svn checkout of the pos_settings files.  Many of the
commandline options are for testing / debugging against previous
representations of the DESI instrument.  A typical use of the this script would
be::

    desi_generate_focalplane \
        --startvalid 2019-10-10T00:00:00 \
        --pos_settings /path/to/svn/fp_settings/pos_settings \
        --collision /path/to/svn/fp_settings/collision_settings/_collision_settings_DEFAULT.conf \
        --exclusion default \
        --petal_id2loc XXXXXXXXXX

This creates a new focalplane model in the desimodel data directory.  It starts
with all positioner devices in a "good" state.

Updating the State of a Model
================================

After a focalplane model is created, one can update the state of an individual
positioner with the following command line tool:

.. include:: update_focalplane.inc


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
polygons for the positioners.  Each named exclusion entry actually has 2
polygons- one for the theta arm and one for the phi arm.  These define the
shape of the polygon at the origin, which is then translated and rotated
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

- We should add an event for CANbus failure.  We only need to mark this once in
  the log and can then map that to multiple devices.  See
  https://github.com/desihub/desimodel/issues/123

- We should build this focalplane model from the online instrument DB.  See
  https://github.com/desihub/desimodel/issues/124
