=================
Fiber Positioners
=================

Files
=====

data/inputs/focalplane/DESI-0530-posloc.txt defines the location of the
positioners on the focal plane, taken from DESI-0530 tab
PositionerAndFiducialLocations cells B40-J589 .

data/inputs/focalplane/cassette_order.txt defines the order of which
cassettes map to which ranges of fibers on the spectrograph slithead.

bin/randomize_fibers combines these while randomizing fibers within a
cassette, with output to data/focalplane/fiberpos.* in
fits and ascii formats.

The square cutouts in the corner of each focal plane wedge are for the
guider/focus cameras.  Additionally there are positioner holes reserved
for unmoving fiducials.

fiberpos.txt is a text file with columns

* fiber [0-4999]
* positioner [1-5000]  (may change to 0-4999 in the future)
* spectro [0-9]
* x, y, z [mm]

fiberpos.fits contains a binary table with the same information.

fiberpos-all.* includes information about non-spectro positioner holes, e.g.
locations of fidicials (postype=FIF) and sky monitor fibers (postype=ETC).

Coordinate Systems
==================

As defined in DESI-0481v1 table 1 (last row "CS5 DESI Focal Plane"),
when the telescope is parked at zenith the +x direction is East (+RA),
+z is pointed toward the ground, and thus +y points south (-dec).

Note that y and z are sign swapped wrt to the other DESI coordinate systems:
most coordinate systems think about +z pointing toward the sky, while the
focal plane thinks about +z as pointing away from the focal plane and thus
away from the sky.
