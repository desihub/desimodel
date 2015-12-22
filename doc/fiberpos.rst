=================
Fiber Positioners
=================

Files
=====

etc/data/inputs/focalplane/pos_on_z1.txt is a text file version of DESI-0403v2,
mapping the location of the fiber positioners in the DESI focal plane.

bin/randomize_fibers randomly assigns fibers to positioners and makes
a mapping for each spectrograph, output to etc/data/focalplane/fiberpos.* in
fits and ascii formats.

The square cutouts in the corner of each focal plane wedge are for the
guider/focus cameras.  Additionally there are random holes within
the plane without a fiber assigned to mimic placement of fiducials.
Hopefully the final DESI system will have more evenly distributed fiducials.

fiberpos.txt is a text file with columns

* fiber [0-4999]
* positioner [1-5000]  (may change to 0-4999 in the future)
* spectro [0-9]
* x, y, z [mm]

fiberpos.fits contains a binary table with the same information.

Coordinate Systems
==================

As defined in DESI-0481v1 table 1 (last row "CS5 DESI Focal Plane"),
when the telescope is parked at zenith the +x direction is East (+RA),
+z is pointed toward the ground, and thus +y points south (-dec).

Note that y and z are sign swapped wrt to the other DESI coordinate systems:
most coordinate systems think about +z pointing toward the sky, while the
focal plane thinks about +z as pointing away from the focal plane and thus
away from the sky.
