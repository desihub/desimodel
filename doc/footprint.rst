=========
Footprint
=========

Files
=====

`desi-tiles.fits` contains a tiling of the sky in 10 overlapping layers,
with column `IN_DESI` indicating whether a particular tile is within the
DESI footprint.  See DESI-0717 for a description of how `desi-tiles.fits`
was generated.

`desi-healpix-weights.fits` contains a nside=256 nested healpix map of weights
for what fraction of those healpix is covered by the footprint.  This was
generated with::

  from desimodel.footprint import pixweight
  newmap = pixweight(256,outfile="desi-healpix-weights.fits")
