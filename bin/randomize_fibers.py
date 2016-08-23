#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sys import exit
#
#
#
def main():
    """
    Randomize fibers from focal plane to spectrograph.

    User pos_on_z1.txt positioner location file from Joe Silber.
    Knock out positions for fiducials if necessary and
    randomly assign fibers to positioners within a sector.
    Write out both a fits file and a text file with the fiber -> (x,y) map.

    Stephen Bailey, LBL
    Fall 2013
    """
    import random
    import numpy as np
    import fitsio
    #- input parameters
    import argparse
    from sys import argv

    parser = argparse.ArgumentParser(prog=argv[0])
    parser.add_argument("-i", "--input", action='store', metavar='FILE',
        default='pos_on_z1.txt', help="positioner location filename")
    parser.add_argument("-o", "--output", action='store', metavar='FILE',
        default='fiberpos.txt', help="fiber position filename")
    parser.add_argument("--debug",   help="debug with ipython prompt at end", action="store_true")

    opts = parser.parse_args()

    #- Random but reproducible
    seed = 2
    random.seed(seed)                  #- for select
    np.random.seed(seed)               #- for shuffle

    #- Basic DESI parameters hardcoded
    nspectro = 10
    nfib_per_spectro = 500
    nfib = nspectro * nfib_per_spectro

    #- Read positioner location file
    columns = [('n', int),
               ('x', float),
               ('y', float),
               ('z', float),
               ('precess', float),
               ('nutation', float),
               ('spin', float),
               ('rotsym', int),
               ('remove', int),
              ]
    data = np.loadtxt(opts.input, dtype=columns).view(np.recarray)

    #- Remove the locations flagged for the guide cameras
    #- (and eventually the fiducials)
    ii = data['remove'] == 0
    data = data[ii]

    #- For now # sectors == # spectrographs.
    #- If that changes in a later design, crash now and fix code
    assert nspectro == len(set(data['rotsym']))

    #- setup fiber position table
    fibercols = [('fiber', int),
                 ('positioner', int),
                 ('spectrograph', int),
                 ('x', float),
                 ('y', float),
                 ('z', float),
                ]
    fiberpos = np.zeros(nfib, dtype=fibercols).view(np.recarray)

    #- Randomize fiber assignment within a spectrograph
    ifiber = np.arange(nfib_per_spectro)
    for i in range(nspectro):
        #- Get positions for this spectrograph
        pos = data[(data.rotsym == i) & (data.remove == 0)].copy()

        #- Mark any extras for removal as fiducials
        nfid = len(pos) - nfib_per_spectro
        assert nfid >= 0
        if nfid > 0:
            #- Check min/max angles to avoid removing edge positioners
            #- Oddly, np.math.atan2 only accepts scalar input
            phi = [np.math.atan2(y, x) for y,x in zip(pos.y, pos.x)]
            phi = (np.array(phi)+2*np.pi)%(2*np.pi)  #- [0,2pi]

            #- interior angles ok to knock out fibers
            ii = np.where((np.min(phi)+0.1 < phi) & (phi < np.max(phi)-0.1))[0]
            remove = random.sample(ii, nfid)
            pos['remove'][random.sample(ii, nfid)] = 2
            pos = pos[pos['remove'] == 0]
            assert len(pos) == nfib_per_spectro

        #- Shuffle positioner -> fiber assignment
        #- Directly shuffling pos doesn't appear to work.  Bug in shuffle ?!?
        ii = np.arange(len(pos))
        np.random.shuffle(ii)
        pos = pos[ii]

        #- Fill fiberpos array
        ii = list(range(i*nfib_per_spectro, (i+1)*nfib_per_spectro))
        fiberpos['fiber'][ii] = ifiber + i*nfib_per_spectro
        fiberpos['positioner'][ii] = pos['n']
        fiberpos['spectrograph'][ii] = i
        fiberpos['x'][ii] = pos['x']
        fiberpos['y'][ii] = pos['y']
        fiberpos['z'][ii] = pos['z']

    if opts.output.endswith('.txt'):
        textout = opts.output
        fitsout = opts.output.replace('.txt', '.fits')
    elif opts.output.endswith('.fits'):
        fitsout = opts.output
        textout = opts.output.replace('.fits', '.txt')

    pngout = textout.replace('.txt', '.png')

    #- Write fits file, then fix up comments (!)
    fitsio.write(fitsout, fiberpos, extname='FIBERPOS', clobber=True)
    fx = fitsio.FITS(fitsout, 'rw')
    hdr = fx['FIBERPOS'].read_header()
    fx['FIBERPOS'].write_key('TTYPE1', hdr['TTYPE1'], 'fiber number [0-4999]')
    fx['FIBERPOS'].write_key('TTYPE2', hdr['TTYPE2'], 'positioner number [0-4999]')
    fx['FIBERPOS'].write_key('TTYPE3', hdr['TTYPE3'], 'spectrograph number [0-9]')
    fx['FIBERPOS'].write_key('TTYPE4', hdr['TTYPE4'], 'positioner x center [mm]')
    fx['FIBERPOS'].write_key('TTYPE5', hdr['TTYPE5'], 'positioner y center [mm]')
    fx['FIBERPOS'].write_key('TTYPE6', hdr['TTYPE6'], 'positioner z location [mm]')
    fx.close()

    #- Write a text file format
    fxlines = [
        "#- Fiber to positioner mapping; x,y,z in mm on focal plane",
        "#- See doc/fiberpos.rst for more details.",
        "#- Coordinates at zenith: +x = East = +RA; +y = South = -dec",
        "",
        '#- fiber positioner spectro  x  y  z']
    for row in fiberpos:
        ### fxlines.append("{fiber:4d}  {positioner:4d}  {spectrograph:2d}  {x:12.6f}  {y:12.6f}  {z:12.6f}".format(**row))
        fxlines.append("{:4d}  {:4d}  {:2d}  {:12.6f}  {:12.6f}  {:12.6f}".format(*row))
    with open(textout,'w') as fx:
        fx.write('\n'.join(fxlines)+'\n')
    fx.close()

    import pylab as P
    P.figure(figsize=(7,7))
    P.scatter(fiberpos.x, fiberpos.y, c=fiberpos.fiber, edgecolor='none')
    P.grid()
    P.xlim(-420,420)
    P.ylim(-420,420)
    P.xlabel('x [mm]')
    P.ylabel('y [mm]')
    P.savefig(pngout, dpi=80)

    if opts.debug:
        #--- DEBUG ---
        P.ion()
        P.show()
        import IPython
        IPython.embed()
        #--- DEBUG ---

    return 0
#
#
#
exit(main())
