#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sys import exit

def write_text_fiberpos(filename, fiberpos):
    '''
    Writes a fiberpos table to filename, maintaining backwards compatibility
    with the original fiberpos.txt format

    Args:
        filename: output file name string
        fiberpos: astropy Table of fiber positions
    '''

    #- Write the old text file format for backwards compatibility
    fxlines = [
        "#- Fiber to positioner mapping; x,y,z in mm on focal plane",
        "#- See doc/fiberpos.rst for more details.",
        "#- Coordinates at zenith: +x = East = +RA; +y = South = -dec",
        "",
        '#- fiber positioner spectro  x  y  z']
    for row in fiberpos:
        fxlines.append("{:4d}  {:4d}  {:2d}  {:12.6f}  {:12.6f}  {:12.6f}".format(
            row['fiber'], row['positioner'], row['spectro'],
            row['x'], row['y'], row['z'],
        ))

    fx = open(filename, 'w')
    fx.write('\n'.join(fxlines)+'\n')
    fx.close()


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
    import argparse
    import os
    import sys

    import numpy as np
    from astropy.table import Table, vstack

    parser = argparse.ArgumentParser(prog=sys.argv[0])
    parser.add_argument("-i", "--input", action='store', metavar='FILE',
        default='DESI-0530-posloc.txt', help="positioner location filename")
    parser.add_argument("-c", "--cassettes", action='store', metavar='FILE',
        default='cassette_order.txt', help="Order of cassettes on slit heads")
    parser.add_argument("-o", "--outdir", action='store', metavar='DIR',
        default='.', help="output directory")
    parser.add_argument("--debug",   help="debug with ipython prompt at end", action="store_true")

    opts = parser.parse_args()

    #- Random but reproducible
    seed = 2
    np.random.seed(seed)

    #- Read positioner locations on the focal plane
    colnames = ['positioner', 'postype', 'x', 'y', 'z', 'preces', 'nutat', 'spin', 'cassette']
    posloc = Table.read(opts.input, format='ascii', names=colnames)
    if 'comments' in posloc.meta.keys():
        del posloc.meta['comments']

    #- Read mapping of cassettes on focal plane to fibers on slithead
    colnames = ['fibermin', 'fibermax', 'sp0', 'sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8', 'sp9']
    cassettes = Table.read(opts.cassettes, format='ascii', names=colnames)

    #- Randomize fibers within a cassette
    petals = list()
    for p in range(10):
        fiberpos = posloc.copy(copy_data=True)
        fiberpos['fiber'] = -1
        fiberpos['spectro'] = p
        iipos = (fiberpos['postype'] == 'POS')
        fiberpos['positioner'] += p*len(fiberpos)
        for c in range(1,11):
            ii = (cassettes['sp'+str(p)] == c)
            assert np.count_nonzero(ii) == 1
            fibermin = p*500 + cassettes['fibermin'][ii][0]
            fibermax = p*500 + cassettes['fibermax'][ii][0]

            jj = iipos & (fiberpos['cassette'] == str(c))  #- includes 'N/A' so column is string not int
            assert np.count_nonzero(jj) == 50
            fiber = list(range(fibermin, fibermax+1))
            np.random.shuffle(fiber)
            fiberpos['fiber'][jj] = fiber

        #- Petal 0 is at the "bottom"; See DESI-0530
        phi = np.radians((7*36 + 36*p)%360)
        x = np.cos(phi)*fiberpos['x'] - np.sin(phi)*fiberpos['y']
        y = np.sin(phi)*fiberpos['x'] + np.cos(phi)*fiberpos['y']
        fiberpos['x'] = x
        fiberpos['y'] = y

        petals.append(fiberpos)

    fiberpos = vstack(petals)
    fiberpos.sort('fiber')
    ii = (fiberpos['postype'] == 'POS')

    #- Sanity checks before writing output
    fp = fiberpos[ii]
    assert len(fp) == 5000
    assert len(np.unique(fp['fiber'])) == 5000
    assert min(fp['fiber']) == 0
    assert max(fp['fiber']) == 4999
    assert len(set(fp['spectro'])) == 10
    assert min(fp['spectro']) == 0
    assert max(fp['spectro']) == 9
    assert len(np.unique(fiberpos['positioner'])) == len(fiberpos)

    #- Pick filenames in output directory
    textout = os.path.join(opts.outdir, 'fiberpos.txt')
    fitsout = os.path.join(opts.outdir, 'fiberpos.fits')
    pngout = os.path.join(opts.outdir, 'fiberpos.png')

    #- Drop some columns we don't need
    fiberpos.remove_columns(['preces', 'nutat', 'spin', 'cassette'])

    #- Update i8 -> i4 for integer columns
    for colname in ['fiber', 'positioner', 'spectro']:
        fiberpos.replace_column(colname, fiberpos[colname].astype('i4'))

    #- Set units and descriptions
    fiberpos['x'].unit = 'mm'
    fiberpos['y'].unit = 'mm'
    fiberpos['z'].unit = 'mm'
    fiberpos['x'].description = 'focal surface location [mm]'
    fiberpos['y'].description = 'focal surface location [mm]'
    fiberpos['z'].description = 'focal surface location [mm]'
    fiberpos['fiber'].description = 'fiber number [0-4999]'
    fiberpos['positioner'].description = 'focal plane positioner hole number'
    fiberpos['spectro'].description = 'spectrograph number [0-9]'
    fiberpos.meta['comments'] = ["Coordinates at zenith: +x = East = +RA; +y = South = -dec"]

    #- Write old text format with just fiber, positioner, spectro, x, y, z
    write_text_fiberpos(textout, fiberpos[ii])
    cols = ['fiber', 'positioner', 'spectro', 'x', 'y', 'z']
    fiberpos[ii][cols].write(fitsout, format='fits', overwrite=True)

    #- Write ecsv and fits format with all columns and rows, including
    #- fiducials (postype='FIF') and sky monitor (postype='ETC')
    fiberpos.sort('positioner')
    fiberpos.write(fitsout.replace('.fits', '-all.fits'), format='fits', overwrite=True)
    fiberpos.write(textout.replace('.txt', '-all.ecsv'), format='ascii.ecsv')

    #- Visualize mapping
    ii = (fiberpos['postype'] == 'POS')
    import pylab as P
    P.figure(figsize=(7,7))
    P.scatter(fiberpos['x'][ii], fiberpos['y'][ii], c=fiberpos['fiber'][ii]%500, edgecolor='none')
    P.grid()
    P.xlim(-420,420)
    P.ylim(-420,420)
    P.xlabel('x [mm]')
    P.ylabel('y [mm]')
    P.savefig(pngout, dpi=80)

    if opts.debug:
        #--- DEBUG ---
        import IPython
        IPython.embed()
        #--- DEBUG ---

    return 0

exit(main())
