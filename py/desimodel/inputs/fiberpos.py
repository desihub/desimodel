'''
Utilities for updating positioner to fiber mapping
'''
import os

import numpy as np
from astropy.table import Table, vstack

from . import docdb
from ..io import datadir

def update(outdir=None, seed=2):
    '''
    TODO: document
    '''
    #- Download input files from DocDB
    cassette_file = docdb.download(2721, 1, 'cassette_order.txt')
    xls_fp_layout = docdb.download(530, 10, 'DESI-0530-v10 (Focal Plane Layout).xlsx')
    
    #- Random but reproducible
    np.random.seed(seed)

    fn = xls_fp_layout                          #- File Name
    sn = 'PositionerAndFiducialLocations'       #- Sheet Name
    posloc = Table()
    posloc['positioner'] = docdb.xls_read_col(fn, sn, 'B', 41, 581, dtype=int)
    posloc['postype'] = docdb.xls_read_col(fn, sn, 'C', 41, 581, dtype=str)
    posloc['x'] = docdb.xls_read_col(fn, sn, 'D', 41, 581, dtype=float)
    posloc['y'] = docdb.xls_read_col(fn, sn, 'E', 41, 581, dtype=float)
    posloc['z'] = docdb.xls_read_col(fn, sn, 'F', 41, 581, dtype=float)
    
    #- Cassette N/A -> -1, and parse string -> float -> int
    c = docdb.xls_read_col(fn, sn, 'J', 41, 581)
    c = [x.replace('N/A', '-1') for x in c]
    posloc['cassette'] = np.array(c, dtype=float).astype(int)

    #- Read mapping of cassettes on focal plane to fibers on slithead
    colnames = ['fibermin', 'fibermax', 'sp0', 'sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8', 'sp9']
    cassettes = Table.read(cassette_file, format='ascii', names=colnames)

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

            jj = iipos & (fiberpos['cassette'] == c)
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
    if outdir is None:
        outdir = os.path.join(datadir(), 'focalplane')

    textout = os.path.join(outdir, 'fiberpos.txt')
    fitsout = os.path.join(outdir, 'fiberpos.fits')
    pngout = os.path.join(outdir, 'fiberpos.png')

    #- Drop some columns we don't need
    fiberpos.remove_column('cassette')

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

    return 0
    
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

    
