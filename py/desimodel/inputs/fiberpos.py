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
    Update positioner to fiber number mapping from DocDB

    Options:
        outdir: string output directory, default $DESIMODEL/data/focalplane
        seed: integer random number seed for randomization within a cartridge

    Writes outdir/fiberpos*
    '''
    from desiutil.log import get_logger
    log = get_logger()

    #- Download input files from DocDB
    cassette_file = docdb.download(2721, 2, 'cassette_order.txt')
    xls_fp_layout = docdb.download(530, 10, 'DESI-0530-v10 (Focal Plane Layout).xlsx')

    #- Random but reproducible
    np.random.seed(seed)

    #- DESI-0530 file name (fn) and sheet name (sn) shortcuts
    fn = xls_fp_layout
    sn = 'PositionerAndFiducialLocations'

    #- Sanity check that columns are still in the same place
    headers = docdb.xls_read_row(fn, sn, 40, 'B', 'S')

    assert headers[0] == 'device_location_id'
    assert headers[1] == 'device_type'
    assert headers[2] == 'X'
    assert headers[3] == 'Y'
    assert headers[4] == 'Z'
    assert headers[8] == 'cassetteID'
    assert headers[15] == 'Q'
    assert headers[17] == 'S'

    #- Read Excel table with device locations
    posloc = Table()
    posloc['device'] = docdb.xls_read_col(fn, sn, 'B', 41, 583, dtype=int)
    posloc['device_type'] = docdb.xls_read_col(fn, sn, 'C', 41, 583, dtype=str)
    posloc['x'] = docdb.xls_read_col(fn, sn, 'D', 41, 583, dtype=float)
    posloc['y'] = docdb.xls_read_col(fn, sn, 'E', 41, 583, dtype=float)
    posloc['z'] = docdb.xls_read_col(fn, sn, 'F', 41, 583, dtype=float)

    #- Q and S not calculated for final to 'GIF' devices in 0530v10,
    #- thus requiring special handling
    q = docdb.xls_read_col(fn, sn, 'Q', 41, 583)
    s = docdb.xls_read_col(fn, sn, 'S', 41, 583)
    assert np.all(posloc['device_type'][-2:] == 'GIF')
    q[-2:] = s[-2:] = '0'
    posloc['q'] = np.array(q, dtype='float')
    posloc['s'] = np.array(s, dtype='float')

    #- Cassette N/A -> -1, and parse string -> float -> int
    c = docdb.xls_read_col(fn, sn, 'J', 41, 583)
    not_spectro_fiber = (c == 'N/A') | (c == '')
    c[not_spectro_fiber] = '-1'
    posloc['cassette'] = np.array(c, dtype=float).astype(int)

    #- Sanity check on values
    assert len(posloc) == 543  #- 543 holes have been drilled
    assert len(np.unique(posloc['device'])) == len(posloc['device'])
    assert set(posloc['device_type']) == set(['POS', 'FIF', 'GIF', 'NON', 'OPT', 'ETC'])
    assert 0 < np.min(posloc['x']) and np.max(posloc['x']) < 410
    assert 0 <= np.min(posloc['q']) and np.max(posloc['q']) < 36.0
    assert 0 <= np.min(posloc['s']) and np.max(posloc['s']) < 412.3
    assert np.all((posloc['s']**2 > posloc['x']**2 + posloc['y']**2 + posloc['z']**2) | (posloc['device_type'] == 'GIF'))
    assert np.min(posloc['cassette']) == -1
    assert np.max(posloc['cassette']) == 10
    assert 0 not in posloc['cassette']

    #- Read mapping of cassettes on focal plane to fibers on slithead
    colnames = ['fibermin', 'fibermax', 'sp0', 'sp1', 'sp2', 'sp3', 'sp4', 'sp5', 'sp6', 'sp7', 'sp8', 'sp9']
    cassettes = Table.read(cassette_file, format='ascii', names=colnames)

    #- Randomize fibers within a cassette
    petals = list()
    for p in range(10):
        fiberpos = posloc.copy(copy_data=True)
        fiberpos['fiber'] = -1
        fiberpos['petal'] = p
        fiberpos['slit'] = p
        fiberpos['spectro'] = p
        iipos = (fiberpos['device_type'] == 'POS')
        fiberpos['device'] += p*len(fiberpos)
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

        #- Additional columns
        fiberpos['slitblock'] = ['B{}'.format(tmp) for tmp in (fiberpos['fiber'] % 500) // 25]
        fiberpos['blockfiber'] = ['F{}'.format(tmp) for tmp in (fiberpos['fiber'] % 500) % 25]

        #- Petal 0 is at the "bottom"; See DESI-0530
        phi = np.radians((7*36 + 36*p)%360)
        x = np.cos(phi)*fiberpos['x'] - np.sin(phi)*fiberpos['y']
        y = np.sin(phi)*fiberpos['x'] + np.cos(phi)*fiberpos['y']
        fiberpos['x'] = x
        fiberpos['y'] = y

        petals.append(fiberpos)

    fiberpos = vstack(petals)
    fiberpos.sort('fiber')
    ii = (fiberpos['device_type'] == 'POS')

    #- devices that don't go to spectrographs don't have slitblock, blockfiber
    fiberpos['slitblock'][~ii] = '--'
    fiberpos['blockfiber'][~ii] = '--'

    #- More sanity checks before writing output
    fp = fiberpos[ii]
    assert len(fp) == 5000
    assert len(np.unique(fp['fiber'])) == 5000
    assert min(fp['fiber']) == 0
    assert max(fp['fiber']) == 4999
    assert len(set(fp['spectro'])) == 10
    assert min(fp['spectro']) == 0
    assert max(fp['spectro']) == 9
    assert len(np.unique(fiberpos['device'])) == len(fiberpos)

    #- Pick filenames in output directory
    if outdir is None:
        outdir = os.path.join(datadir(), 'focalplane')

    textout = os.path.join(outdir, 'fiberpos.txt')
    fitsout = os.path.join(outdir, 'fiberpos.fits')
    pngout = os.path.join(outdir, 'fiberpos.png')

    #- Drop some columns we don't need
    fiberpos.remove_column('cassette')

    #- Update i8 -> i4 for integer columns
    for colname in ['fiber', 'device', 'spectro', 'petal', 'slit']:
        fiberpos.replace_column(colname, fiberpos[colname].astype('i4'))

    #- Set units and descriptions; see DESI-2724
    fiberpos['x'].unit = 'mm'
    fiberpos['y'].unit = 'mm'
    fiberpos['z'].unit = 'mm'
    fiberpos['q'].unit = 'deg'
    fiberpos['s'].unit = 'mm'
    fiberpos['x'].description = 'focal surface location [mm]'
    fiberpos['y'].description = 'focal surface location [mm]'
    fiberpos['z'].description = 'focal surface location [mm]'
    fiberpos['q'].description = 'azimuthal angle on focal surface [deg]'
    fiberpos['s'].description = 'radial distance along focal surface [mm]'
    fiberpos['fiber'].description = 'fiber number [0-4999]'
    fiberpos['device'].description = 'focal plane device number i%543 in [0-542]'
    fiberpos['spectro'].description = 'spectrograph number [0-9]'
    fiberpos['petal'].description = 'focal plane petal number [0-9]'
    fiberpos['slit'].description = 'spectrograph slit number [0-9]'
    fiberpos['slitblock'].description = 'id of the slitblock on the slit [B0-B19]'
    fiberpos['blockfiber'].description = 'id of the fiber on the slitblock [F0-F24]'
    fiberpos.meta['comments'] = ["Coordinates at zenith: +x = East = +RA; +y = South = -dec"]

    #- Write old text format with just fiber, device, spectro, x, y, z
    write_text_fiberpos(textout, fiberpos[ii])
    log.info('Wrote {}'.format(textout))
    cols = ['fiber', 'device', 'spectro', 'x', 'y', 'z']
    fiberpos[ii][cols].write(fitsout, format='fits', overwrite=True)
    log.info('Wrote {}'.format(fitsout))

    #- Write ecsv and fits format with all columns and rows, including
    #- fiducials (device_type='FIF') and sky monitor (device_type='ETC')
    fiberpos.sort('device')
    fitsallout = fitsout.replace('.fits', '-all.fits')
    textallout = textout.replace('.txt', '-all.ecsv')
    fiberpos.write(fitsallout, format='fits', overwrite=True)
    fiberpos.write(textallout, format='ascii.ecsv')
    log.info('Wrote {}'.format(fitsallout))
    log.info('Wrote {}'.format(textallout))

    #- Visualize mapping
    POS = (fiberpos['device_type'] == 'POS')
    FIF = (fiberpos['device_type'] == 'FIF')
    ETC = (fiberpos['device_type'] == 'ETC')
    import pylab as P
    P.jet()     #- With apologies to viridis
    P.figure(figsize=(7,7))
    P.scatter(fiberpos['x'][POS], fiberpos['y'][POS], c=fiberpos['fiber'][POS]%500, edgecolor='none', s=20)
    # P.scatter(fiberpos['x'][FIF], fiberpos['y'][FIF], s=5, color='k')
    # P.plot(fiberpos['x'][ETC], fiberpos['y'][ETC], 'kx', ms=3)
    P.grid(alpha=0.2, color='k')
    P.xlim(-420,420)
    P.ylim(-420,420)
    P.xlabel('x [mm]')
    P.ylabel('y [mm]')
    P.savefig(pngout, dpi=80)
    log.info('Wrote {}'.format(pngout))

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
        "#- THIS FILE IS PROVIDED FOR BACKWARDS COMPATIBILITY",
        "#- Please use fiberpos-all.[ecsv,fits] for additional columns",
        '#- and non-spectrofiber device hole locations.',
        '#-'
        "#- Fiber to focal plane device hole mapping; x,y,z in mm on focal plane",
        "#- See doc/fiberpos.rst and DESI-0530 for more details.",
        "#- Coordinates at zenith: +x = East = +RA; +y = South = -dec",
        "",
        "#- fiber=at spectrograph; device=numbering on focal plane",
        "",
        '#- fiber device spectro  x  y  z']
    for row in fiberpos:
        fxlines.append("{:4d}  {:4d}  {:2d}  {:12.6f}  {:12.6f}  {:12.6f}".format(
            row['fiber'], row['device'], row['spectro'],
            row['x'], row['y'], row['z'],
        ))

    with open(filename, 'w') as fx:
        fx.write('\n'.join(fxlines)+'\n')
