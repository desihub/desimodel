# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.inputs.gfa
====================

Utilities for updating GFA data.
"""
import os.path
import numpy as np
from ..focalplane import get_radius_deg, get_radius_mm, xy2qs
from ..io import datadir

def build_gfa_table(testdir=None):
    '''
    Builds the GFA table given the data from DESI-0530-v14 Excel spreadsheet
    and writes an .ecsv file using the astropy table library.

    The GFA corners come from the P1-P4 "Reference projection of active area"
    in rows 26-29 of columns C-E of DESI-0530-v14.

    The saved table has columns PETAL,CORNER,X,Y,Z,Q,S,RADIUS_DEG.

    Args:
        testdir: If not None, write files here instead of standard locations
            under $DESIMODEL/data/
    '''
    from desiutil.log import get_logger
    log = get_logger()

    if testdir is None:
        outdir = os.path.join(datadir(), 'focalplane')
    else:
        outdir = testdir

    if not os.path.isdir(outdir):
        raise ValueError("Missing directory {}".format(testdir))

    # Uses the reference projection of active area to create data table of GFAs
    from astropy.table import Table
    # Initial x and y coordinates for the GFAs
    # Data obtained from DESI-0530-v14 Excel spreadsheet, which is for petal 3
    x = [318.703, 331.075, 349.121, 336.748]
    y = [225.816, 234.805, 209.944, 200.955]
    z = [-17.053, -18.487, -18.631, -17.198]

    def get_rotatemat(angle):
        '''
        Return a rotation matrix for angle in degrees
        '''
        a = np.radians(angle)
        rotatemat = np.zeros(shape=(2,2))
        rotatemat[0] = [np.cos(a), -np.sin(a)]
        rotatemat[1] = [np.sin(a), np.cos(a)]
        return rotatemat

    # Note: the corners are 0 indexed
    gfatable = Table(names = ('PETAL', 'CORNER', 'X', 'Y', 'Z', 'Q', 'S', 'RADIUS_DEG'),
                 dtype = ('int', 'int', 'float', 'float', 'float', 'float', 'float', 'float'))
    # Sets the units for the GFA table
    gfatable['X'].unit = 'mm'
    gfatable['Y'].unit = 'mm'
    gfatable['Z'].unit = 'mm'
    gfatable['Q'].unit = 'degrees'
    gfatable['S'].unit = 'mm'
    gfatable['RADIUS_DEG'].unit = 'degrees'

    for petal in range(10):
        #- x,y,z from DESI-0530 starts at petal 3
        angle = petal*36 - 3*36
        R = get_rotatemat(angle)
        xx, yy = R.dot([x,y])
        assert len(xx) == len(x)
        assert len(yy) == len(y)

        q, s = xy2qs(xx, yy)
        r_deg = get_radius_deg(xx, yy)

        for i in range(len(xx)):
            gfatable.add_row([petal, i, xx[i], yy[i], z[i], q[i], s[i], r_deg[i]])

    # Saves the table of data as an ecsv file
    outfile = os.path.join(outdir, 'gfa.ecsv')
    gfatable.write(outfile, format='ascii.ecsv')
    log.info('Wrote {}'.format(outfile))
