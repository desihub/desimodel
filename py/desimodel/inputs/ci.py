# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
'''
desimodel.inputs.ci
===================

Utilities for updating Commissioning Instrument corners in
data/focalplane/ci-corners.ecsv.
'''
import os

import numpy as np
from astropy.table import Table

from . import docdb
from ..io import datadir
from ..focalplane import xy2qs

def update():
    '''Update $DESIMODEL/data/focalplane/ci-corners.ecsv from DESI-4633v11'''

    from desiutil.log import get_logger
    log = get_logger()

    #- Get as-build measured CI corners from DESI-4633v11 Corners.txt
    docnum = 4633
    docver = 11
    corners_file = docdb.download(docnum, docver, 'Corners.txt', overwrite=True)

    #- Read file and build dictionary corners[camera][x/y/z] = list len 4
    with open(corners_file) as fx:
        corners = dict()
        for line in fx:
            line = line.strip()
            if line.startswith('#') and len(line)>2 and line[-2] == 'C':
                camera = line[-2:]
                corners[camera] = dict(x=list(), y=list(), z=list())
            if line.startswith('#') or len(line)<1:
                continue
            x, y, z = map(float, line.split(','))
            corners[camera]['x'].append(x)
            corners[camera]['y'].append(y)
            corners[camera]['z'].append(z)

    #- Convert that into the GFA-corners table format
    citable = Table(
        names = ('GFA_LOC', 'CORNER', 'X', 'Y', 'Z', 'Q', 'S', 'RADIUS_DEG'),
        dtype = ('int','int','float','float','float','float','float','float')
        )

    for icam in range(5):
        cam = 'C'+str(icam+1)
        for corner in range(4):
            x = corners[cam]['x'][corner]
            y = corners[cam]['y'][corner]
            z = corners[cam]['z'][corner]
            r = np.sqrt(x**2 + y**2)
            q, s = xy2qs(x, y)
            citable.add_row([icam+1, corner, x, y, z, q, s, r])

    outfile = os.path.join(datadir(), 'focalplane', 'ci-corners.ecsv')
    citable.meta['VERSION'] = 'DESI-{}v{}'.format(docnum, docver)
    citable.write(outfile, overwrite=True)
    log.info('Wrote {}'.format(outfile))
