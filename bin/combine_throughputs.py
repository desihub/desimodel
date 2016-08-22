#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np
import fitsio

#
#
#
def load_throughput(filename):
    """
    Load throughputs copied and pasted from summary throughput document
    DESI-0347v2.  First row is wavelengths in nanometers; subsequent rows
    are throughput components which should be multiplied together and then
    interpolated with a cubic spline.

    Returns InterpolatedUnivariateSpline instance.
    """
    tmp = np.loadtxt(filename).T
    wavelength = tmp[0]*10  #- nm -> Angstroms

    #- Throughput, removing spectrograph*CCD contribution
    throughput = tmp[1] / tmp[2]

    return InterpolatedUnivariateSpline(wavelength, throughput, k=3)
#
#
#
def load_fiberinput(filename):
    """
    Load fiberinput as calculated by fiberloss.py

    Returns InterpolatedUnivariateSpline instance.
    """
    tmp = np.loadtxt(filename).T
    wavelength = tmp[0]  #- nm -> Angstroms
    throughput = tmp[1]

    return InterpolatedUnivariateSpline(wavelength, throughput, k=3)
#
#
#
def load_spec_throughput(filename):
    """
    Spectrograph throughputs from DESI-0334 have wavelength [nm] in the
    first column and total throughput in the last column.  However, the NIR
    throughput has the total in the second-to-the-last column, thus
    the support for the column= keyword.

    Returns InterpolatedUnivariateSpline instance.
    """
    tmp = np.loadtxt(filename)
    wavelength = tmp[:, 0] * 10  #- nm -> Angstroms
    throughput = tmp[:, -1]
    return InterpolatedUnivariateSpline(wavelength, throughput, k=3)
#
#
#
def get_waveminmax(psffile):
    """
    return wmin, wmax as taken from the header of a PSF file
    """
    hdr = fitsio.read_header(psffile)
    return hdr['WAVEMIN'], hdr['WAVEMAX']
#
#-------------------------------------------------------------------------
#
def main():
    """Combine the various throughput parameters and generate a Specter
    compatible throughput fits file.
    """
    from argparse import ArgumentParser
    import yaml
    import os

    parser = ArgumentParser(description="Create a specter-compatible throughput FITS file.",
                            prog=sys.argv[0])
    parser.add_argument("-m", "--modeldir", action='store', metavar='DIR',
                        help="model directory")
    parser.add_argument("-o", "--outdir", action='store', metavar='DIR',
                        help="output directory")
    opts = parser.parse_args()

    if opts.modeldir is None:
        if 'DESIMODEL' in os.environ:
            datadir = os.path.join(os.environ['DESIMODEL'], 'data')
        else:
            print('ERROR: you must give --modeldir or set DESIMODEL',
                  file=sys.stderr)
            return 1
    else:
        datadir = os.path.join(opts.modeldir, 'data')

    if opts.outdir is None:
        opts.outdir = '.'

    #- Load telescope parameters
    fx = open(os.path.join(datadir, "desi.yaml"))
    params = yaml.load(fx)
    fx.close()

    #- Telescope geometric area m^2 -> cm^2
    params['area']['geometric_area'] *= 100**2

    #- Load atmospheric extinction
    d = fitsio.read(os.path.join(datadir, 'inputs', 'throughput',
                                 'ZenithExtinction-KPNO.fits'), 'EXTINCTION')
    extinction = InterpolatedUnivariateSpline(d['WAVELENGTH'], d['EXTINCTION'])

    #- Load pre-spectrograph throughputs
    thru = load_throughput(os.path.join(datadir, 'inputs', 'throughput',
                                        'DESI-0347-throughput.txt'))
    fiberinput = dict()
    for objtype in ['elg', 'lrg', 'sky', 'star']:
        fiberinput[objtype] = load_fiberinput(os.path.join(datadir,
                                                           'throughput',
                                                           'fiberloss-{}.dat'.format(objtype)))

    #- Spectrograph throughputs
    spec = dict()
    spec['b'] = load_spec_throughput(datadir+'/inputs/throughput/DESI-0334-spectro/blue-thru.txt')
    spec['r'] = load_spec_throughput(datadir+'/inputs/throughput/DESI-0334-spectro/red-thru.txt')
    spec['z'] = load_spec_throughput(datadir+'/inputs/throughput/DESI-0334-spectro/nir-thru-500.txt')

    #- Load wavemin / wavemax from PSF files
    wmin = dict()
    wmax = dict()
    wmin['b'], wmax['b'] = get_waveminmax(datadir+'/specpsf/psf-b.fits')
    wmin['r'], wmax['r'] = get_waveminmax(datadir+'/specpsf/psf-r.fits')
    wmin['z'], wmax['z'] = get_waveminmax(datadir+'/specpsf/psf-z.fits')

    for channel in ('b', 'r', 'z'):
        dw = 0.1
        ww = np.arange(wmin[channel], wmax[channel]+dw/2, dw)
        tt = thru(ww) * spec[channel](ww)

        data = dict(wavelength=ww, throughput=tt,
                    extinction=extinction(ww),
                    fiberinput=fiberinput['elg'](ww) )
        hdr = list()
        hdr.append(dict(name='EXPTIME', value=params['exptime'], comment='default exposure time [sec]'))
        hdr.append(dict(name='GEOMAREA', value=params['area']['geometric_area'], comment='geometric area of mirror - obscurations'))
        hdr.append(dict(name='FIBERDIA', value=params['fibers']['diameter_arcsec'], comment='average fiber diameter [arcsec]'))
        hdr.append(dict(name='WAVEMIN', value=wmin[channel], comment='Minimum wavelength [Angstroms]'))
        hdr.append(dict(name='WAVEMAX', value=wmax[channel], comment='Maximum wavelength [Angstroms]'))

        outfile = opts.outdir + '/thru-{0}.fits'.format(channel)
        fitsio.write(outfile, data, header=hdr, clobber=True, extname='THROUGHPUT')

        #- Write another header with fiberinput for multiple object types
        data = np.rec.fromarrays([ww, fiberinput['elg'](ww), fiberinput['lrg'](ww),
                                      fiberinput['star'](ww), fiberinput['sky'](ww)],
                                names='wavelength,elg,lrg,star,sky')
        fitsio.write(outfile, data, extname='FIBERINPUT')

    return 0
#
#
#
if __name__ == '__main__':
    sys.exit(main())
