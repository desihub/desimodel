'''
Utilities for updating throughput model
'''
import os

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.table import Table, vstack
from astropy.io import fits
import yaml

from . import docdb
from ..io import datadir, findfile

def update(outdir=None):
    
    master_thru_file = docdb.download(347, 11, 'DESI-347-v11 Throughput Noise SNR Calcs.xlsx')
    desi_yaml_file = docdb.download(347, 11, 'desi.yaml')
    
    ccd_thru_file = dict()
    ccd_thru_file['b'] = docdb.download(334, 3, 'blue-thru.txt')
    ccd_thru_file['r'] = docdb.download(334, 3, 'red-thru.txt')
    ccd_thru_file['z'] = docdb.download(334, 3, 'nir-thru-250.txt')
    
    with open(desi_yaml_file) as fx:
        params = yaml.load(fx)

    #- Telescope geometric area m^2 -> cm^2
    params['area']['geometric_area'] *= 100**2

    #- Load atmospheric extinction
    d = fits.getdata(findfile('inputs/throughput/ZenithExtinction-KPNO.fits'), 'EXTINCTION')
    extinction = InterpolatedUnivariateSpline(d['WAVELENGTH'], d['EXTINCTION'])

    #- Load pre-spectrograph throughputs
    thru = load_throughput(master_thru_file)
    
    #- Load pre-computed fiberloss for reference objects
    fiberinput = dict()
    for objtype in ['elg', 'lrg', 'sky', 'star']:        
        fiberinput[objtype] = load_fiberinput(
            findfile('throughput/fiberloss-{}.dat'.format(objtype)) )

    #- Spectrograph throughputs
    specthru = dict()

    #- Min/Max wavelength coverage
    wmin = dict()
    wmax = dict()

    for channel in ('b', 'r', 'z'):
        specthru[channel] = load_spec_throughput(ccd_thru_file[channel])
        wmin[channel], wmax[channel] = get_waveminmax(findfile('specpsf/psf-{}.fits'.format(channel)))

        dw = 0.1
        ww = np.arange(wmin[channel], wmax[channel]+dw/2, dw)
        tt = thru(ww) * specthru[channel](ww)

        data = dict(wavelength=ww, throughput=tt,
                    extinction=extinction(ww),
                    fiberinput=fiberinput['elg'](ww) )
        hdr = list()
        hdr.append(dict(name='EXPTIME', value=params['exptime_dark'], comment='default exposure time [sec]'))
        hdr.append(dict(name='GEOMAREA', value=params['area']['geometric_area'], comment='geometric area of mirror - obscurations'))
        hdr.append(dict(name='FIBERDIA', value=params['fibers']['diameter_arcsec'], comment='average fiber diameter [arcsec]'))
        hdr.append(dict(name='WAVEMIN', value=wmin[channel], comment='Minimum wavelength [Angstroms]'))
        hdr.append(dict(name='WAVEMAX', value=wmax[channel], comment='Maximum wavelength [Angstroms]'))

        if outdir is None:
            outdir = os.path.join(datadir(), 'throughput')

        import fitsio
        outfile = outdir + '/thru-{0}.fits'.format(channel)
        fitsio.write(outfile, data, header=hdr, clobber=True, extname='THROUGHPUT')

        #- Write another header with fiberinput for multiple object types
        data = np.rec.fromarrays([ww, fiberinput['elg'](ww), fiberinput['lrg'](ww),
                                      fiberinput['star'](ww), fiberinput['sky'](ww)],
                                names='wavelength,elg,lrg,star,sky')
        fitsio.write(outfile, data, extname='FIBERINPUT')



def load_throughput(filename):
    """
    Load throughputs from DESI-0347, removing the spectrograph contributions
    which will be loaded separately from higher resolution data.

    Returns InterpolatedUnivariateSpline instance of thru vs. wave[Angstroms]
    """
    wave = docdb.xls_read_row(filename, 'Throughput', 3, 'C', 'P')*10
    thru = docdb.xls_read_row(filename, 'Throughput', 112, 'C', 'P')
    specthru = docdb.xls_read_row(filename, 'Throughput', 93, 'C', 'P')
    
    assert len(wave) == 14
    assert len(wave) == len(thru)
    assert len(wave) == len(specthru)

    return InterpolatedUnivariateSpline(wave, thru/specthru, k=3)

def load_fiberinput(filename):
    """
    Load fiberinput as calculated by fiberloss.py

    Returns InterpolatedUnivariateSpline instance.
    """
    tmp = np.loadtxt(filename).T
    wavelength = tmp[0]  #- nm -> Angstroms
    throughput = tmp[1]

    return InterpolatedUnivariateSpline(wavelength, throughput, k=3)

def load_spec_throughput(filename):
    """
    Spectrograph throughputs from DESI-0334 have wavelength [nm] in the
    first column and total throughput in the last column.

    Returns InterpolatedUnivariateSpline instance.
    """
    tmp = np.loadtxt(filename)
    wavelength = tmp[:, 0] * 10  #- nm -> Angstroms
    throughput = tmp[:, -1]
    return InterpolatedUnivariateSpline(wavelength, throughput, k=3)

def get_waveminmax(psffile):
    """
    return wmin, wmax as taken from the header of a PSF file
    """
    hdr = fits.getheader(psffile)
    return hdr['WAVEMIN'], hdr['WAVEMAX']
