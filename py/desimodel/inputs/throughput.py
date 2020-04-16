# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
'''
desimodel.inputs.throughput
===========================

Utilities for updating throughput model.
'''
import os
import shutil

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from astropy.table import Table, vstack
from astropy.io import fits
import yaml

from . import docdb
from ..io import datadir, findfile

def update(testdir=None, desi347_version=16, desi5501_version=3, desi5501_KOSI=True):
    '''
    Update thru-\*.fits from DESI-0347 and DESI-0344

    Args:
        testdir: If not None, write files here instead of standard locations
            under $DESIMODEL/data/
        desi347_version: version of DESI-347 to use
        desi334_version: version of DESI-334 to use
    '''
    from desiutil.log import get_logger
    log = get_logger()

    master_thru_file = docdb.download(
        347, desi347_version,
        'DESI-347-v{}_Throughput-Noise-SNR-Calcs.xlsx'.format(desi347_version))
    desi_yaml_file   = docdb.download(
        347, desi347_version, 'desi_v{}.yaml'.format(desi347_version))

    ccd_thru_file = []
    suffix = 'KOSI' if desi5501_KOSI else ''
    for spectro in range(10):
        ccd_thru_file.append(docdb.download(
            5501, desi5501_version, 'Spectrograph{0}{1}.xlsx'.format(spectro + 1, suffix)))

    with open(desi_yaml_file) as fx:
        params = yaml.safe_load(fx)

    if testdir is None:
        outfile_desiyaml = os.path.join(datadir(), 'desi.yaml')
        thrudir = os.path.join(datadir(), 'throughput')
    elif os.path.isdir(testdir):
        outfile_desiyaml = os.path.join(testdir, 'desi.yaml')
        thrudir = testdir
    else:
        raise ValueError("Missing directory {}".format(testdir))

    shutil.copy(desi_yaml_file, outfile_desiyaml)
    log.info('Wrote {}'.format(outfile_desiyaml))

    #- Telescope geometric area m^2 -> cm^2
    params['area']['geometric_area'] *= 100**2

    #- Load atmospheric extinction
    d = fits.getdata(findfile('inputs/throughput/ZenithExtinction-KPNO.fits'), 'EXTINCTION')
    extinction = InterpolatedUnivariateSpline(d['WAVELENGTH'], d['EXTINCTION'])

    #- Load pre-spectrograph throughputs
    thru, (xlswave, xlstotthru, xlsspecthru) = load_throughput(master_thru_file)

    #- Load pre-computed fiberloss for reference objects
    fiberinput = dict()
    for objtype in ['elg', 'lrg', 'sky', 'star']:
        fiberinput[objtype] = load_fiberinput(
            findfile('throughput/fiberloss-{}.dat'.format(objtype)) )

    #- Read spectrograph throughputs and calculate mean over spectrographs.
    specwave, specthru_data = load_spec_throughputs(ccd_thru_file)
    specthru_mean = np.mean(specthru_data, axis=0)

    #- Min/Max wavelength coverage
    wmin = dict()
    wmax = dict()

    for j, channel in enumerate('brz'):
        # Build an interpolation the DESI-5501 mean spectrograph throughput in Angstrom units.
        specthru = InterpolatedUnivariateSpline(10 * specwave, specthru_mean[j], k=3)
        wmin[channel], wmax[channel] = get_waveminmax(findfile('specpsf/psf-{}.fits'.format(channel)))

        dw = 0.1
        ww = np.arange(wmin[channel], wmax[channel]+dw/2, dw)
        tt = thru(ww) * specthru(ww)

        data = np.rec.fromarrays([ww, tt, extinction(ww), fiberinput['elg'](ww)],
                                names='wavelength,throughput,extinction,fiberinput')

        hdr = fits.Header()
        hdr['EXPTIME']  = (params['exptime_dark'], 'default exposure time [sec]')
        hdr['GEOMAREA'] = (params['area']['geometric_area'], 'geometric area of mirror - obscurations')
        hdr['FIBERDIA'] = (params['fibers']['diameter_arcsec'], 'average fiber diameter [arcsec]')
        hdr['WAVEMIN']  = (wmin[channel], 'Minimum wavelength [Angstroms]')
        hdr['WAVEMAX']  = (wmax[channel], 'Maximum wavelength [Angstroms]')

        fiberinput_data = np.rec.fromarrays([ww, fiberinput['elg'](ww), fiberinput['lrg'](ww),
                                      fiberinput['star'](ww), fiberinput['sky'](ww)],
                                names='wavelength,elg,lrg,star,sky')

        outfile = thrudir + '/thru-{0}.fits'.format(channel)

        hdus = fits.HDUList()
        hdus.append(fits.PrimaryHDU())
        hdus.append(fits.BinTableHDU(data, hdr, name='THROUGHPUT'))
        hdus.append(fits.BinTableHDU(fiberinput_data, name='FIBERINPUT'))
        hdus.writeto(outfile, overwrite=True)
        log.info('Wrote {}'.format(outfile))

def load_throughput(filename, specthru_row=95, thru_row=97):
    """
    Load throughputs from DESI-0347, removing the spectrograph contributions
    which will be loaded separately from higher resolution data.

    Args:
        filename:
            DESI-0347 Excel file location

    Returns (thruspine, xlsdata), where

        thruspline: InterpolatedUnivariateSpline of thru vs. wave[Angstroms]
        xlsdata: tuple of (wave, totalthru, specthru)

    Notes:

      * Alas, DESI-0347 doesn't fill in the final throughput for
        3500 and 9950 Angstroms, even though the inputs are there.
    """
    wave = docdb.xls_read_row(filename, 'Throughput', 3, 'C', 'P')*10

    rowlabel = docdb.xls_read_row(filename, 'Throughput', specthru_row, 'A', 'A')[0]
    assert rowlabel.startswith('Throughput:  Spectrograph'), 'Has the spectrograph throughput row moved?'
    specthru = docdb.xls_read_row(filename, 'Throughput', specthru_row, 'C', 'P')

    rowlabel = docdb.xls_read_row(filename, 'Throughput', thru_row, 'A', 'A')[0]
    assert rowlabel.startswith('sky throughput:'), 'Has the sky throughput row moved?'
    thru = docdb.xls_read_row(filename, 'Throughput', thru_row, 'C', 'P')

    assert len(wave) == 14
    assert len(wave) == len(thru)
    assert len(wave) == len(specthru)
    assert np.all(np.diff(wave)>0)
    assert np.min(wave) == 3600 and np.max(wave) == 9800
    assert 0 <= np.min(thru) and np.max(thru) <= 1
    assert 0 <= np.min(specthru) and np.max(specthru) <= 1

    return InterpolatedUnivariateSpline(wave, thru/specthru, k=3), (wave, thru, specthru)

def load_fiberinput(filename):
    """
    Load fiberinput as calculated by fiberloss.py

    Args:
        filename: fiberloss input file,
            e.g. $DESIMODEL/data/throughput/fiberloss-elg.dat

    Returns InterpolatedUnivariateSpline instance.
    """
    tmp = np.loadtxt(filename).T
    wavelength = tmp[0]  #- nm -> Angstroms
    throughput = tmp[1]

    return InterpolatedUnivariateSpline(wavelength, throughput, k=3)

def load_spec_throughputs(filenames, columns='ABCD', first_row=2, last_row=647):
    """
    Loads spectrograph*CCD throughputs from DESI-5501 excel files.

    Args:
        filenames: list of per-spectrograph filenames.

    Returns arrays of wavelength in nm and throughput per spectrograph.
    """
    wave = None
    nspectro = len(filenames)
    for spectro in range(nspectro):
        fname = filenames[spectro]
        # Sanity check on column headers.
        headers = docdb.xls_read_row(fname, 'Summary', first_row - 1, columns[0], columns[-1])
        assert np.array_equal(headers, ['Wav', 'Blue', 'Red', 'NIR']), 'Unexpected column headers.'
        wave_in = docdb.xls_read_col(fname, 'Summary', columns[0], first_row, last_row, dtype=float)
        if wave is None:
            wave = wave_in
            assert np.allclose(np.diff(wave), 1), 'Unexpected wavelength grid.'
            nwave = len(wave)
            thru = np.zeros((nspectro, 3, nwave))
        else:
            assert np.array_equal(wave, wave_in), 'Wavelength arrays do not match.'
        for k, (camera, column) in enumerate(zip('brz', columns[1:])):
            # Read column values as strings since some cells are empty.
            values = docdb.xls_read_col(fname, 'Summary', column, first_row, last_row)
            # Replace empty cells with zeros and convert non-empty cells to floats.
            values[values == ''] = '0'
            thru[spectro, k] = values.astype(float)

    return wave, thru

def get_waveminmax(psffile):
    """
    return wmin, wmax as taken from the header of a PSF file,
    e.g. $DESIMODEL/data/specpsf/psf-b.fits
    """
    hdr = fits.getheader(psffile)
    return hdr['WAVEMIN'], hdr['WAVEMAX']
