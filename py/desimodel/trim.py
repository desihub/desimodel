# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.trim
==============

Code for trimming desimodel/data into smaller files.
"""
import os
import shutil
import numpy as np
from astropy.table import Table
from astropy.io import fits

from .footprint import pixweight

def trim_data(indir, outdir, overwrite=False):
    '''Trim a ``$DESIMODEL/data`` directory into a lightweight version for
    testing.

    Args:
        indir (str): A ``$DESIMODEL/data`` directory from svn.
        outdir (str): Output data directory location.
        overwrite (bool, optional):  If ``True``, remove `outdir` if it already
            exists.
    '''
    assert os.path.abspath(indir) != os.path.abspath(outdir)
    if os.path.exists(outdir) and overwrite:
        shutil.rmtree(outdir)

    os.makedirs(outdir)

    #- python note:
    #- *inout(indir, outdir, filename) -> indir/filename, outdir/filename

    shutil.copy(*inout(indir, outdir, 'desi.yaml'))
    trim_focalplane(*inout(indir, outdir, 'focalplane'))
    trim_footprint(*inout(indir, outdir, 'footprint'))
    trim_inputs(*inout(indir, outdir, 'inputs'))
    trim_sky(*inout(indir, outdir, 'sky'))
    trim_specpsf(*inout(indir, outdir, 'specpsf'))
    trim_spectra(*inout(indir, outdir, 'spectra'))
    trim_targets(*inout(indir, outdir, 'targets'))
    trim_throughput(*inout(indir, outdir, 'throughput'))
    trim_weather(*inout(indir, outdir, 'weather'))


def inout(indir, outdir, filename):
    '''returns os.path.join(indir, filename) and .join(outdir, filename)'''
    infile = os.path.join(indir, filename)
    outfile = os.path.join(outdir, filename)
    return infile, outfile


def trim_focalplane(indir, outdir):
    '''copy everything in focalplane'''
    assert os.path.basename(indir) == 'focalplane'
    shutil.copytree(indir, outdir)


def trim_footprint(indir, outdir):
    '''Copies subset of desi-tiles.fits and .ecsv but not .par. Also
    creates a corresponding version of desi-healpix-weights.fits.
    '''
    assert os.path.basename(indir) == 'footprint'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    infile, outfile = inout(indir, outdir, 'desi-tiles.fits')
    with fits.open(infile) as hdulist:
        t = Table(hdulist[1].data)
    #- Pick a subset on edge of footprint where healpix testing is done
    ii = (35 < t['RA']) & (t['RA'] < 55) & (-10 < t['DEC']) & (t['DEC'] < 20)
    tx = t[ii]
    tx.write(outfile, format='fits')

    #- Remove multidimensional columns before writing ecsv
    for colname in ['BRIGHTRA', 'BRIGHTDEC', 'BRIGHTVTMAG']:
        if colname in tx.colnames:
            tx.remove_column(colname)

    tx.write(outfile.replace('.fits', '.ecsv'), format='ascii.ecsv')
    infile, outfile = inout(indir, outdir, 'desi-healpix-weights.fits')
    # Use a low precision to speed up the calculation.  Use lower nside
    # to create a smaller file.
    w = pixweight(64, tiles=tx, precision=0.01, outfile=outfile)


def trim_inputs(indir, outdir):
    '''Don't copy any inputs'''
    pass


def trim_sky(indir, outdir):
    '''copy solarspec file as-is'''
    assert os.path.basename(indir) == 'sky'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    infile = os.path.join(indir, 'solarspec.txt')
    outfile = os.path.join(outdir, 'solarspec.txt')
    shutil.copy(infile, outfile)


def trim_specpsf(indir, outdir):
    '''trim specpsf files to be much smaller'''
    assert os.path.basename(indir) == 'specpsf'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    trim_quickpsf(indir, outdir, 'psf-quicksim.fits')
    trim_psf(indir, outdir, 'psf-b.fits')
    trim_psf(indir, outdir, 'psf-r.fits')
    trim_psf(indir, outdir, 'psf-z.fits')


def trim_spectra(indir, outdir):
    '''downsample spectra, and only a few of them'''
    assert os.path.basename(indir) == 'spectra'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for filename in ('spec-ABmag22.0.dat',
                     'spec-elg-o2flux-8e-17-average-line-ratios.dat',
                     'spec-lrg-z0.8-zmag20.38.dat',
                     'spec-qso-z1.5-rmag22.81.dat',
                     'spec-sky.dat',
                     'spec-sky-grey.dat',
                     'spec-sky-bright.dat',
                     'ZenithExtinction-KPNO.dat'):
        i = 0
        with open(os.path.join(indir, filename)) as infx:
            with open(os.path.join(outdir, filename), 'w') as outfx:
                for line in infx:
                    if line.startswith('#'):
                        outfx.write(line)
                    elif i % 20 == 0:
                        outfx.write(line)
                    i += 1


def trim_weather(indir, outdir):
    '''copy everything in weather/'''
    assert os.path.basename(indir) == 'weather'
    shutil.copytree(indir, outdir)


def trim_targets(indir, outdir):
    '''copy everything in targets/'''
    assert os.path.basename(indir) == 'targets'
    shutil.copytree(indir, outdir)


def trim_throughput(indir, outdir):
    '''downsample throughput files'''
    assert os.path.basename(indir) == 'throughput'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for targettype in ('elg', 'lrg', 'perfect', 'qso', 'sky', 'star'):
        filename = 'fiberloss-{}.dat'.format(targettype)
        shutil.copy(os.path.join(indir, filename),
                    os.path.join(outdir, filename))

    for filename in ('thru-b.fits', 'thru-r.fits', 'thru-z.fits'):
        with fits.open(os.path.join(indir, filename)) as fx:
            hdus = fits.HDUList()
            hdus.append(fx[0])
            hdus.append(fits.BinTableHDU(fx[1].data[::20], header=fx[1].header))
            hdus.append(fits.BinTableHDU(fx[2].data[::20], header=fx[2].header))
            hdus.writeto(os.path.join(outdir, filename))

    # galsim-fiber-acceptance.fits is about 230 KB, and it's a fairly
    # complex file, so just copy as-is.
    for filename in ('DESI-0347_blur.ecsv', 'DESI-0347_offset.ecsv',
                     'DESI-0347_static_offset_1.fits',
                     'galsim-fiber-acceptance.fits'):
        shutil.copy(os.path.join(indir, filename), os.path.join(outdir, filename))


#
# Triming PSF files
#
def rebin_image(image, n):
    """
    rebin 2D array pix into bins of size n x n

    New binsize must be evenly divisible into original pix image
    """
    assert image.shape[0] % n == 0
    assert image.shape[1] % n == 0

    s = image.shape[0]//n, n, image.shape[1]//n, n
    return image.reshape(s).sum(-1).sum(1)


def trim_psf(indir, outdir, filename):
    assert os.path.abspath(indir) != os.path.abspath(outdir)
    infile = os.path.join(indir, filename)
    outfile = os.path.join(outdir, filename)

    fx = fits.open(infile)
    hdus = fits.HDUList()

    #- HDU 0 XCOEFF - data unchanged but update keywords for less samples
    xcoeff = fx[0].data
    hdr = fx[0].header
    hdr['NWAVE'] = 3    #- down from 11
    hdr['CRPIX1'] = 23  #- 23=45//2+1, down from 113=225//2+1
    hdr['CRPIX1'] = 23
    hdr['CDELT1'] = 0.005   #- 5mm instead of 1mm
    hdr['CDELT2'] = 0.005   #- 5mm instead of 1mm
    hdr['PIXSIZE'] = 0.005   #- 5mm instead of 1mm

    hdus.append(fits.PrimaryHDU(xcoeff, header=hdr))
    hdus.append(fx['YCOEFF'])

    #- subsample spots
    inspots = fx['SPOTS'].data
    spots = np.zeros((3,3,45,45))
    spots[0,0] = rebin_image(inspots[0,0], 5)
    spots[1,0] = rebin_image(inspots[5,0], 5)
    spots[2,0] = rebin_image(inspots[10,0], 5)
    spots[0,1] = rebin_image(inspots[0,5], 5)
    spots[1,1] = rebin_image(inspots[5,5], 5)
    spots[2,1] = rebin_image(inspots[10,5], 5)
    spots[0,2] = rebin_image(inspots[0,10], 5)
    spots[1,2] = rebin_image(inspots[5,10], 5)
    spots[2,2] = rebin_image(inspots[10,10], 5)
    hdus.append(fits.ImageHDU(spots, header=fx['SPOTS'].header))

    #- subsample spots x,y locations
    dx = fx['SPOTX'].data
    hdus.append(fits.ImageHDU(dx[::5, ::5], header=fx['SPOTX'].header))
    dy = fx['SPOTY'].data
    hdus.append(fits.ImageHDU(dy[::5, ::5], header=fx['SPOTY'].header))

    #- Fiberpos unchanged
    hdus.append(fx['FIBERPOS'])

    #- Subsample SPOTPOS and SPOTWAVE
    d = fx['SPOTPOS'].data
    hdus.append(fits.ImageHDU(d[::5], header=fx['SPOTPOS'].header))
    d = fx['SPOTWAVE'].data
    hdus.append(fits.ImageHDU(d[::5], header=fx['SPOTWAVE'].header))

    hdus.writeto(outfile, overwrite=True)
    fx.close()


def trim_quickpsf(indir, outdir, filename):
    assert os.path.abspath(indir) != os.path.abspath(outdir)
    infile = os.path.join(indir, filename)
    outfile = os.path.join(outdir, filename)
    fx = fits.open(infile)
    hdus = fits.HDUList()
    hdus.append(fx[0])
    for i in [1,2,3]:
        d = fx[i].data
        hdus.append(fits.BinTableHDU(d[::10], header=fx[i].header))
    hdus.writeto(outfile, overwrite=True)
    fx.close()
