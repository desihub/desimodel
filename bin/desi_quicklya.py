#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute the Lyman-alpha S/N for the mean QSO at different g-mag, z,
as a function of wavelength.
It is the python version of the IDL code desi_quicklya.pro by D. Schlegel (LBL),
and it is heavily based on quicksim.py by D. Kirkby (UC Irvine)

Use desi_quicklya.py --help for instructions on running this program.

The DESIMODEL environment variable should be set to the root of your
desimodel package and PYTHONPATH should include $DESIMODEL/py/, e.g.

export DESIMODEL=`pwd`
export PYTHONPATH=$DESIMODEL/py/:$PYTHONPATH

Created 06-May-2015 by Andreu Font-Ribera (afont@lbl.gov)
Modified Jan-2018 by J. Guy
"""

import time
import argparse
import os
import os.path
import numpy as np
import astropy.units as units

#import desimodel.simulate as sim

from desiutil.log import get_logger
from desispec.io.filters import load_filter
import desisim.simexp
import desisim.obs
import desisim.io
import desisim.util
import desitarget
import desispec.io
import desimodel.io    

def sim_spectra(wave, flux, program, obsconditions=None,
                sourcetype=None, expid=0, seed=0, specsim_config_file="desi"):
    
    """
    Simulate spectra from an input set of wavelength and flux and writes a FITS file in the Spectra format that can
    be used as input to the redshift fitter.

    Args:
        wave : 1D np.array of wavelength in Angstrom (in vacuum) in observer frame (i.e. redshifted)
        flux : 1D or 2D np.array. 1D array must have same size as wave, 2D array must have shape[1]=wave.size
               flux has to be in units of 10^-17 ergs/s/cm2/A
        program : dark, lrg, qso, gray, grey, elg, bright, mws, bgs
            ignored if obsconditions is not None
    
    Optional:
        obsconditions : dictionnary of observation conditions with SEEING EXPTIME AIRMASS MOONFRAC MOONALT MOONSEP
        sourcetype : list of string, allowed values are (sky,elg,lrg,qso,bgs,star), type of sources, used for fiber aperture loss , default is star
        expid : this expid number will be saved in the Spectra fibermap
        seed : random seed
        skyerr : fractional sky subtraction error
    """
    
    log = get_logger()
    
    if len(flux.shape)==1 :
        flux=flux.reshape((1,flux.size))
    nspec=flux.shape[0]
    
    log.info("Starting simulation of {} spectra".format(nspec))
    
    if sourcetype is None :        
        sourcetype = np.array(["star" for i in range(nspec)])
    log.debug("sourcetype = {}".format(sourcetype))
    
    tileid  = 0
    telera  = 0
    teledec = 0    
    dateobs = time.gmtime()
    night   = desisim.obs.get_night(utc=dateobs)
    program = program.lower()
        
       
    frame_fibermap = desispec.io.fibermap.empty_fibermap(nspec)    
    frame_fibermap.meta["FLAVOR"]="custom"
    frame_fibermap.meta["NIGHT"]=night
    frame_fibermap.meta["EXPID"]=expid
    
    # add DESI_TARGET 
    tm = desitarget.desi_mask    
    frame_fibermap['DESI_TARGET'][sourcetype=="star"]=tm.STD_FSTAR
    frame_fibermap['DESI_TARGET'][sourcetype=="lrg"]=tm.LRG
    frame_fibermap['DESI_TARGET'][sourcetype=="elg"]=tm.ELG
    frame_fibermap['DESI_TARGET'][sourcetype=="qso"]=tm.QSO
    frame_fibermap['DESI_TARGET'][sourcetype=="sky"]=tm.SKY
    frame_fibermap['DESI_TARGET'][sourcetype=="bgs"]=tm.BGS_ANY
    
    # add dummy TARGETID
    frame_fibermap['TARGETID']=np.arange(nspec).astype(int)
         
    # spectra fibermap has two extra fields : night and expid
    # This would be cleaner if desispec would provide the spectra equivalent
    # of desispec.io.empty_fibermap()
    spectra_fibermap = desispec.io.empty_fibermap(nspec)
    spectra_fibermap = desispec.io.util.add_columns(spectra_fibermap,
                       ['NIGHT', 'EXPID', 'TILEID'],
                       [np.int32(night), np.int32(expid), np.int32(tileid)],
                       )

    for s in range(nspec):
        for tp in frame_fibermap.dtype.fields:
            spectra_fibermap[s][tp] = frame_fibermap[s][tp]
    
    if obsconditions is None:
        if program in ['dark', 'lrg', 'qso']:
            obsconditions = desisim.simexp.reference_conditions['DARK']
        elif program in ['elg', 'gray', 'grey']:
            obsconditions = desisim.simexp.reference_conditions['GRAY']
        elif program in ['mws', 'bgs', 'bright']:
            obsconditions = desisim.simexp.reference_conditions['BRIGHT']
        else:
            raise ValueError('unknown program {}'.format(program))
    elif isinstance(obsconditions, str):
        try:
            obsconditions = desisim.simexp.reference_conditions[obsconditions.upper()]
        except KeyError:
            raise ValueError('obsconditions {} not in {}'.format(
                obsconditions.upper(),
                list(desisim.simexp.reference_conditions.keys())))
    try:
        params = desimodel.io.load_desiparams()
        wavemin = params['ccd']['b']['wavemin']
        wavemax = params['ccd']['z']['wavemax']
    except KeyError:
        wavemin = desimodel.io.load_throughput('b').wavemin
        wavemax = desimodel.io.load_throughput('z').wavemax

    if wave[0] > wavemin:
        log.warning('Minimum input wavelength {}>{}; padding with zeros'.format(
                wave[0], wavemin))
        dwave = wave[1] - wave[0]
        npad = int((wave[0] - wavemin)/dwave + 1)
        wavepad = np.arange(npad) * dwave
        wavepad += wave[0] - dwave - wavepad[-1]
        fluxpad = np.zeros((flux.shape[0], len(wavepad)), dtype=flux.dtype)
        wave = np.concatenate([wavepad, wave])
        flux = np.hstack([fluxpad, flux])
        assert flux.shape[1] == len(wave)
        assert np.allclose(dwave, np.diff(wave))
        assert wave[0] <= wavemin

    if wave[-1] < wavemax:
        log.warning('Maximum input wavelength {}<{}; padding with zeros'.format(
                wave[-1], wavemax))
        dwave = wave[-1] - wave[-2]
        npad = int( (wavemax - wave[-1])/dwave + 1 )
        wavepad = wave[-1] + dwave + np.arange(npad)*dwave
        fluxpad = np.zeros((flux.shape[0], len(wavepad)), dtype=flux.dtype)
        wave = np.concatenate([wave, wavepad])
        flux = np.hstack([flux, fluxpad])
        assert flux.shape[1] == len(wave)
        assert np.allclose(dwave, np.diff(wave))
        assert wavemax <= wave[-1]

    ii = (wavemin <= wave) & (wave <= wavemax)

    flux_unit = 1e-17 * units.erg / (units.Angstrom * units.s * units.cm ** 2 )
    
    wave = wave[ii]*units.Angstrom
    flux = flux[:,ii]*flux_unit
    
    nspec = flux.shape[0]
    
    sim = desisim.simexp.simulate_spectra(wave, flux, fibermap=frame_fibermap,
                                          obsconditions=obsconditions, seed=seed,specsim_config_file = specsim_config_file)
    
    
    
    
    
    
    # full wave array
    wmin=1e12
    wmax=0.
    dwave=0.
    for table in sim.camera_output :
        twave = table['wavelength'].astype(float)
        wmin = min(wmin,np.min(twave))
        wmax = max(wmax,np.max(twave))
        if dwave==0 : 
            dwave=twave[1]-twave[0]
        else :
            assert(np.abs(dwave-(twave[1]-twave[0]))<0.0001)
    
    wave=np.linspace(wmin,wmax,int((wmax-wmin)/dwave)+1)
    log.debug("wmin wmax dwave wave= {} {} {} {}".format(wmin,wmax,dwave,wave))
    
    sivarflux = np.zeros((nspec,wave.size))
    sivar     = np.zeros((nspec,wave.size))
        
    # total signal on all cameras
    for table in sim.camera_output :
        twave = table['wavelength'].astype(float)
        tivar = table['flux_inverse_variance'].T.astype(float)
        tivarflux = table['flux_inverse_variance'].T.astype(float)*table['observed_flux'].T.astype(float)
        for s in range(nspec) :
            sivar[s]     += np.interp(wave,twave,tivar[s],left=0,right=0)
            sivarflux[s] += np.interp(wave,twave,tivarflux[s],left=0,right=0)
    
            
    scale=1e17
    flux = np.zeros(sivar.shape)
    for s in range(nspec) :
        ii=(sivar[s]>0)
        flux[s,ii] = sivarflux[s,ii]/sivar[s,ii] * scale
    ivar = sivar / scale**2
    
    return wave,flux,ivar
    


def main():



    # parse command-line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v','--verbose', action = 'store_true',
        help = 'provide verbose output on progress')
    parser.add_argument('--ab-magnitude', type = str, default = None,
                        help = 'max magnitude to compute, e.g. g=22.0 or r=21.5')
    parser.add_argument('--exptime', type = float, default = 4000,
                        help = 'overrides exposure time specified in the parameter file (secs)')
    parser.add_argument('--config', type = str, default = "desi",
                        help = 'path to specsim configuration file')
    parser.add_argument('--prefix', type = str, default = "sn-spec-lya",
                        help = 'prefix for output S/N files')
    args = parser.parse_args()
    
    log = get_logger()

    obsconditions = desisim.simexp.reference_conditions['DARK']
    obsconditions["EXPTIME"]=args.exptime
    

    # We require that the DESIMODEL environment variable is set.
    if 'DESIMODEL' not in os.environ:
        raise RuntimeError('The environment variable DESIMODEL must be set.')

    # Load the source spectrum to use.
    infile=os.environ['DESIMODEL']+'/data/spectra/spec-lya.dat'
    if not os.path.isfile(infile):
        print('Could not find file: %s' % infile)
        return -1

    x=np.loadtxt(infile).T
    wave = x[0]
    flux = x[1:] # these are mean QSO spectra at different redshifts
    # this is the same grid in the infile file, can't be changed here!
    zqs = np.linspace(2.0, 4.75, 12)
    
    assert(flux.shape[0] == zqs.size)
    
    if args.ab_magnitude is not None :
        # figure out magnitude
        try:
            band = args.ab_magnitude[0]
            abmag = float(args.ab_magnitude[2:])
            assert band in 'ugriz' and args.ab_magnitude[1] == '='
        except(AssertionError,ValueError):
            print('Invalid ab-magnitude parameter. '
                  +'Valid syntax is, e.g. g=22.0 or r=21.5.')
            return -1
        mags=[abmag]
    else :
        # will generate a file for each g magnitude
        band="g"
        min_m=19.25
        dm=0.5
        Nm=np.ceil((abmag-min_m)/dm)
        mags = np.linspace(min_m, min_m+Nm*dm, Nm+1)
        if args.verbose: print('mags', mags)
        
    # compute magnitudes of QSO spectra in input file
    # assuming flux is prop to ergs/s/cm2/A 
    # (norme does not matter because the spectra will be rescaled)
    
    filter_response = load_filter("SDSS_"+band.upper())
    fluxunits       = 1e-17 * units.erg / units.s / units.cm**2 / units.Angstrom
    infile_mags     = np.zeros(flux.shape[0])
    for s in range(flux.shape[0]) :
        infile_mags[s] = filter_response.get_ab_magnitude(flux[s]*fluxunits,wave)
        print(s,infile_mags[s])
    
    for mag in mags:
        if args.verbose: print('generate file for %s = %f' % (band,mag) )

        # scaling fluxes
        scaled_flux = np.zeros(flux.shape)
        for s in range(flux.shape[0]) :
            scaled_flux[s] = 10**(-0.4*(mag-infile_mags[s])) * flux[s]
        
        # simulating
        sim_wave,sim_flux,sim_ivar = sim_spectra(wave, scaled_flux, program="dark", obsconditions=obsconditions, sourcetype="qso", specsim_config_file=args.config)
        sim_snr = np.sqrt(sim_ivar)*sim_flux/np.sqrt(np.gradient(sim_wave)) # S/N per sqrt(A)
        
        # for each wavelength, we will store a SN per redshift
        collect_results = []
        
        #  min_wavelength,max_wavelength,wavelength_step
        for i,zq in enumerate(zqs):
            # store values for each wavelength
            for j in range(sim_wave.size):
                if i==0:
                    new_row = {'wave': sim_wave[j], 'snr_'+str(zq): sim_snr[0,j]}
                    collect_results.append(new_row)
                else:
                    collect_results[j]['snr_'+str(zq)] = sim_snr[i,j]
        
        # Save the results to file
        fname = args.prefix+'-'+band+str(mag)+'-t'+str(int(args.exptime))+'.dat'
        
        log.info('Saving results to %s' % fname)
        # Try opening the requested output file.
        with open(fname,'w') as out:
            print('# Lyman-alpha forest S/N per Ang for mean quasar with '
                                                      +'mean forest',file=out)
            print('# INFILE=',infile,file=out)
            print('# BAND=',band,file=out)
            print('# MAG=',mag,file=out)
            print('# EXPTIME=',obsconditions["EXPTIME"],file=out)
            print('#',file=out)
            z_header = '# Wave'
            for zq in zqs: z_header += ' SN(z='+str(zq)+')'
            print(z_header,file=out)

            for row in collect_results:
                toprint = str('%9.2f' % row['wave'])
                for zq in zqs:
                    toprint += ' %9.2f' % row['snr_'+str(zq)]
                print(toprint, file=out)

if __name__ == '__main__':
    main()
