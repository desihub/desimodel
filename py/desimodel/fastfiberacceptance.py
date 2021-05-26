import os
import astropy.io.fits as pyfits
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d


def gaussian_fwhm(sigma):
    return 2. * np.sqrt(2. * np.log(2.)) * sigma

class FastFiberAcceptance(object):
    """
    This class reads an input fits file generated with specsim.fitgalsim
    ($DESIMODEL/data/throughput/galsim-fiber-acceptance.fits)
    and instanciates RegularGridInterpolator objects for 2D and 3D interpolation
    of the pre-computed galsim fiber acceptance as a function of sigma (atmosphere+telescope blur, in um on focal surface), fiber offset from source (in um on focal surface), and half light radius (in arcsec) from extended source.
    The average and rms interpolation function for POINT,DISK and BULGE profiles
    are loaded.
    """
    def __init__(self,filename=None):
        if filename is None :
            if not "DESIMODEL" in os.environ :
                print("need environment variable DESIMODEL or specify filename in constructor")
                raise RuntimeError("need environment variable DESIMODEL or specify filename in constructor")
            filename=os.path.join(os.environ["DESIMODEL"],"data/throughput/galsim-fiber-acceptance.fits")
        hdulist=pyfits.open(filename)

        sigma=hdulist["SIGMA"].data
        offset=hdulist["OFFSET"].data
        hlradius=hdulist["HLRAD"].data

        self._sigma=sigma
        self._offset=offset
        self._hlradius=hlradius

        self._data = {}

        self.fiber_acceptance_func = {}
        self.fiber_acceptance_rms_func = {}

        self.psf_seeing_func = {}

        for source in ["POINT","DISK","BULGE"] :
            data=hdulist[source].data
            rms=hdulist[source[0]+"RMS"].data
            dim=len(data.shape)

            self._data[source] = data

            if dim == 2 :
                assert source == 'POINT'

                # POINT: zero offset.
                self.psf_seeing_func[source] = interp1d(data[::-1,0], sigma[::-1], kind='linear', copy=True, bounds_error=False, assume_sorted=False, fill_value=(sigma[-1],sigma[0]))

                self.fiber_acceptance_func[source] = RegularGridInterpolator(points=(sigma,offset),values=data,method="linear",bounds_error=False,fill_value=None)
                self.fiber_acceptance_rms_func[source] = RegularGridInterpolator(points=(sigma,offset),values=rms,method="linear",bounds_error=False,fill_value=None)
            elif dim == 3 :
                # Not POINT
                self.fiber_acceptance_func[source] = RegularGridInterpolator(points=(hlradius,sigma,offset),values=data,method="linear",bounds_error=False,fill_value=None)
                self.fiber_acceptance_rms_func[source] = RegularGridInterpolator(points=(hlradius,sigma,offset),values=rms,method="linear",bounds_error=False,fill_value=None)

        hdulist.close()

    def psf_seeing_sigma(self, psf_fiberfrac):
        return  self.psf_seeing_func["POINT"](psf_fiberfrac)

    def psf_seeing_fwhm(self, psf_fiberfrac):
        sigma = self.psf_seeing_func["POINT"](psf_fiberfrac)

        return gaussian_fwhm(sigma)

    def rms(self,source,sigmas,offsets=None,hlradii=None) :
        """
        returns fiber acceptance fraction rms for the given source,sigmas,offsets

        Args:
            source (string) : POINT, DISK or BULGE for point source, exponential profile or De Vaucouleurs profile
            sigmas (np.array) : arbitrary shape, values of sigmas in um for the PSF due to atmosphere and telescope blur

        Optional:
            hlradii (np.array) : same shape as sigmas, half light radius in arcsec for source
            offsets (np.array) : same shape as sigmas, values of offsets on focal surface between fiber and source, in um

        Returns np.array with same shape as input
        """
        was_scalar = np.isscalar(sigmas)
        sigmas = np.atleast_1d(sigmas)
        original_shape = sigmas.shape

        if offsets is None :
            offsets=np.zeros(sigmas.shape)
        else :
            offsets=np.atleast_1d(offsets)
            assert(sigmas.shape==offsets.shape)

        if hlradii is not None :
            hlradii=np.atleast_1d(hlradii)
            assert(hlradii.shape==sigmas.shape)



        res = None
        if source == "POINT" :

            res = self.fiber_acceptance_rms_func[source](np.array([sigmas.ravel(),offsets.ravel()]).T)

        else :

            if hlradii is None :
                if source == "DISK" :
                    hlradii = 0.45 * np.ones(sigmas.shape)
                elif source == "BULGE" :
                    hlradii = 1. * np.ones(sigmas.shape)
            res = self.fiber_acceptance_rms_func[source](np.array([hlradii.ravel(),sigmas.ravel(),offsets.ravel()]).T)

        res[res<0] = 0.
        res[res>1] = 1.

        if was_scalar :
            return float(res[0])
        return res.reshape(original_shape)

    def value(self,source,sigmas,offsets=None,hlradii=None) :
        """
        returns the fiber acceptance for the given source,sigmas,offsets

        Args:
            source (string) : POINT, DISK or BULGE for point source, exponential profile or De Vaucouleurs profile
            sigmas (np.array) : arbitrary shape, values of sigmas in um for the PSF due to atmosphere and telescope blur
            offsets (np.array) : same shape as sigmas, values of offsets on focal surface between fiber and source, in um

        Optional:
            hlradii (np.array) : same shape as sigmas, half light radius in arcsec for source

        Returns np.array with same shape as input
        """

        was_scalar = np.isscalar(sigmas)
        sigmas = np.atleast_1d(sigmas)
        original_shape = sigmas.shape

        if offsets is None :
            offsets=np.zeros(sigmas.shape)
        else :
            offsets=np.atleast_1d(offsets)
            assert(sigmas.shape==offsets.shape)

        if hlradii is not None :
            hlradii=np.atleast_1d(hlradii)
            assert(hlradii.shape==sigmas.shape)

        res = None
        if source == "POINT" :

            res = self.fiber_acceptance_func[source](np.array([sigmas.ravel(),offsets.ravel()]).T)

        else :

            if hlradii is None :
                if source == "DISK" :
                    hlradii = 0.45 * np.ones(sigmas.shape)
                elif source == "BULGE" :
                    hlradii = 1. * np.ones(sigmas.shape)

            res = self.fiber_acceptance_func[source](np.array([hlradii.ravel(),sigmas.ravel(),offsets.ravel()]).T)

        res[res<0] = 0.
        res[res>1] = 1.

        if was_scalar :
            return float(res[0])
        return res.reshape(original_shape)


if __name__ == '__main__':
    import numpy as np
    import pylab as pl

    from fastfiberacceptance import FastFiberAcceptance


    x = FastFiberAcceptance()

    fiberfracs= np.arange(0.0,1.0, 0.01)
    seeings= x.psf_seeing_sigma(fiberfracs)

    avg_platescale = 1.52 / 107. # [''/microns].

    seeings *= avg_platescale

    print(x._sigma[::-1])
    print(x._sigma[::-1] * avg_platescale)
    print(x._data['POINT'][::-1,0])

    pl.figure()

    pl.subplot(121)

    pl.plot(fiberfracs, seeings)
    pl.plot(x._data['POINT'][::-1,0], x._sigma[::-1] * avg_platescale, marker='^', alpha=0.5)
    pl.xlabel('PSF FIBERFRAC')
    pl.ylabel('SEEING SIGMA [ARCSECONDS]')

    pl.subplot(122)

    fwhms= x.psf_seeing_fwhm(fiberfracs)
    fwhms *= avg_platescale

    pl.plot(fiberfracs, fwhms)

    pl.axhline(1.1, c='k', lw=0.5)
    pl.axvline(0.6, c='k', lw=0.5)

    pl.xlabel('PSF FIBERFRAC')
    pl.ylabel('SEEING FWHM [ARCSECONDS]')
    pl.show()
