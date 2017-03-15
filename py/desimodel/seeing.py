# See LICENSE.rst for BSD 3-clause license info
# -*- coding: utf-8 -*-
"""
================
desimodel.seeing
================

Model of the expected DESI zenith atmospheric seeing at 6355A.

The model is based on seeing data collected with the MOSAIC camera on the
Mayall telescope at KPNO, as described and analyzed in this reference:

  A. Dey and F. Valdes, "The Delivered Image Quality with the MOSAIC Cameras
  at the Kitt Peak 4 m Mayall and Cerro Tololo 4 m Blanco Telescopes",
  2014 PASP, 126:296–311, `doi:10.1086/675808
  <https://doi.org/10.1086/675808>`__.

The seeing values modeled here should normally be scaled to the observing
wavelength and airmass.  For DESI simulations, these scalings are applied
in the `specsim package <http://specsim.readthedocs.io/>`__.
"""
from __future__ import print_function, division

import numpy as np
import scipy.special


"""Seeing distribution parameters.

Parameters for equation 1 of Dey & Valdes 2014 that describe the nominal
DESI delivered image quality distribution at 6355A.  These values correspond
to the distribution of all image qualities observed with the "r SDSS k1018"
filter.  Note that the normalization of this parameterization refers to the
number of samples used in the Dey & Valdes 2014 analysis.
"""
parameters = {
    'A': 3013631.5, 'C': 0.52700073, 'B': 0.4630833,
    'E': 30.287271, 'D': 14.632419}


def relative_probability(fwhm_arcsec, blur_arcsec=0.219):
    """Tabulate relative probabilities of FWHM zenith seeing at 6355A.

    Values are derived from equation 1 of Dey & Valdes 2014 with the following
    additional steps: subtract instrumental RMS blur in quadrature
    to convert from delivered image quality to seeing and normalize so that the
    sum of the returned array values is one.  Note that the tabulated values
    are not true probabilities, which would require a numerical
    integration up to infinity, but are close when the input array is equally
    spaced and covers most of the non-zero range ~ [0, 10] arcsec.

    Parameters
    ----------
    fwhm_arcsec : array
        Array of FWHM values in arcseconds where the probability density
        should be tabulated.
    blur_arcsec : float
        RMS instrumental blur to subtract in quadrature.  Set this to zero
        to tabulate the PDF of delivered image quality instead of seeing.

    Returns
    -------
    array
        Array of relative probabilities calculated at each input value and
        normalized to one.
    """
    # Convert from seeing to delivered image quality.
    const = 2 * np.sqrt(2 * np.log(2))
    fwhm_arcsec = const * np.sqrt(
        (fwhm_arcsec / const) ** 2 + blur_arcsec ** 2)
    eqn1 = (parameters['A'] *
            (1 + ((fwhm_arcsec - parameters['B']) /
                  parameters['C']) **2 ) ** (-parameters['D']) *
                  fwhm_arcsec ** parameters['E'])
    return eqn1 / eqn1.sum()


def sample(n_sample, dt_sec=300., max_seeing=5., seed=None,
           psd_tau1=0.025, psd_tau2=5.):
    """Generate random samples of FWHM zenith seeing at 6355A.

    Samples are generated as a time series on a uniform grid of observation
    times, with a 1D probability distribution given by
    :func:`relative_probability` and an autocorrelation power spectrum
    specified by two time constants.

    The default power spectrum parameters are matched to a two-year series of
    DIMM measurements at Cerro Pachon.  The power is parameterized by two
    time constants that yield a spectrum which is initially flat, then
    transitions to 1/f scaling at the smaller time constant and 1/f**2
    at the larger time constant.

    The power spectrum describes "whitened" seeing fluctuations that are then
    transformed to obtain the desired 1D distribution of seeing values. As a
    result, the PSD time constants can be specified independently of the
    seeing distribution.

    The algorithm is fast enough to use even when the time structure is not
    needed.  The time required to simulate the full 5-year survey at 5-minute
    resolution (~175K samples) is < 1 second.  The memory requirements
    and timing are both approximately linear in n_sample.

    Parameters
    ----------
    n_sample : int
        Number of samples to generate. Must be at least 2.
    dt_sec : float
        Elapsed time in seconds between generated samples.
    max_seeing : float
        Maximum FWHM zenith seeing at 6355A to generate. This value
        effectively truncates the 1D distribution that will be sampled.
    seed : int or None
        Random number seed to use for reproducible samples. Use the default
        value of None to generate different samples with each call.
    psd_tau1 : float
        Time constant associated with the autocorrelation power spectral
        density, in units of days.
    psd_tau2 : float
        Time constant associated with the autocorrelation power spectral
        density, in units of days.

    Returns
    -------
    array:
        Array of generated seeing samples on a uniform time grid.  All values
        will be in the range [0, max_seeing] and represent FWHM zenith
        seeing at 6355A.
    """
    if n_sample < 2:
        raise ValueError('n_sample must be at least 2.')

    # Build a grid covering the range of allowed FWHM seeing values.
    seeing_grid = np.linspace(0., max_seeing, 1000)

    # Build a table of cummulative probabilities.
    seeing_cdf = np.cumsum(relative_probability(seeing_grid))

    # Build a linear grid of frequencies present in the Fourier transform
    # of the requested time series.  Frequency units are 1/day.
    dt_day = dt_sec / (24. * 3600.)
    df_day = 1. / (n_sample * dt_day)
    f_grid = np.arange(1 + (n_sample // 2)) * df_day

    # Tabulate the power spectral density at each frequency.  The PSD
    # describes seeing fluctuations that have been "whitened", i.e., mapped
    # via a non-linear monotonic transform to have unit Gaussian probability
    # density.
    omega = 2 * np.pi * f_grid
    psd = 1. / (1. + omega * psd_tau1) / (1. + omega * psd_tau2)
    # Force the mean to zero.
    psd[0] = 0.
    # Force the variance to one.
    psd[1:] /= psd[1:].sum() * df_day ** 2

    # Generate random whitened samples using the specified seed.
    gen = np.random.RandomState(seed)
    n_psd = len(psd)
    x_fft = np.ones(n_psd, dtype=complex)
    x_fft[1:-1].real = gen.normal(size=n_psd - 2)
    x_fft[1:-1].imag = gen.normal(size=n_psd - 2)
    x_fft *= np.sqrt(psd) / (2 * dt_day)
    x_fft[0] *= np.sqrt(2)
    x = np.fft.irfft(x_fft, n_sample)

    # Un-whiten the samples to recover the desired 1D PDF.
    x_cdf = 0.5 * (1 + scipy.special.erf(x / np.sqrt(2)))
    seeing = np.interp(x_cdf, seeing_cdf, seeing_grid)

    return seeing
