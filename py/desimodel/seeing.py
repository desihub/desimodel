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
