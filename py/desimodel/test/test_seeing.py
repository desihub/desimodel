# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.seeing.
"""
from __future__ import print_function, division

import unittest
import numpy as np
from ..seeing import relative_probability, sample


class TestSeeing(unittest.TestCase):
    """Test desimodel.seeing.
    """
    def test_pdf_norm(self):
        """Test that PDF is normalized and has expected mean.
        """
        fwhm = np.linspace(0., 10., 200)
        pdf = relative_probability(fwhm)
        self.assertAlmostEqual(pdf.sum(), 1.)
        # Calculate mean seeing.
        mean = (fwhm * pdf).sum() / pdf.sum()
        self.assertAlmostEqual(mean, 1.726, places=3)


    def test_samples(self):
        """Test that samples have expected PDF.
        """
        # Histogram 1M generated samples.
        samples = sample(1000000, seed=123)
        bin_edges = np.linspace(0., 10., 51)
        bin_prob, _ = np.histogram(samples, bin_edges, density=True)
        bin_prob *= bin_edges[1]
        # Calculate the expected number of samples in each bin.
        fwhm = np.linspace(0., 10., 10 * len(bin_prob) + 1)
        fwhm_midpt = 0.5 * (fwhm[:-1] + fwhm[1:])
        pdf = relative_probability(fwhm_midpt)
        expected = pdf.reshape(len(bin_prob), -1).sum(axis=1)
        # Check for expected bin counts.
        self.assertTrue(np.allclose(bin_prob, expected, rtol=1e-2, atol=1e-2))


def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
