# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.weather.
"""
from __future__ import print_function, division

import unittest
import numpy as np
from desimodel.weather import *


class TestWeather(unittest.TestCase):
    """Test desimodel.weather.
    """
    def test_whiten_identity(self):
        """Check that whitener for samples from a Gaussian is linear.
        """
        mu, sigma = 1.2, 3.4
        gen = np.random.RandomState(seed=123)
        x = mu + sigma * gen.normal(size=100000)
        assert np.allclose(mu, np.mean(x), atol=1e-2)
        assert np.allclose(sigma, np.std(x), atol=1e-2)
        F, G = whiten_transforms(x, data_min=-1e6, data_max=+1e6)
        y = F(x)
        assert np.allclose(0, np.mean(y), atol=1e-2)
        assert np.allclose(1, np.std(y), atol=1e-2)
        ypred = (x - mu) / sigma
        assert np.allclose(ypred, y, atol=0.2)

    def test_whiten_unwhiten_roundtrip(self):
        """Check that F(G(x)) = x for a uniform distribution.
        """
        gen = np.random.RandomState(seed=123)
        lo, hi = -1.2, +3.4
        x = gen.uniform(size=10000)
        F, G = whiten_transforms(x, data_min=lo, data_max=hi)
        assert np.allclose(G(F(x)), x)

    def test_whiten_monotonic(self):
        """Check G(x) is monotonically increasing for a uniform distribution.
        """
        gen = np.random.RandomState(seed=123)
        lo, hi = -1.2, +3.4
        x = gen.uniform(size=10000)
        F, G = whiten_transforms(x, data_min=lo, data_max=hi)
        xs = np.sort(x)
        ys = F(xs)
        assert np.all(np.diff(ys) >= 0)

def test_suite():
    """Allows testing of only this module with the command::
        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
