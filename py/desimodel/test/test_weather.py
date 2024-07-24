# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.weather.
"""
import unittest
import datetime
import os
import numpy as np
import astropy.table
from .. import weather as w


class TestWeather(unittest.TestCase):
    """Test desimodel.weather.
    """
    def setUp(self):
        self.gen = np.random.RandomState(seed=123)

    def test_whiten_transforms_from_cdf(self):
        """Test inputs to whiten_transforms_from_cdf().
        """
        x = np.ones((50,), dtype=np.float32)
        cdf = np.linspace(0, 1, 50)
        with self.assertRaises(ValueError) as e:
            foo = w.whiten_transforms_from_cdf(cdf[::-1], cdf)
        self.assertEqual(str(e.exception), 'Values of x must be non-decreasing.')
        with self.assertRaises(ValueError) as e:
            foo = w.whiten_transforms_from_cdf(x, cdf[::-1])
        self.assertEqual(str(e.exception), 'Values of cdf must be non-decreasing.')
        with self.assertRaises(ValueError) as e:
            foo = w.whiten_transforms_from_cdf(x, cdf[:49])
        self.assertEqual(str(e.exception), 'Input arrays must have same shape.')
        xx = np.ones((50,5), dtype=np.float32)
        with self.assertRaises(ValueError) as e:
            foo = w.whiten_transforms_from_cdf(xx, xx)
        self.assertEqual(str(e.exception), 'Input arrays must be 1D.')

    def test_whiten_transforms_inputs(self):
        """Test inputs to whiten_transforms().
        """
        x = np.linspace(0, 1, 50)
        with self.assertRaises(ValueError) as e:
            F, G = w.whiten_transforms(x, data_min=0.1)
        self.assertEqual(str(e.exception), 'data_min > min(data)')
        with self.assertRaises(ValueError) as e:
            F, G = w.whiten_transforms(x, data_max=0.9)
        self.assertEqual(str(e.exception), 'data_max < max(data)')
        F, G = w.whiten_transforms(x)

    def test_whiten_identity(self):
        """Check that whitener for samples from a Gaussian is linear.
        """
        mu, sigma = 1.2, 3.4
        x = mu + sigma * self.gen.normal(size=100000)
        self.assertTrue(np.allclose(mu, np.mean(x), atol=1e-2))
        self.assertTrue(np.allclose(sigma, np.std(x), atol=1e-2))
        F, G = w.whiten_transforms(x, data_min=-1e6, data_max=+1e6)
        y = F(x)
        self.assertTrue(np.allclose(0, np.mean(y), atol=1e-2))
        self.assertTrue(np.allclose(1, np.std(y), atol=1e-2))
        ypred = (x - mu) / sigma
        self.assertTrue(np.allclose(ypred, y, atol=0.2))

    def test_whiten_unwhiten_roundtrip(self):
        """Check that F(G(x)) = x for a uniform distribution.
        """
        lo, hi = -1.2, +3.4
        x = self.gen.uniform(size=10000)
        F, G = w.whiten_transforms(x, data_min=lo, data_max=hi)
        self.assertTrue(np.allclose(G(F(x)), x))

    def test_whiten_monotonic(self):
        """Check G(x) is monotonically increasing for a uniform distribution.
        """
        lo, hi = -1.2, +3.4
        x = self.gen.uniform(size=10000)
        F, G = w.whiten_transforms(x, data_min=lo, data_max=hi)
        xs = np.sort(x)
        ys = F(xs)
        self.assertTrue(np.all(np.diff(ys) >= 0))

    def test_seeing_pdf_norm(self):
        """Check that seeing PDF is normalized.
        """
        with self.assertRaises(ValueError) as e:
            fwhm, pdf = w.get_seeing_pdf(5.0)
        self.assertEqual(str(e.exception), 'Requested median is outside allowed range.')
        fwhm, pdf = w.get_seeing_pdf()
        norm = np.sum(pdf * np.gradient(fwhm))
        self.assertTrue(np.allclose(norm, 1))

    def test_sample_timeseries(self):
        """Test sampling white noise with uniform 1D PDF.
        """
        x = np.linspace(-1, 1, 500)
        pdf = np.ones_like(x)
        psd = lambda freq: np.ones_like(freq)
        n, nb = 1000000, 10
        gen = np.random.RandomState(1)
        xs = w.sample_timeseries(x, pdf, psd, n, gen=gen)
        bins, _ = np.histogram(xs, range=(-1,1), bins=nb)
        pred = n / float(nb)
        self.assertTrue(np.allclose(bins, pred, atol=5*np.sqrt(pred)))
        xx = np.ones((500,), dtype=pdf.dtype)
        with self.assertRaises(ValueError) as e:
            xs = w.sample_timeseries(xx, pdf, psd, n, gen=gen)
        self.assertEqual(str(e.exception), 'x_grid values are not increasing.')
        with self.assertRaises(ValueError) as e:
            xs = w.sample_timeseries(x[:400], pdf, psd, n, gen=gen)
        self.assertEqual(str(e.exception), 'x_grid and pdf_grid arrays have different shapes.')

    def test_same_seed(self):
        """Same seed should give same samplesself.
        """
        x_grid = np.linspace(-1, 1, 500)
        pdf_grid = np.ones_like(x_grid)
        psd = lambda freq: np.ones_like(freq)
        n_sample = 1000
        gen1 = np.random.RandomState(seed=123)
        x1 = w.sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen1)
        gen2 = np.random.RandomState(seed=123)
        x2 = w.sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen2)
        self.assertTrue(np.all(x1 == x2))

    def test_different_seed(self):
        """Different seeds should give different samples.
        """
        x_grid = np.linspace(-1, 1, 500)
        pdf_grid = np.ones_like(x_grid)
        psd = lambda freq: np.ones_like(freq)
        n_sample = 1000
        gen1 = np.random.RandomState(seed=1)
        x1 = w.sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen1)
        gen2 = np.random.RandomState(seed=2)
        x2 = w.sample_timeseries(x_grid, pdf_grid, psd, n_sample, gen=gen2)
        self.assertTrue(not np.any(x1 == x2))

    def test_seeing_median(self):
        """Check that seeing has expected median.
        """
        n = 100000
        for m in (0.9, 1.0, 1.1, 1.2):
            x = w.sample_seeing(n, median_seeing=m)
            self.assertTrue(np.fabs(np.median(x) - m) < 0.01)

    def test_transp_range(self):
        """Check that transparency has expected range.
        """
        n = 100000
        x = w.sample_transp(n)
        self.assertTrue(np.min(x) >= 0 and np.min(x) < 0.0001)
        self.assertTrue(np.max(x) <= 1 and np.max(x) > 0.9999)

    def test_dome_frac_means(self):
        """Check that mean weather fractions are reasonable.
        """
        first = datetime.date(2019, 12, 1)
        last = datetime.date(2024, 11, 30)
        n = (last - first).days
        for yr in range(2007, 2018):
            probs = w.dome_closed_fractions(first, last, replay='Y{}'.format(yr))
            self.assertTrue(len(probs) == n)
            mean = probs.mean()
            self.assertTrue((mean > 0.25) & (mean < 0.35))

    def test_dome_frac_values(self):
        """Check correct replay of two weather years.
        """
        DESIMODEL = os.getenv('DESIMODEL')
        path = os.path.join(DESIMODEL, 'data', 'weather', 'daily-2007-2017.csv')
        t = astropy.table.Table.read(path)
        probs = w.dome_closed_fractions(
            datetime.date(2021, 1, 1), datetime.date(2023, 1, 1),
            replay='Y2017,Y2007')
        self.assertTrue(np.all(probs[:365] == t['Y2017']))
        self.assertTrue(np.all(probs[365:] == t['Y2007']))
