# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Test desimodel.io.
"""
import os
import sys
import tempfile
import datetime
import numpy as np
from astropy.table import Table
import unittest
from .. import io
#
# Try to import specter.
#
# specter_available = True
# specter_message = "The specter package was not detected."
# try:
#     import specter
# except ImportError:
#     specter_available = False
#
# Check if desimodel data is available.
#
desimodel_available = os.path.isdir(io.datadir())
desimodel_message = "The desimodel data set was not detected."

#
# Try to load the DESI_SURVEYOPS environment variable.
#
surveyops_available = True
surveyops_message = "The DESI_SURVEYOPS directory was not detected."
try:
    _ = os.environ['DESI_SURVEYOPS']
except KeyError:
    surveyops_available = False


class TestIO(unittest.TestCase):
    """Test desimodel.io.
    """
    @classmethod
    def setUpClass(cls):
        # cls.specter_available = specter_available
        cls.desimodel_available = desimodel_available
        cls.surveyops_available = surveyops_available
        cls.tempdir = tempfile.mkdtemp(prefix='testio-')
        cls.trimdir = os.path.join(cls.tempdir, 'trim')
        cls.testfile = os.path.join(cls.tempdir, 'test-abc123.fits')
        cls.testecsvfile = os.path.join(cls.tempdir, 'test-abc123.ecsv')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.tempdir):
            _ = os.system(f"chmod -R +w {cls.tempdir}")
            import shutil
            shutil.rmtree(cls.tempdir)

    def setUp(self):
        """Ensure that any desimodel.io caches are clear before running
        any test.
        """
        io.reset_cache()
        if os.path.exists(self.testfile):
            os.remove(self.testfile)

        self.orig_desimodel = os.getenv('DESIMODEL') # None if not set

    def tearDown(self):
        if self.orig_desimodel is not None:
            os.environ['DESIMODEL'] = self.orig_desimodel
        elif 'DESIMODEL' in os.environ:
            del os.environ['DESIMODEL']

    def test_findfile(self):
        datadir = io.datadir()
        self.assertEqual(io.findfile('blat/foo'), f'{datadir}/blat/foo')

        os.environ['DESIMODEL'] = '/path/to/somewhere'
        datadir = io.datadir()
        self.assertEqual(datadir, '/path/to/somewhere/data')
        self.assertEqual(io.findfile('blat/foo'), f'{datadir}/blat/foo')

        del os.environ['DESIMODEL']
        datadir = io.datadir()
        self.assertEqual(io.findfile('blat/foo'), f'{datadir}/blat/foo')

    def test_reset_cache(self):
        """Test cache reset (two examples at least)
        """
        self.assertTrue(io._fiberpos is None)
        self.assertTrue(isinstance(io._tiles, dict))
        self.assertEqual(len(io._tiles), 0)
        fiberpos = io.load_fiberpos()
        tiles = io.load_tiles()
        self.assertTrue(io._fiberpos is not None)
        self.assertTrue(isinstance(io._tiles, dict))
        self.assertEqual(len(io._tiles), 1)
        io.reset_cache()
        self.assertTrue(io._fiberpos is None)
        self.assertTrue(isinstance(io._tiles, dict))
        self.assertEqual(len(io._tiles), 0)

    #- test requires both specter and desimodel data
    # @unittest.skipUnless(specter_available, specter_message)
    # @unittest.skipUnless(desimodel_available, desimodel_message)
    # def test_load_throughput(self):
    #     """Test loading of throughput files.
    #     """
    #     for channel in ('b', 'r', 'z'):
    #         t = io.load_throughput(channel)
    #         ww = np.arange(t.wavemin, t.wavemax, 0.1)
    #         msg = f'{channel} has negative throughput'
    #         self.assertTrue(np.all(t.fiberinput_throughput(ww)>=0.0), msg)
    #         self.assertTrue(np.all(t.hardware_throughput(ww)>=0.0), msg)
    #         self.assertTrue(np.all(t.thru(ww)>=0.0), msg)
    #         t1 = t.atmospheric_throughput(ww, airmass=1.0)
    #         t2 = t.atmospheric_throughput(ww, airmass=2.0)
    #         self.assertTrue(np.all(t1>=0.0))
    #         self.assertTrue(np.all(t2>=0.0))
    #         self.assertTrue(np.all(t1>=t2))

    #- test requires both specter and desimodel data
    # @unittest.skipUnless(specter_available, specter_message)
    # @unittest.skipUnless(desimodel_available, desimodel_message)
    # def test_load_psf(self):
    #     """Test loading of PSF files.
    #     """
    #     for channel in ('b', 'r', 'z'):
    #         t = io.load_psf(channel)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_desiparams(self):
        """Test loading of basic DESI parameters.
        """
        p = io.load_desiparams()
        self.assertTrue(isinstance(p, dict))

        #- Samity checks on augmented parameters
        for channel in ['b', 'r', 'z']:
            x = p['ccd'][channel]
            self.assertLess(x['wavemin'], x['wavemax'])
            self.assertLess(x['wavemin'], x['wavemin_all'])
            self.assertGreater(x['wavemax'], x['wavemax_all'])

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_gfa(self):
        """Test loading of gfa location data.
            """
        gfa = io.load_gfa()
        # length of the GFA table should be 40 since there are each corner of the 10 GFAs are included
        self.assertEqual(len(gfa), 40)
        for key in ('PETAL', 'CORNER', 'X', 'Y', 'Z', 'Q', 'S', 'RADIUS_DEG'):
            self.assertIn(key, gfa.dtype.names)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_fiberpos(self):
        """Test loading of fiber positioner data.
        """
        fiberpos = io.load_fiberpos()
        self.assertEqual(len(fiberpos), 5000)
        # Check workaround for astropy 1.0.x bug for lower -> upper col names.
        for key in ('FIBER', 'LOCATION', 'SPECTRO', 'X', 'Y', 'Z'):
            self.assertIn(key, fiberpos.dtype.names)
            x = fiberpos[key]

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_platescale(self):
        """Test loading platescale.txt file.
        """
        self.assertTrue(os.path.exists(io.findfile('focalplane/platescale.txt')))
        p1 = io.load_platescale()
        p2 = io.load_platescale()
        self.assertIs(p1, p2)  #- caching worked

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_get_focalplane_dates(self):
        dates = io.get_focalplane_dates()
        self.assertEqual(dates[:5], [
            datetime.datetime(2010,  1,  1, 0, 0, tzinfo=datetime.timezone.utc),
            datetime.datetime(2019,  9, 16, 0, 0, tzinfo=datetime.timezone.utc),
            datetime.datetime(2020,  3,  6, 0, 0, tzinfo=datetime.timezone.utc),
            datetime.datetime(2021,  1, 17, 0, 0, tzinfo=datetime.timezone.utc),
            datetime.datetime(2021,  1, 19, 0, 0, tzinfo=datetime.timezone.utc)])

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_focalplane(self):
        utc = datetime.timezone.utc
        fp,exclude,state,time = io.load_focalplane(
            time=datetime.datetime(2010, 1, 15, tzinfo=utc))
        self.assertEqual(type(fp), Table)
        self.assertEqual(len(fp), 5020)
        self.assertEqual(fp.colnames, ['PETAL', 'DEVICE', 'LOCATION', 'PETAL_ID', 'DEVICE_ID', 'DEVICE_TYPE', 'SLITBLOCK', 'BLOCKFIBER', 'CABLE', 'CONDUIT', 'FIBER', 'FWHM', 'FRD', 'ABS', 'OFFSET_X', 'OFFSET_Y', 'OFFSET_T', 'OFFSET_P', 'LENGTH_R1', 'LENGTH_R2', 'MAX_T', 'MIN_T', 'MAX_P', 'MIN_P'])
        self.assertEqual(type(exclude), dict)
        self.assertEqual(set(exclude.keys()), set(['legacy']))
        self.assertEqual(type(state), Table)
        self.assertEqual(len(state), 5020)
        self.assertEqual(set(state.colnames), set(['PETAL','DEVICE', 'LOCATION',
                                                   'STATE','EXCLUSION']))
        self.assertEqual(time, "2010-01-01T00:00:00")
        fp,exclude,state,time,time_low,time_high = io.load_focalplane(
            time=datetime.datetime(2010, 1, 15, tzinfo=utc),
            get_time_range=True)
        self.assertEqual(time_low,  datetime.datetime(2010, 1,  1, tzinfo=utc))
        self.assertEqual(time_high, datetime.datetime(2019, 9, 16, tzinfo=utc))

        # This one has state changes within the focalplane updates
        fp,exclude,state,time,time_low,time_high = io.load_focalplane(
            time=datetime.datetime(2021, 4, 8, tzinfo=utc),
            get_time_range=True)
        self.assertEqual(time_low,  datetime.datetime(2021,  4,  3, 23, 53, 23, tzinfo=utc))
        self.assertEqual(time_high, datetime.datetime(2021,  4, 10, 20,  0, 39, tzinfo=utc))

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_targets(self):
        """Test loading of tile files.
        """
        data = io.load_target_info()
        #- Test a few keys, but not everything
        self.assertIn('ntarget_lrg', data.keys())
        self.assertIn('nobs_elg', data.keys())
        self.assertIn('success_qso', data.keys())

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_pix_file(self):
        """Test loading of the file of HEALPixel weights.
        """
        nside = 16
        pix = io.load_pixweight(nside)
        npix = len(pix)
        #ADM Test the length of the returned array is appropriate to the passed nside
        self.assertEqual(npix,12*nside*nside)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_tiles_old(self):
        """Test loading of tile files for old default DESIMODEL case.
        """
        # starting clean
        self.assertEqual(io._tiles, {})
        t0 = io.load_tiles(cache=False, surveyops=False)
        self.assertEqual(len(io._tiles), 0)
        # loading tiles fills the cache with one items
        t1 = io.load_tiles(onlydesi=False, surveyops=False)
        self.assertEqual(len(io._tiles), 1)
        tile_cache_id1 = id(list(io._tiles.values())[0])
        # reloading, even with a filter, shouldn't change cache
        t2 = io.load_tiles(onlydesi=True, surveyops=False)
        self.assertEqual(len(io._tiles), 1)
        tile_cache_id2 = id(list(io._tiles.values())[0])
        self.assertEqual(tile_cache_id1, tile_cache_id2)
        #- Temporarily support OBSCONDITIONS as u2 (old) or i4 (new)
        self.assertTrue(np.issubdtype(t1['OBSCONDITIONS'].dtype, np.signedinteger) or
                        np.issubdtype(t1['OBSCONDITIONS'].dtype, np.unsignedinteger) )
        self.assertTrue(np.issubdtype(t2['OBSCONDITIONS'].dtype, np.signedinteger) or
                        np.issubdtype(t2['OBSCONDITIONS'].dtype, np.unsignedinteger) )
        self.assertLess(len(t2), len(t1))
        # All tiles in DESI are also in full set.
        self.assertTrue(np.all(np.isin(t2['TILEID'], t1['TILEID'])))
        # I think this is the exact same test as above, except using set theory.
        self.assertEqual(len(set(t2.TILEID) - set(t1.TILEID)), 0)
        t3 = io.load_tiles(onlydesi=False, surveyops=False)
        tile_cache_id3 = id(list(io._tiles.values())[0])
        self.assertEqual(tile_cache_id1, tile_cache_id3)
        self.assertTrue(np.issubdtype(t3['OBSCONDITIONS'].dtype, np.signedinteger) or
                        np.issubdtype(t3['OBSCONDITIONS'].dtype, np.unsignedinteger) )
        # Check for extra tiles.
        a = io.load_tiles(extra=False, surveyops=False)
        self.assertEqual(np.sum(np.char.startswith(a['PROGRAM'], 'EXTRA')), 0)
        b = io.load_tiles(extra=True, surveyops=False)
        self.assertGreater(np.sum(np.char.startswith(b['PROGRAM'], 'EXTRA')), 0)
        self.assertLess(len(a), len(b))

    @unittest.skipUnless(surveyops_available, surveyops_message)
    def test_load_tiles(self):
        """Test loading of tile files for default DESI_SURVEYOPS case.
        """
        # starting clean
        self.assertEqual(io._tiles, {})
        t0 = io.load_tiles(cache=False)
        self.assertEqual(len(io._tiles), 0)
        # loading tiles fills the cache with one items
        t1 = io.load_tiles(onlydesi=False)
        self.assertEqual(len(io._tiles), 1)
        tile_cache_id1 = id(list(io._tiles.values())[0])
        # reloading, even with a filter, shouldn't change cache
        t2 = io.load_tiles(onlydesi=True)
        self.assertEqual(len(io._tiles), 1)
        tile_cache_id2 = id(list(io._tiles.values())[0])
        self.assertEqual(tile_cache_id1, tile_cache_id2)
        self.assertLess(len(t2), len(t1))
        # All tiles in DESI are also in full set.
        self.assertTrue(np.all(np.isin(t2['TILEID'], t1['TILEID'])))
        t3 = io.load_tiles(onlydesi=False)
        tile_cache_id3 = id(list(io._tiles.values())[0])
        self.assertEqual(tile_cache_id1, tile_cache_id3)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_load_tiles_alt_old(self):
        """Test alternative tile-load for old default DESIMODEL case
        """
        # starting clean
        self.assertEqual(io._tiles, {})
        t1 = Table(io.load_tiles(surveyops=False))
        t1.write(self.testfile)
        #- no path; should fail since that file isn't in $DESIMODEL/data/footprint/
        if sys.version_info.major == 2:
            with self.assertRaises(IOError):
                t2 = io.load_tiles(tilesfile=os.path.basename(self.testfile),
                                   surveyops=False)
        else:
            with self.assertRaises(FileNotFoundError):
                t2 = io.load_tiles(tilesfile=os.path.basename(self.testfile),
                                   surveyops=False)

        #- with path, should work:
        t2 = Table(io.load_tiles(tilesfile=self.testfile, surveyops=False))
        self.assertTrue(np.all(t1 == t2))

        #- cache should have two items now
        self.assertEqual(len(io._tiles), 2)

    @unittest.skipUnless(surveyops_available, surveyops_message)
    def test_load_tiles_alt(self):
        """Test alternative tile-load for default DESI_SURVEYOPS case.
        """
        # starting clean
        self.assertEqual(io._tiles, {})
        t1 = Table(io.load_tiles())
        t1.write(self.testecsvfile)
        #- no path; should fail since that file isn't in $DESIMODEL/data/footprint/
        with self.assertRaises(FileNotFoundError):
            t2 = io.load_tiles(tilesfile=os.path.basename(self.testecsvfile))

        #- with path, should work:
        t2 = Table(io.load_tiles(tilesfile=self.testecsvfile))
        self.assertTrue(np.all(t1 == t2))

        #- cache should have two items now
        self.assertEqual(len(io._tiles), 2)

    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_tiles_consistency_old(self):
        """Test consistency of tile files.

        - Validate different tile loading schemes.
        - Test numerical consistency between FITS & ECSV versions.
        - Assure that PROGRAM column has no trailing whitespace.
        """
        fitstiles = io.findfile('footprint/desi-tiles.fits')
        ecsvtiles = io.findfile('footprint/desi-tiles.ecsv')
        # tf=TilesFits  tt=TilesTable  te=TilesEcsv
        tf = io.load_tiles(onlydesi=False, extra=True, surveyops=False)
        tt = Table.read(fitstiles)
        te = Table.read(ecsvtiles, format='ascii.ecsv')

        self.assertEqual(sorted(tf.dtype.names), sorted(tt.colnames))

        #- ECSV is a subset of the columns
        missing = set(te.colnames) - set(tf.dtype.names)
        self.assertEqual(len(missing), 0)

        for col in tt.colnames:
            self.assertTrue(np.all(tf[col]==tt[col]), 'fits[{col}] != table[{col}]'.format(col=col))
            if col in te.colnames:
                if np.issubdtype(tf[col].dtype, np.floating):
                    self.assertTrue(np.allclose(tf[col], te[col], atol=1e-4, rtol=1e-4), 'fits[{col}] != ecsv[{col}]'.format(col=col))
                else:
                    self.assertTrue(np.all(tf[col]==te[col]), 'fits[{col}] != ecsv[{col}]'.format(col=col))

        for program in set(tf['PROGRAM']):
            self.assertTrue((program[-1] != ' ') and (program[-1] != b' '))
        for program in set(tt['PROGRAM']):
            self.assertTrue((program[-1] != ' ') and (program[-1] != b' '))
        for program in set(te['PROGRAM']):
            self.assertTrue((program[-1] != ' ') and (program[-1] != b' '))


    @unittest.skipUnless(surveyops_available, surveyops_message)
    def validate_tile_load(self):
        """Validate different tile loading schemes.

        Also assure that PROGRAM column has no trailing whitespace.
        """
        ecsvtiles = io.findfile("tiles-main.ecsv", surveyops=True)
        tt = io.load_tiles(onlydesi=False, extra=True)
        te = Table.read(ecsvtiles, format='ascii.ecsv')

        self.assertEqual(sorted(tt.dtype.names), sorted(te.colnames))

        for program in set(tt['PROGRAM']):
            self.assertTrue((program[-1] != ' ') and (program[-1] != b' '))
        for program in set(te['PROGRAM']):
            self.assertTrue((program[-1] != ' ') and (program[-1] != b' '))


    @unittest.skipUnless(desimodel_available, desimodel_message)
    def test_trim_data(self):
        '''Test trimming data files for lightweight tests'''
        psffile = io.findfile('specpsf/psf-b.fits')
        if os.path.getsize(psffile) > 20e6:
            from .. import trim
            indir = io.datadir()
            trim.trim_data(indir, self.trimdir)
            self.assertTrue(os.path.isdir(self.trimdir))
            self.assertGreater(len(list(os.walk(self.trimdir))), 1)
