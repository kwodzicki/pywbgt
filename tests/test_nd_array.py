import unittest

import pandas
import numpy as np
import xarray as xr

import pywbgt


class TestNDArray(unittest.TestCase):
    """
    Test of multidimension/xarray support

    """

    def setUp(self):

        dims = ('t', 'y', 'x')
        dates = pandas.date_range(
            start='2000-06-01T12',
            end='2000-06-02T12',
            freq='1D',
        )
        lats = np.asarray((33, 43, 59))
        lons = np.asarray((-84, 22, 50, 120))

        shape = (len(dates), len(lats), len(lons))

        self.data = xr.Dataset(
            data_vars=dict(
                solar=(
                    dims,
                    np.random.random(shape) * 1000,
                    {'units': 'watt/meter**2'},
                ),
                pres=(
                    dims,
                    np.random.random(shape) * 40 + 985,
                    {'units': 'hPa'},
                ),
                temp_air=(
                    dims,
                    np.random.random(shape) * 10 + 25,
                    {'units': 'degree_Celsius'},
                ),
                temp_dew=(
                    dims,
                    np.random.random(shape) * 10 + 15,
                    {'units': 'degree_Celsius'},
                ),
                speed=(
                    dims,
                    np.random.random(shape) * 4 + 1,
                    {'units': 'mile/hour'},
                ),
            ),
            coords=dict(
                t=(('t'), dates, {'axis': 'T'}),
                y=(('y'), lats, {'axis': 'Y'}),
                x=(('x'), lons, {'axis': 'X'}),
            ),
        )

        self.new_dim = 'stacked'
        self.ref_ds = self.data.metpy.quantify().stack({self.new_dim: dims})

    def check_method(self, method):

        ref = pywbgt.wbgt(
            self.ref_ds.t.data,
            self.ref_ds.y.data,
            self.ref_ds.x.data,
            self.ref_ds.solar.data,
            self.ref_ds.pres.data,
            self.ref_ds.temp_air.data,
            self.ref_ds.temp_dew.data,
            self.ref_ds.speed.data,
            method=method,
        )

        ref = self.ref_ds.assign(
            {
                key: (self.new_dim, val)
                for key, val in zip(pywbgt.OUTPUT_ORDER, ref)
                if val.shape == ref[0].shape
            }
        ).unstack().metpy.dequantify()

        res = pywbgt.wbgt(self.data, method=method).metpy.dequantify()

        for var in res.data_vars:
            np.testing.assert_equal(res[var].values, ref[var].values)

    def test_liljegren(self):
        """Test Liljegren method for multidim"""
        self.check_method('liljegren')

    def test_dimiceli(self):
        """Test Dimiceli method for multidim"""
        self.check_method('dimiceli')

    def test_dimiceli_nws(self):
        """Test Dimiceli NWS method for multidim"""
        self.check_method('dimiceli_nws')

    def test_bernard(self):
        """Test Bernard method for multidim"""
        self.check_method('bernard')

    def test_map_blocks(self):
        """Test mapping blocks to function"""

        method = 'liljegren'

        ref = pywbgt.wbgt(
            self.ref_ds.t.data,
            self.ref_ds.y.data,
            self.ref_ds.x.data,
            self.ref_ds.solar.data,
            self.ref_ds.pres.data,
            self.ref_ds.temp_air.data,
            self.ref_ds.temp_dew.data,
            self.ref_ds.speed.data,
            method=method,
        )

        ref = self.ref_ds.assign(
            {
                key: (self.new_dim, val)
                for key, val in zip(pywbgt.OUTPUT_ORDER, ref)
                if val.shape == ref[0].shape
            }
        ).unstack().metpy.dequantify()

        res = self.data.map_blocks(
            pywbgt.wbgt,
            kwargs={'method': method},
        ).metpy.dequantify()

        for var in res.data_vars:
            np.testing.assert_equal(res[var].values, ref[var].values)
