import unittest

import pandas
import numpy as np
import xarray as xr

from pvlib import spa
import pywbgt


class TestSolar(unittest.TestCase):
    """
    Test of multidimension/xarray support

    """

    def setUp(self):

        dims = ('time', 'latitude', 'longitude')
        dates = pandas.date_range(
            start='2000-01-15T00',
            end='2000-12-15T00',
            freq='30D',
        )
        lats = np.linspace(-90, 90, 15)
        lons = np.linspace(-180, 180, 15)

        shape = (len(dates), len(lats), len(lons))

        self.elev = 0.0
        self.pressure = 101325.0
        self.temp = 12.0
        self.delta_t = 0
        self.atmos_refract = 0

        tmp = xr.Dataset(
            data_vars={
                'solar': (dims, np.full(shape, 1400)),
            },
            coords={
                'time': (dims[0], dates, {'axis': 'T'}),
                'latitude': (dims[1], lats, {'axis': 'Y'}),
                'longitude': (dims[2], lons, {'axis': 'X'}),
            },
        )

        self.new_dim = 'stacked'
        self.ref_ds = tmp.stack({self.new_dim: dims[1:]})

    def test_cza_by_location(self):
        """Test cosine of zenith looping over locations"""

        subset = self.ref_ds
        for i in range(subset[self.new_dim].size):
            tmp = subset.isel({self.new_dim: i})
            unixtime = tmp.time.astype('datetime64[s]').astype('float').data
            lat = tmp.latitude.data.astype(float)
            lon = tmp.longitude.data.astype(float)
            solar = tmp.solar.data
            theta, theta0, e, e0, phi, eot = spa.solar_position(
                unixtime,
                lat,
                lon,
                self.elev,
                self.pressure,
                self.temp,
                self.delta_t,
                self.atmos_refract,
            )
            cza_ref = np.cos(np.deg2rad(90.0 - e))

            solar, cza, fdir = pywbgt.solar.solar_parameters(
                tmp.time,
                lat,
                lon,
                solar,
                avg=0,
                elev=np.asarray(self.elev),
                pressure=np.asarray(self.pressure),
                temp=np.asarray(self.temp),
            )
            np.testing.assert_almost_equal(cza, cza_ref, decimal=12)

    def test_cza_array(self):
        """Test cosine of zenith array operation"""

        subset = self.ref_ds
        cza_ref = []
        for i in range(subset[self.new_dim].size):
            tmp = subset.isel({self.new_dim: i})
            unixtime = tmp.time.astype('datetime64[s]').astype('float').data
            lat = tmp.latitude.data
            lon = tmp.longitude.data
            solar = tmp.solar
            theta, theta0, e, e0, phi, eot = spa.solar_position(
                unixtime,
                lat,
                lon,
                self.elev,
                self.pressure,
                self.temp,
                self.delta_t,
                self.atmos_refract,
            )
            cza_ref.append(
                np.cos(np.deg2rad(90.0 - e))
            )

        sizes = {val: key for key, val in subset.sizes.items()}
        cza_ref = np.asarray(cza_ref)
        cza_ref = xr.DataArray(
            data=cza_ref,
            dims=[sizes.get(size) for size in cza_ref.shape],
            coords=subset.coords,
        ).unstack()

        subset = self.ref_ds.unstack()
        # Ensure transposed data are contiguous with copy
        time = subset.time.broadcast_like(subset['solar']).copy()
        lat = subset.latitude.broadcast_like(subset['solar']).copy()
        lon = subset.longitude.broadcast_like(subset['solar']).copy()

        solar, cza, fdir = pywbgt.solar.solar_parameters(
            time.data.ravel(),
            lat.data.ravel(),
            lon.data.ravel(),
            subset['solar'].data.ravel(),
            avg=0,
            elev=np.asarray(self.elev),
            pressure=np.asarray(self.pressure),
            temp=np.asarray(self.temp),
        )

        cza = cza.reshape(time.shape)
        np.testing.assert_almost_equal(cza, cza_ref, decimal=12)
