import unittest

import pandas
import numpy
from metpy.units import units

from pywbgt import dimiceli


def degMinSec2Frac(degree, minute, second):

    return degree + (minute + second / 60.0) / 60.0


class TestDimiceli(unittest.TestCase):

    def setUp(self):

        lats = (33, 43, 59)
        lons = (-84, 22, 59)

        self.dates = pandas.date_range(
            '20000101T16',
            '20010101T16',
            freq='MS',
            inclusive='left'
        )

        self.solar = units.Quantity([500.0, 805.0], 'watt/meter**2')
        self.pres = units.Quantity([985.0, 1013.0], 'hPa')
        self.Tair = units.Quantity([25.0, 35.0], 'degree_Celsius')
        self.Tdew = units.Quantity([15.0, 25.0], 'degree_Celsius')
        self.speed = units.Quantity([1.0, 5.0], 'mile/hour')
        self.zspeed = units.Quantity(2.0, 'meters')
        self.lats = numpy.full(self.solar.size, degMinSec2Frac(*lats))
        self.lons = numpy.full(self.solar.size, degMinSec2Frac(*lons))

    def compute_wbgt(self):

        return dimiceli.wetbulb_globe(
            self.dates.values,
            numpy.resize(self.lats, self.dates.size),
            numpy.resize(self.lons, self.dates.size),
            numpy.resize(self.solar, self.dates.size),
            numpy.resize(self.pres, self.dates.size),
            numpy.resize(self.Tair, self.dates.size),
            numpy.resize(self.Tdew, self.dates.size),
            numpy.resize(self.speed, self.dates.size),
            zspeed=self.zspeed,
            natural_wetbulb='malchaire',
            min_speed=dimiceli.DIMICELI_MIN_SPEED,
        )

    def test_conv_heat_flow_coeff(self):

        numpy.testing.assert_equal(
            dimiceli.conv_heat_flow_coeff(),
            0.315,
        )

    def test_atmos_vapor_pres(self):

        ref_vals = [16.053150020118238, 29.254910801080698]
        test_vals = dimiceli.atmospheric_vapor_pressure(
            self.Tair.to('degC').magnitude,
            self.Tdew.to('degC').magnitude,
            self.pres.to('hPa').magnitude,
        )

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_thermal_emiss(self):

        ref_vals = [0.8548516210648657, 0.9313755074245822]
        test_vals = dimiceli.thermal_emissivity(
            self.Tair.to('degC').magnitude,
            self.Tdew.to('degC').magnitude,
            self.pres.to('hPa').magnitude,
        )

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_factor_b(self):

        ref_vals = [
            5.8044407892149286e+09,
            7.5804950605705595e+09,
            7.2534252044102097e+09,
            7.0166825107848835e+09,
            8.7643116674216957e+09,
            8.0659939244838324e+09,
            8.8269380062377396e+09,
            7.6432888060895243e+09,
            8.3862067629410048e+09,
            6.3811488840500526e+09,
            6.7730918917926073e+09,
            7.9294556381097021e+09,
        ]

        solar, cosz, f_db = dimiceli.solar_parameters(
            self.dates.values,
            numpy.resize(self.lats, self.dates.size),
            numpy.resize(self.lons, self.dates.size),
            numpy.resize(
                self.solar.to('watt/m**2').magnitude,
                self.dates.size,
            ),
        )

        test_vals = dimiceli.factor_b(
            numpy.resize(self.Tair.to('degC').magnitude, self.dates.size),
            numpy.resize(self.Tdew.to('degC').magnitude, self.dates.size),
            numpy.resize(self.pres.to('hPa').magnitude, self.dates.size),
            numpy.resize(
                self.solar.to('watt/m**2').magnitude,
                self.dates.size,
            ),
            f_db,
            cosz,
        )

        # We decrease precision here as numbers are very large making
        # decimal values less important
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=5)

    def test_factor_c(self):

        ref_vals = [435690588.0876065, 1077116913.6427102]
        ss, _ = dimiceli.adjust_speed_2m(
            self.speed,
            zspeed=self.zspeed,
            min_speed=dimiceli.DIMICELI_MIN_SPEED,
        )

        test_vals = dimiceli.factor_c(
            ss.to('meter/hour').magnitude
        )

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_psychrometric_wetbulb(self):

        # Vals for metpy>1.6.3
        # ref_vals = [18.60753537424834, 27.48209036683773]

        # Vals for metpy==1.6.3
        ref_vals = [18.59627755506145, 27.44575252795573]
        test_vals = dimiceli.psychrometric_wetbulb(self.Tair, self.Tdew)

        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_speed(self):

        ref_vals = (
            numpy
            .resize(self.speed, self.dates.size)
            .to('meter per second')
        )
        ref_vals = numpy.clip(
            ref_vals,
            dimiceli.DIMICELI_MIN_SPEED,
            None,
        ).magnitude

        test_vals = (
            self.compute_wbgt()
            ['speed']
            .to('meter per second')
            .magnitude
        )
        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_globe_temp(self):

        ref_vals = [
            38.31750482251332, 40.68555668641199, 41.64127074886598,
            41.51158240749257, 45.10703124406651, 42.48553618005505,
            45.25068723433534, 42.09318817033472, 44.23971190997252,
            40.92169044094403, 40.53945386179109, 40.60748956330349,
        ]

        test_vals = self.compute_wbgt()['Tg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=13)

    def test_natural_wetbulb(self):

        # Vals for metpy>1.6.3
        # ref_vals = [
        #     24.16624303085751, 29.96360332651852, 25.32004951539023,
        #     30.21727416354683, 26.52314775449753, 30.51637340756263,
        #     26.57301625209539, 30.39588412315881, 26.22206798714569,
        #     30.03611952744798, 24.93756666393694, 29.9396290700784,
        # ]

        # Vals for metpy==.16.3
        ref_vals = [
            24.15500178036731, 29.92728728234052, 25.30881127133950,
            30.18095980416511, 26.51191264532434, 30.48006103469735,
            26.56178127286306, 30.35957095004429, 26.21083209345783,
            29.99980396489816, 24.92632742326195, 29.90331286667146,
        ]
        test_vals = self.compute_wbgt()['Tnwb'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_wetbulb_globe(self):

        # Vals for metpy>1.6.3
        # ref_vals = [
        #     27.07987108610292, 32.61163366584536, 28.55228881054636,
        #     32.95440839598130, 30.08760967696157, 33.35856862130485,
        #     30.15124882333384, 33.19575652027811, 29.70338997299648,
        #     32.70962175740239, 28.06418743711408, 32.57923826171557,
        # ]

        # Vals for metpy==1.6.3
        ref_vals = [
            27.07200221075978, 32.58621243492076, 28.54442203971085,
            32.92898834441409, 30.07974510054034, 33.33314996029915,
            30.14338433787120, 33.17033729909794, 29.69552484741498,
            32.68420086361751, 28.05631996864158, 32.55381691933071,
        ]

        test_vals = self.compute_wbgt()['Twbg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)
