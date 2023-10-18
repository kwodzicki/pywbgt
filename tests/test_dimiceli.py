import unittest
from itertools import starmap

import pandas
import numpy
from metpy.units import units

from pywbgt import dimiceli

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestDimiceli(unittest.TestCase):

    def setUp( self ):

        lats = ( 33, 43, 59)
        lons = (-84, 22, 59)

        self.dates   = pandas.date_range(
            '20000101T16', 
            '20010101T16', 
            freq      = 'MS',
            inclusive = 'left'
        )

        self.solar  = units.Quantity( [500.0,  805.0], 'watt/meter**2')
        self.pres   = units.Quantity( [985.0, 1013.0], 'hPa')
        self.Tair   = units.Quantity( [ 25.0,   35.0], 'degree_Celsius')
        self.Tdew   = units.Quantity( [ 15.0,   25.0], 'degree_Celsius')
        self.speed  = units.Quantity( [  1.0,    5.0], 'mile/hour')
        self.zspeed = units.Quantity( 2.0, 'meters' )
        self.lats   = numpy.full( self.solar.size, degMinSec2Frac(*lats) )
        self.lons   = numpy.full( self.solar.size, degMinSec2Frac(*lons) )

    def compute_wbgt(self):

        return dimiceli.wetbulb_globe(
            self.dates,
            numpy.resize(self.lats,  self.dates.size),
            numpy.resize(self.lons,  self.dates.size),
            numpy.resize(self.solar, self.dates.size),
            numpy.resize(self.pres,  self.dates.size),
            numpy.resize(self.Tair,  self.dates.size),
            numpy.resize(self.Tdew,  self.dates.size),
            numpy.resize(self.speed, self.dates.size),
            zspeed = self.zspeed,
            natural_wetbulb = 'malchaire',
        )

    def test_conv_heat_flow_coeff(self):

        numpy.testing.assert_equal(
            dimiceli.conv_heat_flow_coeff(),
            0.315,
        )

    def test_atmos_vapor_pres(self):

        ref_vals  = [16.053150020118238, 29.254910801080698]
        test_vals = dimiceli.atmospheric_vapor_pressure(
            self.Tair.to('degC').magnitude,
            self.Tdew.to('degC').magnitude,
            self.pres.to('hPa').magnitude,
        )

        numpy.testing.assert_equal(test_vals, ref_vals)
        
    def test_thermal_emiss(self):

        ref_vals  = [0.8548516210648657, 0.9313755074245822]
        test_vals = dimiceli.thermal_emissivity(
            self.Tair.to('degC').magnitude,
            self.Tdew.to('degC').magnitude,
            self.pres.to('hPa').magnitude,
        )

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_factor_b(self):

        ref_vals = [
            5804600426.414478, 7580461375.370829, 7253912926.414478,
            7017295450.370829, 8764512926.414478, 8066641125.370829,
            8827311926.414478, 7644227035.370829, 8386488926.414478,
            6381388115.370829, 6772964426.414478, 7929655080.370829,
        ]
 
        solar, cosz, f_db = dimiceli.solar_parameters(
            self.dates,
            numpy.resize(self.lats, self.dates.size),
            numpy.resize(self.lons, self.dates.size),
            numpy.resize(self.solar.to('watt/m**2').magnitude, self.dates.size),
        )

        test_vals = dimiceli.factor_b(
            numpy.resize(self.Tair.to('degC').magnitude, self.dates.size),
            numpy.resize(self.Tdew.to('degC').magnitude, self.dates.size),
            numpy.resize(self.pres.to('hPa').magnitude, self.dates.size),
            numpy.resize(self.solar.to('watt/m**2').magnitude, self.dates.size),
            f_db,
            cosz,
        )

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_factor_c(self):

        ref_vals  = [435690588.0876065, 1077116913.6427102]
        test_vals = dimiceli.factor_c(self.speed.to('meter/hour').magnitude)

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_psychrometric_wetbulb(self):

        ref_vals = [18.59627755506145, 27.445752527955726]

        test_vals = dimiceli.psychrometric_wetbulb(
            self.Tair.to('degC').magnitude,
            self.Tdew.to('degC').magnitude,
        )

        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_speed(self):

        ref_vals  = (
            numpy
            .resize(self.speed, self.dates.size)
            .to('meter per second')
        )
        ref_vals  = numpy.clip(
            ref_vals,
            dimiceli.MIN_SPEED,
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
            38.31786982410732,  40.68589666196336,  41.64238894549278,
            41.5121511377606,   45.1074916009044,   42.48613705731751,
            45.25154565182722,  42.09405894708196,  44.24035917154351,
            40.921912545396204, 40.539161538415684, 40.607525521450356, 
        ] 

        test_vals = self.compute_wbgt()['Tg'].magnitude
        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_natural_wetbulb(self):

        ref_vals = [
            24.155128486733737, 29.927391688841762, 25.30919944123756,
            30.181134461246167, 26.512072453235177, 30.4802455641213,
            26.562079263236427, 30.359838365607484, 26.211056783330807,
            29.99987217318158,  24.926225946355192, 29.903323909419353,
        ] 

        test_vals = self.compute_wbgt()['Tnwb'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_wetbulb_globe(self):

        ref_vals = [
            27.072163905535078, 32.5863535145819,   28.54491739796485,
            32.929224350424434, 30.079949037445502, 33.33339930634841,
            30.143764614630943, 33.17069864534163,  29.695811582640268,
            32.68429303030634,  28.056190470131774, 32.55383184088362,
        ] 

        test_vals = self.compute_wbgt()['Twbg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)
        
