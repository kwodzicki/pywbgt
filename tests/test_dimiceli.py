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
            zspeed          = self.zspeed,
            natural_wetbulb = 'malchaire',
            min_speed       = dimiceli.DIMICELI_MIN_SPEED,
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
            5.804440926414478e+09, 7.580495185370829e+09, 7.253425926414478e+09,
            7.016682845370829e+09, 8.764311926414478e+09, 8.065993905370829e+09,
            8.826937926414478e+09, 7.643289210370829e+09, 8.386207926414478e+09,
            6.381148627870829e+09, 6.773091926414478e+09, 7.929456245370829e+09,
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
        ss, _ = dimiceli.adjust_speed_2m(
            self.speed,
            zspeed    = self.zspeed,
            min_speed = dimiceli.DIMICELI_MIN_SPEED,
        )
 
        test_vals = dimiceli.factor_c(
            ss.to('meter/hour').magnitude
        )

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
            38.31750513722974, 40.68555674157649, 41.64127240504196,
            41.51158271804979, 45.10703183815943, 42.48553616231467,
            45.25068705123203, 42.09318854558204, 44.2397145788166,
            40.92169020316264, 40.53945394120878, 40.60748977555843, 
        ] 

        test_vals = self.compute_wbgt()['Tg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_natural_wetbulb(self):

        ref_vals = [
            24.15500188961772, 29.92728729928154, 25.30881184626314,
            30.18095989953724, 26.51191285155729, 30.48006102924928,
            26.56178120930073, 30.35957106528275, 26.21083301991828,
            29.99980389187549, 24.92632745083094, 29.90331293185496,
        ] 

        test_vals = self.compute_wbgt()['Tnwb'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_wetbulb_globe(self):

        ref_vals = [
            27.07200235017836, 32.58621245781237, 28.54442277339258,
            32.92898847328603, 30.07974536372199, 33.33314995293743,
            30.14338425675692, 33.17033745481433, 29.69552602970612,
            32.68420076494537, 28.05632000382341, 32.55381700741015,            
        ] 

        test_vals = self.compute_wbgt()['Twbg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)
        
