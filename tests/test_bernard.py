import unittest
from itertools import starmap

import pandas
import numpy
from metpy.units import units

from pywbgt import bernard 

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestBernard(unittest.TestCase):

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

        return bernard.wetbulb_globe(
            self.dates,
            numpy.resize(self.lats,  self.dates.size),
            numpy.resize(self.lons,  self.dates.size),
            numpy.resize(self.solar, self.dates.size),
            numpy.resize(self.pres,  self.dates.size),
            numpy.resize(self.Tair,  self.dates.size),
            numpy.resize(self.Tdew,  self.dates.size),
            numpy.resize(self.speed, self.dates.size),
            zspeed = self.zspeed,
        )

    def test_conv_heat_flow_coeff(self):

        ref_vals = [
             7.335649258049876, 17.239835857124724,  7.297282852987763,
            17.237462152087303,  7.301638680294556, 17.23687275664054,
             7.30179322137431,  17.236670700514786,  7.298272284705485,
            17.23822423913286,   7.300158403564784, 17.23957048141723,
        ]

        temp_g = self.compute_wbgt()['Tg']
        test_vals = bernard.conv_heat_trans_coeff(
            temp_g.to('degC').magnitude,
            numpy.resize(self.Tair, self.dates.size).to('degC').magnitude,
            numpy.resize(self.speed, self.dates.size).to('m/s').magnitude,
        )

        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_factor_e(self):

        speeds   = numpy.asarray(
            [0.09, 0.1, 0.5, 0.99, 1.0, 1.1]
        )
        ref_vals = [
            1.1, 1.0589254117941675, 0.014354692507258626, -0.0988883294140563, -0.1, -0.1,
        ]
 
        test_vals = bernard.factor_e(speeds)

        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_factor_c(self):
        
        speeds    = numpy.asarray(
            [0.029, 0.03, 1.5, 3.0, 3.1]
        ) 
        ref_vals  = [0.850, 0.8549213665756566, 0.972150296874842, 0.9929213665756567, 1.0]
        test_vals = bernard.factor_c(speeds)

        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_psychrometric_wetbulb(self):
        
        ref_vals = [17.96148884743604, 27.140425028129574] 

        test_vals = bernard.psychrometric_wetbulb(
            self.Tair.to('degC').magnitude,
            temp_dew=self.Tdew.to('degC').magnitude,
        )
        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_globe_temp(self):
        
        ref_vals = [
            43.37438061016138, 47.25171782800197, 40.871948729574626,
            46.48048629244579, 41.14928878477076, 46.291353905818596,
            41.15916066936984, 46.22673501041123, 40.9347933837588,
            46.72643804797883, 41.05484196199592, 47.16474553354357,
        ]

        test_vals = self.compute_wbgt()['Tg'].magnitude
        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_natural_wetbulb(self):

        ref_vals = [
            22.597532207220826, 30.103354485130065, 21.971924237074138,
            29.91054660124102, 22.04125925087317, 29.86326350458422,
            22.04372722202294, 29.84710878073238, 21.98763540062018,
            29.97203454012428, 22.01764754517946, 30.081611411515464,
        ] 

        test_vals = self.compute_wbgt()['Tnwb'].magnitude
        numpy.testing.assert_equal(test_vals, ref_vals)

    def test_wetbulb_globe(self):
        
        ref_vals = [
            26.993148667086857, 34.02269170519144, 26.054736711866823,
            33.73347987935787, 26.15873923256537, 33.662555234372675,
            26.16244118929003, 33.638323148594914, 26.07830345718589,
            33.82571178768276, 26.123321674024808, 33.99007709476954,
        ] 

        test_vals = self.compute_wbgt()['Twbg'].magnitude
        numpy.testing.assert_equal(test_vals, ref_vals)
