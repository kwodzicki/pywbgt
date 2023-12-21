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
            min_speed = units.Quantity(0, 'm/s'),
        )

    def test_conv_heat_flow_coeff(self):

        ref_vals = [
             7.33563906422347, 17.23983369440773,  7.29728047905653,
            17.23746167679371,  7.30163605847636, 17.23687164492864,
             7.30178758577696, 17.23666888682254,  7.29826760007174,
            17.23822295266038,  7.30015437019847, 17.23957011909292,
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

#        test_vals = bernard.psychrometric_wetbulb(
#            self.Tair.to('degC').magnitude,
#            temp_dew=self.Tdew.to('degC').magnitude,
#        )
        test_vals = bernard.psychrometric_wetbulb(
            self.Tair.to('degC').magnitude,
            temp_dew=self.Tdew.to('degC'),
        )
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=14)

    def test_globe_temp(self):
        
        ref_vals = [
            43.373697999106, 47.251008278667, 40.871798055794, 46.480333391746,
            41.149121325362, 46.290998065833, 41.158800636094, 46.226155488287,
            40.934495622227, 46.726021531819, 41.054584894706, 47.164626916066,
        ]

        test_vals = self.compute_wbgt()['Tg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=12)

    def test_natural_wetbulb(self):

        ref_vals = [
            22.597361554457 , 30.1031770977963, 21.9718865686289,
            29.910508376066 , 22.0412173860209, 29.8631745445878,
            22.043637213704 , 29.8469639002014, 21.9875609602373,
            29.9719304110844, 22.0175832783571, 30.081581757146,
        ] 

        test_vals = self.compute_wbgt()['Tnwb'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=13)

    def test_wetbulb_globe(self):
        
        ref_vals = [
            26.9928926879412, 34.0224256241908, 26.054680209199 ,
            33.7334225415954, 26.158676435287 , 33.6624217943781,
            26.1623061768116, 33.6381058277984, 26.0781917966116,
            33.8255555941229, 26.1232252737912, 33.9900326132153,
        ] 

        test_vals = self.compute_wbgt()['Twbg'].magnitude
        numpy.testing.assert_almost_equal(test_vals, ref_vals, decimal=13)
