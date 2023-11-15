import unittest
from itertools import starmap

import pandas
import numpy
from metpy.units import units

from pywbgt import wbgt
from pywbgt.liljegren import solar_parameters

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestSolarParams( unittest.TestCase ):

  def setUp( self ):

    self.lats = ( 33, 43, 59)
    self.lons = (-84, 22, 59)

    self.dates   = pandas.date_range(
      '20000101T16', 
      '20010101T16', 
      freq      = 'MS',
      inclusive = 'left'
    )

    self.ref_COSZ = numpy.array(
      [0.479295 , 0.5436556, 0.6659838, 0.7999455, 0.8850405, 0.9173493, 
       0.9125813, 0.8857513, 0.8315094, 0.7408212, 0.617371 , 0.5131574]
    )

  def test_cosz(self):

    lats = numpy.full( self.dates.size, degMinSec2Frac( *self.lats ) )
    lons = numpy.full( self.dates.size, degMinSec2Frac( *self.lons ) )

    solar = numpy.full( self.dates.size, 1000 )

    res = solar_parameters(
        self.dates, lats, lons,
        solar, 
        use_spa=False,
    )
    numpy.testing.assert_almost_equal( res[1,:], self.ref_COSZ )
  
  def test_compare_spa(self):
    dates = pandas.date_range('2000-06-01', '2000-07-01', freq='1H')

    lats  = numpy.full( dates.size, degMinSec2Frac( *self.lats ) )
    lons  = numpy.full( dates.size, degMinSec2Frac( *self.lons ) )

    solar = numpy.full( dates.size, 1000 )

    expAct = []
    for use_spa in [False, True]:
      expAct.append(
        solar_parameters(
        dates,
        lats, lons,
        solar, 
        use_spa=use_spa,
      )
    )

    numpy.testing.assert_almost_equal(
      expAct[0][1,:], expAct[1][1,:], decimal=4
    )

class TestWBGT( unittest.TestCase ):

  def setUp( self ):

    dates       = pandas.to_datetime( ['2000-06-01T16:00:00'] ) 
    lats        = ( 33, 43, 59)  
    lons        = (-84, 22, 59)

    self.solar  = units.Quantity( [500.0, 1000.0], 'watt/meter**2')
    self.solar  = units.Quantity( [500.0,  805.0], 'watt/meter**2')
    self.pres   = units.Quantity( [985.0, 1013.0], 'hPa')
    self.Tair   = units.Quantity( [ 25.0,   35.0], 'degree_Celsius')
    self.Tdew   = units.Quantity( [ 15.0,   25.0], 'degree_Celsius')
    self.speed  = units.Quantity( [  1.0,    5.0], 'mile/hour')
    self.zspeed = units.Quantity( 2.0, 'meters' )
    self.lats   = numpy.full( self.solar.size, degMinSec2Frac(*lats) )
    self.lons   = numpy.full( self.solar.size, degMinSec2Frac(*lons) )

    self.res    = wbgt('liljegren', 
      dates.repeat( self.solar.size ), self.lats, self.lons,
      self.solar,
      self.pres,
      self.Tair,
      self.Tdew,
      self.speed,
      avg = 1.0,
      zspeed = self.zspeed,
      min_speed = units.Quantity(0.13, 'm/s'),
     )

  def test_Tg(self):
    
    Tg   = [41.970542908, 48.868949890] 
    numpy.testing.assert_almost_equal(Tg, self.res['Tg'].magnitude)

  def test_Tpsy(self):
    
    Tpsy = [18.278467178, 27.251123428]
    numpy.testing.assert_almost_equal(Tpsy, self.res['Tpsy'].magnitude)

  def test_Tnwb(self):

    Tnwb = [23.320886612, 29.015771866]
    numpy.testing.assert_almost_equal(Tnwb, self.res['Tnwb'].magnitude)

  def test_Twbg(self):

    Twbg = [27.218729019, 33.584831238]  
    numpy.testing.assert_almost_equal(Twbg, self.res['Twbg'].magnitude)
