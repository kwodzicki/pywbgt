import unittest
from itertools import starmap

import pandas
import numpy
from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt
from seus_hvi_wbgt.wbgt.liljegren import solar_parameters

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestSolarParams( unittest.TestCase ):

  def setUp( self ):

    self.lats = ( 33, 43, 59)
    self.lons = (-84, 22, 59)

    self.months  = numpy.arange( 12 ) + 1
    self.years   = numpy.full( self.months.size, 2000 )
    self.days    = numpy.full( self.months.size,  1 )
    self.hours   = numpy.full( self.months.size, 16 )
    self.minutes = numpy.full( self.months.size,  0 )
    self.seconds = numpy.full( self.months.size,  0 )

    self.ref_COSZ = numpy.array(
      [0.479295 , 0.5436556, 0.6659838, 0.7999455, 0.8850405, 0.9173493, 
       0.9125813, 0.8857513, 0.8315094, 0.7408212, 0.617371 , 0.5131574]
    )

  def test_cosz(self):

    lats = numpy.full( self.years.size, degMinSec2Frac( *self.lats ) )
    lons = numpy.full( self.years.size, degMinSec2Frac( *self.lons ) )

    solar = numpy.full( self.years.size, 1000 )

    res = solar_parameters(
        lats, lons,
        self.years, 
        self.months,
        self.days,
        self.hours,
        self.minutes,
        self.seconds,
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
        lats, lons,
        dates.year.values, 
        dates.month.values,
        dates.day.values,
        dates.hour.values,
        dates.minute.values,
        dates.second.values,
        solar, 
        use_spa=use_spa,
      )
    )

    numpy.testing.assert_almost_equal(
      expAct[0][1,:], expAct[1][1,:], decimal=4
    )

class TestWBGT( unittest.TestCase ):

  def setUp( self ):

    self.year   = 2000
    self.month  =    6
    self.day    =    1
    self.hour   =   16
    self.minute =    0
    self.second =    0

    lats       = ( 33, 43, 59)  
    lons       = (-84, 22, 59)

    self.solar = units.Quantity( [500.0, 1000.0], 'watt/meter**2')
    self.solar = units.Quantity( [500.0,  805.0], 'watt/meter**2')
    self.pres  = units.Quantity( [985.0, 1013.0], 'hPa')
    self.Tair  = units.Quantity( [ 25.0,   35.0], 'degree_Celsius')
    self.Tdew  = units.Quantity( [ 15.0,   25.0], 'degree_Celsius')
    self.speed = units.Quantity( [  1.0,    5.0], 'mile/hour')

    self.lats  = numpy.full( self.solar.size, degMinSec2Frac(*lats) )
    self.lons  = numpy.full( self.solar.size, degMinSec2Frac(*lons) )

    self.res   = wbgt('liljegren', self.lats, self.lons,
      numpy.full( self.solar.size, self.year ),
      numpy.full( self.solar.size, self.month),
      numpy.full( self.solar.size, self.day  ),
      numpy.full( self.solar.size, self.hour ),
      numpy.full( self.solar.size, self.minute ),
      numpy.full( self.solar.size, self.second ),
      self.solar,
      self.pres,
      self.Tair,
      self.Tdew,
      self.speed,
      avg = 1.0
     )

  def test_Tg(self):
    
    Tg   = [41.970542908, 48.868949890] 
    numpy.testing.assert_almost_equal(Tg, self.res['Tg'])

  def test_Tpsy(self):
    
    Tpsy = [18.278467178, 27.251123428]
    numpy.testing.assert_almost_equal(Tpsy, self.res['Tpsy'])

  def test_Tnwb(self):

    Tnwb = [23.320886612, 29.015771866]
    numpy.testing.assert_almost_equal(Tnwb, self.res['Tnwb'])

  def test_Twbg(self):

    Twbg = [27.218729019, 33.584831238]  
    numpy.testing.assert_almost_equal(Twbg, self.res['Twbg'])
