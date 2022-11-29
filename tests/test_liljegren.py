import unittest
from itertools import starmap

import numpy
from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt
from seus_hvi_wbgt.wbgt.liljegren import solar_parameters

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestSolarParams( unittest.TestCase ):

  def test_cosz(self):

    lats = [
      ( 33, 43, 59),
    ]
    lons = [ 
      (-84, 22, 59),
    ]
    lats = numpy.asarray( list( starmap( degMinSec2Frac, lats ) ) )
    lons = numpy.asarray( list( starmap( degMinSec2Frac, lons ) ) )

    years   = numpy.asarray( [2000] )
    months  = numpy.asarray( [   6] )
    days    = numpy.full( years.size,  1 )
    hours   = numpy.full( years.size, 16 )
    minutes = numpy.full( years.size,  0 )

    solar   = numpy.asarray( [1000] )

    for use_spa in [False, True]:
      tmp     = solar_parameters(
        lats, lons,
        years, months, days, hours, minutes,
        solar, use_spa=use_spa
      )
      print( tmp )

    #numpy.testing.assert_almost_equal(tmp[1], 9.1734928e-01)

"""
class TestWBGT( unittest.TestCase ):

  def setUp( self ):

    self.year   = 2000
    self.month  =    6
    self.day    =    1
    self.hour   =   16
    self.minute =    0

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
      self.solar,
      self.pres,
      self.Tair,
      self.Tdew,
      self.speed,
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

"""

