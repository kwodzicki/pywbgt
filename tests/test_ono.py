import unittest
from itertools import starmap

import pandas
import numpy
from metpy.units import units

from seus_hvi_wbgt.wbgt.ono import ono 

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestWBGT( unittest.TestCase ):

  def setUp( self ):

    dates      = pandas.to_datetime( ['2000-06-01T16:00:00'] ) 
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

    self.res   = ono(self.lats, self.lons,
      dates.repeat( self.solar.size ),
      self.solar,
      self.pres,
      self.Tair,
      self.Tdew,
      self.speed,
      avg = 1.0
     )

  def test_Twbg(self):

    Twbg = [27.218729019, 33.584831238]  
    numpy.testing.assert_almost_equal(Twbg, self.res['Twbg'])
