import unittest
from itertools import starmap
from datetime import datetime, timedelta

import numpy
from metpy.units import units

from pywbgt import spa

def degMinSec2Frac( degree, minute, second ):

  return degree + (minute + second/60.0)/60.0

class TestSPA( unittest.TestCase ):

  def test_cosz(self):

    lat      =   39.742476
    lon      = -105.178600

    year     = 2003
    month    = 10
    day      = 17
    hour     = 12
    minute   = 30
    second   = 30
    timezone = -7.0

    delta_ut1     = numpy.asarray( [    0.0000] )#, 0.0, 1.0, 2.0000] );
    delta_t       = numpy.asarray( [   67.0000] )#, 0.0, 1.0, 2.0000] );
    elevation     = numpy.asarray( [ 1830.1400] )#, 0.0, 1.0, 2.0000] );
    pressure      = numpy.asarray( [  820.0000] )#, 0.0, 1.0, 2.0000] );
    temperature   = numpy.asarray( [   11.0000] )#, 0.0, 1.0, 2.0000] );
    slope         = numpy.asarray( [   30.0000] )#, 0.0, 1.0, 2.0000] );
    azm_rotation  = numpy.asarray( [  -10.0000] )#, 0.0, 1.0, 2.0000] );
    atmos_refract = numpy.asarray( [    0.5667] )#, 0.0, 1.0, 2.0000] );

    lat           = numpy.full(  delta_ut1.shape, lat )
    lon           = numpy.full(  delta_ut1.shape, lon )

    years         = numpy.full(  delta_ut1.shape, year )
    months        = numpy.full(  delta_ut1.shape, month )
    days          = numpy.full(  delta_ut1.shape, day )
    hours         = numpy.full(  delta_ut1.shape, hour )
    minutes       = numpy.full(  delta_ut1.shape, minute )
    seconds       = numpy.full(  delta_ut1.shape, second )
    timezone      = numpy.full(  delta_ut1.shape, timezone )
    solar         = numpy.zeros( delta_ut1.shape, dtype=numpy.float32 ) 
     
    spaRes = spa.solar_parameters(
      lat, lon,
      years, months, days, hours, minutes, seconds,
      function      = spa.SPA_ALL, 
      delta_ut1     = delta_ut1, 
      delta_t       = delta_t,   
      elevation     = elevation,
      pressure      = pressure,
      temperature   = temperature,
      slope         = slope,
      azm_rotation  = azm_rotation,
      atmos_refract = atmos_refract,
      timezone      = timezone
    )

    numpy.testing.assert_almost_equal(spaRes['zenith'],     50.111622, decimal=6)
    numpy.testing.assert_almost_equal(spaRes['azimuth'],   194.340241, decimal=6) 
    numpy.testing.assert_almost_equal(spaRes['incidence'],  25.187000, decimal=6)

