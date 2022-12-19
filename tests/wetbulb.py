import unittest

import numpy as np
from metpy.units import units

from seus_hvi_wbgt.wbgt.wetbulb import wetBulb

class TestWetBulb( unittest.TestCase ):

  def setUp( self ):

    self.Ta   = units.Quantity( np.arange( 10, 26 ), 'degC' )
    self.Td   = units.Quantity( np.full( self.Ta.shape, 8), self.Ta.units )
    self.p    = units.Quantity( np.full( self.Ta.shape, 1013.25 ), 'hPa' )

  def test_wbt(self):

    stull = wetBulb( self.Ta.m, self.Td.m, method='stull' )
    dimi  = wetBulb( self.Ta.m, self.Td.m, method='dimiceli' )
    irib  = wetBulb( self.Ta.m, self.Td.m, self.p.m, method='iribarne' )

    np.testing.assert_almost_equal( stull, dimi, decimal=2 )
  
