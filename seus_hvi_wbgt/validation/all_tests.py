import sys, os
from datetime import datetime, timedelta

import numpy
from pandas import to_datetime, date_range, read_pickle

from matplotlib import pyplot as plt

from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt

from ncsuCLOUDS import io
from ncsuCLOUDS.api import CLOUDS

#def addUnits( df, vName, newUnit, dt = None, dtype = numpy.float32 ):
def addUnits( df, params, vName ):
  """To wrap values in DataFrame into pint.Quantity"""

  unit = params[ params['requested'] == vName ]['unit'].values
  unit = '' if unit.size != 1 else unit[0]
  return units.Quantity( df.loc[:,vName,:]['value'].values, unit )

def main( *root ):

  if len(root) == 0: root = ('/Users/kwodzicki/Data/', )
  data, meta  = io.read( *root )
  locs        = meta['Location'].set_index( 'Location ID' )
  params      = meta['Parameter']
  raw_labels  = ['solar', 'pressure', 'T', 'RH', 'u']
  wbgt_labels = ['est. speed', 'globe', 'natural', 'psychrometric', 'WBGT']

  labels      = ['Pyschrometric Wet Bulb', 'Globe', 'Wet Bulb Globe']
  xLabel      = [f"ECONet {l}" for l in labels]

  fig, axes   = plt.subplots(3, 3)
  xyrange     = [float('inf'), -float('inf')]

  solar    = addUnits( data, params, 'rad2m_total')    # Add units to values, convert to appropriate units, and reindex
  pres     = addUnits( data, params, 'airpres|kPa')
  Tair     = addUnits( data, params, 'airtemp2m|K')
  Tdew     = addUnits( data, params, 'dewtemp2m|K')
  speed    = addUnits( data, params, 'windspeed2m|mph')

  Twbt     = addUnits( data, params, 'wetbulbtemp|K'     ).to( 'degree_Celsius' ).magnitude
  Tglobe   = addUnits( data, params, 'blackglobetemp2m|K').to( 'degree_Celsius' ).magnitude
  Twbgt    = addUnits( data, params, 'wbgt2m|K'          ).to( 'degree_Celsius' ).magnitude
  #TwbgtB   = addUnits( data, params, 'wbgt2m_bernard|K'  ).to( 'degree_Celsius' ).magnitude

  lon      = 'Longitude [degrees East]'
  lat      = 'Latitude [degrees North]'
  stations = data.index.get_level_values( data.index.names.index( 'location' ) )
  lon      = stations.map( locs[ lon ] ).values.astype( numpy.float32 ) 
  lat      = stations.map( locs[ lat ] ).values.astype( numpy.float32 )

  tmpVar   = data.index[0][1]
  dt       = data.loc[:,tmpVar,:]
  dt       = dt.index.get_level_values( dt.index.names.index( 'datetime' ) )

  for row, method in enumerate( ['liljegren' , 'dimiceli', 'bernard'] ):
    res = wbgt(method,
      lat, lon, dt,
      solar, pres, Tair, Tdew, speed, 
    )

    for col, (xx, yy) in enumerate( zip( [Twbt, Tglobe, Twbgt], [res['Tpsy'], res['Tg'], res['Twbg']] ) ):
      xyrange = [min( [numpy.nanmin(xx), numpy.nanmax(yy), xyrange[0]] ), 
                 max( [numpy.nanmin(xx), numpy.nanmax(yy), xyrange[1]] ) ]
      axes[row,col].scatter( xx, yy, label=method, alpha=0.6 )
      axes[row,col].set_title( labels[col] )
      axes[row,col].set_ylabel( method.title() )

  for col, ax in enumerate(axes[-1,:]): 
    ax.set_xlabel( 'ECONet' )

  xyrange = numpy.asarray( xyrange ) + [-2, 2]
  for ax in axes.flatten():
    ax.set_aspect( 'equal' ) 
    ax.set_xlim( xyrange )
    ax.set_ylim( xyrange )
    ax.plot( xyrange, xyrange, '-k' ) 
  
  plt.show()

if __name__ == "__main__":
  if len(sys.argv) >= 1: 
    main( *sys.argv[1:] )
  else:
    main()
