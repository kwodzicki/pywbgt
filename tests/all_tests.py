from datetime import datetime, timedelta

import numpy
from pandas import to_datetime, date_range, read_pickle

from matplotlib import pyplot as plt

from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt

from ncsuCLOUDS.api import CLOUDS

c       = CLOUDS()
c._var  = "airpres,airtemp2m,rh2m,dewtemp2m,wetbulbtemp2m,blackglobetemp2m,windspeed2m,rad2m_total,wbgt2m,wbgt2m_bernard"
c.setAttributes('location', 'var', 'value', 'unit', 'obtime', 'timezone') 
c.setLoc( network = 'ECONET' )
c.end   = datetime.utcnow()
c.start = c.end - timedelta(days=5)

dfile  = '/Users/kwodzicki/Data/clouds_api_wbgt_data.pic'
mdfile = dfile.split('_')[:-1]
mdfile.append( 'metadata.pic' )
mdfile = '_'.join( mdfile )

def loadData():

  return read_pickle( dfile ), read_pickle( mdfile )

def solar_parameters_test( n = 1000, sd = datetime(1950, 1, 1, 0)):

  ed      = sd + timedelta( minutes= n )
  dates   = date_range( start = sd, end = ed, freq = '1min', inclusive = 'left')
  lon     = numpy.full( n, -80.0, dtype = numpy.float32 )
  lat     = numpy.full( n,  30.0, dtype = numpy.float32 )

  solar   = numpy.full( n, 0, dtype = numpy.float32 )

  return solar_parameters( 
    dates.year.values.astype(numpy.int32), 
    dates.month.values.astype(numpy.int32), 
    dates.day.values.astype(numpy.int32), 
    dates.hour.values.astype(numpy.int32), 
    dates.minute.values.astype(numpy.int32), 
    dates.second.values.astype(numpy.int32), lat, lon, solar )
  
#def addUnits( df, vName, newUnit, dt = None, dtype = numpy.float32 ):
def addUnits( df, vName, dt = None ):
  """To wrap values in DataFrame into pint.Quantity"""

  #unit = df['unit'].unique()
  #if unit.size != 1:
  #  raise Exception( "More then one unit for array!" )
  #unit = unit[0]

  dd = df[ df['var'] == vName ]
  if dt is not None:
    dd = dd.reindex( dt )

  unit = dd['unit'].values.astype( str )
  unit = unit[ unit != 'nan' ][0]
  if unit == 'mb': 
    unit = 'mbar'
  elif unit == 'F':
    unit = 'degree_Fahrenheit'
  
  return units.Quantity( dd['value'].values, unit )

def main( *args ):
  if len( args ) != 2:
    data, meta = c.download()
  else:
    data, meta = args

  raw_labels  = ['solar', 'pressure', 'T', 'RH', 'u']
  wbgt_labels = ['est. speed', 'globe', 'natural', 'psychrometric', 'WBGT']
  fig, axes   = plt.subplots(1, 3)
  xyrange     = [float('inf'), -float('inf')]

  for loc, df in data.groupby('location'): 

    #df     = df.sort_values( 'obtime' ).set_index('obtime')
    df       = df.set_index('obtime')                                           # Set index for dataframe         
    df.index = df.index.tz_convert(None)                                        # Convert index times to UTC
    dt       = to_datetime( df.index.unique() )                                 # Get only the unique times

    solar    = addUnits( df, 'rad2m_total', dt )    # Add units to values, convert to appropriate units, and reindex
    pres     = addUnits( df, 'airpres',     dt )
    Tair     = addUnits( df, 'airtemp2m',   dt )
    Tdew     = addUnits( df, 'dewtemp2m',   dt )
    relhum   = addUnits( df, 'rh2m',        dt )
    speed    = addUnits( df, 'windspeed2m', dt )

    Tglobe   = addUnits( df, 'blackglobetemp2m', dt ).to( 'degree_Celsius' ).magnitude
    Twbgt    = addUnits( df, 'wbgt2m',           dt ).to( 'degree_Celsius' ).magnitude
    TwbgtB   = addUnits( df, 'wbgt2m_bernard',   dt ).to( 'degree_Celsius' ).magnitude

    metaLoc  = meta[ meta['Location ID'] == loc ] 
    lon      = metaLoc['Longitude [degrees East]' ].values.astype( numpy.float32 )
    lat      = metaLoc['Latitude [degrees North]' ].values.astype( numpy.float32 )

    print( f"Location : {loc}; lat : {lat}; lon : {lon}" )    
    nn = len( solar )

    #line  =  ax[0].plot( dt, solar, 'r')
    #ax[0].yaxis.label.set_color('r')
    #for i, (var, color) in enumerate( zip( [pres, Tair, relhum, speed], ['g', 'b', 'm', 'k']) ):
    #  subax = ax[0].twinx()
    #  line  = subax.plot( dt, var, color )
    #  subax.yaxis.label.set_color( color )
    #  if i > 0: subax.spines.right.set_position( ("axes", 1.2 * i) )
    #  
    ##lines = ax[0].plot( dt, solar/10, 'r', dt, pres/10, 'g', dt, Tair, 'b', dt, relhum, 'm', dt, speed, 'k' )
    ##for line, label in zip(lines, raw_labels): line.set_label( label )
    ##ax[0].legend() 
    #plt.show()
    for mm, method in enumerate( ['liljegren' , 'dimiceli', 'bernard'] ):
      res = wbgt(method,
        lat, lon, 
        dt.year.values,
        dt.month.values,
        dt.day.values,
        dt.hour.values,
        dt.minute.values,
        solar, pres, Tair, Tdew, speed, 
      )

      for ii, (xx, yy) in enumerate( zip( [Tglobe, Twbgt, TwbgtB], [res['Tg'], res['Twbg'], res['Twbg']] ) ):
        xyrange = [min( [xx.min(), yy.min(), xyrange[0]] ), 
                   max( [xx.max(), yy.max(), xyrange[1]] ) ]
        print( xyrange )
        axes[ii].scatter( xx, yy, label=loc, alpha=0.6 )

  xyrange = numpy.asarray( xyrange ) + [-2, 2]
  axes[-1].legend()
  for ax in axes:
    ax.set_aspect( 'equal' ) 
    ax.set_xlim( xyrange )
    ax.set_ylim( xyrange )
    ax.plot( xyrange, xyrange, '-k' ) 
  
  plt.show()

if __name__ == "__main__":
  
  main()
