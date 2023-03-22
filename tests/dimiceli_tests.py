from datetime import datetime, timedelta

import numpy
from pandas import to_datetime, date_range, read_pickle

from matplotlib import pyplot as plt

from metpy.units import units

from seus_hvi_wbgt.wbgt.dimiceli import dimiceli, globeTemperature 

from ncsuCLOUDS.api import CLOUDS

c       = CLOUDS()
c._var  = "airpres,airtemp2m,rh2m,dewtemp2m,blackglobetemp2m,windspeed2m,rad2m_total,wbgt2m,wbgt2m_bernard"
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

def addUnits( df, vName, dt = None ):
  """To wrap values in DataFrame into pint.Quantity"""

  #unit = df['unit'].unique()
  #if unit.size != 1:
  #  raise Exception( "More then one unit for array!" )
  #unit = unit[0]

  dd = df[ df['var'] == vName ]
  if dt is not None:
    dd = dd.reindex( dt )

  unit = dd['unit'].values[0]
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
  wbgt_labels = ['globe', 'natural', 'WBGT']

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

    Ta = Tair.to('kelvin').magnitude
    u  = speed.to('meter/second').magnitude
    P  = pres.to('hPa').magnitude
    S  = solar.to('watt/meter**2').magnitude

    #Tg = globeTemperature( Ta[0], u[0], P[0], S[0] )

    #return Ta, u, P, S, Tg, Tglobe

    fig, ax = plt.subplots(2, 1)
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
    Tg, Tnwb, Twbg = dimiceli(
        lat, lon, dt,
        solar, pres, Tair, Tdew, speed
    )

    lines = ax[1].plot( dt, Tg, 'g', dt, Tnwb, 'b', dt, Twbg, 'k' )
    for line, label in zip(lines, wbgt_labels): line.set_label( label )
    ax[1].legend() 

    fig, ax = plt.subplots(2, 2)
    for i, (xx, yy) in enumerate( zip( [Tg, Twbg, Twbg], [Tglobe, Twbgt, TwbgtB] ) ):
      xyrange = [min( [numpy.nanmin(xx), numpy.nanmin(yy)] )-2.0, 
                 max( [numpy.nanmax(xx), numpy.nanmax(yy)] )+2.0 ]
      ax[i//2, i%2].scatter( xx, yy )
      ax[i//2, i%2].set_aspect( 'equal' ) 
      ax[i//2, i%2].set_xlim( xyrange )
      ax[i//2, i%2].set_ylim( xyrange ) 
   
    plt.show()
    res = input('Type Y to quit : ')
    if res == 'Y': return 
