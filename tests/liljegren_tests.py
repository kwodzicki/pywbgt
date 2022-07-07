from datetime import datetime, timedelta

import numpy
from pandas import date_range

from matplotlib import pyplot as plt

from metpy.units import units

from seus_hvi_wbgt.wbgt.liljegren import liljegren, solar_parameters

from ncsuCLOUDS.api import CLOUDS

c       = CLOUDS()
c._var  = "airpres,airtemp2m,rh2m,dewtemp2m,blackglobetemp2m,windspeed2m,rad2m_total"
c.setAttributes('location', 'var', 'value', 'unit', 'obtime', 'timezone') 
c.setLoc( network = 'ECONET' )
c.end   = datetime.utcnow()
c.start = c.end - timedelta(days=1)

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
  

def main( *args ):
  if len( args ) != 2:
    data, meta = c.download()
  else:
    data, meta = args

  raw_labels  = ['solar', 'pressure', 'T', 'RH', 'u']
  wbgt_labels = ['solar', 'globe', 'natural', 'psychrometric', 'WBGT']

  for loc, df in data.groupby('location'): 

    df      = df.sort_values( 'obtime' )
    solar   = df[ df['var'] == 'rad2m_total' ]['value'].values.astype( numpy.float32 )
    pres    = df[ df['var'] == 'airpres'     ]['value'].values.astype( numpy.float32 )
    Tair    = df[ df['var'] == 'airtemp2m'   ]['value'].values.astype( numpy.float32 )
    relhum  = df[ df['var'] == 'rh2m'        ]['value'].values.astype( numpy.float32 )
    speed   = df[ df['var'] == 'windspeed2m' ]['value'].values.astype( numpy.float32 )

    metaLoc = meta[ meta['Location ID'] == loc ] 
    lon     = metaLoc['Longitude [degrees East]' ].values.astype( numpy.float32 )
    lat     = metaLoc['Latitude [degrees North]' ].values.astype( numpy.float32 )
    
    nn = len( solar )
    dt = df['obtime'].dt.tz_convert(None)

    xx = dt.unique()
    lines = plt.plot( xx, solar/1000, 'r', xx, pres/1000, 'g', xx, Tair, 'b', xx, relhum, 'm', xx, speed, 'k' )
    for line, label in zip(lines, raw_labels): line.set_label( label )
    plt.legend() 
    plt.show()

    sp, Tg, Tnwb, Tpsy, Twbg = liljegren(
        lat.repeat( nn ), lon.repeat( nn ), 
        numpy.ones( nn, dtype=numpy.int32 ), 
        dt.dt.year.values.astype(numpy.int32), 
        dt.dt.month.values.astype(numpy.int32), 
        dt.dt.day.values.astype(numpy.int32), 
        dt.dt.hour.values.astype(numpy.int32), 
        dt.dt.minute.values.astype(numpy.int32), 
        numpy.zeros( nn, dtype=numpy.int32), 
        numpy.ones( nn, dtype=numpy.int32 ),
        solar, pres, Tair, relhum, speed, 
        numpy.full(nn, 2, dtype=numpy.float32), 
        numpy.full(nn, -1, dtype=numpy.float32) ) 

    lines = plt.plot( dt, sp, 'r', dt, Tg, 'g', dt, Tnwb, 'b', dt, Tpsy, 'm', dt, Twbg, 'k' )
    for line, label in zip(lines, wbgt_labels): line.set_label( label )
    plt.legend() 
    plt.show()
    res = input('Type Y to quit')
    if res == 'Y': return 
