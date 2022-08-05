from datetime import datetime
import pytz

from matplotlib import pyplot as plt

import numpy
from pandas import date_range
from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt

lat   = numpy.full( 1,  36.1540 )
lon   = numpy.full( 1, -95.9928 )
#lat   = numpy.full( 1,  0.0 )
#lon   = numpy.full( 1,  0.0)
print( lat )
print( lon ) 
date0 = datetime(2020, 6, 21, 12, 0, 0, tzinfo=pytz.timezone('US/Central'))
date0 = date0.astimezone( pytz.timezone('UTC'))

dates = date_range( '2020-01-01T00:00:00', 
                    '2021-01-01T00:00:00', 
                    freq      = '1H',
                    inclusive = 'left'
)

solar0 = units.Quantity( [1000], 'watt/meter**2')
solar  = numpy.arange(0, 1000, 10) * units('watt/meter**2')

pres0  = units.Quantity( [1013], 'hPa')
pres   = numpy.arange( 950, 1050, 10 ) * units('hPa')

Tair0  = units.Quantity( [20.0], 'degree_Celsius')
Tair   = numpy.arange( 5, 45, 1 ) * units('degree_Celsius')

Tdew0  = units.Quantity( [5.0], 'degree_Celsius')
Tdew   = numpy.arange( 5, 45, 1 ) * units('degree_Celsius')

speed0 = units.Quantity( [5.0], 'mile/hour')
speed  = numpy.logspace(-2, numpy.log10( 50 ) ) * units('mile/hour')
speed  = numpy.concatenate( (units.Quantity([0], 'mile/hour'), speed) )

methods  = ['liljegren', 'dimiceli', 'bernard']
methods  = ['dimiceli']
lineFmt  = ['-r',        '-k',       '-g']
lineFmt  = {m:l for m, l in zip( methods, lineFmt )}

keyTitle = ( ('Tg',   'Globe Temperature'),
             ('Tpsy', 'Psychrometric Wet Bulb Temperature'),
             ('Tnwb', 'Natural Wet Bulb Temperature'),
             ('Twbg', 'Wet Bulb Globe Temperature')
)

def legend_without_duplicate_labels( fig ):
  handles, labels = plt.gca().get_legend_handles_labels()
  by_label = dict(zip(labels, handles))
  fig.legend(by_label.values(), by_label.keys(), loc='lower right')

def plotter( xx, data, figAxes = None, nn=-1, suptitle=None, **kwargs):

  if figAxes is None:
    fig, axes = plt.subplots( 4, 1 )
  else:
    fig, axes = figAxes

  kwargs.update( {'linewidth' : 3.0 } )
  if 'alpha' not in kwargs:
    kwargs['alpha'] = 0.6

  for i, (key, title) in enumerate( keyTitle ):
    if key in data:
      axes[i].plot( xx[:nn], data[key][:nn], 
        lineFmt.get( kwargs.get('label', ''), None ),
        **kwargs ) 
      axes[i].set_title( title ) 
  
  fig.subplots_adjust( left=0.1, right=0.95, bottom=0.08, top=0.9, hspace = 1.0 )
  if isinstance(suptitle, str):
    fig.suptitle( suptitle )

  return fig, axes


#figAxes = None 
#for method in methods: 
#  # Influence of datetime
#  tmp = wbgt(method, lat, lon, 
#    dates.year.values, 
#    dates.month.values, 
#    dates.day.values, 
#    dates.hour.values, 
#    dates.minute.values,
#    solar0.repeat( dates.size ), 
#    pres0.repeat(  dates.size ),
#    Tair0.repeat(  dates.size ),
#    Tdew0.repeat(  dates.size ),
#    speed0.repeat( dates.size )
#   )
#  
#  figAxes = plotter( dates, tmp, suptitle='Datetime', figAxes=figAxes, label=method )
#legend_without_duplicate_labels( figAxes[0] )


f_db  = numpy.linspace(0, 1, solar.size)
f_dif = 1.0 - f_db
tmp   = wbgt('dimiceli', lat, lon, 
  numpy.full( solar.size, date0.year ),
  numpy.full( solar.size, date0.month),
  numpy.full( solar.size, date0.day  ),
  numpy.full( solar.size, date0.hour ),
  numpy.full( solar.size, date0.minute ),
  solar0.repeat( solar.size ), 
  pres0.repeat(  solar.size ),
  Tair0.repeat(  solar.size ),
  Tdew0.repeat(  solar.size ),
  speed0.repeat( solar.size ),
  f_db = f_db,
  cosz = numpy.full( solar.size, 0.01 )
 )
figAxes = plotter( solar, tmp, suptitle = 'Solar')


figAxes = None
for method in methods: 
  # influence of solar irradiance
  tmp = wbgt(method, lat, lon, 
    numpy.full( solar.size, date0.year ),
    numpy.full( solar.size, date0.month),
    numpy.full( solar.size, date0.day  ),
    numpy.full( solar.size, date0.hour ),
    numpy.full( solar.size, date0.minute ),
    solar, 
    pres0.repeat(  solar.size ),
    Tair0.repeat(  solar.size ),
    Tdew0.repeat(  solar.size ),
    speed0.repeat( solar.size )
   )
  figAxes = plotter( solar, tmp, suptitle = 'Solar', figAxes=figAxes, label=method)
legend_without_duplicate_labels( figAxes[0] )
#exit()  
  
figAxes = None
for method in methods: 
  # influence of atmospheric pressure
  tmp = wbgt(method, lat, lon, 
    numpy.full( pres.size, date0.year ),
    numpy.full( pres.size, date0.month),
    numpy.full( pres.size, date0.day  ),
    numpy.full( pres.size, date0.hour ),
    numpy.full( pres.size, date0.minute ),
    solar0.repeat( pres.size), 
    pres,
    Tair0.repeat(  pres.size ),
    Tdew0.repeat(  pres.size ),
    speed0.repeat( pres.size )
   )
  
  figAxes = plotter( pres, tmp, suptitle = 'Pressure', figAxes=figAxes, label=method )
legend_without_duplicate_labels( figAxes[0] )

figAxes = None
for method in methods: 
  # influence of air temperature 
  tmp = wbgt(method, lat, lon, 
    numpy.full( Tair.size, date0.year ),
    numpy.full( Tair.size, date0.month),
    numpy.full( Tair.size, date0.day  ),
    numpy.full( Tair.size, date0.hour ),
    numpy.full( Tair.size, date0.minute ),
    solar0.repeat( Tair.size), 
    pres0.repeat(  Tair.size ),
    Tair,
    Tdew0.repeat(  Tair.size ),
    speed0.repeat( Tair.size )
   )
  
  figAxes = plotter( Tair, tmp, suptitle = 'Air Temperature' , figAxes=figAxes, label=method )
legend_without_duplicate_labels( figAxes[0] )

figAxes = None
for method in methods: 
  # influence of dew point temperature 
  tmp = wbgt(method, lat, lon, 
    numpy.full( Tdew.size, date0.year ),
    numpy.full( Tdew.size, date0.month),
    numpy.full( Tdew.size, date0.day  ),
    numpy.full( Tdew.size, date0.hour ),
    numpy.full( Tdew.size, date0.minute ),
    solar0.repeat( Tdew.size), 
    pres0.repeat(  Tdew.size ),
    numpy.full( Tdew.size, Tdew.max() ) * Tair0.units,
    Tdew,
    speed0.repeat( Tdew.size )
   )
  
  figAxes = plotter( Tdew, tmp, suptitle = 'Dew Point Temperature' , figAxes=figAxes, label=method )
legend_without_duplicate_labels( figAxes[0] )
  
figAxes = None
for method in methods: 
  # influence of wind speed 
  tmp = wbgt(method, lat, lon, 
    numpy.full( speed.size, date0.year ),
    numpy.full( speed.size, date0.month),
    numpy.full( speed.size, date0.day  ),
    numpy.full( speed.size, date0.hour ),
    numpy.full( speed.size, date0.minute ),
    solar0.repeat( speed.size), 
    pres0.repeat(  speed.size ),
    Tair0.repeat(  speed.size ),
    Tdew0.repeat(  speed.size ),
    speed
   )
  
  figAxes = plotter( speed, tmp, suptitle='Wind speed', figAxes=figAxes, label=method )
legend_without_duplicate_labels( figAxes[0] )

plt.show()

