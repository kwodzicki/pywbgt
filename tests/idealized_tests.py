import sys, os
from datetime import datetime
import pytz

from matplotlib import pyplot as plt

import numpy
from pandas import date_range
from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt
from seus_hvi_wbgt.wbgt.utils import relative_humidity

figKWArgs = {'dpi' : 300}


# Tulsa, OK
lat   = numpy.full( 1,  36.1540 )
lon   = numpy.full( 1, -95.9928 )
date0 = datetime(2020, 6, 21, 12, 0, 0, tzinfo=pytz.timezone('US/Central'))

# Charlotte NC
lat   = numpy.full( 1,  35.2271 )
lon   = numpy.full( 1, -80.8431 )
date0 = datetime(2020, 6, 21, 12, 0, 0, tzinfo=pytz.timezone('US/Eastern'))

#lat   = numpy.full( 1,  0.0 )
#lon   = numpy.full( 1,  0.0)

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

speed0 = units.Quantity( [1.0], 'mile/hour')
speed  = numpy.logspace(-2, numpy.log10( 25 ) ) * units('mile/hour')
speed  = numpy.concatenate( (units.Quantity([0], 'mile/hour'), speed) )

methods  = ['liljegren', 'dimiceli', 'bernard']
#methods  = ['dimiceli']
lineFmt  = ['-r',        '-k',       '-g']
lineFmt  = {m:l for m, l in zip( methods, lineFmt )}

keyTitle = ( ('Tpsy', 'Psychrometric Wet Bulb'),
             ('Tnwb', 'Natural Wet Bulb'),
             ('Tg',   'Globe'),
             ('Twbg', 'Wet Bulb Globe')
)

def legend_without_duplicate_labels( fig ):
  handles, labels = plt.gca().get_legend_handles_labels()
  by_label = dict(zip(labels, handles))
  fig.legend(by_label.values(), by_label.keys(), 
    ncol = 3,
    loc='lower left', 
    fontsize='small'
  )

def plotter( xx, data, data2=None, data2Label='data2', figAxes = None, nn=None, suptitle=None, **kwargs):

  position = [0.125, 0.15, 0.95, 0.93]
  wspace   = None
  hspace   = 0.35

  if figAxes is None:
    fig, axes = plt.subplots( 2, 2 )
    axes      = axes.flatten()
  else:
    fig, axes = figAxes

  kwargs.update( {'linewidth' : 3.0 } )
  if 'alpha' not in kwargs:
    kwargs['alpha'] = 0.9

  unit = None
  for i, (key, title) in enumerate( keyTitle ):
    if key in data:
      ee = len(data[key]) if nn is None else nn
      axes[i].plot( xx[:ee], data[key][:ee], 
        lineFmt.get( kwargs.get('label', ''), '' ),
        **kwargs ) 
      axes[i].set_title( title ) 
      axes[i].set_xlabel('')
      if unit is None: unit = xx.units
    if data2 is not None:
      ax2 = axes[i].twinx()
      ax2.plot( xx[:ee], data2[:ee], label=data2Label )

  fig.text( 0.01, 0.5, 'Temperature [degree_Celsius]', 
    fontsize            = 'large',
    rotation            = 'vertical',
    horizontalalignment = 'left',
    verticalalignment   = 'center'
  ) 
  #fig.subplots_adjust( left=0.15, right=0.95, bottom=0.15, top=0.9, hspace = 0.5 )
  fig.subplots_adjust( *position, wspace, hspace ) 
  if isinstance(suptitle, str):
    #fig.suptitle( suptitle )
    if unit is not None:
      suptitle = f"{suptitle} [{unit}]"
    fig.text( sum(position[::2])/2.0, 0.07, suptitle, 
      horizontalalignment = 'center',
      fontsize            = 'large'
    )

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

def main( *args, fbase = 'Idealized WBGT' ):
  path = args[1] if len(args) > 1 else None

  for method in methods:
    f_db  = numpy.linspace(0, 1, solar.size)
    tmp   = wbgt(method, lat, lon, 
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
      cosz = numpy.full( solar.size, 0.8 )
     )
    figAxes = plotter( solar, tmp, suptitle = 'Direct Beam Solar Radiation')
  
  figAxes  = None
  supTitle = 'Direct Bean Solar Radiation'
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
    figAxes = plotter( solar, tmp, suptitle=supTitle, figAxes=figAxes, label=method)
    #if method == 'bernard':
    #  print( solar[ (tmp['Tg']-Tair0.repeat( solar.size ).to('degree_Celsius').magnitude) > 4 ] )
  legend_without_duplicate_labels( figAxes[0] )
  if path is not None:
    figAxes[0].savefig( os.path.join( path, f"{fbase} {supTitle}.png" ), **figKWArgs)
  
  figAxes  = None
  supTitle = 'Pressure'
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
    
    figAxes = plotter( pres, tmp, suptitle=supTitle, figAxes=figAxes, label=method )
  legend_without_duplicate_labels( figAxes[0] )
  if path is not None:
    figAxes[0].savefig( os.path.join( path, f"{fbase} {supTitle}.png" ) , **figKWArgs)
  
  figAxes  = None
  supTitle = 'Air Temperature'
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
    
    figAxes = plotter( Tair, tmp, suptitle=supTitle, figAxes=figAxes, label=method )
  legend_without_duplicate_labels( figAxes[0] )
  if path is not None:
    figAxes[0].savefig( os.path.join( path, f"{fbase} {supTitle}.png" ) , **figKWArgs)
  
  figAxes  = None
  supTitle = 'Dew Point Temperature'
  Tair_tmp = numpy.full( Tdew.size, Tdew.max() ) * Tair0.units
  rh       = relative_humidity( Tair_tmp.to('degree_Celsius').m, Tdew.to('degree_Celsius').m )
  rh       = None
  
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
      Tair_tmp,
      Tdew,
      speed0.repeat( Tdew.size )
    )
    #print( method, ' : ', tmp['Tg'] ) 
    figAxes = plotter( Tdew, tmp, rh, data2Label='Relative Humidity', suptitle=supTitle, figAxes=figAxes, label=method )
  legend_without_duplicate_labels( figAxes[0] )
  if path is not None:
    figAxes[0].savefig( os.path.join( path, f"{fbase} {supTitle}.png" ) , **figKWArgs)
    
  figAxes  = None
  supTitle = 'Wind Speed'
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
    
    figAxes = plotter( speed, tmp, suptitle=supTitle, figAxes=figAxes, label=method )
  legend_without_duplicate_labels( figAxes[0] )
  if path is not None:
    figAxes[0].savefig( os.path.join( path, f"{fbase} {supTitle}.png" ) , **figKWArgs)
  else:  
    plt.show()

if __name__ == "__main__":
  main( *sys.argv )
