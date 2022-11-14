import sys, os
from datetime import datetime
import pytz

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import numpy
from pandas import date_range
from metpy.units import units

from ncsuCLOUDS import io
from ncsuCLOUDS.units import addUnits

from seus_hvi_wbgt.wbgt import wbgt
from seus_hvi_wbgt.wbgt.utils import relative_humidity
from seus_hvi_wbgt.wbgt.liljegren import solar_parameters

EXT       = 'pdf'
EXT       = 'png'
figKWArgs = {'dpi' : 300}
pltKWArgs = {
  'linewidth' : 3.0, 
  'alpha'     : 0.9
}
legend_kwargs = {
  'loc'      : 'center', 
}

 
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

Tair0  = units.Quantity( [30.0], 'degree_Celsius')
Tair   = numpy.arange( 20, 35, 1 ) * units('degree_Celsius')

Tdew0  = units.Quantity( [20.0], 'degree_Celsius')
Tdew   = numpy.arange( 20, 35, 1 ) * units('degree_Celsius')

speed0 = units.Quantity( [1.0], 'mile/hour')
speed  = numpy.logspace(-2, numpy.log10( 25 ) ) * units('mile/hour')
speed  = numpy.concatenate( (units.Quantity([0], 'mile/hour'), speed) )

METHODS    = ['liljegren', 'dimiceli', 'bernard']
colors     = ['black',     'red',      'green']
linestyle  = '-'
marker     = '.'
markersize = 2

METHODS    = [m.title() for m in METHODS]
lineFmt    = {
    method : {
        'color'     : colors[i],
        'linestyle' : linestyle
    }  for i, method in enumerate( METHODS )
}

keyTitle = ( ('Tpsy', 'Psychrometric Wet Bulb'),
             ('Tnwb', 'Natural Wet Bulb'),
             ('Tg',   'Globe'),
             ('Twbg', 'Wet Bulb Globe')
)


def mae(y_true, predictions):
  """Compute mean absolute error"""

  return numpy.mean( numpy.abs( y_true - predictions ) )

def uplot( ax, xx, data, 
    nn            = None, 
    label         = None, 
    color         = 'gray', 
    ylim          = None,
    ylabel        = None, 
    labelsize     = 'small',
    labelrotation = 90,
    **kwargs ):

  ee = len(data) if nn is None else nn
  ax2 = ax.twinx()
  ax2.tick_params('y',
      labelrotation = labelrotation,
      colors        = color,
      labelsize     = labelsize 
  )
  ax2.plot( xx[:ee], data[:ee], 
    color     = color, 
    linestyle = '-', 
    label     = label,
    **kwargs
  )
  if ylim   is not None: ax2.set_ylim( ylim )
  if ylabel is not None:
    ax2.set_ylabel( ylabel )
    ax2.yaxis.label.set_color( color )

  ax.set_zorder( ax2.get_zorder()+1 )
  ax.patch.set_visible( False )

def plotter( ax, xx, data, data2=None, data2Label=None, figAxes = None, nn=None, suptitle=None, label='', **kwargs):

  position = [0.125, 0.15, 0.95, 0.93]
  wspace   = 0.25
  hspace   = 0.35

  kwargs.update( {'linewidth' : 3.0 } )
  if 'alpha' not in kwargs:
    kwargs['alpha'] = 0.9

  unit = None
  for i, (key, title) in enumerate( keyTitle ):
    if key in data:
      ee  = len(data[key]) if nn is None else nn
      fmt = lineFmt.get( label, '' )
      ax.plot( xx[:ee], data[key][:ee], fmt, 
        label = label if i == 0 else None,
        **kwargs ) 
      ax.set_title( title ) 
      ax.set_xlabel('')
      ax.xaxis.set_minor_locator(AutoMinorLocator())
      ax.yaxis.set_minor_locator(AutoMinorLocator())

      if unit is None: unit = xx.units
  return
  fig.text( 0.01, 0.5, 'Temperature [degree_Celsius]', 
    fontsize            = 'large',
    rotation            = 'vertical',
    horizontalalignment = 'left',
    verticalalignment   = 'center'
  ) 
  #fig.subplots_adjust( left=0.15, right=0.95, bottom=0.15, top=0.9, hspace = 0.5 )
  return 
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

def idealized( *args, 
        methods  = METHODS,
        fKey     = 'Tg', 
        label    = '' ,
        position = [0.10, 0.13, 0.99, 0.95],
        hspace   = 0.55,
        wspace   = 0.55,
        **kwargs):

  
  variables = (
    ('Air Temperature',       '$^{\circ}$C'),
    ('Dew Point Temperature', '$^{\circ}$C'),
    ('Pressure',              'hPa'),
    ('Solar Irradiance',      'W/m$^{2}$'),
    ('Wind Speed',            'mph'),
  )

  x_label = [ f"{l} [{u}]" for l, u in variables ]

  path = args[1] if len(args) > 1 else None

 
  ss, cosz, f_db = solar_parameters(
    lat, lon,
      numpy.full( solar.size, date0.year ),
      numpy.full( solar.size, date0.month),
      numpy.full( solar.size, date0.day  ),
      numpy.full( solar.size, date0.hour ),
      numpy.full( solar.size, date0.minute ),
      solar.to('watt/m**2').magnitude, 
  ) 
  #############################################################################  
  fig, ax = plt.subplots( 2, 3, figsize = (6,3.5))
  ax      = ax.flatten()
  ax[-1].axis('off')

  pp = 0
  #############################################################################  
  # influence of air temperature 
  rh       = relative_humidity( Tair.to('degree_Celsius').m, Tdew0.to('degree_Celsius').m )
  for method in methods: 
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
      speed0.repeat( Tair.size ),
      **kwargs
     )
    ax[pp].plot( Tair.to('degree_Celsius'), tmp[fKey], 
        label = method.title(), 
        **lineFmt[method], 
        **pltKWArgs
    )
  uplot( ax[pp], Tair, rh, ylabel = 'RH', ylim=[0, 1] )
 
  pp += 1

  #############################################################################  
  # influence of dew point temperature 
  Tair_tmp = numpy.full( Tdew.size, Tdew.max() ) * Tair0.units
  rh       = relative_humidity( Tair_tmp.to('degree_Celsius').m, Tdew.to('degree_Celsius').m )
  for j, method in enumerate(methods): 
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
      speed0.repeat( Tdew.size ),
      **kwargs
    )
    ax[pp].plot( Tdew.to('degree_Celsius'), tmp[fKey], 
        label = method.title(), 
        **lineFmt[method], 
        **pltKWArgs
    )
  uplot( ax[pp], Tdew, rh, ylabel = 'RH', ylim=[0, 1] )

  pp += 1

  #############################################################################  
  # influence of atmospheric pressure
  for method in methods: 
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
      speed0.repeat( pres.size ),
      **kwargs
     )
    
    ax[pp].plot( pres.to('hPa'), tmp[fKey], 
        label = method.title(), 
        **lineFmt[method], 
        **pltKWArgs
    )

  pp += 1

  #############################################################################  
  # influence of solar irradiance
  for method in methods: 
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
      speed0.repeat( solar.size ),
      **kwargs
    )
    ax[pp].plot( solar.to('watt/m**2'), tmp[fKey], 
        label = method.title(), 
        **lineFmt[method], 
        **pltKWArgs
    )
    #if method == 'bernard':
    #  print( solar[ (tmp['Tg']-Tair0.repeat( solar.size ).to('degree_Celsius').magnitude) > 4 ] )
  uplot( ax[pp], solar, f_db, ylabel = 'F$_{direct}$', ylim=[0, 1] )

  pp += 1

  #############################################################################  
  # influence of wind speed 
  for method in methods:
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
      speed,
      **kwargs
     )
    ax[pp].plot( speed.to('mile/hour'), tmp[fKey], 
        label = method.title(), 
        **lineFmt[method], 
        **pltKWArgs
    )

  pp += 1

  for l, a in zip( ['a', 'b', 'c', 'd', 'e'], ax ):
    a.text(0.01, 1.0, f"({l})", transform=a.transAxes,# + trans,
            fontsize            = 'small', 
            verticalalignment   = 'bottom',
            horizontalalignment = 'right' )
            #bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

  for i, l in enumerate( x_label ):
    ax[i].set_xlabel( l )

  fig.text( 0.02, sum(position[1::2])/2.0, f"{label} [$^{{\circ}}$C]",
    horizontalalignment='center',
    verticalalignment='center',
    rotation=90,
    fontsize='medium'
  )


  fig.subplots_adjust( *position, wspace, hspace ) 
  ax[-1].legend( *ax[0].get_legend_handles_labels(), **legend_kwargs )

  if path is not None:
    fig.savefig( os.path.join( path, f"Idealized {fKey}.{EXT}" ) , **figKWArgs)
  else:  
    plt.show()

def observations( obs, meta, *args,
        methods  = METHODS,
        color_by = {'' :()},
        obsKey   = 'blackglobetemp2m|K', 
        fKey     = 'Tg', 
        label    = '',
        position = [0.07, 0.08, 0.87, 0.95],
        hspace   = 0.10,
        wspace   = 0.05,
        **kwargs ):

  if obsKey == '' or obsKey not in obs: return
 
  path = args[1] if len(args) > 1 else None

  obs       = obs.loc[:, '2020-06-01':, :]
  dates     = obs.index.get_level_values( obs.index.names.index( 'datetime' ) )

  nRowCol   = (len(color_by), len(methods)+1 )
  fig, ax   = plt.subplots( 
                *nRowCol, 
                figsize = (6, 1.75*len(color_by)),
                gridspec_kw = {'width_ratios' : [*[10]*3, 1]}
  )
  ax        = ax.flatten()

  minmax = [float('inf'), -float('inf')]

  axID   = 0
  for j, (cb_label, cb_info) in enumerate( color_by.items() ):
    for i, method in enumerate( methods ):
      tmp = wbgt( method, lat, lon,
        dates.year.values, 
        dates.month.values, 
        dates.day.values, 
        dates.hour.values, 
        dates.minute.values,
        addUnits( 'rad2m_total',     obs, meta ),
        addUnits( 'airpres|kPa',     obs, meta ),
        addUnits( 'airtemp2m|K',     obs, meta ),
        addUnits( 'dewtemp2m|K',     obs, meta ),
        addUnits( 'windspeed2m|mph', obs, meta )
      )

      x = addUnits( obsKey, obs, meta ).to('degree_Celsius').magnitude
      y = tmp[fKey]

      if j == 0:
        idx  = numpy.where( numpy.isfinite( x ) & numpy.isfinite( y ) )
        r    = numpy.corrcoef( x[idx], y[idx] )[0,1]
        text = [
          f"r   = {r:0.3f}",
          f"MAE = {mae( x[idx], y[idx] ):0.3f}",
          f"n   = {idx[0].size}",
        ]
            
        ax[axID].text( 0.95, 0.05, os.linesep.join( text ),
          family              = 'monospace',
          fontsize            = 'small',
          transform           = ax[axID].transAxes,
          horizontalalignment = 'right',
        )

      minmax = [
          min( [minmax[0], *x, *y] ),
          max( [minmax[1], *x, *y] )
      ]

      color_opts = {
          'c'      : 'black',
          'marker' : marker
      }

      if len( cb_info ) == 4:
        color_opts.update(
          { 
              'c'    : addUnits( cb_info[0], obs, meta ),
              'vmin' : cb_info[1],
              'vmax' : cb_info[2],
          }
        )
        if cb_info[-1] != '':
          color_opts[ 'c'] = color_opts['c'].to( cb_info[-1] )
        color_opts[ 'c'] = color_opts['c'].magnitude
        

      ss = ax[axID].scatter( x, y, markersize, **color_opts )

      if axID < len( methods ):
        ax[axID].set_title( method )

      if i != 0:
        ax[axID].yaxis.set_ticklabels( [] )

      axID += 1

    if len( color_by ) > 0:
      fig.colorbar( ss, 
        label    = cb_label, 
        cax      = ax[axID],
        extend   = 'both', 
        fraction = 0.05,
      )

    axID += 1

  ax = ax.reshape( nRowCol )
  for i in range( len(color_by)-1 ):
    for j in range( len(methods) ):
      ax[i,j].xaxis.set_ticklabels( [] )

  fig.text( (position[0]+position[2])/2.0, 0.01, 
    f"ECONet {label} [$^{{\circ}}$C]",
    horizontalalignment = 'center',
    verticalalignment='bottom' 
  )

  minmax = (numpy.asarray( minmax )//10 + [0,1]) * 10
  for i in range( len( color_by ) ):
    for j in range( len(methods ) ):
      ax[i,j].set( xlim=minmax, ylim=minmax )
      ax[i,j].plot( minmax, minmax, color='gray', linestyle='-', zorder=0 )
      ax[i,j].xaxis.set_ticks( numpy.arange( *minmax, 10 ) )
      ax[i,j].xaxis.set_minor_locator( AutoMinorLocator() ) 
      ax[i,j].yaxis.set_ticks( numpy.arange( *minmax, 10 ) )
      ax[i,j].yaxis.set_minor_locator( AutoMinorLocator() ) 
      ax[i,j].set_aspect('equal', 'box')

  fig.text( 0.01, 0.5, f"{label} [$^{{\circ}}$C]",
    rotation=90,
    horizontalalignment = 'left',
    verticalalignment   = 'center'
  )
  
  fig.subplots_adjust( *position, wspace, hspace )

  if path is not None:
    fig.savefig( os.path.join( path, f"Observed {fKey}.{EXT}" ) , **figKWArgs)
  else:  
    plt.show()



if __name__ == "__main__":

  color_by  = {
     'Air Temperature\n[$^{\circ}$C]'       : ('airtemp2m|K',     5,   30, 'degC'),
     'Dew Point Temperature\n[$^{\circ}$C]' : ('dewtemp2m|K',     5,   25, 'degC'),
     'Wind Speed\n[mile/hour]'              : ('windspeed2m|mph', 0,    8, 'mph'),
     'Solar Irradiance\n[W/m$^{2}$]'        : ('rad2m_total',     0, 1000, ''),
  }

  variables = [
    {'fKey'     : 'Twbg', 
     'obsKey'   : 'wbgt2m|K',
     'label'    : 'Wet Bulb Globe Temperature',
     'color_by' : color_by
    }
  ]

  data, meta = io.read( '/Users/kwodzicki/Data' )

  for info in variables:
    idealized( *sys.argv, **info)
    observations( data, meta, *sys.argv, **info )
