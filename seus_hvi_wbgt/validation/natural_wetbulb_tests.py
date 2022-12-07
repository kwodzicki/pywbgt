import sys, os
from datetime import datetime
import pytz

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import numpy
from pandas import date_range, to_datetime
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

date0 = to_datetime( [date0.astimezone( pytz.timezone('UTC'))] )

dates = date_range( '2020-01-01T00:00:00', 
                    '2021-01-01T00:00:00', 
                    freq      = '1H',
                    inclusive = 'left'
)

class Values( object ):

  def __init__(self):
    super().__init__()

    self.solar = {
        'scalar' : 1000,
        'array'  : numpy.arange(0, 1000, 10),
        'units'  : 'watt/meter**2'
    }
    self.pres = {
        'scalar' : 1013,
        'array'  : numpy.arange( 950, 1050, 10),
        'units'  : 'hPa'
    }
    self.Tair =  {
        'scalar' : 30.0,
        'array'  : numpy.arange( 20, 35, 1 ),
        'units'  : 'degree_Celsius'
    }
    self.Tdew =  {
        'scalar' : 20.0,
        'array'  : numpy.arange( 20, 35, 1 ),
        'units'  : 'degree_Celsius'
    }
    self.speed =  {
        'scalar' : 1.0,
        'array'  : numpy.concatenate( [ [0], numpy.logspace(-2, numpy.log10( 25 ) ) ] ),
        'units'  : 'mile/hour'
    }

  def __setitem__(self, key, item):
    self.__dict__[key] = item

  def __getitem__(self, key):
    return self.__dict__[key]


  def keys(self):
    return self.__dict__.keys()

  def units( self, var ):

    return units( self[var]['units'] )

  def scalar(self, var, repeat = None):

    tmp = self[var]
    tmp = units.Quantity( tmp['scalar'], tmp['units'] )
    if repeat is None:
      return tmp
    return tmp.repeat( repeat )

  def array(self, var):

    tmp = self[var]
    return units.Quantity( tmp['array'], tmp['units'] )

  def wbgt( self, method, refVar, **kwargs ):

    arr = self.array( refVar )
    args = [
      arr if k == refVar else self.scalar(k, arr.size) for k in self.keys() 
    ]

    return wbgt(method, lat, lon, date0.repeat(arr.size),
          *args,
          **kwargs
    )
 

methods    = ['liljegren', 'dimiceli', 'bernard']
colors     = ['black',     'red',      'green']
linestyle  = ['-', '--']
marker     = '.'
markersize = 2

methods    = [m.title() for m in methods]
lineFmt    = {
    method : {
        'color'     : colors[i],
        'linestyle' : linestyle
    }  for i, method in enumerate( methods )
}

keyTitle = ( ('Tpsy', 'Psychrometric Wet Bulb'),
             ('Tnwb', 'Natural Wet Bulb'),
             ('Tg',   'Globe'),
             ('Twbg', 'Wet Bulb Globe')
)



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

def idealized( *args, 
        fKey     = 'Tg', 
        label    = '' ,
        position = [0.10, 0.13, 0.99, 0.95],
        hspace   = 0.55,
        wspace   = 0.55,
        **kwargs):


  vals = Values()
 
  variables = (
    ('Tair',  'Air Temperature',       '$^{\circ}$C'),
    ('Tdew',  'Dew Point Temperature', '$^{\circ}$C'),
    ('pres',  'Pressure',              'hPa'),
    ('solar', 'Solar Irradiance',      'W/m$^{2}$'),
    ('speed', 'Wind Speed',            'mph'),
  )

  #x_label = [ f"{l} [{u}]" for l, u in variables ]

  path = args[1] if len(args) > 1 else None


  solar = vals.array('solar') 
  ss, cosz, f_db = solar_parameters(
    lat, lon, date0.repeat(solar.size),
      solar.to( 'watt/m**2' ).m
  ) 
  #############################################################################  
  fig, ax = plt.subplots( 2, 3, figsize = (6,3.5))
  ax      = ax.flatten()
  ax[-1].axis('off')

  natural_wetbulb = kwargs.pop('natural_wetbulb', (None) )

  
  #############################################################################  
  # influence of air temperature 
  for pp, (key, name, unit) in enumerate( variables ):
    ax[pp].set_xlabel( f"{name} [{unit}]" )
    for i, method in enumerate(methods):
      nn = natural_wetbulb if method.upper() == 'DIMICELI' else (None,) 

      xx = vals.array( key )
      for ii, n in enumerate(nn):
        tmp = vals.wbgt( method, key, natural_wetbulb = n )
        ax[pp].plot( xx, tmp[fKey], 
            label     = n.title() if n is not None else method.title(), 
            color     = colors[i],
            linestyle = linestyle[ii],
            **pltKWArgs
        )



  Tair = vals.array(  'Tair' )
  Tdew = vals.scalar( 'Tdew' )
  rh   = relative_humidity( Tair.to('degC').m, Tdew.to('degC').m )
  uplot( ax[0], Tair, rh, ylabel = 'RH', ylim=[0, 1] )

  Tdew = vals.array( 'Tdew' )
  Tair = Tdew.max().repeat( Tdew.size )
  rh   = relative_humidity( Tair.to('degC').m, Tdew.to('degC').m )
  uplot( ax[1], Tdew, rh, ylabel = 'RH', ylim=[0, 1] )

  solar = vals.array( 'solar' )
  uplot( ax[3], solar, f_db, ylabel = 'F$_{direct}$', ylim=[0, 1] )

  for l, a in zip( ['a', 'b', 'c', 'd', 'e'], ax ):
    a.text(0.01, 1.0, f"({l})", transform=a.transAxes,# + trans,
            fontsize            = 'small', 
            verticalalignment   = 'bottom',
            horizontalalignment = 'right' )
            #bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

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

if __name__ == "__main__":

  infos = [
    #{'fKey'   : 'Twbg', 
    # 'obsKey' : '', 
    # 'label'  : 'Wet Bulb Globe Temperature',
    # 'natural_wetbulb' : ['HUNTER_MINYARD', 'MALCHAIRE']
    #},

    {'fKey'   : 'Tnwb', 
     'obsKey' : '', 
     'label'  : 'Natural Wet Bulb Temperature',
     'natural_wetbulb' : ['HUNTER_MINYARD', 'MALCHAIRE']
    },

  ]

  for info in infos:
    idealized( *sys.argv, **info )
