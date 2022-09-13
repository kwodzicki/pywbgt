import sys, os
from datetime import datetime
import pytz

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

#plt.rcParams['text.usetex'] = True

import numpy
from pandas import date_range
from metpy.units import units

from seus_hvi_wbgt.wbgt.bernard import psychrometricWetBulb as bernard_wetbulb
from seus_hvi_wbgt.wbgt.liljegren import psychrometricWetBulb as liljegren_wetbulb#, solar_parameters
from seus_hvi_wbgt.wbgt.wetbulb import wetBulb
from seus_hvi_wbgt.wbgt.utils import relative_humidity

from ncsuCLOUDS import io
from ncsuCLOUDS.units import addUnits

figKWArgs  = {'dpi' : 300}
markersize = 2


from idealized_tests import uplot

# Charlotte NC
lat   = numpy.full( 1,  35.2271 )
lon   = numpy.full( 1, -80.8431 )
date0 = datetime(2020, 6, 21, 12, 0, 0, tzinfo=pytz.timezone('US/Eastern'))

date0 = date0.astimezone( pytz.timezone('UTC'))

dates = date_range( '2020-01-01T00:00:00', 
                    '2021-01-01T00:00:00', 
                    freq      = '1H',
                    inclusive = 'left'
)

Tair   = numpy.arange( 20, 35, 1 ) #'degree_Celsius')
Tdew   = numpy.arange( 20, 35, 1 ) #'degree_Celsius')

Tair0  = numpy.full( Tair.size,   35.0)# 'degree_Celsius')
Tdew0  = numpy.full( Tair.size,   20.0)# 'degree_Celsius')
solar0 = numpy.full( Tair.size, 1000.0)# 'watt/meter**2')
pres0  = numpy.full( Tair.size, 1013.0)# 'hPa')
speed0 = numpy.full( Tair.size,    0.5)# 'meter/second')


methods  = ['liljegren', 'dimiceli', 'bernard', 'stull']
#methods  = ['dimiceli']
lineFmt  = ['-r',        '-k',       '-g']
lineFmt  = {m:l for m, l in zip( methods, lineFmt )}

year   = numpy.full( solar0.size, date0.year )
month  = numpy.full( solar0.size, date0.month)
day    = numpy.full( solar0.size, date0.day  )
hour   = numpy.full( solar0.size, date0.hour )
minute = numpy.full( solar0.size, date0.minute )
 
def calcRange( *args, ref = None ):

  lims = (
    min( [numpy.nanmin(arg) for arg in args] ), 
    max( [numpy.nanmax(arg) for arg in args] )
  )
  if ref is not None:
    lims = [ min( [ref[0], lims[0]] ), max( [ref[1], lims[1]] ) ]

  return numpy.asarray( lims )

def idealized( outpath = None ):
  kwargs = {
      'linewidth'     : 3.0,
  }
  xyrange = [15, 40]
  
  fig, axes = plt.subplots( 1, 2, figsize=(6,3) )
  
  #solar0, cosz, f_db = solar_parameters( lat, lon, year, month, day, hour, minute, solar0 )
  
  rad = numpy.full( Tair.size, 0, dtype = numpy.int32 )
  
  rh  = relative_humidity( Tair, Tdew0 )
  
  liljegren = liljegren_wetbulb(Tair, Tdew0, pres0)#, speed0, solar0, f_db, cosz, rad )
  dimiceli  = wetBulb( Tair, Tdew0, 'dimiceli')
  stull     = wetBulb( Tair, Tdew0, 'stull')
  bernard   = bernard_wetbulb( Tair, Td = Tdew0 )
  
  axes[0].plot( 
      Tair, liljegren, '-k', 
      Tair, dimiceli,  '-r', 
      Tair, bernard,   '-g', 
      Tair, stull,     '-b',
      **kwargs
  )
  uplot( axes[:1], Tair, rh, 
      labelsize     = 'medium', 
      labelrotation = 0,
      **kwargs 
  )
  rh  = relative_humidity( Tair0, Tdew )
  
  liljegren = liljegren_wetbulb(Tair0, Tdew, pres0)#, speed0, solar0, f_db, cosz, rad )
  dimiceli  = wetBulb( Tair0, Tdew, 'dimiceli')
  stull     = wetBulb( Tair0, Tdew, 'stull')
  bernard   = bernard_wetbulb( Tair0, Td = Tdew )
  axes[1].plot( Tdew, liljegren, '-k', 
      label = 'Liljegren', 
      **kwargs 
  ) 
  axes[1].plot( Tdew, dimiceli,  '-r', label = 'Dimiceli',  **kwargs ) 
  axes[1].plot( Tdew, bernard,   '-g', label = 'Bernard',   **kwargs ) 
  axes[1].plot( Tdew, stull,     '-b', label = 'Stull',     **kwargs )
  uplot( axes[1:], Tdew, rh, 
      label         = 'RH', 
      labelsize     = 'medium', 
      ylabel        = 'Relative Humidity', 
      labelrotation = 0,
      **kwargs
  )
  
  axes[0].set_ylabel( 'Psychrometric Wet Bulb [$^{\circ}$C]' )
  axes[0].set_xlabel( 'Air Temperature [$^{\circ}$C]' )
  axes[1].set_xlabel( 'Dew Point Temperature [$^{\circ}$C]' )
  
  #  ax.set_zorder( ax2.get_zorder()+1 )
  #  ax.patch.set_visible( False )
  
  labels = ['a', 'b', 'c', 'd', 'e']
  for i, ax in enumerate(axes): 
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.text( 0.05, 0.95, 
      f"({labels[i]})",
      fontsize            = 'small',
      transform           = ax.transAxes,
      horizontalalignment = 'left',
    )

  #  ax.set_aspect('equal')
  #  ax.set_ylim( xyrange )
  #  ax.set_xlim( xyrange )
  
  fig.subplots_adjust( left=0.1, bottom = 0.25, top=0.98, wspace=0.45 )
  fig.legend(ncol = 5, loc='lower center', fontsize='small')
  if outpath is None:
    plt.show()
  else:
    fname = os.path.join( outpath, 'Idealized WBT compare.png' )
    fig.savefig( fname, **figKWArgs )

########################################
def realData( dataDir, outpath = None ):

  refDate    = pytz.timezone('UTC').localize( datetime(2020, 2, 1, 0) )
  lims       = [float('inf'), -float('inf')]

  data, meta = io.read( dataDir ) 
  dates      = data.index.get_level_values( data.index.names.index( 'datetime' ) )


  fig, axes  = plt.subplots( 2, 4, figsize=(6,4) )

  for i in range( 2 ):
    idx  = dates < refDate if (i == 0) else dates > refDate
    tmp  = data[idx]
    Tair = addUnits(   'airtemp2m|K', tmp, meta ).to('degree_Celsius').m 
    Tdew = addUnits(   'dewtemp2m|K', tmp, meta ).to('degree_Celsius').m
    pres = addUnits(   'airpres|kPa', tmp, meta ).to('hPa').m
    Twb  = addUnits( 'wetbulbtemp|K', tmp, meta ).to('degree_Celsius').m

    bernard   =   bernard_wetbulb( Tair, Td = Tdew)
    liljegren = liljegren_wetbulb( Tair, Tdew, pres)
    dimiceli  =           wetBulb( Tair, Tdew, 'dimiceli')
    stull     =           wetBulb( Tair, Tdew, 'stull')

    lims = calcRange( bernard, liljegren, dimiceli, stull, ref=lims )

    for j, y in enumerate( [liljegren, dimiceli, bernard, stull] ):
      idx = numpy.where( numpy.isfinite( Twb ) & numpy.isfinite( y ) )
      r   = numpy.corrcoef( Twb[idx], y[idx] )[0,1]
      axes[i,j].plot( Twb, y, markersize,
            marker='.', 
            color='black' 
      )
      axes[i,j].text( 0.95, 0.05, 
        f"r = {r:0.3f}\nn = {idx[0].size}",
        fontsize            = 'small',
        transform           = axes[i,j].transAxes,
        horizontalalignment = 'right',
      )

  lims = lims + [-2, 2]

  #axes[0,0].set_title('Jan. 2020')
  #axes[0,1].set_title('Jun. 2020')
    
  axes[0,0].set_title('Liljegren')
  axes[0,1].set_title('Dimiceli')
  axes[0,2].set_title('Bernard')
  axes[0,3].set_title('Stull')

  for i in range(2):
    axes[i,0].set_ylabel( "$T_{w}\, [^{\circ}\\text{C}]$" )
  for i in range(4):
    axes[-1,i].set_xlabel( "ECONet $T_{w}\, [^{\circ}\\text{C}]$" ) 
    if i > 0:
      axes[0,i].get_yaxis().set_ticklabels( [] )
      axes[1,i].get_yaxis().set_ticklabels( [] )

  for ax in axes.flatten():
    ax.plot( lims, lims, linestyle='-', color='gray', zorder=0 )
    ax.set_xlim( lims )
    ax.set_ylim( lims )
    ax.set_aspect('equal')

  figLoc = {
    'left'   : 0.085, 
    'right'  : 0.99, 
    'bottom' : 0.085, 
    'top'    : 0.92, 
    'wspace' : 0.1, 
    'hspace' : 0.2
  }

  xx = (figLoc['left'] + figLoc['right']) / 2.0
  fig.text( xx, 0.96, 'January 2020', fontsize='large', horizontalalignment='center' ) 
  fig.text( xx, 0.46, 'June 2020',    fontsize='large', horizontalalignment='center' ) 
  fig.subplots_adjust( **figLoc )

  if outpath is None:
    plt.show()
  else:
    fname = os.path.join( outpath, 'WetBulb compare.png' )
    fig.savefig( fname, **figKWArgs )

if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument( 'datadir', type=str, help='Top-level directory of data files to read in' )
  parser.add_argument( '--outpath', type=str, help='Top-level output directory to save PNG files to' )

  args = parser.parse_args()

  idealized( args.outpath )
  realData( args.datadir, args.outpath )
