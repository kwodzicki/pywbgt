from datetime import datetime, timedelta

import numpy
from pandas import to_datetime, date_range

import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D

from metpy.units import units

from seus_hvi_wbgt.wbgt import wbgt

from .utils import *

matplotlib.rcParams.update({'font.size': 8})
ALPHA = 0.5

def main( *args, methods = None, nCol = 7, **kwargs ):
  if len( args ) != 2:
    data, meta = c.download()
  else:
    data, meta = args

  if methods is None:
    methods     = ['liljegren', 'dimiceli', 'bernard']                            # Set names for different WBGT methods

  if len(methods) > 3:
    raise Exception( 'Can only test three (3) methods at a time!' )

  markers     = ['o', 'v', '^']                                                 # Set marker types for different WBGT methods
  variables   = ['Tg', 'Tpsy', 'Tnwb', 'Twbg']
  titles      = ['Globe Temperature', 'Psychrometric Wet Bulb', 'Natural Wet Bulb', 'Wet Bulb Globe Temperature']
  ylabels     = ['Globe Temperature (C)', 'Psychrometric WB', 'Natural WB', 'WBGT (C)', 'WBGT (C)']
  xlabels     = ['Globe Temperature (C)', 'Psychrometric WB', 'Psychrometric WB', 'WBGT (C)', 'WBGT (C; Bernard)']

  nVars       = len(variables)                                                  # Number of variables being plotted
  fig, axes   = plt.subplots(2, 2, figsize=(6,5))                                          # Create figure and subplots
  xyrange     = [ [float('inf'), -float('inf')] ] * nVars                       # Initialize ranges from each of the subplots

  cmap        = cm.get_cmap( 'jet' )                                            # Initialize colormap for station names
  nLocations  = data['location'].unique().size - 1                              # Maximum value number of stations; zero (0) based
  locations   = []                                                              # Empty list for station location names
  clines      = []                                                              # Empty list for lines used in custom legend

  RMSE_data   = { v : {'x' : [], 'y' : []} for v in variables }
  RMSE_data   = { m : RMSE_data for m in methods }
  RMSE_data   = {}
  for m in methods:
    for v in variables:
      RMSE_data[ f'{m}:{v}' ] = {'x' : [], 'y' : [] }

  xvals       = [ [] ] * nVars
  yvals       = [ [] ] * nVars

  for color, (location, df) in enumerate( data.groupby('location') ):           # Iterate over all stations within the dataframe

    locations.append( location )                                                # Append the locaiton name to locations list
    clines.append( Line2D([0], [0], color=cmap(color/nLocations), alpha=ALPHA, lw=4) )       # Append line for the station in the legend

    df       = df.set_index('obtime')                                           # Set index for dataframe         
    df.index = df.index.tz_convert(None)                                        # Convert index times to UTC
    df       = df.drop_duplicates()
    dt       = to_datetime( df.index.unique() )                                 # Get only the unique times

    solar    = addUnits( df, 'rad2m_total', dt )                                # Add units to values, convert to appropriate units, and reindex
    pres     = addUnits( df, 'airpres',     dt )
    Tair     = addUnits( df, 'airtemp2m',   dt )
    Tdew     = addUnits( df, 'dewtemp2m',   dt )
    speed    = addUnits( df, 'windspeed2m', dt )

    Tglobe   = addUnits( df, 'blackglobetemp2m', dt )
    Tpsy     = addUnits( df, 'wetbulbtemp2m',    dt )
    Tnwb     = addUnits( df, 'wetbulbtemp2m',    dt )
    Twbgt    = addUnits( df, 'wbgt2m',           dt )

    Tglobe   = Tglobe.to( 'degree_Celsius' ).magnitude                          # Globe temperature as measured at station
    Tpsy     = Tpsy.to(   'degree_Celsius' ).magnitude
    Tnwb     = Tnwb.to(   'degree_Celsius' ).magnitude
    Twbgt    = Twbgt.to(  'degree_Celsius' ).magnitude                          # Wet bulb globe as measured at station

    metaLoc  = meta[ meta['Location ID'] == location ] 
    lon      = metaLoc['Longitude [degrees East]' ].values.astype( numpy.float32 )
    lat      = metaLoc['Latitude [degrees North]' ].values.astype( numpy.float32 )

    print( f"Location : {location}; lat : {lat}; lon : {lon}" )    

    for midx, method in enumerate( methods ):

      print( f'  {method}' )

      res = wbgt(method,
        lat, lon, dt, 
        solar, pres, Tair, Tdew, speed, 
      )

      #for vidx, (xx, yy) in enumerate( zip( [Tglobe, Twbgt, TwbgtB], [res['Tg'], res['Twbg'], res['Twbg']] ) ):
      for vidx, (var, xx) in enumerate( zip( variables, [Tglobe, Tpsy, Tnwb, Twbgt] ) ):
        yy = res.get( var, None )
        if yy is None: continue

        xyrange[vidx] = [min( [numpy.nanmin(xx), numpy.nanmax(yy), xyrange[vidx][0]] ), 
                         max( [numpy.nanmin(xx), numpy.nanmax(yy), xyrange[vidx][1]] ) ]

        print( f'    Var index : {vidx};  Range: {xyrange[vidx]}' )

        key = f'{method}:{var}'
        RMSE_data[key]['x'].append( xx )
        RMSE_data[key]['y'].append( yy )
  
        axes[vidx//2][vidx%2].scatter( xx, yy, 
          marker = markers[midx], 
          c      = [cmap(color/nLocations)]*len(xx), 
          label  = location, 
          alpha  = ALPHA )

  rmse = {}
  for key, vals in RMSE_data.items():
    try:
      xx   = numpy.concatenate( vals['x'] )
      yy   = numpy.concatenate( vals['y'] )
    except:
      continue

    N    = numpy.sum( numpy.isfinite(xx) & numpy.isfinite(yy) )
    rmse =  numpy.sqrt( numpy.nansum( (xx-yy)**2 ) / N )
    print( f'{key} - {rmse}' )

  xyrange = numpy.asarray( xyrange ) + [-2, 2]
  ii = 0
  for i, axi in enumerate(axes):
    for j, ax in enumerate( axi ):
      ax.set_aspect( 'equal' )
      ax.set_xlabel( xlabels[ii] )
      ax.set_ylabel( ylabels[ii] )
      ax.set_xlim( xyrange[ii] )
      ax.set_ylim( xyrange[ii] )
      ax.set_title( titles[ii] )
      ax.plot( xyrange[ii], xyrange[ii], '-k' )
      ii += 1

  fig.legend( flip(clines, nCol), flip(locations, nCol),
         ncol     = nCol, 
         loc      = [0.02, 0.02],
         fontsize = 'x-small'
  )
 

  fig.legend( 
    [Line2D([0], [0], marker=m, markersize=5, linewidth=0.0) for m in markers], 
    [m.title() for m in methods], 
    loc = [0.85, 0.02],
  ) 

  plt.subplots_adjust(
    left   = 0.1,
    bottom = 0.25, 
    right  = 0.98, 
    top    = 0.95, 
    wspace = 0.05, 
    hspace = 0.45
  )

  fig.savefig( '/Users/kwodzicki/Documents/NCICS/NOAA_Dello_113020_SE-HVI-WBGT/Figures/WBGT_by_station_method.png',
     dpi = 300 )
