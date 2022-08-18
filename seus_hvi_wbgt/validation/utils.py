import os
from datetime import datetime, timedelta
import itertools

import pandas
from metpy.units import units

from ncsuCLOUDS import STREAM
from ncsuCLOUDS.api import CLOUDS, io

STREAM.setLevel( 0 )

DFILE   = '/Users/kwodzicki/Data/clouds_api_wbgt_data'
DATEFMT = '%Y%m%dT%H%M'
def fnameDate( fname, sdate, edate ):

  fname = os.path.splitext( fname )
  return fname[0] + '_' + sdate.strftime(DATEFMT)+'-'+edate.strftime(DATEFMT) + fname[1]

def dateWindow( dt, sdate, edate=None ):

  second = timedelta( seconds = 1 )
  if edate is None:
    return sdate, sdate + dt - second
  return dateWindow( dt, edate + second )

def download( sdate = None, edate = None, ndays = 5, 
    dt        = timedelta(hours=6), 
    fileBase  = DFILE,
    **kwargs ):
  """
  Download data from CLOUDS API and write to pickle file

  """

  if edate is None:
    edate = datetime.utcnow()

  if sdate is None:
    sdate = edate - timedelta(days = ndays)

  if sdate > edate:
    raise Exception( 'Starting date cannot be after ending date' )

  if 'variables' not in kwargs:
    kwargs['variables'] = [
        'rad2m_total',
        'airpres|kPa',
        'windspeed2m|mph',
        'airtemp2m|K',
        'dewtemp2m|K',
        'blackglobetemp2m|K',
        'wetbulbtemp|K'
        'wbgt2m|K'
      ] 

  if (edate-sdate) > dt:
    ss, ee = dateWindow( dt, sdate )

    print( f'Downloading : {ss} - {ee}' )

    fname = fnameDate( fileBase,  ss, ee )
    download( ss, ee, fileBase=fname, **kwargs )

    ss, ee = dateWindow( dt, ss, ee )
    if ee > edate: ee = edate

    while ss < edate:
      print( f'Downloading : {ss} - {ee}' )
      fname = fnameDate( fileBase,  ss, ee )
      download( ss, ee, fileBase=fname, **kwargs)

      ss, ee = dateWindow( dt, ss, ee )
      if ee > edate: ee = edate

  else:
    print( fileBase )
    c       = CLOUDS()
    c.setVariables( *kwargs['variables'])
    c.setLoc( network = 'ECONET' )
    c.end   = edate
    c.start = sdate 
    c.setInterval( hour=1 )

    if not kwargs.get('dryrun', False): _ = c.download( fileBase )

  return fileBase 

def combineDataFiles( files, outfile ):

  files = sorted( files )
  data  = pandas.concat( [pandas.read_pickle( f ) for f in files], ignore_index=True )

  data.to_pickle( outfile )

def combineMetadataFiles( files, outfile ):

  files = sorted( files )
  data  = pandas.read_pickle( files[0] )
  for f in files[1:]:
    data = pandas.merge( data, pandas.read_pickle( f ) )

  data.to_pickle( outfile )


def loadData( dfile=DFILE ):
  """
  Load data and metadata from pickled files

  """

  data, meta = io.read( dfile )

  for key in ['Longitude [degrees East]', 'Latitude [degrees North]']:
    data[key] = data['location'].map(meta.set_index('Location ID')[key])

  return data, meta
 
def flip(items, ncol):
  """Reorder data for column spanning in legend"""

  return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def addUnits( vName, df, meta ):
  """To wrap values in DataFrame into pint.Quantity"""

  dd   = df.loc[:, vName, :] 
  unit = meta['Parameter']
  unit = unit[ unit['requested'] == vName ]['unit'][0] 
  
  return units.Quantity( dd['value'].values, unit )
