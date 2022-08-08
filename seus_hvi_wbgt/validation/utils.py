import os
from datetime import datetime, timedelta
import itertools

import pandas
from metpy.units import units

from ncsuCLOUDS import STREAM
from ncsuCLOUDS.api import CLOUDS

STREAM.setLevel( 0 )

DFILE  = '/Users/kwodzicki/Data/clouds_api_wbgt_data.pic'
MDFILE = DFILE.split('_')[:-1]
MDFILE.append( 'metadata.pic' )
MDFILE = '_'.join( MDFILE )

DATEFMT = '%Y%m%dT%H%M'
def fnameDate( fname, sdate, edate ):

  fname = os.path.splitext( fname )
  return fname[0] + '_' + sdate.strftime(DATEFMT)+'-'+edate.strftime(DATEFMT) + fname[1]

def download( sdate = None, edate = None, ndays = 5, dfile = DFILE, mdfile = MDFILE, dryrun = False):
  """
  Download data from CLOUDS API and write to pickle file

  """

  if edate is None:
    edate = datetime.utcnow()

  if sdate is None:
    sdate = edate - timedelta(days = ndays)

  if sdate > edate:
    raise Exception( 'Starting date cannot be after ending date' )

  maxDays = 10 
  if (edate-sdate).days > maxDays:
    ss = sdate
    ee = ss + timedelta( days=maxDays )

    print( f'Downloading : {ss} - {ee}' )

    ddfile  = fnameDate( dfile,  ss, ee )
    mmdfile = fnameDate( mdfile, ss, ee )
    download( ss, ee, dfile=ddfile, mdfile=mmdfile, dryrun=dryrun)

    ss = ee
    ee = ss + timedelta( days=maxDays )
    if ee > edate: ee = edate

    while ss < edate:
      print( f'Downloading : {ss} - {ee}' )
      ddfile  = fnameDate( dfile,  ss, ee )
      mmdfile = fnameDate( mdfile, ss, ee )
      download( ss, ee, dfile=ddfile, mdfile=mmdfile, dryrun=dryrun)

      ss = ee
      ee = ss + timedelta( days=maxDays )
      if ee > edate: ee = edate

  else:
    print( dfile )
    print( mdfile )

    if os.path.isfile( dfile ) and os.path.isfile( mdfile ):
      print( 'Files already downloaded')
      return dfile, mdfile

    c       = CLOUDS()
    c._var  = "airpres,airtemp2m,dewtemp2m,blackglobetemp2m,windspeed2m,rad2m_total,wbgt2m"
    c.setAttributes('location', 'var', 'value', 'unit', 'obtime', 'timezone') 
    c.setLoc( network = 'ECONET' )
    c.end   = edate
    c.start = sdate 
    #c.setInterval( hour = 1 )

    if not dryrun: data, meta = c.download()

  if not dryrun:
    data.to_pickle( dfile )
    meta.to_pickle( mdfile )

  return dfile, mdfile

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


def loadData( dfile=DFILE, mdfile=MDFILE ):
  """
  Load data and metadata from pickled files

  """

  data = pandas.read_pickle( dfile )
  meta = pandas.read_pickle( mdfile )

  for key in ['Longitude [degrees East]', 'Latitude [degrees North]']:
    data[key] = data['location'].map(meta.set_index('Location ID')[key])

  return data, meta
 
def flip(items, ncol):
  """Reorder data for column spanning in legend"""

  return itertools.chain(*[items[i::ncol] for i in range(ncol)])

def addUnits( df, vName, dt = None ):
  """To wrap values in DataFrame into pint.Quantity"""

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
