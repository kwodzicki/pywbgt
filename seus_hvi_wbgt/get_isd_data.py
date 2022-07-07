import os
from urllib.request import urlopen
import gzip
import re

from bs4 import BeautifulSoup as BS
from metpy.calc import relative_humidity_from_dewpoint as RH

from .readers.noaa_isd import parseRecord

URL = 'https://www.ncei.noaa.gov/pub/data/noaa/'

def saveData( uri, data, outdir = '/Data/NCICS/Southeash_HVI/ASOS_ISD', **kwargs):
  if not os.path.isdir( outdir ):
    os.makedirs( outdir )

  fname = uri.split('/')[-1]

  if fname.endswith( '.gz' ):
    fname = os.path.splitext(fname)[0] + '.dat'
  fpath   = os.path.join( outdir, fname )

  with open( fpath, 'w' ) as oid:
    oid.write( f"Latitude Longitude Year GMT_offset Averaging_Period wind_height urban{os.linesep}" )
    records = data.splitlines()
    info    = parseRecord( records[0] )
    oid.write( f"{info[2]:f} {info[3]:f} {info[4].year:d}, 0, 1, 10, 1{os.linesep}" )
    oid.write( f"Jul HrMn    u_30m    u_10m     u_2m    solar     Pres       RH        T   dT30_2   dT10_2{os.linesep}" )
    for record in records:
      info = parseRecord( record )
      out  = "{date:%j} {time:%H%M} {u_30m:>8.2f} {u_10m:>8.2f} {u_2m:>8.2f} {solar:>8.2f} {pres:>8.2f} {RH:>8.2f} {T:>8.2f} {dT30_2:>8.2f} {dT10_2:>8.2f}"
      out  = out.format( date   = info[4],
                         time   = info[4], 
	  										 u_30m  = info[8].speed, 
                         u_10m  = info[8].speed, 
	  										 u_2m   = info[8].speed, 
                         solar  = info[9], 
                         pres   = info[5],
	  									   T      = info[6].m,
	  										 RH     = RH(info[6],  info[7]).m, 
                         dT30_2 = 0,   
                         dT10_2 = 0)
      oid.write( f"{out}{os.linesep}" ) 
    
def download( uri, **kwargs ):
  """Actually download the data"""

  with urlopen( uri ) as resp:
    print( f"Getting data for : {uri} ")
    data = gzip.decompress( resp.read() ).decode()
    saveData( uri, data, **kwargs )
    #if re.search('GM1\d{30}', data):
    #  print( f'   Solar data found : {uri}' )
    #  saveData( uri, data, **kwargs )
 
def getYear(year, usaf = None, wban = None, **kwargs ):
  url = f'{URL}{year}'

  if (usaf is not None) and (wban is not None):
    download( f"{url}/{usaf}-{wban}-{year}.gz", **kwargs )
    return

  else:
    with urlopen( url ) as resp:
      data = resp.read()
    soup = BS(data, 'lxml')
    for i, link in enumerate( soup.find_all( 'a', href=True ) ):
      if link['href'].endswith('.gz'):
        print( f"{url}/{link['href']}" )
        exit()
        download( f"{url}/{link['href']}", **kwargs )
      if i > 10: break

if __name__ == "__main__":
  import sys
  getYear( *sys.argv[1:] ) 
