import re
from datetime import datetime

import numpy
from metpy.units import units

class Record( object ):

  def __init__(self, line=None):

    self._date = None
    self._temp = None

class Wind( object ):
  def __init__(self, info):
    self.units            = "m/s"
    self.direction        = float(info[ 0: 3])
    self.directionQuality = info[ 3] 
    self.type             = info[ 4] 
    self.speed            = float(info[ 5: 9])/10.0
    self.speedQuality     = info[ 9] 

  def __repr__(self):

    return f"<Wind Dir : {self.direction}; Wind Speed : {self.speed} {self.units}>"

  @property
  def u(self):
    return self.speed * numpy.cos( (270.0-self.direction) * numpy.pi/180.0 )
    
  @property
  def v(self):
    return self.speed * numpy.sin( (270.0-self.direction) * numpy.pi/180.0 )

def parseSolar( line ):

  try:
    data = re.findall('GM1(\d{30})', line )[0]
  except: 
    return numpy.nan 
  else:
    dur = int( data[ :4] )
    val = int( data[4:8] )
    if val == 9999:
      return numpy.nan
    return val/1.0 
    
def parseRecord( line ):
  totChars  = int(line[   :  4]) + 105
  usaf      = line[  4: 10]
  wban      = line[ 10: 15]
  dateTime  = datetime.strptime(line[15:27], "%Y%m%d%H%M")
 
  data      = line[ 27]

  lat       = float(line[ 28: 34])/1000.0
  lon       = float(line[ 34: 41])/1000.0
  tcode     = line[ 41: 46]

  elev      = line[ 46: 51]
  letterid  = line[ 51: 56]
  qcName    = line[ 56: 60]

  wind = Wind( line[60:70] )

  ceil      = float(line[ 70: 75])
  ceilMeth  = line[ 75] 
  ceilQC    = line[ 76] 
  CAVOK     = line[ 77]
 
  visDist   = float(line[ 78: 84])
  visQC     = line[ 84] 
  visVar    = line[ 85] 
  visVarQC  = line[ 86] 

  temp      = float(line[ 87: 92])/10.0 * units.degC
  tempQC    = line[ 92] 
  dtemp     = float(line[ 93: 98])/10.0 * units.degC
  dtempQC   = line[ 98] 

  pres      = int(line[ 99: 104])
  pres      = numpy.nan if pres == 99999 else pres / 10.0
  presQC    = line[104]

  if line[105:108] == 'ADD':
    solar = parseSolar( line[108:] )
  else:
    print( 'No additional data block found!' )
  print( usaf, wban, lat, lon, dateTime, pres, temp, dtemp, wind, solar )
  return usaf, wban, lat, lon, dateTime, pres, temp, dtemp, wind, solar

def readISD( fpath ):
  if isinstance(fpath, str):
    if fpath.endswith('.gz'):
      iid = gzip.open( fpath, 'rt' )
    else:
      iid = open(fpath, 'r')

  for record in iid.readlines():
    print( parseRecord( record ) )
    return

  iid.close()
if __name__ == "__main__":
  import sys
  readISD( sys.argv[1] )
