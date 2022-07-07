import logging
from datetime import datetime

from ncei.mshr import getEMSHR
from ncei.ISD import utils, reader

def checkConsecutive( data, syear, eyear, ndays ):
  """
  
  Arguments:
    data (DataFrame) : Station data
    syear (int) : Starting year for consecutive check
    eyear (int) : Ending year for consecutive check
    ndays (int)  : Number of days in a give year that must have data
      for year to be considered

  Returns:
    tuple : First element is list of years for and second is
      number of consecutive years that meet criteria. Not that
      this list 'resets' if a year does not meet the number of
      days criteria.

  """

  consecutive = [0]
  years       = []
  for i, year in enumerate( range(syear, eyear+1) ):
    years.append( year )
    dd = sum(data['datetime'].dt.year == year) / 24														# hourly data so divide by 24 for number of days of data
    if dd > ndays:
      consecutive.append( consecutive[-1] + 1 )
    else:
      consecutive.append( 0 )

  return years, consecutive[1:]

def checkStation( data, syear, eyear = None, ndays = 300, nyears = 10 ):
  """
 
  Arguments:
    data (DataFrame) : Data from station to check if meets time criteria
    syear (int) : Starting year to check
 
  Keyword arguments:
    ndays (int)  : Number of days in a give year that must have data
      for year to be considered
    nyears (int) : Threshold for how many consecutive years must meet
      the days per year threshold


  """

  if eyear is None: eyear = syear

  years, consecutive = checkConsecutive( data, syear, eyear, ndays )

  return any( [c >= nyears for c in consecutive] )


def main(wban, syear, eyear):
  variables = ['date', 'time', 'air temperature', 'dew point temperature'] 

  urlGen    = utils.genRemoteList( range(syear, eyear+1), wban=wban )
  data      = reader.read( urlGen, variables = variables )
