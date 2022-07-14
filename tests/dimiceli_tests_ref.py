#!/usr/bin/env python3

from numpy import deg2rad

from metpy.calc import wet_bulb_temperature
from metpy.units import units
from seus_hvi_wbgt.wbgt.dimiceli import *
from seus_hvi_wbgt.wbgt.utils import *


prelim_tests = {
	'20100909' : {
		'f_db' 		: 0.67,
		'f_dif'		: 0.33,
  	'Tair'			: units.Quantity(30.0,  'degree_Celsius'),
		'Tdew'			: units.Quantity(20.56, 'degree_Celsius'),
		'pres'				: units.Quantity(30.08, 'inHg'),
		'solar'				: units.Quantity(336,   'watt per meter**2'),
		'z'				: units.Quantity(38.44, 'degree'),
		'speed'				: units.Quantity(6,     'meters per second'),
		'RH'			: 67.5,
		'Tw'			: 26.77550,
		'Ea'			: 22.64868,
		'Epsa'		: 0.897936,
		'B'				: 3614652151,
		'C'				: 1197110170,
		'Tg'			: 33.02,
		'E'				: 28.54383,
		'Tg act'	: 32.78,
    'WBGT'		: 28.35,
		'Actual'	: 28.30
	},
	'20100910' : {
		'f_db' 		: 0.75,
		'f_dif' 	: 0.25,
  	'Tair'			: units.Quantity(33.89, 'degree_Celsius'),
		'Tdew'			: units.Quantity(24.44, 'degree_Celsius'),
		'pres'				: units.Quantity(29.75, 'inHg'),
		'solar'				: units.Quantity(754,   'watt per meter**2'),
		'z'				: units.Quantity(36.65, 'degree'),
		'speed'				: units.Quantity(7,     'meter per second'),
		'RH'			: 54.27,
		'Tw'			: 26.575561,
		'Ea'			: 28.487086,
		'Epsa'		: 0.9278434,
		'B'				: 7098450550,
		'C'				: 1309071170,
		'Tg'			: 39.31,
		'E'				: 28.584421,
		'Tg act'	: 39.44,
    'WBGT'		: 29.85,
		'Actual'	: 29.88
	},
	'20100917' : {
		'f_db' 		: 0.75,
		'f_dif' 	: 0.25,
  	'Tair'			: units.Quantity(34.44, 'degree_Celsius'),
		'Tdew'			: units.Quantity(21.11, 'degree_Celsius'),
		'pres'				: units.Quantity(30.05, 'inHg'),
		'solar'				: units.Quantity(579,   'watt per meter**2'),
		'z'				: units.Quantity(41.41, 'degree'),
		'speed'				: units.Quantity(3.7,   'meter per second'),
		'RH'			: 52.0,
		'Tw'			: 26.34527,
		'Ea'			: 22.48619,
		'Epsa'		: 0.897013,
		'B'				: 56176782,
		'C'				: 90440592,
		'Tg'			: 40.65,
		'E'				: 28.24707,
		'Tg act'	: 40.56,
    'WBGT'		: 30.02,
		'Actual'	: 30.00
	}
}

def relDiff( true, exp ):

  return (exp-true)/true * 100.0

def printError( key, true, exp ):
  
  print( f"{key:<12}  : {true:>20.4f} vs. {exp:>20.4f}  ---> {relDiff(true, exp):10.4f}%" )


for key, val in prelim_tests.items():
  #val['pres'] *=   33.8638864
  #val['speed'] *= 1609.344
  #val['z']  = deg2rad( val['z'] ) 

  val['pres'] = val['pres'].to('inHg') 
  val['speed'] = val['speed'].to('meters per minute')
  val['z'] = val['z'].to('radians')

  ea   = atmosVaporPres( val['Tair'].m, val['Tdew'].m, val['pres'].m )
  eps  = thermalEmissivity( val['Tair'].m, val['Tdew'].m, val['pres'].m )
  B    = factorB( val['Tair'].m, val['Tdew'].m, val['pres'].m, val['solar'].m, val['f_db'], val['f_dif'], val['z'].m  )
  C    = factorC( val['speed'].m ) 
  RH   = relativeHumidity( val['Tair'].m, val['Tdew'].m )
  Tw   = wetBulb( val['Tair'].m, val['Tdew'].m )
  Tg   = globeTemperature( val['speed'].m, val['Tair'].m, val['Tdew'].m, val['pres'].m, val['solar'].m, val['f_db'], 
          val['f_dif'], val['z'].m )
  WBGT = dimiceli( **val )#**{k : v.m if isinstance(v, units.Quantity) else v for k, v in val.items()} ) 

  Tw2  = wet_bulb_temperature( units.Quantity(val['pres'], 'hPa'), units.Quantity(val['Tair'], 'degree_Celsius'), units.Quantity(val['Tdew'], 'degree_Celsius') ) 

  print( key )
  printError( "Ea",    			val['Ea'],   ea )
  printError( "Eps",   			val['Epsa'], eps )
  printError( "B",     			val['B'],    B )
  printError( "C",     			val['C'],    C )
  printError( "RH",    			val['RH'],   RH )
  printError( "Tw (Td)", 		val['Tw'],   Tw )
  printError( "Tw (RH)",    val['Tw'],   wetBulbRH( val['Tair'].m, val['RH'] ) )
  printError( "Tw (metpy)",	val['Tw'],   Tw2.m)
  printError( "Tg",    			val['Tg'],   Tg )
  printError( "WBGT",  			val['WBGT'], WBGT.m )
  #print( dimiceli( **val ) )

  print()
