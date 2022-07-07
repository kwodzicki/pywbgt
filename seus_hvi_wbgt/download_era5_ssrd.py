import os

import cdsapi

months = list( range( 1, 13 ) )
days   = list( range( 1, 32 ) )
times  = [ f"{i:02d}:00" for i in range(0, 24) ]

months = [1]
days   = [1]
times  = [ f"{i:02d}:00" for i in range(0, 24) ]
outdir = '/Data/ERA5/SE_HVI'

os.makedirs( outdir, exist_ok = True )

if __name__ == "__main__":
  c = cdsapi.Client()

  info = {
    'product_type': 'reanalysis',
    'format'      : 'netcdf',
    'variable'    : 'surface_solar_radiation_downwards',
    'month'       : months,
    'day'         : days,
    'time'        : times,
  }

  years = []
  for year in range(1979, 1981):
    if (year % 5) == 0:
      info['year'] = years
      c.retrieve(
          'reanalysis-era5-single-levels',
          info,
          os.path.join( outdir, f'ERA5_SSRD_{min(years)}-{max(years)}.nc' ) )
      years = [year]
    else:
      years.append( year )

