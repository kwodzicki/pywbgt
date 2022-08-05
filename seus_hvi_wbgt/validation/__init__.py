from seus_hvi_wbgt.validation import wbgt_by_wind_method
from seus_hvi_wbgt.validation import wbgt_by_station_method
from seus_hvi_wbgt.validation import wbgt_by_solar_method
from seus_hvi_wbgt.validation import wbgt_by_rh_method


def main():
  data, meta = wbgt_by_wind_method.loadData()

  xx = wbgt_by_wind_method.main(data, meta)
  xx = wbgt_by_station_method.main(data, meta)
  xx = wbgt_by_solar_method.main(data, meta)
  xx = wbgt_by_rh_method.main(data, meta)
  
