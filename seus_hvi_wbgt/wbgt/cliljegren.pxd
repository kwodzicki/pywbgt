cdef extern from "src/liljegren.c" nogil:
  int calc_wbgt(int year, int month, int day, int hour, int minute, int second, int gmt, int avg,
    float lat, float lon, float solar, float pres, float Tair, float relhum,
	  float speed, float zspeed, float dT, int urban, int use_spa, 
    float* est_speed, float* Tg, float* Tnwb, float* Tpsy, float* Twbg)

  float h_sphere_in_air(float diameter, float Tair, float Pair, float speed)

  int calc_solar_parameters(int year, int month, double day, 
    float lat, float lon, int use_spa, float* solar, float* cza, float* fdir) 

  float Tglobe( float Tair, float rh, float Pair, float speed,
    float solar, float fdir, float cza )

  float Twb( float Tair, float rh, float Pair, float speed,
    float solar, float fdir, float cza, int rad )
