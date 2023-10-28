"""
Header definitions for C imports

To make the liljegren C-code available to
cython, muse create hearder information

"""

cdef extern from "src/liljegren_c.c" nogil:
    # Expose define to cython
    float _D_GLOBE   "D_GLOBE"
    float _MIN_SPEED "MIN_SPEED"

    # Expose functions to cython
    int calc_wbgt(
        int year,
        int month,
        int day,
        int hour,
        int minute,
        int second,
        int gmt,
        int avg,
        float lat,
        float lon,
        float solar,
        float pres,
        float Tair,
        float relhum,
        float speed,
        float zspeed,
        float dT,
        int urban,
        int use_spa,
        float min_speed,
        float d_globe,
        float* est_speed,
        float* solar_adj,
        float* Tg,
        float* Tnwb,
        float* Tpsy,
        float* Twbg,
    )

    float h_sphere_in_air(
        float diameter,
        float Tair,
        float Pair,
        float speed,
    )

    int calc_solar_parameters(
        int year,
        int month,
        double day, 
        float lat,
        float lon,
        int use_spa,
        float* solar,
        float* cza,
        float* fdir,
    ) 

    float Tglobe(
        float Tair,
        float rh,
        float Pair,
        float speed,
        float solar,
        float fdir,
        float cza,
        float d_globe,
    )

    float Twb(
        float Tair,
        float rh,
        float Pair,
        float speed,
        float solar,
        float fdir,
        float cza,
        int rad,
    )
