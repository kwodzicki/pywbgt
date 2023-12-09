"""
Header definitions for C imports

To make the liljegren C-code available to
cython, muse create hearder information

"""

cdef extern from "src/liljegren_c.c" nogil:
    # Expose define to cython
    float _D_GLOBE       "D_GLOBE"
    float _MIN_SPEED     "MIN_SPEED"
    float _CZA_MIN       "CZA_MIN"
    float _NORMSOLAR_MAX "NORMSOLAR_MAX"
    float _SOLAR_CONST   "SOLAR_CONST"
    float _REF_HEIGHT    "REF_HEIGHT"

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
        float* solar,
        float* cza,
        float* fdir,
    )

    int stab_srdt(
        int daytime,
        float speed,
        float solar,
        float dT,
    )

    float est_wind_speed(
        float speed,
        float zspeed,
        int stability_class,
        int urban,
        float min_speed,
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
