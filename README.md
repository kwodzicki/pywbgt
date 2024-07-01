This python package contains three (3) algorithms/models for estimating WBGT from standard, widely available meteorological variables.
The algorithms/models come from the following papers:

  - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. (2008). Modeling the wet bulb globe temperature using standard meteorological measurements. Journal of occupational and environmental hygiene, 5(10), 645-655. 
  - Dimiceli, V. E., Piltz, S. F., & Amburn, S. A. (2013). Black globe temperature estimate for the WBGT index. In IAENG Transactions on Engineering Technologies (pp. 323-334). Springer, Dordrecht.
  - Bernard, T. E., & Pourmoghani, M. (1999). Prediction of workplace wet bulb global temperature. Applied occupational and environmental hygiene, 14(2), 126–134. https://doi.org/10.1080/104732299303296

# Table of Contents
[TOC]

# Installation
This package depends on some C/Cython modules; namely the Liljegren and Bernard methods and the National Renewable Energy Laboratory (NREL) Solar Position Algorithm (SPA; Reda and Andreas 2003).
These codes are bundled with this package and *should* compile during installation.

For these to be compiled/installed, a C compiler **must** be installed on your system.
The Cython wrappers for the C-codes are designed to take advantage of multi-core CPUs.
For this to work, the OpenMP library is required and can be installed using your preferred package manager.

To actually install the package, `cd` into this directory and run the following command

    pip install ./

If you must specify the path to your C-compiler and/or OpenMP installations, you can use the `CC`, `CFLAGS`, and `LDFLAGS` environmental variables to specify paths.
For example

    CC=/path/to/gcc CFLAGS=-I/path/to/include LDFLAGS=-L/path/to/lib pip install ./

Or, these can be exported before running pip
    
    export CC=/path/to/gcc
    export CFLAGS=-I/path/to/include
    export LDFLAGS=-L/path/to/lib
    pip install ./

## Advanced Option
As some of the code in this package is dependant on Cython, some program files must be 'cythonized' into C-code before they can be fully compiled.
For ease of use, both the source (`.pyx`) and cythonized (`.c`) versions of the files are provided in the repo/distribution, with the cythonized files being used during installation.
It is possible to re-cythonize the source files using the following command:

    python3 setup.py build_ext

After this, the standard `pip install ./` command can be run to install the package using the new, re-cythonized files.

Note that this re-cythonizing step is necessary when making changes to any of the source Cython files because the install script only compiles the `.c` files.
So, if the `.pyx` files are not run through Cython after making changes, none of the changes will appear in the installed package.

Also note that Cython 3+ must be installed for this to work, as it is not installed as part of this package.

# Using the package
This package includes codes a few algorithms for estimating wetbulb temperature and natural wetbulb temperature that are not associated with the three main WBGT algorithms.
Codes for adjusting wind speeds to given heights above ground level are also included.
The NREL SPA code is also provided as part of this package.

Then, there are the three modules that provide the WBGT algorithms, which are provided as a mix of pure python and Cython wrappers for C codes.
These modules include the main WBGT functions along with various other helper functions associated with the algorithms.
A main `wbgt` function provides a simple API for calling any of the three WBGT algorithms by specifying the algorithm via string and providing the required meteorological parameters.

## Example

    from pywbgt import wbgt
    vals = wbgt(
        'liljegren',
        dates,
        latitudes,
        longitudes,
        solar,
        pres,
        temp_air,
        temp_dew,
        speed,
    )
 
## Input Variables and Unit Handling
There are a large number of input parameters that are required for running the various WBGT algorithms.
These are outlined in the below table:

| Variable | Units | Description |
| - | - | - |
| datetime | various | Datetime of observation as pandas DatetimeIndex object |
| Latitude | degrees North | Latitude of observation; positive North |
| Longitude | degrees East | Longitude of observation; positive East |
| solar | various | Solar irradiance; Quantity |
| pres | various | Barometric pressure; Quantity |
| temp\_air | various | Ambient (dry bulb) temperature; Quantity |
| temp\_dew | various | Dew point temperature; Quantity |
| speed | various | Wind speed; Quantity |

Variables with 'various' units must be `pint.Quantity` objects for seamless unit conversion.
As the WBGT algorithms require different units for input arguments, making these argument unit aware takes some burden off the end-user and moves unit conversions into the functions.
Units can be specified using the `metpy.units` sub-package, which is installed as a dependency of this package.
For example, to tag air temperature values with units of Kelvin, one could do:

    from metpy.units import units
    temp_air = units.Quantity([283, 293], 'K')

Or, using a numpy array:

    import numpy as np
    from metpy.units import units
    temp_air = np.asarray([283, 293]) * units('K')

All the user needs to know is the units of their data and all conversions for the algorithms are done by the algorithms.
Values returned from the algorithms are also `pint.Quanity` objects so that the user knows the units for the values.

## Keyword Arguments:
There are various keywords associated with three main WBGT algorithms that control different aspects of the algorithms.
The most relevant for most users will be the `zspeed` keyword, which sets the height at which the wind measurement was taken as a unit aware (i.e., `pint.Quantity`) value.
This is important because the WBGT algorithms estimate the 2 meter wind speed for use in their WBGT estimates.
As most weather stations measure wind at 10 meters, this is the default value for the keyword.
However, if the station makes the measurement at 3 feet, then setting the keyword will ensure that the wind speed is adjusted properly to the 2 meter height.
For example:

    from metpy.units import units
    from pywbgt import wbgt
    vals = wbgt('liljegren', datetime, lat, lon, ..., zspeed=units.Quantity(3, 'ft'))

Another useful keyword argument is `min_speed`, wherein the minimum speed allowed for the 2m-adjusted wind speeds is set.
After adjusting wind speeds to 2m height, this value is use to clip the wind speeds so that none are below this value.
By default, `min_speed = Quantity(2.0, 'knot')` as ASOS stations report any wind speed of <= 2 knots as calm.
This value can be overridden to the extent possible by the various algorithms.
For example, there is a hard minimum speed of 1690 m/hr (~1 mph) for the Dimiceli method and 0.13 m/s for the Liljegren method.
If a user inputs a value for `min_speed` less than either of these values, the user input is overridden and the hard limit used.
Data dictionaries returned by the methods include a `min_speed` key/value pair indicating the value used in the algorithm.
Note that this value must be a unit-aware object.

For more information about other keyword arguments, please refer to the function docstrings.

## Xarray Support
Xarray DataArray objects are 'supported' for the main `wbgt` function; however, there is some work that the user will likely have to do.
First, the DataArray objects MUST be unit aware objects; i.e., `metpy` integration is enabled/working correctly.
Next, all dimensions of the data must be stacked as the algorithms currently only support 1-D arrays as input.
Then, the DataArray (variables) can be passed to the function.
The output data from the function can then be merged into the stacked Dataset and unstacked.

An example of this process is outlined below:

    import xarray as xr
    from metpy.calc import wind_speed

    stackvar = 'stacked'
    dataset  = xr.open_dataset('/path/to/file.nc')
    dataset  = dataset.stack(
        {stackvar : list(dataset.dims)}
    )

    wetbulb_data = wbgt(
        'dimiceli',
        dataset.time,
        dataset.latitude.values,
        dataset.longitude.values,
        dataset.ssrd,
        dataset.sp,
        dataset.t2m,
        dataset.d2m,
        wind_speed(dataset.u10, dataset.v10),
    )

    dataset = dataset.assign(
        {key : (stackvar, val) for key, val in wetbulb_data.items()}
    ).unstack() 

One major point to note is that the wind speed will likely need to be calculated from u- and v-components.
To do this, use the `metpy.calc.wind_speed()` function to maintain unit information.

If working with accumlated fields (e.g., downward solar radation from a model or reanalysis), some 'unit trickery' may be required.
For example, in the ERA5-Land dataset, solar radiation is accumlated over a given interval as denoted by the `step` variable in the grib files.
This means the data are in units of `Joules / m**2`.
To get to `Watt/m**2` units we can do the following:

    ssrd = (
        dataset.ssrd.metpy.quantify() /
        (dataset.step.dt.seconds*units('second'))
    )

This ensures that the `ssrd` DataArray is explicitly tagged with units using the `.metpy.quantify()` method and then is divided by the accumulation time in seconds.
It is important to note that this will give the average radiation over the entire accumulation period NOT the instanteous value measured at the given model/reanalysis time step.

# Solar position calculations
The Liljegren code provides an algorithm for calculating solar position parameters; however, the algorithm is only valid from 1950 to 2050.
To get around this limiation, the Python pvlib package is used to calculate the solar position using their implementation of the National Renewable Energy Laboratory Solar Postition Algorithm (SPA).
The Python implementation of the SPA code is combined with the solar parameters code of the Liljegren algorithm to create a hybrid function for calculating the adjusted solar irradiance, cosine of solar zenith angle, and fraction of direct beam radiation.
This new `solar_parameters()` function is used in all the included algorithms to compute the parameters requried to estimate WBGT.

# Notes on the algorithms

### Liljegren et al.

The most robust algorithm with explicit calculation of many aspects of the WBGT, this has been the de facto standard for WBGT computation.
One limitation of this algorithm is that the included code for estimating solar position is most accurate for dates between 1950 to 2050.
As previously mentioned, this is overriden by the pvlib SPA implemenation and augmented code for computing the solar parameters.

### Dimiceli et al.

In the Dimiceli paper, they provide an equation for calculating the convective heat transfer coefficient; however, the constants a, b, and c are not provided.
In another paper by Dimiceli and Piltz a constant value of `h = 0.315` is recommended.
Liljegren et al. (2008) does provide a formula/function to estimate this value that *could* be used in the Dimiceli method; this is an option I should code into the Dimiceli method (use the constant, use Dimiceli function if can get their formula working, or use the Liljegren method).

#### Limitations

The Dimiceli has a limitation around wind speed: when wind speeds are below 1 mile per hour, globe temperatures become exponentially large.
To get around this limitation, wind speeds are clamped to be at least 1 mile per hour.
As the code needs wind speeds in meters per hour, the minimum value (after conversion) for wind speed is 1690.0 meters/hour.

A second limitation is that no formula for the natural wet bulb temperature is provided. 
Two algorithms for computing the natural wet bulb temperature are included to address this limitation: Malchaire (1976) and Hunter and Minyard (1999)
Either of these algorithms can be selected when running the Dimiceli method.

After some testing, it was discovered that the psychrometric wet bulb algorithm provided by the Dimiceli method was not the most accurate.
This has been addressed by providing the option to use the Stull (2011) algorithm instead of the Dimiceli wet bulb algorithm.

By default the Dimiceli wet bulb and the Malchaire natural wet bulb algorithms are used.
 
### Bernard and Pourmoghani

The Bernard method focused on indoor WBGT, so there is no discussion of globe temperature estimation methods.
Thus, the method outlined in the article has been modified to include solar radiation as a source in the calculations.
This is done using an approach similar to that of Liljegren et al. (2008), using an iterative approach to determine black globe temperature.

# References
  - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. (2008). Modeling the wet bulb globe temperature using standard meteorological measurements. Journal of occupational and environmental hygiene, 5(10), 645-655. 
  - Dimiceli, V. E., Piltz, S. F., & Amburn, S. A. (2013). Black globe temperature estimate for the WBGT index. In IAENG Transactions on Engineering Technologies (pp. 323-334). Springer, Dordrecht.
  - Dimiceli, V. E., & Piltz, S. F., Estimation of Black Globe Temperature for Calculation of the WBGT Index. https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  - Bernard, T. E., & Pourmoghani, M. (1999). Prediction of workplace wet bulb global temperature. Applied occupational and environmental hygiene, 14(2), 126–134. https://doi.org/10.1080/104732299303296
  - Malchaire, J. B., (1976) EVALUATION OF NATURAL WET BULB AND WET GLOBE THERMOMETERS. The Annals of Occupational Hygiene, Volume 19, Issue 3-4, December 1976, Pages 251–258, https://doi.org/10.1093/annhyg/19.3-4.251
  - Hunter, Charles H., and C. Olivia Minyard. (1999) Estimating wet bulb globe temperature using standard meteorological measurements." Proceedings of the Conference: 2nd Conference on Environmental Applications, Long Beach, CA, USA. Vol. 18. 1999.
  - Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and Air Temperature, Journal of Applied Meteorology and Climatology, 50(11), 2267-2269. Retrieved Jul 20, 2022, from https://journals.ametsoc.org/view/journals/apme/50/11/jamc-d-11-0143.1.xml
  - Reda, I. and Andreas, A. (2003). Solar Position Algorithm for Solar Radiation Applications. 55 pp.; NREL Report No. TP-560-34302, Revised January 2008.

NCICS/CICSNC K. R. Wodzicki 2023
