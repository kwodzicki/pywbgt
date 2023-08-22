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

# JupyterLab and Docker
A Dockerfile and `run-jupyter.sh` shell script are provided to create and run a consistent environment for the package.
This requires the user to have the Docker containerization program installed to work.
Once the use has Docker installed, they can run the shell script using the following:

    bash run-jupyter.sh

This will build the Docker image including installing the `pywbgt` package in development mode.
When the building of the image is complete, a JupyterLab instance will be reachable at `http://localhost:8080`.
The passcode for JupyterLab is displayed in the command line output as part of the server URL.

From here the user can run any Notebooks or scripts they wish with the confidence that their work is in a fully reproducible environment.

# Using the package
As this package is intended to be a complete processing and analysis package, various sub-packages are included for different stages in the analysis pipeline.

## `pywbgt` Modules 
The main modules in the package are used for computing WBGT estimates.
This package includes codes a few algorithms for estimating wetbulb temperature and natural wetbulb temperature that are not associated with the three main WBGT algorithms.
Codes for adjusting wind speeds to given heights above ground level are also included.
The NREL SPA code is also provided as part of this package.

Then, there are the three modules that provide the WBGT algorithms, which are provided as a mix of pure python and Cython wrappers for C codes.
These modules include the main WBGT functions along with various other helper functions associated with the algorithms.
A main `wbgt` function provides a simple API for calling any of the three WBGT algorithms by specifying the algorithm via string and providing the required meteorological parameters.

### Example

    from pywbgt import wbgt
    vals = wbgt(
        'liljegren',
        latitudes,
        longitudes,
        dates,
        solar,
        pres,
        temp_air,
        temp_dew,
        speed,
    )
 
### Input Variables and Unit Handling
There are a large number of input parameters that are required for running the various WBGT algorithms.
These are outlined in the below table:

| Variable | Units | Description |
| - | - | - |
| Latitude | degrees North | Latitude of observation; positive North |
| Longitude | degrees East | Longitude of observation; positive East |
| datetime | various | Datetime of observation as pandas Datetime object |
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

### Keyword Arguments:
There are various keywords associated with three main WBGT algorithms that control different aspects of the algorithms.
The most relevant for most users will be the `zspeed` keyword, which sets the height at which the wind measurement was taken as a unit aware (i.e., `pint.Quantity`) value.
This is important because the WBGT algorithms estimate the 2 meter wind speed for use in their WBGT estimates.
As most weather stations measure wind at 10 meters, this is the default value for the keyword.
However, if the station makes the measurement at 3 feet, then setting the keyword will ensure that the wind speed is adjusted properly to the 2 meter height.
For example:

    from metpy.units import units
    from pywbgt import wbgt
    vals = wbgt('liljegren', lat, lon, ..., zspeed=units.Quantity(3, 'ft'))

For more information about other keyword arguments, please refer to the function docstrings.

# Notes on the algorithms

### Liljegren et al.

The most robust algorithm with explicit calculation of many aspects of the WBGT, this has been the de facto standard for WBGT computation.
One limitation of this algorithm is that the included code for estimating solar position is most accurate for dates between 1950 to 2050.
To compensate for this limitation, the NREL SPA algorithm has been integrated to enable highly accurate solar position calculations.

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
