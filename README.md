# Southeast US Heat Vulnerability Index (HVI) using Wet Bulb Globe Temperature (WBGT)

[TOC]

This python packge contains three (3) algorithms/models for calculating WBGT from standard, widely available meteorological variables.
The algorthims/models come from the following papers:

  - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. (2008). Modeling the wet bulb globe temperature using standard meteorological measurements. Journal of occupational and environmental hygiene, 5(10), 645-655. 
  - Dimiceli, V. E., Piltz, S. F., & Amburn, S. A. (2013). Black globe temperature estimate for the WBGT index. In IAENG Transactions on Engineering Technologies (pp. 323-334). Springer, Dordrecht.
  - Bernard, T. E., & Pourmoghani, M. (1999). Prediction of workplace wet bulb global temperature. Applied occupational and environmental hygiene, 14(2), 126–134. https://doi.org/10.1080/104732299303296


## Notes on the algorithms

### Liljegren et al.

The most robust algorithm with explicit calculation of many aspects of the WBGT, this has been the de facto standard for WBGT computation.

### Dimiceli et al.

In the Dimiceli paper, they provide an equation for calculating the convective heat transfer coefficient; however, the constants a, b, and c are not provided.
In another paper by Dimiceli and Piltz a constant value of `h = 0.315` is recommended.
Liljegren et al. (2008) does provide a formula/function to estimate this value that *could* be used in the Dimiceli method; this is an option I should code into the Dimiceli method (use the constant, use Dimiceli function if can get their formula working, or use the Liljegren method).

#### Limitations

The Dimiceli has a limitation around wind speed: when wind speeds are below 1 mile per hour, globe temperatures become exponentially large.
To get around this limitation, wind speeds are clamped to be at least 1 mile per hour.
As the code needs wind speeds in meters per hour, the minimum value (after conversion) for wind speed is 1690.0 meters/hour.

A second limitation is that no formula for the natural wet bulb temperature is provided. 
Two algorithms for computing the natrual wet bulb temperature are included to address this limitation: Malchaire (1976) and Hunter and Minyard (1999)
Either of these algorithms can be selected when running the Dimilceli method.

After some testing, it was discovered that the psychrometric wet bulb algorithm provided by the Dimiceli method was not the most accurate.
This has been addressed by providing the option to use the Stull (2011) algorithm instead of the Dimiceli wet bulb algorithm.

By default the Dimiceli wet bulb and the Malcharie natural wet bulb algorithms are used.
 
### Bernard and Pourmaghani

For the Bernard method, they were looking at indoor WBGT, so there is no discussion of globe temperature estimation methods.
From discussion with colleague, we could use the globe temperature estimation methods from either Liljegren or Dimiceli in the Bernard equations.
The main problem with this method is that it does not take into account solar radiation as it was developed for indoor use.
This could be a major draw back of this method as it will not produce accurate values for outdoor/sunny condition as both the globe temperature and natural wet bulb temperature are influence by solar radiation.s

## Road Map

Future versions of code will include:
  - Analysis code is also included used to preform trend analyses, spatial statistcs, etc.
  - Development of HVI
  - Production level code for real-time computation and dicemination of HVI

## References
  - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. (2008). Modeling the wet bulb globe temperature using standard meteorological measurements. Journal of occupational and environmental hygiene, 5(10), 645-655. 
  - Dimiceli, V. E., Piltz, S. F., & Amburn, S. A. (2013). Black globe temperature estimate for the WBGT index. In IAENG Transactions on Engineering Technologies (pp. 323-334). Springer, Dordrecht.
  - Dimiceli, V. E., & Piltz, S. F., Estimation of Black Globe Temperature for Cacluation of the WBGT Index. https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
  - Bernard, T. E., & Pourmoghani, M. (1999). Prediction of workplace wet bulb global temperature. Applied occupational and environmental hygiene, 14(2), 126–134. https://doi.org/10.1080/104732299303296
  - Malchaire, J. B., (1976) EVALUATION OF NATURAL WET BULB AND WET GLOBE THERMOMETERS. The Annals of Occupational Hygiene, Volume 19, Issue 3-4, December 1976, Pages 251–258, https://doi.org/10.1093/annhyg/19.3-4.251
  - Hunter, Charles H., and C. Olivia Minyard. (1999) Estimating wet bulb globe temperature using standard meteorological measurements." Proceedings of the Conference: 2nd Conference on Environmental Applications, Long Beach, CA, USA. Vol. 18. 1999.
  - Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and Air Temperature, Journal of Applied Meteorology and Climatology, 50(11), 2267-2269. Retrieved Jul 20, 2022, from https://journals.ametsoc.org/view/journals/apme/50/11/jamc-d-11-0143.1.xml

NCICS/CICSNC K. R. Wodzicki 2022
