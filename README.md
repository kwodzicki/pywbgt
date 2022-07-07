# Southeast US Heat Vulnerability Index (HVI) using Wet Bulb Globe Temperature (WBGT)

This python packge contains three (3) algorithms/models for calculating WBGT from standard, widely available meteorological variables.
The algorthims/models come from the following papers:

  - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S., & Sharp, R. (2008). Modeling the wet bulb globe temperature using standard meteorological measurements. Journal of occupational and environmental hygiene, 5(10), 645-655. 
  - Dimiceli, V. E., Piltz, S. F., & Amburn, S. A. (2013). Black globe temperature estimate for the WBGT index. In IAENG Transactions on Engineering Technologies (pp. 323-334). Springer, Dordrecht.
  - Bernard, T. E., & Pourmoghani, M. (1999). Prediction of workplace wet bulb global temperature. Applied occupational and environmental hygiene, 14(2), 126–134. https://doi.org/10.1080/104732299303296

## Notes on the algorithms

### Dimiceli

In the Dimiceli paper, they provide an equation for calculating the convective heat transfer coefficient; however, the constants a, b, and c are not provided.
In another paper by [Dimiceli][Dimiceli_wbgt2] a constant value of `h = 0.315` is recommended.
Liljegren et al. (2008) does provide a formula/function to estimate this values that *could* be used in the Dimiceli method; this is an option I should code into the Dimiceli method (use the constant, use Dimiceli function if can get their formula working, or use the Liljegren method).

### Bernard

For the Bernard method, they were looking at indoor WBGT, so there is no discussion of globe temperature estimation methods.
From discussion with colleague, we could use the globe temperature estimation methods from either Liljegren or Dimiceli in the Bernard equations.
The main problem with this method is that it does not take into account solar radiation as it was developed for indoor use.
This could be a major draw back of this method as it will not produce accurate values for outdoor/sunny condition as both the globe temperature and natural wet bulb temperature are influence by solar radiation.s

## Road Map

Future versions of code will include:
  - Analysis code is also included used to preform trend analyses, spatial statistcs, etc.
  - Development of HVI
  - Production level code for real-time computation and dicemination of HVI

NCICS/CICSNC K. R. Wodzicki 2022

[Dimiceli_wbgt2]: https://www.weather.gov/media/tsa/pdf/WBGTpaper2.pdf
