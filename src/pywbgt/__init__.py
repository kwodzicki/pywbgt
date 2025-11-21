"""
Package for estimating wetbulb globe temperature

Various algorithms for estimating wetbulb globe
temperature from standard meteorological variables.

"""

import xarray as xr

from .constants import METHODS
from .liljegren import wetbulb_globe as liljegrenWBGT
from .bernard import wetbulb_globe as bernardWBGT
from .dimiceli import wetbulb_globe as dimiceliWBGT
from .dimiceli_nws import wetbulb_globe as dimiceli_nwsWBGT

# Attributes for Dataset output
ATTRS = {
    'Tg': {
        'long_name': 'black_globe_temperature',
        'description': 'Estimated black globe temperature',
    },
    'Tpsy': {
        'long_name': 'psychrometric_wetbulb_temperature',
        'description': 'Estimated psychrometric wetbulb temperature',
    },
    'Tnwb': {
        'long_name': 'natural_wetbulb_temperature',
        'description': 'Estimated natural wetbulb temperature',
    },
    'Twbg': {
        'long_name': 'wetbulb_globe_temperature',
        'description': 'Estimated wetbulb globe temperature',
    },
    'solar': {
        'long_name': 'surface_net_solar_radiation',
        'description': (
            'Adjust surface net solar radiation used in estimation '
            'of temperatures'
        ),
    },
    'speed': {
        'long_name': 'wind_speed',
        'description': (
            'Height adjusted wind speed used in estimation '
            'of termpatures'
        ),
    },
    'min_speed': {
        'description': 'Minimum allowable wind speed threshold',
    },
}


def wbgt(method: str, *args, **kwargs):
    """
    Estimate wet bulb globe temperature

    Wrapper for the various WBGT algorithms provided
    by this package. Set the method to use for
    estimating WBGT and you're off

    Note that a single Xarray Dataset can be input into this function, but
    it must contain data variables with names that match the Arguments
    listed below. For the datetime, lat, and lon arguments, values can
    be coordinates within the Dataset, but the MUST have CF-compliant
    axis attributes to indicate which axis the data represents; e.g.,
    longitude should have {'axis': 'X'}. The naming for these coordinates
    does not matter because ordering is determined from the axis attribute.

    Arguments:
        method (str) : name of the method to use.
        datetime (pandas.DatetimeIndex) : Datetime(s) corresponding to data
        lat (ndarray) : Latitude corresponding to data values (decimal).
            Can be one (1) element array; will be expanded to match dates/data
        lon (ndarray) : Longitude correspondning to data values (decimal).
            Can be one (1) element array; will be expanded to match dates/data
        solar (Quantity) : solar irradiance; units of any power over area
        pres (Qantity) : barometric pressure; units of pressure
        temp_air (Quantity) : air (dry bulb) temperature; units of temperature
        temp_dew (Quantity) : Dew point temperature; units of temperature
        speed (Quatity) : wind speed; units of speed

    Keyword arguments:
        f_db (float) : Direct beam radiation from the sun. Valid for
            the Dimiceli and Bernard methods. Type: fraction
        cosz (float) : Cosine of solar zenith angle. Valid for
            Dimiceli and Bernard methods
        zspeed (Quantity) : Height of the wind speed measurment.
            Default is 10 meters
        wetbulb (str) : Name of wet bulb algorithm to use in the Dimiceli
            algorithm. Valid options are:
            {dimiceli, stull} DEFAULT = dimiceli
        natural_wetbulb (str) : Name of the natural wet bulb algorithm to use
            in the Dimiceli algorithm. Valid options are:
            {malchaire, hunter_minyard} DEFAULT = malchaire
        urban (ndarray) : Boolean flag indicating if "urban" (1) or
            "rural" (0) for wind speed power law exponent. Valid for the
            Liljegren algorithm.
            Can be one (1) element array; will be expanded to match dates/data
        gmt (ndarray) LST-GMT difference  (hours; negative in USA). Valid for
            the Liljegren algorithm.
        avg (ndarray) : averaging time of the meteorological inputs (minutes).
            Valid for the Liljegren algorithm.
        dT (Quantity) : Vertical temperature difference; upper minus lower;
            unit of temperature. Valid for the Liljegren algorithm.
        use_spa (bool) : If set, use the National Renewable Energy
            Laboratory (NREL) Solar Position Algorithm (SPA) to determine
            sun position. Default is to use the build-it, low precision model.
        d_globe (Quantity) : Diameter of the black globe thermometer
            unit of distance. Valid for the Liljegren algorithm.

    Returns:
        dict :
            - Tg : Globe temperatures as Quantity
            - Tpsy : psychrometric wet bulb temperatures as Quantity
            - Tnwb : Natural wet bulb temperatures as Quantity
            - Twbg : Wet bulb-globe temperatures as Quantity
            - solar : Adjusted solar irradiance as Quantity.
            - speed : Estimated 2m wind speed as Quantity:
                will be same as input if already 2m t

    """

    method = method.lower()
    if method not in METHODS:
        raise Exception(
            f'Unsupported WBGT method : {method}! Must be one of {METHODS}'
        )

    if method == 'liljegren':
        func = liljegrenWBGT
    elif method == 'bernard':
        func = bernardWBGT
    elif method == 'dimiceli':
        func = dimiceliWBGT
    elif method == 'dimiceli_nws':
        func = dimiceli_nwsWBGT
    else:
        raise Exception(
            f'Unsupported WBGT method : {method}! Must be one of {METHODS}'
        )

    # Variables for testing of object types and tracking dims/coords
    is_dataset = False
    is_dataarray = False
    coords = dims = None
    ds_coords = ds_dims = None

    # If only one argument input and is a Dataset
    if len(args) == 1 and isinstance(args[0], xr.Dataset):
        is_dataset = True
        *args, ds_dims, ds_coords = parse_dataset(args[0])

    args = list(args)  # Ensure args is a list
    ndim = 0  # Tracker for maximum number of dimensions
    for i, arg in enumerate(args):  # Iterate over all arguments
        update = False  # Track if number of dims was update
        if arg.ndim > ndim:  # If ndim of arg greater than tracker
            update = True  # We updated ndim
            ndim = arg.ndim  # Update ndim
            shape = arg.shape  # Update shape

        # If has a metpy attribute, then quantify the values
        # We do NOT quantify datetime, lat, or lon; skip first 3 loops
        if i > 2 and hasattr(arg, 'metpy'):
            arg = arg.metpy.quantify()

        # If argument is a DataArray, then get the data out of the object
        if isinstance(args[i], xr.DataArray):
            is_dataarray = True
            # If the input was NOT a Dataset (already have coords/dims)
            # and we updated number of dimensions, we now update
            # coords and dims
            if ds_coords is None and update:
                coords = arg.coords
                dims = arg.dims

            # Force a load of the variable and get the data
            arg = arg.load().data

        # Get a 1-D reference to the data and update in args list
        args[i] = arg.ravel()

    # Run the WBGT function
    res = func(*args, **kwargs)

    # If got coords from Dataset, then set coords/dims to vals from Dataset
    if ds_coords is not None:
        coords = ds_coords
        dims = ds_dims

    # Iterate over all items in the resultant dictionary
    for key, val in res.items():
        # Try to reshape the data to the shape of input args
        try:
            val = val.reshape(shape)
        except Exception:
            continue

        # If any input args were DataArray, then try to convert to DataArray
        if is_dataarray:
            try:
                val = xr.DataArray(
                    data=val,
                    dims=dims,
                    attrs=ATTRS.get(key, None),
                )
            except Exception:
                continue

        # Update value in the dictionary
        res[key] = val

    # If input was NOT a Dataset, then just return
    if not is_dataset:
        return res

    # Iterate over all keys again
    for key in tuple(res.keys()):
        # If value IS a DataArray, ignore it
        if isinstance(res[key], xr.DataArray):
            continue

        # When NOT a DataArray, pop off the value from the dict and add it
        # to the coords object
        val = res.pop(key)
        attrs = {
            'units': str(val.units),
            **ATTRS.get(key, {}),
        }
        coords = coords.assign(
            {key: ([], val.magnitude, attrs)}
        )

    # Return a Dataset
    return xr.Dataset(
        data_vars=res,
        coords=coords,
        attrs={
            'method': method,
        }
    )


def parse_dataset(ds):
    """
    Parse Xarray Dataset for input into algorithms

    Expand the T, Y, and X axes of a Dataset to match the dimension
    of data variables for input into the various WBGT algorithms.
    Coordinates MUST have CF-compliante attributes specifying the
    axis they represent; e.g., longitude should have {'axis': 'X'}.

    This function assumes that the data variables are named to match
    argument names for the functions:
        solar: Solar radiation
        pres: Surface pressure
        temp_air: Surface air temperature
        temp_dew: Surface dew point temperature
        speed: Surface wind speed

    """

    coords = {}
    for cname, cval in ds.coords.items():
        axis = cval.attrs.get('axis', '')
        if axis == 'X':
            coords['X'] = cval
        elif axis == 'Y':
            coords['Y'] = cval
        elif axis == 'T':
            coords['T'] = cval

    if len(coords) != 3:
        raise ValueError("Failed to find coordinate(s)!")

    # Expand coord dimensions to match the data
    for dname, dsize in ds['solar'].sizes.items():
        for cname, cval in coords.items():
            if dname in cval.dims:
                continue
            coords[cname] = cval.expand_dims({dname: dsize})

    # Ensure are same shape
    for cname, cval in coords.items():
        coords[cname] = cval.transpose(*ds['solar'].dims)

    return (
        coords['T'],
        coords['Y'],
        coords['X'],
        ds['solar'],
        ds['pres'],
        ds['temp_air'],
        ds['temp_dew'],
        ds['speed'],
        ds['solar'].dims,
        ds['solar'].coords,
    )
