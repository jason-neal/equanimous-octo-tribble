"""Code to obtain and find row in Baraffe tables."""
import numpy as np
from typing import Tuple, List, Dict


def age_table(age: float, model: str="2003") -> Tuple[Dict[str, List[float]], List[str]]:
    """Determine the correct Baraffe table to load.

    Parameters
    ----------
    age: float
        Stellar age (Gyr).
    model: str
        Baraffe mdel version to use. options=[2003, 2015].

    Returns
    -------
    model_data: numpy.ndarray
        The correct model table data.
    column_names: list of str
        List of the colunms in the the table.

    """
    if not isinstance(model, str):
        raise ValueError("Model is not the valid type 'str'.")
    elif model not in ["2003", "03", "2015", "2003"]:
        raise ValueError("Model value '{}' is not valid".format(model))

    if model in '2003':
        modelages = ["0.001", "0.005", "0.010", "0.050", "0.100", "0.120",
                     "0.500", "1.000", "5.000", "10.000"]
        base_name = "./Baraffe2003/BaraffeCOND2003-"
        skiprows = 18
        cols = ["M/Ms", "Teff", "L/Ls", "g", "R", "Mv",
                "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    else:
        print("Using 2015 models.")
        modelages = ["0.0005", "0.001", "0.003", "0.004", "0.005", "0.008",
                     "0.010", "0.015", "0.020", "0.025", "0.030", "0.040",
                     "0.050", "0.080", "0.100", "0.120", "0.200", "0.300",
                     "0.400", "0.500", "0.625", "0.800", "1.000", "2.000",
                     "3.000", "4.000", "5.000", "8.000", "10.000"]
        base_name = "./Baraffe2015/BaraffeBHAC15-"
        skiprows = 22
        cols = ["M/Ms", "Teff", "L/Ls", "g", "R/Rs", "Li/Li0", "Mv", "Mr",
                "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]

    # Find closest model age.
    model_age = min(modelages, key=lambda x: abs(float(x) - age))  # Closest one
    model_id = "p".join(str(model_age).split("."))   # Replace . with p in number str
    model_name = base_name + model_id + "Gyr.dat"

    model_data = np.loadtxt(model_name, skiprows=skiprows, unpack=False)
    model_data = model_data.T

    # Turn into Dict of values
    data_dict = {}
    for i, col in enumerate(cols):
        data_dict[col] = model_data[i]

    return data_dict, cols


def mass_table_search(companion_mass: float, age: float, model: str="2003") -> Dict[str, float]:
    """Search Baraffe tables to find the companion entry given a mass value.

    Parameters
    ----------
    companion_mass: float
        Companion Mass (Msun)
    age: float
        Age of star/system (Gyr).
    model: int
       Year of Barraffe model to use [2003 (default), 2015].

    Returns
    -------
    companion_parameters: list
        Companion parameters from barraffe table, interpolated to the provided mass.

    """
    model_data, cols = age_table(age, model=model)

    ref_val = companion_mass
    ref_col = "M/Ms"
    companion_parameters = table_interpolation(model_data, ref_val, ref_col, age, model)
    return companion_parameters  # as a dictionary


def magnitude_table_search(magnitudes: Dict[str, float], age: float, band: str="K",
                           model: str="2003") -> Dict[str, float]:
    """Search Baraffe tables to find the companion entry given a band magnitude value.

    Parameters
    ----------
    magnitudes: dict
        Dictionary of (band: magnitude) pairs.
    age: float
        Age of star?system (Gyr).
    band: str
        Wavelength band to use.
    model: int
       Year of Barraffe model to use [2003 (default), 2015].

    Returns
    -------
    companion_parameters: list
        Companion parameters from barraffe table, interpolated between the
        rows to the provided magnitude.

    """
    if not isinstance(band, str):
        raise ValueError('Band {0} was given, when onlnot given as a single string.'.format(band))

    model_data, cols = age_table(age, model=model)

    if "M{}".format(band.lower()) not in cols:
        raise ValueError("Band '{0!s}' is not in Baraffe tables.".format(band))

    if band not in magnitudes.keys():
        print("magnitudes", magnitudes)
        raise ValueError("The band '{0!s}' given is not in the given magnitudes".format(band))

    ref_col = "M{}".format(band.lower())

    ref_val = magnitudes[band]

    companion_parameters = table_interpolation(model_data, ref_val, ref_col, age, model)

    return companion_parameters  # as a dictionary


def table_interpolation(data: Dict[str, List[float]], ref_value: float, ref_col: str,
                        age: float, model: str) -> Dict[str, float]:
    """Interpolate table data from dictionary to the reference value.

    Parameters
    ----------
    data: dict
        Dictionary of table data. keys are the colunm headers.
    ref_value: float
        Value of reference parameter interpolating to.
    ref_col: str
        Colunm name string.
    age: float
        Age of star?system (Gyr).
    model: int
       Year of Barraffe model to use [2003 (default), 2015].

    Returns
    -------
    result_parameters: dict(str, float)
        Result from interpolation of each dict item to the reference.

    """
    result_parameters = {}
    for key in data.keys():
        x_data = data[ref_col][::-1]
        y_data = data[key][::-1]

        if x_data[-1] < x_data[0]:
            # Reverse data if not increaseing.
            x_data = x_data[::-1]
            y_data = y_data[::-1]

        result_parameters[key] = np.interp(ref_value, x_data, y_data)

        if isinstance(result_parameters[key], (np.ndarray, list)):
            result_parameters[key] = result_parameters[key][0]

    return result_parameters
