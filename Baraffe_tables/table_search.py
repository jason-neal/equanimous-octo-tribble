
import numpy as np
from astropy.constants import M_jup, M_sun


def age_table(age, model="2003"):
    """Determine the correct Baraffe table to load.

    Parameters
    ----------
    age: float
        Stellar age (Gyr).
    model: str
        Baraffe mdel version to use. options=[2003, 2015].

    Returns
    -------
    table_name: str

    column_names: list of str

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

    return model_data.T, cols


def mass_table_search(companion_mass, age, model="2003"):
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
    # mass_solar = companion_mass / 1047.56   # covert to solar mass
    # mass_solar = companion_mass * (M_jup / M_sun).value   # covert to solar mass

    model_data, cols = age_table(age, model=model)

    mass_index = cols.index("M/Ms")  # Column of mass

    companion_parameters = dict()
    for col, data in zip(cols, model_data):
        # Interpolate columns to mass of companion
        # The xp column needs to be sorted, increasingly.
        # assert np.all(model_data[mass_index].sort() == model_data[mass_index])
        companion_parameters[col] = np.interp(mass_solar, model_data[mass_index], data)

    return companion_parameters  # as a dictionary


def magnitude_table_search(magnitudes, age, band="K", model="2003"):
    """ Search Baraffe tables to find the companion entry given a band magnitude value.

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
    if isinstance(band, list):
        if len(band) != 1:
            raise ValueError('More than one band was given, when only one required.')
        else:
            band = band[0]

    # mass_solar = companion_mass / 1047.56   # covert to solar mass
    companion_parameters = dict()

    model_data, cols = age_table(age, model=model)

    if "M{}".format(band.lower()) not in cols:
        raise ValueError("Band '{0!s}' is not in Baraffe tables.".format(band))

    if band not in magnitudes.keys():
        print("magnitudes", magnitudes)
        raise ValueError("The band '{0!s}' given is not in the given magnitudes".format(band))

    band_index = cols.index("M{}".format(band.lower()))

    print("Band mag", magnitudes[band])
    # print("Mk", model_data)
    print("model band data", model_data[band_index])

    for col, data in zip(cols, model_data):
        # Interpolate columns to magnitude of companion
        # print("\nThis col =", col, "this col data ", data)
        # print("Band mag", magnitudes[band])
        # print("model data", model_data[band_index])

        # Needs reversing to be increasing in x.
        x_data = model_data[band_index][::-1]
        y_data = data[::-1]
        # print("x_data", x_data, "y_data", y_data)
        companion_parameters[col] = np.interp(magnitudes[band], x_data, y_data)
        if isinstance(companion_parameters[col], (np.ndarray, list)):
            companion_parameters[col] = companion_parameters[col][0]
        # print("{0!s} \tcol interp value = {1:.4f}".format(col, companion_parameters[col]))
    print("comp params", companion_parameters)
    return companion_parameters  # as a dictionary
