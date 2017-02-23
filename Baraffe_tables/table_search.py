import numpy as np


def mass_table_search(companion_mass, age, model="2003"):
    """Search Baraffe tables to find the companion entry given a mass value.

    Parameters
    ----------
    companion_mass: float
        Companion Mass (Mjup)
    age: float
        Age of star/system (Gyr).
    model: int
       Year of Barraffe model to use [2003 (default), 2015].

    Returns
    -------
    companion_parameters: list
        Companion parameters from barraffe table, interpolated to the provided mass.

    """
    mass_solar = companion_mass / 1047.56   # covert to solar mass
    companion_parameters = dict()

    if model in '2003':
        # Find closest age model
        modelages = ["0.001", "0.005", "0.010", "0.050", "0.100", "0.120",
                     "0.500", "1.000", "5.000", "10.000"]
        # Find closest model age.
        model_age = min(modelages, key=lambda x: abs(float(x) - age))
        model_id = "p".join(str(model_age).split("."))   # Replace . with p in number str
        model_name = "./Baraffe2003/BaraffeCOND2003-" + model_id + "Gyr.dat"

        model_data = np.loadtxt(model_name, skiprows=18, unpack=False)

        cols = ["M/Ms", "Teff", "L/Ls", "g", "R", "Mv",
                "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    else:
        print("Using 2015 models.")
        modelages = ["0.0005", "0.001", "0.003", "0.004", "0.005", "0.008",
                     "0.010", "0.015", "0.020", "0.025", "0.030", "0.040",
                     "0.050", "0.080", "0.100", "0.120", "0.200", "0.300",
                     "0.400", "0.500", "0.625", "0.800", "1.000", "2.000",
                     "3.000", "4.000", "5.000", "8.000", "10.000"]
        model_age = min(modelages, key=lambda x: abs(float(x) - age))  # Closest one
        model_id = "p".join(str(model_age).split("."))

        model_name = "./Baraffe2015/BaraffeBHAC15-" + model_id + "Gyr.dat"
        print(model_name)
        model_data = np.loadtxt(model_name, skiprows=22, unpack=False)

        cols = ["M/Ms", "Teff", "L/Ls", "g", "R/Rs", "Li/Li0", "Mv", "Mr",
                "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]

    model_data = model_data.T

    for col, data in zip(cols, model_data):
        # Interpolate columns to mass of companion
        companion_parameters[col] = np.interp(mass_solar, model_data[0], data)

    return companion_parameters  # as a dictionary


def magnitude_table_search(magnitudes, age, band="K", model=2003):
    """ Search Baraffe tables to find the companion entry given a band magnitude value.

    Parameters
    ----------
    magnitudes: dict
        Dictionary of (band: magnitude value) pairs.
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
    # mass_solar = companion_mass / 1047.56   # covert to solar mass
    companion_parameters = dict()

    if model in '2003':
        # Find closest age model
        modelages = ["0.001", "0.005", "0.010", "0.050", "0.100", "0.120",
                     "0.500", "1.000", "5.000", "10.000"]
        model_age = min(modelages, key=lambda x: abs(float(x) - age))  # Closest one
        model_id = "p".join(str(model_age).split("."))

        model_name = "./Baraffe2003/BaraffeCOND2003-" + model_id + "Gyr.dat"

        model_data = np.loadtxt(model_name, skiprows=18, unpack=False)

        cols = ["M/Ms", "Teff", "L/Ls", "g", "R", "Mv",
                "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]

    else:
        print("Using 2015 models")
        modelages = ["0.0005", "0.001", "0.003", "0.004", "0.005", "0.008",
                     "0.010", "0.015", "0.020", "0.025", "0.030", "0.040",
                     "0.050", "0.080", "0.100", "0.120", "0.200", "0.300",
                     "0.400", "0.500", "0.625", "0.800", "1.000", "2.000",
                     "3.000", "4.000", "5.000", "8.000", "10.000"]
        model_age = min(modelages, key=lambda x: abs(float(x) - age))  # Closest one
        model_id = "p".join(str(model_age).split("."))

        model_name = "./Baraffe2015/BaraffeBHAC15-" + model_id + "Gyr.dat"
        print(model_name)
        model_data = np.loadtxt(model_name, skiprows=22, unpack=False)

        cols = ["M/Ms", "Teff", "L/Ls", "g", "R/Rs", "Li/Li0", "Mv", "Mr",
                "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]

    model_data = model_data.T

    if "M" + band.lower() not in cols:
        raise ValueError("Band {} is not in Baraffe tables.".format(band))
    if band not in magnitudes.keys():
        raise ValueError("The band {} given is not in the given magnitudes".format(band))

    band_index = cols.index("M" + band.lower())

    print("Band mag", magnitudes[band])
    print(model_data)
    print("model data", model_data[band_index])
    print("all data")
    for col, data in zip(cols, model_data):
        # Interpolate columns to magnitude of companion
        print("\nThis col =", col)
        print("Band mag", magnitudes[band])
        # print(model_data)
        print("model data", model_data[band_index])
        print("this col data ", data)

        companion_parameters[col] = np.interp(magnitudes[band], model_data[band_index], data)
        print("This col interp value =", companion_parameters[col])

    return companion_parameters  # as a dictionary
