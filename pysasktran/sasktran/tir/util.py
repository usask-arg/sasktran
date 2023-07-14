import numpy as np
from pathlib import Path
from sasktran.tir.constants import Planck, speed_of_light


def convert_radiance_wavelength_to_wavenumber(radiance_wavelength: np.ndarray, wavelengths_nm: np.array):
    """
    Convert a radiance array from photon/s/cm2/sr/nm to nW/cm2/sr/cm-1. Radiance can be multidimensional (e.g.
    wavelength and line of sight) BUT the 1st dimension must be wavelength in order of increasing wavelength.
    Returned radiance is in order of increasing wavenumber.

    Examples
    --------
    >>> import numpy as np
    >>> from sasktran.tir.util import convert_radiance_wavelength_to_wavenumber
    >>> wavelength = 1.0e7 / np.linspace(775.0, 775.003, 4)[::-1]
    >>> rad_wavelength = np.array([[8.73341577e+12, 1.09506519e+13, 1.32697534e+13, 4.41957608e+12],\
                                   [1.47059359e+13, 1.67048571e+13, 2.02328109e+13, 1.51375632e+13],\
                                   [1.06787128e+13, 1.28978370e+13, 1.59440769e+13, 7.41160317e+12],\
                                   [5.22908338e+12, 7.13817250e+12, 7.17447283e+12, 9.33290088e+11]])
    >>> rad_wavenumber = convert_radiance_wavelength_to_wavenumber()
    >>> print(rad_wavenumber)
    [[1340.29561508 1829.62492777 1838.92926844  239.21680372]
     [2737.11708775 3305.91249235 4086.71027574 1899.70702127]
     [3769.3509751  4281.70432858 5185.9715705  3879.98349759]
     [2238.50198968 2806.81198651 3401.23156514 1132.80188521]]
    """
    wavelen = wavelengths_nm
    # make wavelen same number of dimensions as radiance_wavelength so that they can be multiplied
    for i in range(radiance_wavelength.ndim - 1):
        wavelen = wavelen[:, None]
    radiance_wavenumber = radiance_wavelength * Planck * speed_of_light * wavelen * 1e11
    return radiance_wavenumber[::-1]


def convert_radiance_wavenumber_to_wavelength(radiance_wavenumber: np.ndarray, wavenumbers_cm: np.array):
    """
    Convert a radiance array from nW/cm2/sr/cm-1 to photon/s/cm2/sr/nm. Radiance can be multidimensional (e.g.
    wavenumber and line of sight) BUT the 1st dimension must be wavenumber in order of increasing wavenumber.
    Returned radiance is in order of increasing wavelength.

    Examples
    --------
    >>> import numpy as np
    >>> from sasktran.tir.util import convert_radiance_wavenumber_to_wavelength
    >>> wavenumber = np.linspace(775.0, 775.003, 4)
    >>> rad_wavenumber = np.array([[1340.29561508, 1829.62492777, 1838.92926844,  239.21680372],\
                                   [2737.11708775, 3305.91249235, 4086.71027574, 1899.70702127],\
                                   [3769.3509751 , 4281.70432858, 5185.9715705 , 3879.98349759],\
                                   [2238.50198968, 2806.81198651, 3401.23156514, 1132.80188521]])
    >>> rad_wavelength = convert_radiance_wavenumber_to_wavelength()
    >>> print(rad_wavelength)
    [[8.73341577e+12 1.09506519e+13 1.32697534e+13 4.41957608e+12]
     [1.47059359e+13 1.67048571e+13 2.02328109e+13 1.51375632e+13]
     [1.06787128e+13 1.28978370e+13 1.59440769e+13 7.41160317e+12]
     [5.22908338e+12 7.13817250e+12 7.17447283e+12 9.33290088e+11]]
    """
    wavelen = 1.0e7 / wavenumbers_cm
    for i in range(radiance_wavenumber.ndim - 1):
        wavelen = wavelen[:, None]
    radiance_wavelength = radiance_wavenumber / (Planck * speed_of_light * wavelen * 1e11)
    return radiance_wavelength[::-1]


def _atm_file_path(folder_name: str, file_name: str):
    """
    Determine the location of a specified '.atm' file.

    Parameters
    ----------
    folder_name : str
        Name of folder with the file, e.g. 'fascode', 'mipas_1998', or 'mipas_2001'. Do not include slashes
    file_name : str
        Name of the file, e.g. 'std.atm', 'tro.atm', ...

    Returns
    -------
    Path
        Absolute path to the file.
    """
    base_path = Path(__file__).parent                           # absolute path to this python file
    path_to_file = "data/atm/" + folder_name + '/' + file_name  # relative path to atm file from this python file
    file_path = (base_path / path_to_file).resolve()            # absolute path to atm file
    return file_path


def _atm_reader(atm_file: str):
    """
    Reads in a file with the extension '.atm' which contains information about atmospheric profiles and returns a
    dictionary with the atmospheric parameters as keys and the profiles as values. Refer to the input atm file for
    units of each profile but typically heights are in km, temperatures in K, pressures in mb, and species populations
    in ppmv.

    Parameters
    ----------
    atm_file : str

    Returns
    -------
    dict of [species: str, profile: np.array]
    """
    profiles = dict()  # dictionary where keys are the species tags and values are the profiles
    with open(atm_file) as f:
        num_levels_set = False
        num_levels = 0
        cur_profile = ''
        for line in f:
            if line[0] == '!':  # line is a comment
                pass
            elif line[0] == '*':  # line is a profile header
                cur_profile = ''
                for letter in line[1::]:
                    if letter == ' ':
                        break  # end of species name
                    cur_profile += letter
                if cur_profile == 'END\n':
                    break  # end of file
                else:
                    profiles[cur_profile] = np.array([])
            elif not num_levels_set:  # first uncommented line gives number of levels in each profile
                num_levels = int(list(filter(None, line.split(' ')))[0])
                num_levels_set = True
            else:  # line contains profile values
                no_space = list(filter(None, line.split(' ')))
                no_space_or_comma = []
                for item in no_space:
                    no_space_or_comma.append(list(filter(None, item.split(',')))[0])
                for item in no_space_or_comma:
                    try:
                        profiles[cur_profile] = np.append(profiles[cur_profile], float(item))
                    except ValueError:  # occurs when item == '\n'
                        pass
    return profiles
