import numpy as np
import sasktran as sk
from sasktran.tir.constants import Boltzmann
from sasktran.tir.util import _atm_reader, _atm_file_path


class ClimatologyAtmosphericState(sk.ClimatologyUserDefined):
    """
    Creates a custom climatology, containing temperature and pressure, for use as an atmospheric state.
    The profiles are from http://eodg.atm.ox.ac.uk/RFM/atm/ and are divided into three folders:
    'FASCODE/ICRCCM Model Atmospheres', 'MIPAS Model Atmospheres (1998)', and 'MIPAS Model Atmospheres (2001)'.
    Each folder contains several '.atm' files, where each file contains a set of profiles describing various
    latitudinal regions and times.

    NOTE: These are one-dimensional atmospheric profiles. Latitude, longitude, and mjd will not affect the returned values.
    """

    def __init__(self, dataset: str = 'fascode', climatology: str = 'std', interp: str = 'linear',
                 spline: bool = False):
        """
        Create an atmospheric state from data files. By default uses the US Standard Atmosphere.

        Parameters
        ----------
        dataset : str, optional
            Name of dataset to use. Valid choices are 'fascode', 'mipas_1998', and 'mipas_2001'. Default is 'fascode'.
        climatology : str, optional
            Name of the data file within the chosen dataset to load. Default is 'std'.
        interp : str, optional
            One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'.
        spline: bool, optional
            One of True or False.  If True then a quadratic spline will be used to interpolate, if False the piecewise
            linear interpolation is done.  Default is False.

        Examples
        --------
        >>> from sasktran.tir.climatology import ClimatologyAtmosphericState
        >>> us_standard_atmosphere = ClimatologyAtmosphericState()  # US Standard Atmosphere
        >>> us_standard_atmosphere.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', latitude=0, longitude=0,\
                                                 altitudes=[10000, 11000], mjd=54372)
        array([223.3, 216.8])
        >>> mipas_2001_polar_summer = ClimatologyAtmosphericState(dataset='mipas_2001', climatology='sum')
        >>> mipas_2001_polar_summer.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', latitude=0, longitude=0,\
                                                  altitudes=[10000, 11000], mjd=54372)
        array([24397. , 20921.5])
        """

        # get climatology data
        if climatology.endswith('.atm'):
            atm_file = climatology
        else:
            atm_file = climatology + '.atm'
        file_path = _atm_file_path(dataset, atm_file)
        data = _atm_reader(file_path)

        super().__init__(data['HGT'] * 1000.0,
                         {'SKCLIMATOLOGY_TEMPERATURE_K': data['TEM'], 'SKCLIMATOLOGY_PRESSURE_PA': data['PRE'] * 100.0},
                         interp=interp, spline=spline)


class ClimatologySpecies(sk.ClimatologyUserDefined):
    """
    Create a custom climatology using a species profile from data files. By default uses the US Standard Atmosphere.
    The CLIMATOLOGY_HANDLE of the species will be SKCLIMATOLOGY_<species>_CM3. Note that despite the files
    providing profiles of species in PPMV (parts per million by volume) the climatology created will be a profile
    with units of molecules per cm^3. This is because the VMR climatologies are currently not supported in the
    SASKTRAN cross section calculation. The species PPMV is converted to using the temperature, :math:`T`, and
    pressure, :math:`p`, in the specified '.atm' file and using the Ideal Gas Law to compute a neutral density,
    :math:`n`:

    .. math::

        n = p / (k T)

    where k is the Boltzmann constant.  The species VMR is multiplied by the neutral density to obtain the species
    density.

    NOTE: These are one-dimensional atmospheric profiles. Latitude, longitude, and MJD will not affect the returned
    values.

    Parameters
    ----------
    species : str
        Chemical formula of the species to use. E.g. 'O3', 'CO2', etc.
    dataset : str, optional
        Name of dataset to use. Valid choices are 'fascode', 'mipas_1998', and 'mipas_2001'. Default is 'fascode'.
    climatology : str, optional
        Name of the data file within the chosen dataset to load. Default is 'std'.
    interp : str, optional
        One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'.
    spline: bool, optional
        One of True or False.  If True then a quadratic spline will be used to interpolate, if False the piecewise
        linear interpolation is done.  Default is False.

    Examples
    --------
    >>> from sasktran.tir.climatology import ClimatologySpecies
    >>> mipas_2001_polar_summer_co2 = ClimatologySpecies('CO2', dataset='mipas_2001', climatology='sum')
    >>> mipas_2001_polar_summer_co2.get_parameter('SKCLIMATOLOGY_CO2_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([2.94644516e+15, 2.50292282e+15])


    The following table lists the species supported by each climatology.
    Each dataset contains several different climatologies. Each climatology contains unique profiles of temperature,
    pressure, and VMRs of the major species. The VMRs for the minor species are the same for all climatologies in each
    dataset. However, the number density of each species is computed using the pressure and temperature of the selected
    climatology.

    +------------+---------------------------------+---------------------------------+---------------------------------+
    | Dataset    | Climatologies                   | Major Species                   | Minor Species                   |
    +============+=================================+=================================+=================================+
    | fascode    | tro (Tropical),                 | H2O, CO2, O3, N2O, CO, CH4, O2  | NO, SO2, NO2, NH3, HNO3, OH,    |
    |            | mls (Mid-Latitude Summer),      |                                 | HF, HCl, HBr, HI, ClO, OCS,     |
    |            | mlw (Mid-Latitude Winter),      |                                 | H2CO, HOCl, N2, HCN, CH3Cl,     |
    |            | sas (Sub-Arctic Summer),        |                                 | H2O2, C2H2, C2H6, PH3, COF2,    |
    |            | saw (Sub-Arctic Winter),        |                                 | SF6, H2S, CFCl3, CF2Cl2, CClF3, |
    |            | std (US Standard Atmosphere),   |                                 | CF4, CHCl2F, CHClF2, C2Cl3F3,   |
    |            |                                 |                                 | C2Cl2F4, C2ClF5, CCl4, ClONO2,  |
    |            |                                 |                                 | N2O5, HNO4                      |
    +------------+---------------------------------+---------------------------------+---------------------------------+
    | mipas_1998 | day_imk (Mid-Latitude Day)      | N2, O2, O3P, CO2, O3, H2O, CH4, | HCN, H2O2, F12, F14, F22, COF2, |
    |            | ngt_imk (Mid-Latitude Night)    | N2O, HNO3, CO, NO2, N2O5, ClO,  | OCS, NH3, SO2, CFCl3, C2H2,     |
    |            | win_imk (Polar Winter)          | HOCl, ClONO2, NO                | C2H6, CCl4, SF6, HNO4, CH3Cl,   |
    |            | sum_imk (Polar Summer)          |                                 | CClF3, CHCl2F, C2Cl3F3, C2Cl2F4 |
    +------------+---------------------------------+---------------------------------+---------------------------------+
    | mipas_2001 | day (Mid-Latitude Day)          | N2, O2, CO2, O3, H2O, CH4, N2O, | CClF3, CHCl2F, C2Cl3F3,         |
    |            | ngt (Mid-Latitude Night)        | HNO3, CO, NO2, N2O5, ClO, HOCl, | C2Cl2F4, C2ClF5, CH3Cl, H2S     |
    |            | win (Polar Winter)              | ClONO2, NO, HNO4, HCN, NH3,     |                                 |
    |            | sum (Polar Summer)              | F11, F12, F14, F22, CCl4, COF2, |                                 |
    |            | equ (Equatorial Day)            | H2O2, C2H2, C2H6, OCS, SO2, SF6 |                                 |
    +------------+---------------------------------+---------------------------------+---------------------------------+
    """

    def __init__(self, species: str, dataset: str = 'fascode', climatology: str = 'std', interp: str = 'linear',
                 spline: bool = False):
        # get climatology data
        if climatology.endswith('.atm'):
            atm_file = climatology
        else:
            atm_file = climatology + '.atm'
        file_path = _atm_file_path(dataset, atm_file)
        data = _atm_reader(file_path)

        if species in data:
            # Calculate neutral density from Ideal Gas Law, n = p / (k T) then convert from m^-3 to cm^-3
            air_nd = data['PRE'] * 100.0 / Boltzmann / data['TEM'] / 1000000.0
            species_nd = data[species] / 1.0e6 * air_nd
            species_heights = data['HGT'] * 1000.0
        else:
            # This is a minor species, need to load the minor species data
            if dataset == 'mipas_2001':
                file_path = _atm_file_path(dataset, 'extra.atm')
            elif dataset == 'mipas_1998':
                file_path = _atm_file_path(dataset, 'extra_imk.atm')
            elif dataset == 'fascode':
                file_path = _atm_file_path(dataset, 'minor.atm')
            else:
                raise ValueError('{} is not a valid dataset'.format(dataset))
            data_minor = _atm_reader(file_path)
            pre_minor = np.interp(data_minor['HGT'], data['HGT'], data['PRE'])  # get pressure at minor species heights
            tem_minor = np.interp(data_minor['HGT'], data['HGT'], data['TEM'])
            air_nd_minor = pre_minor * 100.0 / Boltzmann / tem_minor / 1000000.0
            species_nd = data_minor[species] / 1.0e6 * air_nd_minor
            species_heights = data_minor['HGT'] * 1000.0

        super().__init__(species_heights, {'SKCLIMATOLOGY_' + species + '_CM3': species_nd},
                         interp=interp, spline=spline)


class ClimatologySpeciesCustomPT(sk.ClimatologyUserDefined):
    """
    Similar to ClimatologySpecies but the pressure and temperature are given by a user specified climatology.
    Designed to be used when the atmospheric state of an engine is set to some climatology such as MSIS90 or ECMWF,
    so that the conversion from VMR to number density uses the same pressure and temperature.

    The climatology created by this class is still a 1D vertical profile, so the latitude, longitude, and MJD must be
    passed as inputs along with the climatology containing the desired pressure and temperature.

    Examples
    --------
    Use MSIS90 pressure and temperature at 0 latitude, 0 longitude to convert the CH4 VMR given by the US Standard
    atmosphere to number density.

    >>> from sasktran.tir.climatology import ClimatologySpeciesCustomPT
    >>> from sasktran import MSIS90
    >>> clim = ClimatologySpeciesCustomPT('CH4', MSIS90(), latitude=0, longitude=0, mjd=54372)
    >>> clim.get_parameter('SKCLIMATOLOGY_CH4_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([1.45841301e+13, 1.30346981e+13])

    Parameters
    ----------
    species : str
        Chemical formula of the species to use. E.g. 'O3', 'CO2', etc.
    atmospheric_state : sk.Climatology
        Climatology which must define pressure and temperature.
    latitude :float
        Latitude where the pressure and temperature profile will be taken from.
    longitude : float
        Longitude where the pressure and temperature profile will be taken from.
    mjd : float
        Time where the pressure and temperature profile will be taken from.
    dataset : str, optional
        Name of dataset to use. Valid choices are 'fascode', 'mipas_1998', and 'mipas_2001'. Default is 'fascode'.
    climatology : str, optional
        Name of the data file within the chosen dataset to load. Default is 'std'.
    interp : str, optional
        One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'.
    spline: bool, optional
        One of True or False.  If True then a quadratic spline will be used to interpolate, if False the piecewise
        linear interpolation is done.  Default is False.
    """

    def __init__(self, species: str, atmospheric_state: sk.Climatology, latitude: float, longitude: float, mjd: float,
                 dataset: str = 'fascode', climatology: str = 'std', interp: str = 'linear', spline: bool = False):
        # get climatology data
        if climatology.endswith('.atm'):
            atm_file = climatology
        else:
            atm_file = climatology + '.atm'
        file_path = _atm_file_path(dataset, atm_file)
        data = _atm_reader(file_path)

        # get altitudes and species ppmv from file
        if species in data:
            altitudes = data['HGT'] * 1000.0
            species_ppmv = data[species]
        else:
            # This is a minor species, need to load the minor species data
            if dataset == 'mipas_2001':
                file_path = _atm_file_path(dataset, 'extra.atm')
            elif dataset == 'mipas_1998':
                file_path = _atm_file_path(dataset, 'extra_imk.atm')
            elif dataset == 'fascode':
                file_path = _atm_file_path(dataset, 'minor.atm')
            else:
                raise ValueError('{} is not a valid dataset'.format(dataset))
            data_minor = _atm_reader(file_path)
            altitudes = data_minor['HGT'] * 1000.0
            species_ppmv = data_minor[species]

        # set air density from climatology if it exists
        try:
            air_nd = atmospheric_state.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', latitude=latitude,
                                                     longitude=longitude, altitudes=altitudes, mjd=mjd)
        except sk.SasktranError:
            pressure = atmospheric_state.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', latitude=latitude,
                                                       longitude=longitude, altitudes=altitudes, mjd=mjd)
            temperature = atmospheric_state.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', latitude=latitude,
                                                          longitude=longitude, altitudes=altitudes, mjd=mjd)
            # Calculate neutral density from Ideal Gas Law, n = p / (k T) then convert from m^-3 to cm^-3
            air_nd = pressure / Boltzmann / temperature / 1.0e6

        species_nd = species_ppmv / 1.0e6 * air_nd  # convert from PPMV to number density per cm^3

        super().__init__(altitudes, {'SKCLIMATOLOGY_' + species + '_CM3': species_nd}, interp=interp, spline=spline)


class ClimatologyFull(sk.ClimatologyUserDefined):
    """
    Creates a custom climatology from data files. Similar to ``ClimatologySpecies`` and ``ClimatologyAtmosphericState``
    but includes ALL profiles defined in the data file. See ``ClimatologySpecies`` for the species contained in each
    file.
    If minor species are included, they are interpolated to the same altitude grid as the major species.

    Examples
    --------
    >>> from sasktran_tir.climatology import ClimatologyFull
    >>> clim = ClimatologyFull()  # US Standard Atmosphere
    >>> clim.get_parameter('SKCLIMATOLOGY_TEMPERATURE_K', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([223.3, 216.8])
    >>> clim.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([26500., 22700.])
    >>> clim.get_parameter('SKCLIMATOLOGY_O3_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([1.12859623e+12, 1.62974521e+12])

    Parameters
    ----------
    dataset : str, optional
        Name of dataset to use. Valid choices are 'fascode', 'mipas_1998', and 'mipas_2001'. Default is 'fascode'.
    climatology : str, optional
        Name of the data file within the chosen dataset to load. Default is 'std'.
    include_minor : bool, optional
        If True, minor species are included in the climatology. Default is True.
    interp : str, optional
        One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'.
    spline: bool, optional
        One of True or False.  If True then a quadratic spline will be used to interpolate, if False the piecewise
        linear interpolation is done.  Default is False.
    """

    def __init__(self, dataset: str = 'fascode', climatology: str = 'std', include_minor=True, interp: str = 'linear',
                 spline: bool = False):

        # get climatology data
        if climatology.endswith('.atm'):
            atm_file = climatology
        else:
            atm_file = climatology + '.atm'
        file_path = _atm_file_path(dataset, atm_file)
        data = _atm_reader(file_path)

        # Calculate neutral density from Ideal Gas Law, n = p / (k T) then convert from m^-3 to cm^-3
        air_nd = data['PRE'] * 100.0 / Boltzmann / data['TEM'] / 1000000.0

        values = dict()
        for species in data.keys():
            if species == 'HGT':
                continue
            elif species == 'TEM':
                values['SKCLIMATOLOGY_TEMPERATURE_K'] = data['TEM']
            elif species == 'PRE':
                values['SKCLIMATOLOGY_PRESSURE_PA'] = data['PRE'] * 100.0
            else:
                values['SKCLIMATOLOGY_' + species + '_CM3'] = data[species] / 1.0e6 * air_nd

        if include_minor:
            if dataset == 'mipas_2001':
                file_path = _atm_file_path(dataset, 'extra.atm')
            elif dataset == 'mipas_1998':
                file_path = _atm_file_path(dataset, 'extra_imk.atm')
            elif dataset == 'fascode':
                file_path = _atm_file_path(dataset, 'minor.atm')
            else:
                raise ValueError('{} is not a valid dataset'.format(dataset))
            data_minor = _atm_reader(file_path)
            pre_minor = np.interp(data_minor['HGT'], data['HGT'], data['PRE'])  # get pressure at minor species heights
            tem_minor = np.interp(data_minor['HGT'], data['HGT'], data['TEM'])
            air_nd_minor = pre_minor * 100.0 / Boltzmann / tem_minor / 1000000.0

            for species in data_minor.keys():
                if species == 'HGT':
                    continue
                else:
                    values_on_minor_grid = data_minor[species] / 1.0e6 * air_nd_minor
                    values['SKCLIMATOLOGY_' + species + '_CM3'] = np.interp(data['HGT'], data_minor['HGT'],
                                                                            values_on_minor_grid)

        super().__init__(data['HGT'] * 1000.0, values, interp=interp, spline=spline)
