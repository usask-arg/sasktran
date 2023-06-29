import sasktranif.sasktranif as skif
import numpy as np
from sasktran.handles import standard_handles
from sasktran.exceptions import wrap_skif_functionfail, SasktranError
from abc import ABCMeta, abstractmethod
from sasktran.util import ArrayWithCallback
import logging
from sasktran import config


class ClimatologyBase(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def skif_object(self, **kwargs) -> skif.ISKClimatology:
        pass


class Climatology(ClimatologyBase):
    """
    A climatology quantifying something in/about the atmosphere. For example climatologies could be used
    to define profiles for temperature, pressure, gas number densities, etc.

    Parameters
    ----------
    name : str
        Name of the climatology to create

    Examples
    --------
    >>> import sasktran as sk
    >>> msis = sk.Climatology('msis90')
    >>> print(msis)
    SasktranIF Climatology: msis90
    Supported Species: ['SKCLIMATOLOGY_PRESSURE_PA', 'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 'SKCLIMATOLOGY_TEMPERATURE_K', 'SKCLIMATOLOGY_O2_CM3']
    """
    @wrap_skif_functionfail
    def __init__(self, name: str):
        self._iskclimatology = skif.ISKClimatology(name)
        self._name = name

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_iskclimatology'}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._iskclimatology = skif.ISKClimatology(state['_name'])

    def __repr__(self):
        supported_species = self.supported_species()

        representation = ('SasktranIF Climatology: {}\n'
                          'Supported Species: {}'.format(self._name, supported_species))

        return representation

    def skif_object(self, **kwargs) -> skif.ISKClimatology:
        """

        Returns
        -------
        skif.ISKClimatology
            Underlying ISKClimatology object which represents the climatology
        """
        return self._iskclimatology

    def supported_species(self):
        """
        Get a list of :ref:`species` identifiers that are supported by this climatology.

        Returns
        -------
        list
            Species identifiers that are supported by this climatology. Identifiers correspond to
            :py:attr:`sasktran.Species.species`
        """
        supported_species = []

        for handle in standard_handles():
            try:
                self.get_parameter(handle, 0, 0, 1000, 54372)
                supported_species.append(handle)
            except SasktranError:
                pass

        return supported_species

    @wrap_skif_functionfail
    def get_parameter(self, species: str, latitude: float, longitude: float, altitudes: np.ndarray, mjd: float):
        """
        Get a profile from the climatology.

        Parameters
        ----------
        species : str
            Species identifier of the species to query. Corresponds to :py:attr:`sasktran.Species.species`.
        latitude : float
            Latitude in degrees.
        longitude : float
            Longitude in degrees.
        altitudes : np.ndarray
            Profile altitude grid in meters.
        mjd : float
            Modified julian date.

        Returns
        -------
        np.ndarray
            Climatology values at the specified altitudes.
        """
        location = [latitude, longitude, 0, mjd]
        self._iskclimatology.UpdateCache(location)
        return self._iskclimatology.GetHeightProfile(species, location, np.asarray(altitudes))[1]


class ClimatologyUserDefined(Climatology):
    """
    A special climatology which handles user defined values.

    Parameters
    ----------
    altitudes: np.ndarray
        Array of altitudes in meters that the profile is to be specified at.  Same shape as values.
    values: dict
        dictionary of species, value pairs. Values should be the same length as altitudes.
    interp: str, optional
        One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'
    spline: bool, optional
        One of True of False.  If True then a quadratic spline will be used to interpolate, if False then piecewise
        linear interpolation is done.  Default is False.

    Examples
    --------
    >>> from sasktran.climatology import ClimatologyUserDefined
    >>> import numpy as np
    >>> altitudes = np.linspace(500, 99500, 100)
    >>> values = np.ones_like(altitudes)
    >>> climatology = ClimatologyUserDefined(altitudes, {'SKCLIMATOLOGY_O3_CM3': values})
    >>> climatology.get_parameter('SKCLIMATOLOGY_O3_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([ 1.,  1.])
    >>> climatology['SKCLIMATOLOGY_O3_CM3'] *= 2
    >>> climatology.get_parameter('SKCLIMATOLOGY_O3_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([ 2.,  2.])
    """
    def __init__(self, altitudes: np.ndarray, values: dict, interp: str='linear', spline: bool=False, out_of_bounds_value=None):
        super().__init__('USERDEFINED_PROFILE')

        self._altitudes = altitudes
        self._values = {key: ArrayWithCallback(val) for key, val in values.items()}

        for val in self._values.values():
            val.add_modified_callback(self._update_internal_clim)

        self._interp = interp
        self._spline = spline
        self._out_of_bounds_value = out_of_bounds_value

        self._climatology_created = False
        self._create_climatology()

    def __setstate__(self, state):
        super().__setstate__(state)
        self._create_climatology()

    def supported_species(self):
        return list(self._values.keys())

    @wrap_skif_functionfail
    def _create_climatology(self):
        if self._interp.lower() == 'linear':
            self._iskclimatology.SetProperty('DoLogInterpolation', 0)
        elif self._interp.lower() == 'log':
            self._iskclimatology.SetProperty('DoLogInterpolation', 1)
        else:
            raise ValueError('interp must be one of linear or log')

        if self._spline:
            self._iskclimatology.SetProperty('DoPiecewiseLinear', 0)
        else:
            self._iskclimatology.SetProperty('DoPiecewiseLinear', 1)

        if self._out_of_bounds_value is not None:
            self._iskclimatology.SetProperty('badvalue', self._out_of_bounds_value)

        self._update_internal_clim()

        self._climatology_created = True

    @wrap_skif_functionfail
    def _update_internal_clim(self):
        self._iskclimatology.SetProperty('Heights', self._altitudes)
        for species, values in self._values.items():
            if hasattr(skif, 'AddGeneratedGlobalClimatologyHandleIfNotExists'):
                skif.AddGeneratedGlobalClimatologyHandleIfNotExists(species)
            else:
                logging.warning('SASKTRANIF Version is out of date')
            self._iskclimatology.SetPropertyUserDefined(species, values)

    def __getitem__(self, item):
        return self._values[item]

    @wrap_skif_functionfail
    def __setitem__(self, item, value):
        self._values[item] = ArrayWithCallback(value)
        self._iskclimatology.SetPropertyUserDefined(item, value)

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if hasattr(self, '_climatology_created') and self._climatology_created:
            self._update_internal_clim()

    @property
    def altitudes(self):
        return self._altitudes


class ClimatologyUserDefined2D(Climatology):
    """
    A two dimensional user defined climatology in altitude and angle within a user specified plane.

    Parameters
    ----------
    angles_deg: np.ndarray
        One dimensional array of the angular grid in degrees
    alts_m: np.ndarray
        One dimensional array of the height grid in meters
    values: dict
        dictionary of species, value pairs. Values should be two dimensional arrays of shape (len(angles_deg), len(alts_m))
    reference_vector: np.ndarray
        Length 3 unit vector defining the location where the angle is 0, i.e., the x-axis of the plane
    normal_vector: np.ndarray
        Length 3 unit vector defining the normal vector of the plane.
    """
    def __init__(self, angles_deg: np.ndarray, alts_m: np.ndarray, values: dict, reference_vector: np.ndarray,
                 normal_vector: np.ndarray):
        self._angular_grid = angles_deg
        self._alt_grid = alts_m
        self._reference = reference_vector
        self._normal = normal_vector
        self._values = {key: ArrayWithCallback(val) for key, val in values.items()}

        super().__init__('USERDEFINED_PROFILE_PLANE')

        self._climatology_created = False
        self._create_climatology()

    def supported_species(self):
        return list(self._values.keys())

    @wrap_skif_functionfail
    def _create_climatology(self):
        self._iskclimatology.SetProperty('normalandreference', np.concatenate((self._normal, self._reference)))

        self._update_internal_clim()

        self._climatology_created = True

    @wrap_skif_functionfail
    def _update_internal_clim(self):
        self._iskclimatology.SetProperty('Heights', self._alt_grid)
        self._iskclimatology.SetProperty('Angles', self._angular_grid)
        for species, values in self._values.items():
            if hasattr(skif, 'AddGeneratedGlobalClimatologyHandleIfNotExists'):
                skif.AddGeneratedGlobalClimatologyHandleIfNotExists(species)
            else:
                logging.warning('SASKTRANIF Version is out of date')
            self._iskclimatology.SetPropertyUserDefined(species, values.flatten())

    def __getitem__(self, item):
        return self._values[item]

    @wrap_skif_functionfail
    def __setitem__(self, item, value):
        self._values[item] = ArrayWithCallback(value)
        self._iskclimatology.SetPropertyUserDefined(item, value.flatten())

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if hasattr(self, '_climatology_created') and self._climatology_created:
            self._update_internal_clim()


class ClimatologyUserDefined3D(Climatology):
    """ A user defined climatology class for profiles which have a geographic dependence.

    For more information see `sk.ClimatologyUserDefined`.

    Parameters
    ----------
    lats : np.ndarray
        The latitude grid (1D) for the climatology. The latitudes for the region will be queried must be defined.
        Latitudes don't need to be evenly spaced.
    lons : np.ndarray
        The longitude grid (1D) for the climatology. Must span 0 to 360 longitude (but doesn't need even spacing).
    alts : np.ndarray
        The altitude grid for the climatology in meters. Altitudes don't need to be evenly spaced. Altitudes at the
        top of that atmosphere (TOA) and the ground (elevation) must be defined.
    values:
        A 3D array of values with dimensions: latitude, longitude, altitude.
    interp: str, optional
        One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'
    spline: bool, optional
        One of True of False.  If True then a quadratic spline will be used to interpolate, if False then piecewise
        linear interpolation is done.  Default is False.
    """
    def __init__(self, lats: np.ndarray, lons: np.ndarray, alts: np.ndarray, values: dict, interp: str='linear', spline: bool=False):
        super().__init__('USERDEFINED_PROFILE3D_LATLONHEIGHT')
        self._latitudes = lats
        self._longitudes = lons
        self._altitudes = alts

        self._values = {key: ArrayWithCallback(val) for key, val in values.items()}

        for val in self._values.values():
            val.add_modified_callback(self._update_internal_clim)

        self._interp = interp
        self._spline = spline

        self._climatology_created = False
        self._create_climatology()

    def __setstate__(self, state):
        super().__setstate__(state)
        self._create_climatology()

    @wrap_skif_functionfail
    def _create_climatology(self):
        self._update_internal_clim()
        self._climatology_created = True

    @wrap_skif_functionfail
    def _update_internal_clim(self):
        self._iskclimatology.SetProperty('Heights', self._altitudes)
        self._iskclimatology.SetProperty('Latitudes', self._latitudes)
        self._iskclimatology.SetProperty('Longitudes', self._longitudes)
        for species, values in self._values.items():
                self._iskclimatology.SetPropertyUserDefined(species, values.flatten())

    def __getitem__(self, item):
        return self._values[item]

    @wrap_skif_functionfail
    def __setitem__(self, item, value):
        self._values[item] = ArrayWithCallback(value)
        self._iskclimatology.SetPropertyUserDefined(item, value)

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if hasattr(self, '_climatology_created') and self._climatology_created:
            self._update_internal_clim()

    def supported_species(self):
        return [species for species in self._values.keys()]


class Labow(Climatology):
    """
    The Labow volume mixing ratio climatology of ozone. This is a climatology of the volume mixing ratio from
    0 km to 60 km in 1 km increments for 18 latitudes in steps of 10 degrees from -85 to +85 for each month of the year.
    The climatology seems to be unpublished as all references refer to an unpublished piece of work by McPeters and
    Labow in 2002 or 2003. We have extended the climatology above 60 km using the differences in scale height between
    the neutrals (7 km scale height ) and ozone (4.5 km scale height) to calculate a new scale height to extrapolate
    the signal at 60 km. We have also extended the volume mixing ratio so it can be converted to an ozone density.
    This requires an atmospheric number density climatology for which the MSIS90 climatology is used:
    """
    def __init__(self):
        super().__init__('O3LABOW')


class Pratmo(Climatology):
    """
    A climatology for the NO2 number density based upon the Prather photo-chemical box model. This climatological
    model is based upon monthly mean solutions to the box model for a range of latitudes, local solar times and
    altitudes.
    """
    def __init__(self):
        super().__init__('NO2PRATMO')

    def __repr__(self):
        supported_species = ['SKCLIMATOLOGY_NO2_CM3']

        representation = ('SasktranIF Climatology: {}\n'
                          'Supported Species: {}'.format(self._name, supported_species))

        return representation

    def supported_species(self):
        return ['SKCLIMATOLOGY_NO2_CM3']


class MSIS90(Climatology):
    """
    A climatology that implements the MSIS-90 atmospheric model, (Hedin 1991).
    This is typically used in Sasktran applications as a quick and robust background atmospheric state.
    The MSIS-90 model was built by the ionospheric/thermospheric community using mass spectrometer and incohorent
    radar scatter data to provide background atmospheric state in the thermosphere under varying geomagnetic conditions.
    The model was coupled to CIRA-86 (Chandra 1990 and Fleming 1990) to provide atmospheric state for altitudes between
    the ground and ~120 km.

    Most sasktran applications only require atmospheric state below 100 km and only utilize the CIRA-86 part of the
    MSIS model. Thus we have configured the default implementation of the MSIS-90 model to only provide the 3 basic
    atmospheric state parameters, pressure, temperature and number density.

    Users may configure the MSIS-90 climatology using the objects properties outlined below to fetch seven other species
    over a larger height range with specific geomagnetic and solar flux conditions. Users are referred to Hedin’s 1991
    publication for details on how to configure the parameters. The current sasktran MSIS model does not provide any
    method to access the TSELEC function describes in the fortran code.

    Parameters
    ----------
    includes: List[str], optional
        Extra species to include in the model, supported species are 'O2', 'O', 'He', 'N2', 'Ar', 'N', and 'H'.
        This option is not case sensitive. Default None.

    f107_avg_flux: float, optional
        The three day average of the F10.7 solar flux in sfu.  Default 150.0.

    f107_prev_day_flux: float, optional
        The F10.7 solar flux for the previous day in sfu.  Default 150.0.

    ap_index: arraylike, optional
        The Ap index that defines the prevailing geomagnetic conditions used by the MSIS model.  The elements
        are

        1. Daily Ap.
        2. 3 hour Ap index for the current time.
        3. 3 hour Ap index for 3 hours before the current time.
        4. 3 hour Ap index for 6 hours before the current time.
        5. 3 hour Ap index for 9 hours before the current time.
        6. Average of eight 3 hour Ap indicies from 12 to 33 hours prior to current time.
        7. Average of eight 3 hour Ap indicies from 36 to 59 hours prior to current time.

        The default value is 4.0 for each of these.

    max_altitude_km: float, optional
        The nominal maximum altitude that will be calculated by the model in km.  Default is 120.

    height_spacing_km: float, optional
        The spacing in km between sample points internally stored by the model.  Default is 1.

    Examples
    --------
    Using MSIS below 120 km, the default settings are typically adequate.

    >>> import sasktran as sk
    >>> msis = sk.MSIS90()
    >>> pressure = msis.get_parameter('SKCLIMATOLOGY_PRESSURE_PA', latitude=0, longitude=0, altitudes=[25000], mjd=54372)
    >>> print(pressure)
    [2541.21983765]

    Above 120 km extra options should be set

    >>> import sasktran as sk
    >>> msis = sk.MSIS90(max_altitude_km=200, includes=['O2', 'N2'])
    >>> N2 = msis.get_parameter('SKCLIMATOLOGY_N2_CM3', latitude=0, longitude=0, altitudes=[150000], mjd=54372)
    >>> O2 = msis.get_parameter('SKCLIMATOLOGY_O2_CM3', latitude=0, longitude=0, altitudes=[150000], mjd=54372)
    >>> print(N2)
    [2.89394292e+10]
    >>> print(O2)
    [2.04250319e+09]


    References
    ----------
    Fleming, E.L., S. Chandra, J.J. Barnett and M. Corney (1990), Zonal mean temperature, pressure, zonal wind, and
    geopotential height as functions of latitude, COSPAR International Reference Atmosphere: 1986,
    Part II: Middle Atmosphere Models, Adv. Space Res., 10, 12, 11-59, doi:10.1016/0273-1177(90)90386-E.

    Chandra, S., E.L. Fleming, M.R. Schoeberl, J.J. Barnett, (1990), Monthly mean global climatology of temperature,
    wind, geopotential height and pressure for 0–120 km, Advances in Space Research, 10, 6, 3-12, doi.org/10.1016/0273-1177(90)90230-W.

    Hedin, A. E. (1991), Extension of the MSIS Thermosphere Model into the middle and lower atmosphere,
    J. Geophys. Res., 96 ( A2), 1159– 1172, doi:10.1029/90JA02125
    """
    def __init__(self, include_o2=False, includes=None, f107_avg_flux=150.0, f107_prev_day_flux=150.0,
                 ap_index=(4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0), max_altitude_km=120,
                 height_spacing_km=1.0):
        super().__init__('MSIS90')

        self._includes = []
        if includes:
            for include in includes:
                self._includes.append(include.lower())

        if include_o2:
            if 'o2' not in self._includes:
                self._includes.append('o2')
            raise DeprecationWarning('Usage of include_o2=True in sk.MSIS90() is deprecated, use includes=[\'O2\'] instead')

        self._f107_avg_flux = f107_avg_flux
        self._f107_prev_day_flux = f107_prev_day_flux
        self._ap_index = ap_index
        self._max_altitude_km = max_altitude_km
        self._height_spacing_km = height_spacing_km

        self._include_map = {'o': 'SKCLIMATOLOGY_O_CM3',
                             'he': 'SKCLIMATOLOGY_He_CM3',
                             'n2': 'SKCLIMATOLOGY_N2_CM3',
                             'ar': 'SKCLIMATOLOGY_AR_CM3',
                             'h': 'SKCLIMATOLOGY_H_CM3',
                             'n': 'SKCLIMATOLOGY_N_CM3'}

        self._update_internal_clim()

    @wrap_skif_functionfail
    def _update_internal_clim(self):
        self.skif_object().SetProperty('F10.7Avg', self._f107_avg_flux)
        self.skif_object().SetProperty('F10.7', self._f107_prev_day_flux)
        self.skif_object().SetProperty('Ap', np.asarray(self._ap_index))
        self.skif_object().SetProperty('MaxHeightKMS', self._max_altitude_km)
        self.skif_object().SetProperty('HeightSpacingKMS', self._height_spacing_km)

        for include in self._includes:
            if include != 'o2':
                self.skif_object().SetProperty('AddSpecies', self._include_map[include])

    def supported_species(self):
        species = super().supported_species()

        if 'o2' not in self._includes:
            species.remove('SKCLIMATOLOGY_O2_CM3')

        return species


class ECMWF(Climatology):
    """
    A climatology of number density, pressure, and temperature based on the ERA interim reanalysis.  This climatology
    requires additional configuration, see :ref:`configuration` for more information.
    """
    def __init__(self):
        super().__init__('ECMWF')

    def __repr__(self):
        representation = ('SasktranIF Climatology: {}\n'
                          'Supported Species: {}'.format(self._name, ['SKCLIMATOLOGY_PRESSURE_PA',
                                                                      'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3',
                                                                      'SKCLIMATOLOGY_TEMPERATURE_K']))
        return representation

    @staticmethod
    def _validate_registry():
        registry_config = config.load_sasktran_registry()

        ecmwf_folder = registry_config['software']['usask-arg']['climatology']['ecmwf']['storage']['erainterimdir']

        if ecmwf_folder.lower() == 'undefined':
            user_config = config.user_config_file_location()
            raise IOError('ECMWF Folder is not set, please add ecmwf_directory: xxxx to the file {}'.format(user_config))

    def get_parameter(self, species: str, latitude: float, longitude: float, altitudes: np.ndarray, mjd: float):
        longitude -= 360 * np.floor(longitude / 360)  # convert longitude to range [0, 360)
        return super().get_parameter(species, latitude, longitude, altitudes, mjd)


class GloSSAC(ClimatologyUserDefined):
    """
    The global stratospheric aerosol extinction climatology at 8 wavelengths from 386-3400nm that spans from 1979-2016
    and incorporates data from numerous satellite datasets. Interpolations in latitude and time use nearest neigbour.

    Note about wavelengths:
        386, 452, 525 and 1020 are the merged aerosol extinction datasets, although 386 and 452nm are not generally
        recommended. 525 and 1020nm channels provide the most robust extinction values so for radiative transfer
        calculations it is recommended to use one of these two channels.

        Other channels (750, 1257 and 3400) are from individual instruments.

    As the GloSSAC climatology is extinction based the species_name is 'SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM', although
    the skif_object is converted to number density for use in sasktran.

    Parameters
    ----------
    wavelength_nm: float
        wavelength of the returned extinction. Possible values are: 386, 452, 525, 750, 1020, 1257, 3400
    extend: bool
        If true the extinction profile is extended above and below the valid range to avoid sharp distcontinuities.
        Profiles are assumed to exponenetially decay outside the valid range. Default is False where all non-valid
        regions are set to zero. The class variables decay_below and decay_above can be used to adjust the decay rates
        below and above the valid range respectively
    """
    def __init__(self, wavelength_nm: int=525, extend=False):

        from netCDF4 import Dataset
        from configparser import ConfigParser, NoSectionError
        import os

        try:
            settings = ConfigParser()
            settings.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'config.ini'))
            self._filename = os.path.normpath(settings.get('Aerosol', 'GloSSACFile'))
            self.dataset = Dataset(self._filename, 'r')
        except (OSError, NoSectionError):
            # Try the global sasktran config file instead
            from sasktran.config import load_user_config, user_config_file_location

            try:
                self._filename = load_user_config()['glossac_file']
                self.dataset = Dataset(self._filename, 'r')
            except KeyError:
                raise IOError('error reading the GloSSAC file, ' + self._filename + ' you may need to update the '
                              'location of the GloSSAC climatology in {}'.format(user_config_file_location()))
            except OSError:
                raise IOError('error reading the GloSSAC file, ' + self._filename + ' you may need to update the '
                              'key glossac_file in {}'.format(user_config_file_location()))

        altitudes = np.arange(0, 100000, 500)
        self._species_name = 'SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM'
        self.wavelength_nm = wavelength_nm
        self.extend = extend
        self.decay_below = 1e-4
        self.decay_above = 5e-4

        super().__init__(altitudes=altitudes, values={self._species_name: altitudes * 0}, interp='linear', spline=False)

    def get_parameter(self, species: str, latitude: float, longitude: float, altitudes: np.ndarray, mjd: float):

        if species not in self._values.keys():
            raise SasktranError('species not supported')

        # convert the time in the nc file to mjds
        year = np.floor(np.array(self.dataset['time'][:])).astype(int)
        month = np.round(np.array(self.dataset['time'][:] % 1 * 12)).astype(int) + 1
        time = np.array([(np.datetime64(str(y + 1979) + '-' + str(m).zfill(2) + '-15') - np.datetime64('1858-11-17')) /
                         np.timedelta64(1, 'D') for y, m in zip(year, month)])

        lat_idx = np.argmin(np.abs(latitude - self.dataset['lat'][:]))
        date_idx = np.argmin(np.abs(mjd - time))

        # in the glossac dataset 2 525nm channels exist (2 is the merged, 3 is OSIRIS only so use 2)
        if np.abs(self.wavelength_nm - 525) <= \
                np.min(np.abs(self.wavelength_nm - self.dataset['measurement_wavelengths'][:] * 1000)):
            wavel_idx = 2
        else:
            wavel_idx = np.argmin(np.abs(self.wavelength_nm - self.dataset['measurement_wavelengths'][:] * 1000))

        ext = self.dataset['Measurements_extinction'][:, lat_idx, date_idx, wavel_idx]
        if type(ext) is np.ma.core.MaskedArray:
            ext = ext.data
        good = (ext != self.dataset['Measurements_extinction'].getncattr('missing_value')) & ~np.isnan(ext)
        if sum(good) == 0:
            raise ValueError('Extinction profile at this wavelength is entirely nans')
        ext = np.interp(altitudes, self.dataset['altitude1'][good] * 1000, ext[good])
        if hasattr(altitudes, '__len__'):
            if self.extend:
                # decay the aerosol above the top point faster than the background
                msis = MSIS90()
                high_alts = altitudes > self.dataset['altitude1'][good][-1] * 1000
                air = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', latitude, longitude,
                                         altitudes[high_alts], mjd)
                air_top = msis.get_parameter('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', latitude, longitude,
                                             self.dataset['altitude1'][good][-1] * 1000, mjd)
                last = np.where(~high_alts)[0][-1]
                ext[high_alts] = ext[high_alts][0] * air / air_top * \
                    np.exp(self.decay_above * -(altitudes[high_alts] - altitudes[last]))

                # exponentially decay the profile below the lower bound
                low_alts = altitudes < self.dataset['altitude1'][good][0] * 1000
                first = np.where(~low_alts)[0][0]
                ext[low_alts] = ext[first] * np.exp(self.decay_below * (altitudes[low_alts] - altitudes[first]))
            else:
                ext[altitudes > self.dataset['altitude1'][-1] * 1000] = 0.0
                ext[altitudes < self.dataset['altitude1'][0] * 1000] = 0.0

        else:
            if altitudes > self.dataset['altitude1'][-1] * 1000:
                ext = 0.0
            if altitudes < self.dataset['altitude1'][0] * 1000:
                ext = 0.0

        return ext

    def skif_object(self, **kwargs):

        clim_values = dict()

        engine = kwargs['engine']
        reference_point = engine.model_parameters['referencepoint']
        latitude = reference_point[0]
        longitude = reference_point[1]
        mjd = reference_point[3]

        opt_prop = kwargs['opt_prop']
        xsec = np.array([opt_prop.calculate_cross_sections(MSIS90(), latitude, longitude, a, mjd,
                                                           self.wavelength_nm).scattering for a in self._altitudes]).flatten()
        clim_values['SKCLIMATOLOGY_AEROSOL_CM3'] = self.get_parameter(self._species_name, latitude, longitude,
                                                                      self._altitudes, mjd) / xsec / 1e5
        userdef_clim = ClimatologyUserDefined(self._altitudes, clim_values)

        return userdef_clim.skif_object()