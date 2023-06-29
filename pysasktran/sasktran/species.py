import sasktran as sk
from sasktran.climatology import ClimatologyUserDefined
import numpy as np
import logging


class Species(object):
    """
    An object representing some atmospheric constituent. A species is an object that associates a constituents
    :ref:`optical_property` and :ref:`climatology`.
    
    Examples
    --------
    >>> import sasktran as sk
    >>> no2 = sk.Species(optical_property=sk.NO2Vandaele1998(), climatology=sk.Pratmo())
    >>> print(no2)
    SKCLIMATOLOGY_NO2_CM3
    SasktranIF Climatology: NO2PRATMO
    Supported Species: ['SKCLIMATOLOGY_NO2_CM3']
    SasktranIF Optical Property: NO2_VANDAELE1998
    """
    def __init__(self, optical_property: sk.OpticalProperty, climatology: sk.Climatology, species: str=None):
        if isinstance(optical_property, str):
            logging.warning('Deprecation Warning: The signature Species(species, optical_property, climatology) has been'
                            ' deprecated in favor of Species(optical_property, climatology, species) where species'
                            ' is now optional.')
            self.__init__(species, optical_property, climatology)
        else:
            if species is None:
                if len(climatology.supported_species()) == 1:
                    self._species = climatology.supported_species()[0]
                else:
                    # Some species might have both _VMR and _CM3 handles, like O3Labow, but only the _CM3 one works
                    # in SASKTRAN
                    supported_species = [species for species in climatology.supported_species()
                                         if self.is_sasktran_supported_handle(species)]
                    if len(supported_species) == 1:
                        self._species = supported_species[0]
                    else:
                        raise ValueError('Could not automatically infer the species key, use species= one of {}'.format(climatology.supported_species()))
            else:
                self._species = species
            self._climatology = climatology
            self._optical_property = optical_property

    def __repr__(self):
        return '{}\n{}\n{}'.format(self._species, self._climatology, self._optical_property)

    @property
    def species(self):
        """
        str : The species identifier. If the ``species`` argument was given during construction, this property will be 
        that value, otherwise an identifier will be generated. There is not setter for this property.
        """
        return self._species

    @property
    def skif_species(self):
        """
        str : The species identifier that is passed to the sasktran engine.  In most cases this is identical to the
        property species, however in some specialized cases this can be different.  skif_species returns the species
        handled used by the climatology returned by climatology.skif_object(), while species returns the species handle
        used by the climatology object.
        """
        return self._species

    @property
    def climatology(self):
        """
        sasktran.Climatology : The climatology associated with this species. There is not setter for this property.
        """
        return self._climatology

    @property
    def optical_property(self):
        """
        sasktran.OpticalProperty : The optical property associated with this species. There is not setter for this 
        property.
        """
        return self._optical_property

    @staticmethod
    def is_sasktran_supported_handle(handle: str):
        return '_CM3' in handle.upper()


class SpeciesAerosol(Species):
    """    
    Special climatology which handles aerosol distributions. Contains a particle size Climatology and
    MieAerosol OpticalProperty internally to handle changes to size parameters.

    Parameters
    ----------
    altitudes: np.ndarray
        Array of altitudes in meters that the profile is to be specified at.  Same shape as values.
    aerosol_values: dict
        Dictionary of species names and values containing the aerosol number density Values should be the same size
        as altitudes.
    particle_size_values: dict
        Dictionary of species names and values containing the particle size profile. Values should be the same size
        as altitudes.
    species: str (optional, default='H2SO4')
        Molecule to use, one of ['H2SO4', 'ICE', 'WATER']
    interp: str, optional
        One of 'linear' or 'log'. Defines the interpolation space. Default is 'linear'
    spline: bool, optional
        One of True of False.  If True then a quadratic spline will be used to interpolate, if False then piecewise
        linear interpolation is done.  Default is False.

    Examples
    --------
    >>> from sasktran.species import SpeciesAerosol
    >>> import numpy as np
    >>> altitudes = np.linspace(500, 99500, 100)
    >>> values = np.ones_like(altitudes)
    >>> species = SpeciesAerosol(altitudes, {'SKCLIMATOLOGY_AEROSOL_CM3': values},\
                                            {'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': values * 0.08,\
                                             'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': values * 1.6}, 'H2SO4')
    >>> species.get_parameter('SKCLIMATOLOGY_AEROSOL_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([ 1.,  1.])
    >>> species['SKCLIMATOLOGY_AEROSOL_CM3'] *= 2
    >>> species.get_parameter('SKCLIMATOLOGY_AEROSOL_CM3', latitude=0, longitude=0, altitudes=[10000, 11000], mjd=54372)
    array([ 2.,  2.])
    >>> species.get_parameter('SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', latitude=0, longitude=0,\
                              altitudes=[10000, 11000], mjd=54372)
    array([ 0.08,  0.08])
    >>> species['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'] *= 2
    >>> species.get_parameter('SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', latitude=0, longitude=0,\
                              altitudes=[10000, 11000], mjd=54372)
    array([ 0.16,  0.16])
    """
    def __init__(self, altitudes: np.ndarray, aerosol_values: dict, particle_size_values: dict,
                 species: str='H2SO4', interp: str='linear', spline: bool=False):
        from sasktran import MieAerosol

        self._particlesize_climatology = ClimatologyUserDefined(altitudes, particle_size_values, interp, spline)
        self._optical_property = MieAerosol(self._particlesize_climatology, species)
        self._climatology = ClimatologyUserDefined(altitudes, aerosol_values, interp, spline)
        self._aerosol_guids = aerosol_values.keys()

        if len(self._aerosol_guids) > 1:
            raise ValueError('Error, trying to add more than one aerosol profile to SpeciesAerosol')

        super().__init__(self._optical_property, self._climatology, list(self._aerosol_guids)[0])

    def get_parameter(self, species: str, latitude: float, longitude: float, altitudes: np.ndarray, mjd: float)\
            -> np.ndarray:
        if species in self._aerosol_guids:
            return self._climatology.get_parameter(species, latitude, longitude, altitudes, mjd)
        else:
            return self._particlesize_climatology.get_parameter(species, latitude, longitude, altitudes, mjd)

    def __getitem__(self, item):
        if item in self._aerosol_guids:
            return self._climatology[item]
        else:
            return self._particlesize_climatology[item]

    def __setitem__(self, item, value):
        if item in self._aerosol_guids:
            self._climatology[item] = value
        else:
            self._particlesize_climatology[item] = value
            self._optical_property.particlesize_climatology = self._particlesize_climatology


class SpeciesAerosolGloSSAC(SpeciesAerosol):
    """
    __init__(self, altitudes: np.ndarray=np.arange(0, 100001, 500), particle_size_values: dict=None, \
    species: str='H2SO4', interp: str='linear', spline: bool=False, extend=False)
    
    Special climatology which handles aerosol distributions. Similar to :py:class:`sasktran.SpeciesAerosol`, however 
    rather than a user defined aerosol profile the :py:class:`sasktran.GloSSAC` climatology is used.

    A dictionary of particle size values can be provided. If not, a lognormal distribution with median radius=80nm and
    width of 1.6 is used.

    Notes
    -----
    Although the user side uses the handle 'SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM', aerosol properties are
    passed to sasktran using the 'SKCLIMATOLOGY_AEROSOL_CM3' handle

    Examples
    --------
    >>> from sasktran.species import SpeciesAerosolGloSSAC
    >>> import numpy as np
    >>> altitudes = np.linspace(500, 99500, 100)
    >>> values = np.ones_like(altitudes)
    >>> species = SpeciesAerosolGloSSAC()
    >>> species.get_parameter('SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM', latitude=0, longitude=0, \
                              altitudes=[22000, 23000], mjd=54372)
    array([ 0.00049428,  0.0004945 ])
    >>> species.get_parameter('SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM', latitude=20, longitude=0, \
                              altitudes=[22000, 23000], mjd=54372)
    array([ 0.00048611,  0.0004257 ])
    >>> species.get_parameter('SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', latitude=0, longitude=0, \
                               altitudes=[10000, 11000], mjd=54372)
    array([ 0.08,  0.08])
    >>> species['SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS'] *= 2
    >>> species.get_parameter('SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', latitude=0, longitude=0, \
                               altitudes=[10000, 11000], mjd=54372)
    array([ 0.16,  0.16])
    """
    def __init__(self, altitudes: np.ndarray=np.arange(0, 100001, 500),
                 particle_size_values: dict=None, species: str='H2SO4', interp: str='linear', spline: bool=False,
                 extend=False):
        from .climatology import GloSSAC

        if particle_size_values is None:
            particle_size_values = {'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS': np.ones_like(altitudes) * 0.08,
                                    'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH': np.ones_like(altitudes) * 1.6}

        super().__init__(altitudes, {'SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM': np.ones_like(altitudes) * 0.0},
                         particle_size_values, species=species, interp=interp, spline=spline)

        self._climatology = GloSSAC(extend=extend)

    @property
    def skif_species(self):
        # Note that skif object for this species returns back a climatology using aerosol CM3, rather than
        # extinction
        return 'SKCLIMATOLOGY_AEROSOL_CM3'


class SpeciesBaumIceCloud(Species):
    """
    Creates a species of Baum ice crystals with a Gaussian cloud shape.

    Parameters
    ----------
    particlesize_microns : float
        Particle size used for the :py:class:`sasktran.BaumIceCrystal` optical property
    cloud_top_m : float
        Cloud top altitude in m.  Defined such that 95% of the cloud is below this height. [m]
    cloud_width_fwhm_m : float
        Full width at half maximum of the Gaussian cloud. [m]
    vertical_optical_depth : float
        Vertical optical depth of the cloud specified at vertical_optical_depth_wavel_nm
    vertical_optical_depth_wavel_nm : float
        Wavelength where the cloud has the optical depth specified by vertical_optical_depth. [nm]
    altitude_resolution_m : float, optional
        Internal resolution to calculate the cloud number density at.  Default 10 m.
    num_sigma : int, optional
        Number of standard deviations to include away from the center of the cloud.  Default 5.
    name : str, optional
        Name of the species, only necessary to change if you want a different name or want to add multiple ice cloud
        layers.
    """
    def __init__(self, particlesize_microns: float, cloud_top_m: float,
                 cloud_width_fwhm_m: float, vertical_optical_depth: float,
                 vertical_optical_depth_wavel_nm: float, altitude_resolution_m: float=10, num_sigma: int=5,
                 name='icecloud'):
        optical_property = sk.BaumIceCrystal(particlesize_microns)

        cloud_sigma = cloud_width_fwhm_m / (2 * np.sqrt(2 * np.log(2)))

        # Cloud top is defined as middle + 2 sigma
        cloud_middle_m = cloud_top_m - 2.0 * cloud_sigma
        cloud_heights = np.arange(cloud_middle_m - num_sigma * cloud_sigma,
                                  cloud_middle_m + num_sigma * cloud_sigma + altitude_resolution_m,
                                  altitude_resolution_m)

        # Unnormalized gaussian since we will normalize to vertical optical depth anyways
        cloud_numden = np.exp(-(cloud_heights - cloud_middle_m)**2 / (2 * cloud_sigma**2))

        xs_at_wavel = optical_property.calculate_cross_sections(sk.MSIS90(), 0, 0, 0, 54372,
                                                                vertical_optical_depth_wavel_nm).total

        current_vert_optical_depth = np.sum(cloud_numden * xs_at_wavel) * altitude_resolution_m * 100  # convert /cm to /m

        cloud_numden *= vertical_optical_depth / current_vert_optical_depth

        climatology = sk.ClimatologyUserDefined(cloud_heights, {name: cloud_numden},
                                                out_of_bounds_value=0.0)

        super().__init__(optical_property, climatology)
        self._cloud_heights = cloud_heights

    @property
    def altitude_grid(self):
        return self._cloud_heights
