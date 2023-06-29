import sasktran as sk
from sasktran.climatology import Climatology
from sasktran.opticalproperty import OpticalProperty
from sasktran.util import DictWithCallback
import logging


class Atmosphere(object):
    """
    The Atmosphere class is a generic container that defines the atmospheric state for a radiative transfer calculation.
    The first component of this is defining a set of atmospheric constituents, which is controlled by the 
    :py:attr:`sasktran.Atmosphere.species` property.

    The second component is defining the ground reflectance. Sasktran support :ref:`BRDF` surfaces. The 
    :ref:`BRDF` can be set with the :py:attr:`sasktran.Atmosphere.brdf` property.

    Emissions can also be added to the atmosphere using the :py:attr:`sasktran.Atmosphere.brdf` property.  Note that
    currently the only radiative transfer model that supports emissions is the HR engine.
    
    Additionally, the atmosphere object can be used to set the temperature and pressure :ref:`climatology` (see 
    :py:attr:`sasktran.Atmosphere.atmospheric_state`), as well as species to calculate weighting functions for (see
    :py:attr:`sasktran.Atmosphere.wf_species`).
    """
    def __init__(self):
        self._species = DictWithCallback()
        self._species.add_modified_callback(self._set_has_changed)
        self._emissions = DictWithCallback()
        self._emissions.add_modified_callback(self._set_has_changed)
        self._atmospheric_state = None
        self._brdf = sk.brdf.Lambertian(0.3)
        self._wf_species = None
        self._has_changed = False

    @property
    def has_changed(self):
        """
        bool: True if the atmosphere has changed since being added to an :py:class:`sasktran.Engine` object.
        """
        return self._has_changed

    def _set_has_changed(self):
        """
        Sets the :py:attr:`sasktran.Atmosphere.has_changed` property to True.  Called whenever an internal member
        of this class is modified.
        """
        self._has_changed = True

    def _added_to_engine(self):
        """
        Sets the :py:attr:`sasktran.Atmosphere.has_changed` property to False.  Called when the atmosphere is added
        to an engine.
        """
        self._has_changed = False

    def add_species(self, species: str, climatology: Climatology, optical_property: OpticalProperty, name: str=None):
        """
        .. note:: Deprecated in Sasktran 0.1. The :py:attr:`sasktran.Atmosphere.species` property should be used 
            directly instead.
        
        Adds a species to the atmosphere container

        Parameters
        ----------
        species : str
            GUID Identifier of the species
        climatology : Climatology
            Climatology object for the species
        optical_property : OpticalProperty
            Optical property object for the species
        name : str, optional
            Internal name used in the atmosphere to refer to the species.  If set to None (default) then the species
            string is used as the name string

        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> atmosphere.add_species('SKCLIMATOLOGY_O3_CM3', sk.Climatology('o3labow'),\
                                   sk.OpticalProperty('o3_osirisres'))
        >>> print('SKCLIMATOLOGY_O3_CM3' in atmosphere.species)
        True

        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> atmosphere.add_species('SKCLIMATOLOGY_O3_CM3', sk.Climatology('o3labow'),\
                                   sk.OpticalProperty('o3_osirisres'), name='o3')
        >>> print(atmosphere['o3'])
        SKCLIMATOLOGY_O3_CM3
        SasktranIF Climatology: o3labow
        Supported Species: ['SKCLIMATOLOGY_O3_CM3', 'SKCLIMATOLOGY_O3_VMR']
        SasktranIF Optical Property: o3_osirisres
        """
        logging.warning('Deprecation Warning: Atmosphere.add_species will be removed in the future, use Atmosphere[name] = sk.Species()'
                        ' instead')
        if name is None:
            name = species
        self._species[name] = sk.Species(optical_property, climatology, species)
        self._set_has_changed()

    def remove_species(self, name: str):
        """
        Removes a species with a given identifier from the atmosphere.

        Parameters
        ----------
        name : str
            Identifier of the species to remove.

        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> atmosphere.species['o3'] = sk.Species(sk.O3DBM(), sk.Labow())            
        >>> print(atmosphere.species['o3'])             
        SKCLIMATOLOGY_O3_CM3
        SasktranIF Climatology: O3LABOW
        Supported Species: ['SKCLIMATOLOGY_O3_CM3', 'SKCLIMATOLOGY_O3_VMR']
        SasktranIF Optical Property: O3_DBM
        >>> atmosphere.remove_species('o3')             # remove the species
        >>> print('o3' in atmosphere.species)
        False
        """
        del self._species[name]
        self._set_has_changed()

    @property
    def species(self):
        """
        dict: A dictionary of :ref:`species` objects with the keys being the species identifying strings. 

        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> ozone = sk.Species(sk.O3DBM(), sk.Labow())  # create a species
        >>> atmosphere.species['o3'] = ozone            # add the species
        >>> print(atmosphere.species['o3'])             # get the species
        SKCLIMATOLOGY_O3_CM3
        SasktranIF Climatology: O3LABOW
        Supported Species: ['SKCLIMATOLOGY_O3_CM3', 'SKCLIMATOLOGY_O3_VMR']
        SasktranIF Optical Property: O3_DBM

        Notes
        -----
        The :py:meth:`sasktran.Atmosphere.__getitem___` and :py:meth:`sasktran.Atmosphere.__setitem___` are
        aliases to this dictionary. For example this means that
        
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> atmosphere.species['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
        
        is equivalent to 
        
        >>> atmosphere['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
        """
        return self._species

    @property
    def emissions(self):
        """
        dict: A dictionary of :ref:`Emission` objects with keys being the identifying emission strings

        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> emission = sk.EmissionThermal()
        >>> atmosphere.emissions['thermal'] = emission
        >>> print(atmosphere.emissions['thermal'])
        Emission of type: THERMAL
        """
        return self._emissions

    @property
    def atmospheric_state(self):
        """ 
        sasktran.Climatology: A special :ref:`climatology` used to set the atmosphere's temperature, and pressure. This 
        climatology must support both 'SKCLIMATOLOGY_TEMPERATURE_K' and 'SKCLIMATOLOGY_PRESSURE_PA'. The default 
        atmospheric state is :py:class:`sasktran.MSIS90`.
        
        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> msis90 = sk.MSIS90()   # MSIS90 is a climatology that supports the temperature and pressure identifiers
        >>> print('SKCLIMATOLOGY_TEMPERATURE_K' in msis90.supported_species())
        True
        >>> print('SKCLIMATOLOGY_PRESSURE_PA' in msis90.supported_species())
        True
        >>> atmosphere.atmospheric_state = msis90
        >>> print(atmosphere)
        sasktran.Atmosphere:
        Species: dict_keys([])
        BRDF: BRDF Object: lambertian
        Atmospheric State: SasktranIF Climatology: MSIS90
        Supported Species: ['SKCLIMATOLOGY_PRESSURE_PA', 'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 'SKCLIMATOLOGY_TEMPERATURE_K']
        WF Species: None
        """
        return self._atmospheric_state

    @atmospheric_state.setter
    def atmospheric_state(self, value: Climatology):
        self._atmospheric_state = value
        self._set_has_changed()

    @property
    def brdf(self):
        """
        sk.BRDF : The surface :ref:`brdf` object. If set to a scalar, a Lambertian albedo is assumed. 

        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> atmosphere.brdf = 1                 # Equivalent to sk.Lambertian(1.0)
        >>> print(atmosphere.brdf)
        BRDF Object: lambertian
        >>> atmosphere.brdf = sk.Kokhanovsky()  # Non Lambertian surface
        >>> print(atmosphere.brdf)
        BRDF Object: SNOW_KOKHANOVSKY2012
        """
        return self._brdf

    @brdf.setter
    def brdf(self, value):
        if isinstance(value, sk.BRDF):
            self._brdf = value
        else:
            self._brdf = sk.brdf.Lambertian(value)
        self._set_has_changed()

    @property
    def wf_species(self):
        """
        str : :ref:`species` identifier to calculate weighing functions for. Weighting functions will be calculated if 
        this property is set. Note that engines that support multi-species weighting functions can handle this 
        property being a list of species identifiers. If this property is None, weighting functions won't be calculated.
        
        Examples
        --------
        >>> import sasktran as sk
        >>> atmosphere = sk.Atmosphere()
        >>> atmosphere.species['o3'] = sk.Species(sk.O3DBM(), sk.Labow())
        >>> print(atmosphere)
        sasktran.Atmosphere:
        Species: dict_keys(['o3'])
        BRDF: BRDF Object: lambertian
        Atmospheric State: None
        WF Species: None
        >>> atmosphere.wf_species = 'o3'
        >>> print(atmosphere)
        sasktran.Atmosphere:
        Species: dict_keys(['o3'])
        BRDF: BRDF Object: lambertian
        Atmospheric State: None
        WF Species: o3
        """
        return self._wf_species

    @wf_species.setter
    def wf_species(self, value):
        self._wf_species = value
        self._set_has_changed()

    def __getitem__(self, item):
        return self._species.get(item)

    def __setitem__(self, item, value):
        self._species[item] = value
        self._set_has_changed()

    def __repr__(self):
        return "sasktran.Atmosphere:\nSpecies: {}\nBRDF: {}\nAtmospheric State: {}\nWF Species: {}".format(
            self.species.keys(), self.brdf, self.atmospheric_state, self.wf_species
        )
