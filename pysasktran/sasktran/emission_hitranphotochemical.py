from typing import Tuple
from copy import copy
import numpy as np
import sasktranif.sasktranif as skif
from .exceptions import wrap_skif_functionfail
from .climatology import Climatology
from .emission import Emission
from .opticalproperty import HITRANChemical


# ------------------------------------------------------------------------------
#           Class HITRANPhotoChemical
# ------------------------------------------------------------------------------
class HITRANPhotoChemical(Emission):
    """
    Calculates the photo-chemical emission from an excited molecular upper state. This work was originally developed
    to model the molecular oxygen atmospheric emissions due to A-band and singlet delta etc. but can also be used for OH and other
    molecular bands in the hitran database. This class wraps the SasktranIF Emission object, ``HITRAN_PHOTOCHEMICAL`` and
    users are encouraged to refer to the available documentation within that package.

    The HITRANPhotoChemical object calculates the photo-chemical emission of a molecular band (assumed vibrational) using the Einstein-A coefficients
    given in the Hitran database. The code calculates the spectral emission of a single electron statistically distributed
    across the vibrational band. All spectral lines are modelled using a Voigt profile. The user must explicitly select
    an appropriate wavenumber range, and upper and lower quanta so only lines from the target vibrational band are selected
    from the Hitran database. The user can specify the climatology used to calculate pressure and temperature used in the Voigt
    profile claculation. By default the MSIS90 model is used.

    The user must provide a climatology object that returns the number of electrons in the excited upper state. The climatology object
    can be as simple as a user-defined profile or as complex as a photochemical model. The HITRANPhotoChemical object thermalizes the
    distribution of these electrons across the upper band using a Boltzmann distribution.

    The HITRANPhotoChemical object can account for self-broadening of the line-shapes but requires the user to provide
    a climatology object for the number density of the molecule of interest. If no climatology is provided then self-braodening
    is ignored.

    Parameters
    ----------
    chemical_name : str
        Chemical abbreviation of the molecule of interest. Typically O2, OH etc.

    wavenumber_range: Tuple[float,float]
        A two element tuple specifying the range of wavenumbers in cm-1 used to select spectral lines from the Hitran database.
        The selection range should span the entire wavenumber range of the target vibrational band. The range is used
        along with the ``lower_state_global_quant_filter`` and ``upper_state_global_quanta_filter`` to select all the lines
        associated with the target vibrational band.  The first element should be the lower wavenumber. The second element
        should be the upper wavenumber. Both are in cm-1.

    lower_state_global_quanta_filter: str, optional
        The quanta filter used to select spectral lines in the given wavenumber that match the lower quantum state.
        These quanta are given in the Hitran data records as the V′ and V′′ fields and are strings similar to “X 0”.
        The code will only accept spectral lines that match the lower state filter.
        The code ignores spaces but is case sensitive and order sensitive. i.e ‘X 0’ is not the same as ‘0 X’. This
        option requires the user to be familiar with the Hitran 160 character records and is comfortable reading the
        Hitran files to ascertain the notation used in the database files. The default is an empty string,"", which
        disables the filter and matches all spectral lines regardless of lower state quanta, this is usually not the
        desired behaviour and the string must be set.

    upper_state_global_quanta_filter: str, optional
        Same as lower_global_quanta_filter except for the upper state. The default is the empty string, "" which is
        not correct for isolating upper state vibrational levels.

    atmospheric_state_climatology(object climatology)
        The climatology object used to calculate pressure and temperature for the Voigt profile calculations.

    excited_upper_state_climatology: Tuple[ object, str]
        A two element tuple providing the climatology object to be used for calculation of the total number of excited
        upper state molecules in molecules/cm3. The first element is the :class:`Climatology` object that  captures the
        photo-chemistry of the atmosphere. It can be as simple as a user defined height profile or a more complicted
        photo-chemical model. The climatology object should provide the number of excited molecules consistent
        with the selected spectral lines, i.e. they are referring to the same vibrational band etc. The second element
        is the climatology handle used by the excited upper-state climatology to acquire the number of excited molecules
        per cm3. A typical value for the oxygen molecule is ‘SKEMISSION_PHOTOCHEMICAL_O2’. A full list of valid
        climatology handles is available in the SasktranIF documentation.

    atmospheric_state_climatology(object climatology)
        The climatology object used to calculate pressure and temperature for the Voigt profile calculations.

    self_broadening_climatology: Tuple[ object, str]
        Optional. A two element tuple providing the climatology object to be used for self-broadening calculations.
        The first element is a :class:`~Climatology` object that must provide the number density of the target molecule
        in molecules/cm3. The second element is the climatology handle used by the self-broadening climatology to acquire
        the number of molecules/cm3 of the target molecule. A typical value for the oxygen molecule is
        ‘SKCLIMATOLOGY_O2_CM3’.  A full list of climatology handles is available in the SasktranIF documentation. This
        climatology is used along with the pressure and temperature values from the atmospheric_state_climatology to
        calculate the partial pressure of the molecule which is used in the self broadening calculation. If no
        climatology is provided then self broadening is ignored. For reference, ignoring this item leads to an error on
        the order of 1% for O2 in the A-band.

    isotope_filter : int, optional
        Optional Allows the HITRAN object to load in just one isotope of the requested molecule. The value set must match one of
        the isotope labels used for the given molecule in the HITRAN database file, molparam.txt. Note that the code
        does not adjust the line strength but uses the line strength value as written in the HITRAN database. This
        means you may have to account for and/or remove the abundance automatically built into the HITRAN database line
        strength values. A value of 0 will load all isotopes in the Hitran database. Default is 0. Most users will not
        change this value.

    use_cache: bool, optional
        Optional. Species that this emission object will create an internal cache of signals in wavelength and location.
        This feature is provided for RT engines, such as HR, SO, Disco and MC, that internally need to calculate the Voigt profile at
        many locations for one wavelength while the Voigt line profile code is much more efficient if performed
        across all wavenumbers at each location. This feature, if active, is only used inside method skif_object when the
        input arguments have an 'engine' keyword. The default is True.
    """
    @wrap_skif_functionfail
    def __init__(self,
                 chemical_name: str,
                 wavenumber_range: Tuple[float, float],
                 atmospheric_state_climatology: Climatology,
                 excited_upper_state_climatology: Tuple[Climatology, str],
                 lower_state_global_quanta_filter: str,
                 upper_state_global_quanta_filter: str,
                 self_broadening_climatology: Tuple[Climatology, str] = (None, ''),
                 isotope_filter: int=0,
                 use_cache: bool = True):

        HITRANChemical._validate_registry()
        super().__init__("HITRAN_PHOTOCHEMICAL")
        self._chemical_name: str = str(chemical_name)
        self._wavenumber_range: Tuple[float, float] = wavenumber_range
        self._excited_upper_state_climatology: Tuple[Climatology, str] = excited_upper_state_climatology
        self._lower_state_global_quanta_filter: str = str(lower_state_global_quanta_filter)           # string describing the lowe state global quanta
        self._upper_state_global_quanta_filter: str = str(upper_state_global_quanta_filter)
        self._self_broadening_climatology: Tuple[Climatology, str] = self_broadening_climatology
        self._atmospheric_state_climatology: Climatology = atmospheric_state_climatology
        self._isotope_filter: int = int(isotope_filter)
        self._use_cache: bool = use_cache
        self._cached_wavelengths: np.ndarray = None
        self._update_properties()

    # ------------------------------------------------------------------------------
    #           __setstate__
    # ------------------------------------------------------------------------------

    def __setstate__(self, state):
        super().__setstate__(state)
        self._update_properties()

    # ------------------------------------------------------------------------------
    #           _update_opticalproperty
    # ------------------------------------------------------------------------------

    def _update_properties(self):
        self._iskemission.SetProperty('set_chemical_name', self._chemical_name)
        self._iskemission.SetProperty('set_wavenumber_range', self._wavenumber_range)
        self._iskemission.SetProperty('set_lower_state_global_quanta_filter', self._lower_state_global_quanta_filter)
        self._iskemission.SetProperty('set_upper_state_global_quanta_filter', self._upper_state_global_quanta_filter)
        self._iskemission.SetProperty('set_excited_upper_state_climatology', self._excited_upper_state_climatology[0].skif_object())
        self._iskemission.SetProperty('set_excited_upper_state_climatology_handle', self._excited_upper_state_climatology[1])
        self._iskemission.SetProperty('setisotopefilter', self._isotope_filter)

        if self._self_broadening_climatology is not None:
            self._iskemission.SetProperty('set_self_broadening_climatology', self._self_broadening_climatology[0].skif_object())
            self._iskemission.SetProperty('set_self_broadening_climatology_handle', self._self_broadening_climatology[1])

        if self._atmospheric_state_climatology:
            self._iskemission.SetProperty('set_atmospheric_state_climatology', self._atmospheric_state_climatology.skif_object())

        if self._use_cache and self._cached_wavelengths is not None:
            self._iskemission.SetProperty('enable_cached_emissions', 1e7 / self._cached_wavelengths)

    # ------------------------------------------------------------------------------
    #           _reset_skif_object
    # ------------------------------------------------------------------------------

    def _reset_skif_object(self):
        self._iskemission = skif.ISKEmission("HITRAN_PHOTOCHEMICAL")
        self._update_properties()

    # ------------------------------------------------------------------------------
    #           skif_object
    # ------------------------------------------------------------------------------

    @wrap_skif_functionfail
    def skif_object(self, **kwargs):

        if 'engine' in kwargs:                                                          # Check to see if we are being called by an engine prior
            self._reset_skif_object()
            if self._use_cache:                                                        # If we are not using the cache, just check if the wavenumber range has changed and if so reset the object
                wavelengths = np.asarray(kwargs['engine'].wavelengths)                  # If we are then fetch the wavelengths the engine plans to use
                self._cached_wavelengths = copy(wavelengths)
                self._reset_skif_object()
        return self._iskemission


# ------------------------------------------------------------------------------
#           HITRANPhotoChemical_O2_ABand
# ------------------------------------------------------------------------------
class HITRANPhotoChemical_O2_ABand(HITRANPhotoChemical):
    """
    Calculates the photo-chemical emission from the O2 A-Band around 762 nm. This is a small wrapper around the
    :class:`~HitranPhotoChemical` class. The O2 A-band spectral lines cover the wavelength range 748 nnm to 780 nm,  correspondiong
    to 12820.0 to 13369.0 cm-1.

    Parameters
    ----------
    atmospheric_state_climatology(object climatology)
        The climatology object used to calculate pressure and temperature for the Voigt profile calculations.

    excited_upper_state_climatology: Tuple[ object, str]
        A two element tuple providing the climatology object to be used for calculation of the total number of O2 A-band
        molecules in molecules/cm3. The first element is the :class:`Climatology` object that  captures the
        photo-chemistry of the atmosphere. It can be as simple as a user defined height profile or a more complicated
        photo-chemical model. The climatology object should provide the number of excited molecules in the oxygen A-band
        The second element is the climatology handle used by the excited upper-state climatology to acquire the number of excited molecules
        per cm3. A typical value for the oxygen molecule is ‘SKEMISSION_PHOTOCHEMICAL_O2’. A full list of valid
        climatology handles is available in the SasktranIF documentation.

    atmospheric_state_climatology(object climatology)
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    self_broadening_climatology: Tuple[ object, str]
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    isotope_filter : int, optional
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    use_cache: boo, optional
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    Examples
    --------
    >>>    import sasktran as sk
    >>>    import numpy as np
    >>>
    >>>    msis = sk.MSIS90()
    >>>    constant = sk.ClimatologyUserDefined(np.array([0.0, 100000.0]), {'SKEMISSION_PHOTOCHEMICAL_O2': [1.0E6, 1.0E6]})
    >>>    emission = sk.HITRANPhotoChemical_O2_ABand(msis, (constant, 'SKEMISSION_PHOTOCHEMICAL_O2'), self_broadening_climatology=(msis, 'SKCLIMATOLOGY_O2_CM3'))
    >>>    wavelen  = np.arange(749.0, 779.0, 0.0001)
    >>>    signal = emission.isotropic_emission(52.0, -106.0, 80000.0, 57005.8, wavelen, False)

    """

    @wrap_skif_functionfail
    def __init__(self,
                 excited_upper_state_climatology: Tuple[Climatology, str],
                 atmospheric_state_climatology: Climatology,
                 self_broadening_climatology: Tuple[Climatology, str] = (None, ''),
                 isotope_filter: int=0,
                 use_cache: bool = True):

        super().__init__('O2',
                         (12820.0, 13369.0),
                         atmospheric_state_climatology,
                         excited_upper_state_climatology,
                         'X 0', 'b 0',
                         self_broadening_climatology=self_broadening_climatology,
                         isotope_filter=isotope_filter,
                         use_cache=use_cache)


# ------------------------------------------------------------------------------
#           class HITRANPhotoChemical_O2_SingletDelta( HITRANPhotoChemical ):
# ------------------------------------------------------------------------------
class HITRANPhotoChemical_O2_SingletDelta(HITRANPhotoChemical):
    """
    Calculates the photo-chemical emission from the O2 Singlet Delta around 1.27 microns. This is a small wrapper around the
    :class:`~HitranPhotoChemical` class. The O2 A-band spectral lines cover the wavelength range 1.223 microns to 1.321 microns,
    correspondiong to 7570 to 8172 cm-1.

    Parameters
    ----------
    atmospheric_state_climatology(object climatology)
        The climatology object used to calculate pressure and temperature for the Voigt profile calculations.

    excited_upper_state_climatology: Tuple[ object, str]
        A two element tuple providing the climatology object to be used for calculation of the total number of O2 A-band
        molecules in molecules/cm3. The first element is the :class:`Climatology` object that  captures the
        photo-chemistry of the atmosphere. It can be as simple as a user defined height profile or a more complicated
        photo-chemical model. The climatology object should provide the number of excited molecules in the oxygen A-band
        The second element is the climatology handle used by the excited upper-state climatology to acquire the number of excited molecules
        per cm3. A typical value for the oxygen molecule is ‘SKEMISSION_PHOTOCHEMICAL_O2’. A full list of valid
        climatology handles is available in the SasktranIF documentation.

    self_broadening_climatology: Tuple[ object, str]
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    isotope_filter : int, optional
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    use_cache: boo, optional
        Optional. See :class:`~HitranPhotoChemical` for a description of this parameter.

    Examples
    --------
    >>>    import sasktran as sk
    >>>    import numpy as np
    >>>
    >>>    msis = sk.MSIS90()
    >>>    constant = sk.ClimatologyUserDefined(np.array([0.0, 100000.0]), {'SKEMISSION_PHOTOCHEMICAL_O2': [1.0E6, 1.0E6]})
    >>>    emission = sk.HITRANPhotoChemical_O2_SingletDelta( msis, (constant, 'SKEMISSION_PHOTOCHEMICAL_O2'), self_broadening_climatology=(msis, 'SKCLIMATOLOGY_O2_CM3'))
    >>>    wavelen  = np.arange(1223.0, 1321.0, 0.0001)
    >>>    signal = emission.isotropic_emission(52.0, -106.0, 80000.0, 57005.8, wavelen, False)

    """

    @wrap_skif_functionfail
    def __init__(self,
                 atmospheric_state_climatology: Climatology,
                 excited_upper_state_climatology: Tuple[Climatology, str],
                 self_broadening_climatology: Tuple[Climatology, str] = (None, ''),
                 isotope_filter: int=0,
                 use_cache: bool = True):

        super().__init__('O2',
                         (7570.0, 8172.0),
                         atmospheric_state_climatology,
                         excited_upper_state_climatology,
                         'X 0', 'a 0',
                         self_broadening_climatology=self_broadening_climatology,
                         isotope_filter=isotope_filter,
                         use_cache=use_cache)
