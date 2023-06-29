import sasktran as sk
from sasktran.exceptions import wrap_skif_functionfail
import sasktranif.sasktranif as skif
from abc import ABCMeta, abstractmethod
from typing import Tuple
import numpy as np
from sasktran.util import to_iter
from .climatology import Climatology


class EmissionBase(object):
    """
    Defines the interface required to create an emission that can be used with a radiative transfer model.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def skif_object(self, **kwargs) -> skif.ISKEmission:
        pass


class Emission(EmissionBase):
    """
    Class which represents the emission of a species in the atmosphere.  Most of the time it is not recommended to use
    this class and more specific classes, such as EmissionTable, should be used instead.

    Parameters
    ----------
    name : str
        Name handle of the emission to create.
    """
    @wrap_skif_functionfail
    def __init__(self, name: str):
        self._iskemission = skif.ISKEmission(name)
        self._name = name

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_iskemission'}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._iskemission = skif.ISKEmission(self._name)

    def __repr__(self):
        return "Emission of type: {}".format(self._name)

    @wrap_skif_functionfail
    def isotropic_emission(self, latitude: float, longitude: float,
                           altitude: float, mjd: float, wavelengths: np.array, isground=False):
        """
        Calculates the isotropic emission rate.

        Parameters
        ----------
        atmospheric_state : sasktran.Climatology
            Atmospheric state (sometime called background climatology, corresponds to
            :py:class:`sasktran.Atmosphere.atmospheric_state`). Typically the background climatology must
            support temperature and pressure, e.g. msis90.  Some optical properties may not need a background
            climatology, for example Mie aerosols, however one must still be passed in for the function to operate.
        latitude : float
            Latitude in degrees.
        longitude : float
            Longitude in degrees.
        altitude : float
            Altitude in meters.
        mjd : float
            Modified Julian Date.
        wavelengths : np.array
            Wavelengths in nm.
        isground : bool, Optional
            Set to true if this is a ground emission, default is False.

        Returns
        -------
        emission : np.ndarray
            Array of emissions in [photons/sec/nm/ster/cm2], same length at wavelengths
        """
        self._iskemission.UpdateLocation([latitude, longitude, altitude, mjd], isground)

        wavenumbers = np.asarray(1e7 / np.asarray(to_iter(wavelengths)))
        sort_index = np.argsort(wavenumbers)
        raw_result = self._iskemission.IsotropicEmission(wavenumbers[sort_index])

        # applying the sorted sorting indices restores the original order
        reverse_index = np.argsort(sort_index)
        emission = raw_result[1][reverse_index]

        return emission

    def skif_object(self, **kwargs) -> skif.ISKEmission:
        return self._iskemission

    def skif_species_override(self):
        # Some emission species require a specific GUID string to be used, if this is not None then that guid will
        # be forced
        return None


class EmissionTable(Emission):
    """
    Defines an emission in the atmosphere represented with a table of volume emission rates on an altitude/wavelength
    grid.

    Parameters
    ----------
    altitude_m : np.ndarray
        Array of altitudes in [m] defining the altitude dimension of the emission grid

    wavel_nm : np.ndarray
        Array of wavelengths in [nm] defining the wavelength dimension of the altitude grid

    volumeemission : np.ndarray
        Two dimensional array of volume emission rates in [photons/nm/sec/ster] with dimension (len(alts), len(wavel))

    Examples
    --------
    >>> import sasktran as sk
    >>> alts = np.arange(500, 100500, 1000)
    >>> wavel = np.linspace(300, 350, 50)
    >>> ver = np.ones((len(alts), len(wavel)))
    >>> emission = sk.EmissionTable(alts, wavel, ver)
    >>> print(emission)
    Emission of type: USERDEFINED_WAVELENGTHHEIGHT
    """
    def __init__(self, altitude_m: np.ndarray, wavel_nm: np.ndarray, volumeemission: np.ndarray):
        super().__init__('USERDEFINED_WAVELENGTHHEIGHT')

        self._altitude_m = altitude_m
        self._wavel_nm = wavel_nm
        self._volumeemission = volumeemission

        self._update_internal()

    def __setstate__(self, state):
        super().__setstate__(state)
        self._update_internal()

    @wrap_skif_functionfail
    def _update_internal(self):
        self._iskemission.SetProperty('Heights', self._altitude_m)
        self._iskemission.SetProperty('Wavelengths', self._wavel_nm)
        self._iskemission.SetProperty('EmissionTable', self._volumeemission)


class EmissionThermal(Emission):
    """
    Class used to enable thermal emissions within the HR radiative transfer model.  Currently only works as part
    of a full radiative transfer calculation.

    Examples
    --------
    >>> import sasktran as sk
    >>> emission = sk.EmissionThermal()
    >>> geometry = sk.VerticalImage()
    >>> geometry.from_sza_saa(sza=60, saa=60, lat=0, lon=0, tanalts_km=[10, 20, 30, 40], mjd=54372, locallook=0, \
        satalt_km=600, refalt_km=20)
    >>> atmosphere = sk.Atmosphere()
    >>> atmosphere['ozone'] = sk.Species(sk.O3OSIRISRes(), sk.Labow())
    >>> atmosphere['air'] = sk.Species(sk.Rayleigh(), sk.MSIS90())
    >>> atmosphere.emissions['thermal'] = emission
    >>> engine = sk.EngineHR(geometry, atmosphere, wavelengths=[1000])
    >>> engine.num_orders_of_scatter = 1
    >>> engine.include_emissions = True
    >>> engine.calculate_radiance()
    array([[ 0.01218222,  0.00270503,  0.00055809,  0.00013423]])
    """
    def __init__(self, emissivity=1.0):
        super().__init__('THERMAL')

        self._emissivity = emissivity

        self._update_internal()

    @wrap_skif_functionfail
    def _update_internal(self):
        self._iskemission.SetProperty('Emissivity', self._emissivity)

    @wrap_skif_functionfail
    def set_atmospheric_state(self, atmospheric_state_climatology: Climatology):
        self._iskemission.SetProperty('atmospheric_state', atmospheric_state_climatology.skif_object())

    @property
    def emissivity(self):
        return self._emissivity

    @emissivity.setter
    def emissivity(self, value):
        self._emissivity = value
        self._update_internal()

    def skif_species_override(self):
        return 'SKEMISSION_THERMAL'
