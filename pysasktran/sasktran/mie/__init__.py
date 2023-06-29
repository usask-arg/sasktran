import sasktranif.sasktranif as skif
import numpy as np
from sasktran.exceptions import wrap_skif_functionfail
import abc


class Mie(abc.ABC):
    def __init__(self, name: str):
        """
        An interface to the Mie code module used internally within SASKTRAN.  This is a low level interface directly
        to the Mie code which returns back low level parameters (scattering efficiencies, S1, S2, etc.) for a single
        particle size.


        Parameters
        ----------
        name : str
            Internal Mie code to use, currently only 'wiscombe' is supported
        """
        self._mie = skif.ISKMie(name.upper())
        self._name = name

    @wrap_skif_functionfail
    def calculate(self, wavelength: float,
                  radius: float,
                  refractive_index):
        """
        Sets the parameters for the Mie calculation.  No calculation is actually performed until one of the accessor
        functions is called, this merely sets the parameters.  The output of the Mie code is accessed through
        separate outputs.

        Parameters
        ----------
        wavelength : float
            Wavelength in the same units as radius
        radius : float
            Particle radius in the same units as wavelength
        refractive_index : complex
            Complex refractive index, the convention that the imaginary part of the refractive index is negative is
            used here
        """
        try:
            import xarray as xr
        except ImportError:
            raise EnvironmentError('xarray is required to use Mie.calculate()')

        self._mie.Calculate(wavelength, radius, np.real(refractive_index), np.imag(refractive_index))

    @wrap_skif_functionfail
    def Qext(self):
        """
        Returns
        -------
        float
            The Mie extinction efficiency parameter
        """
        return self._mie.Qext()

    @wrap_skif_functionfail
    def Qsca(self):
        """
        Returns
        -------
        float
            The Mie scattering efficiency parameter
        """
        return self._mie.Qsca()

    @wrap_skif_functionfail
    def Qabs(self):
        """
        Returns
        -------
        float
            The Mie absorption efficiency parameter
        """
        return self._mie.Qabs()

    @wrap_skif_functionfail
    def Cext(self):
        """
        Returns
        -------
        float
            The total extinction in units of radius**2
        """
        return self._mie.Cext()

    @wrap_skif_functionfail
    def Cabs(self):
        """
        Returns
        -------
        float
            The absorption extinction in units of radius**2
        """
        return self._mie.Cabs()

    @wrap_skif_functionfail
    def Csca(self):
        """
        Returns
        -------
        float
            The scattering extinction in units of radius**2
        """
        return self._mie.Csca()

    @wrap_skif_functionfail
    def S1(self):
        """
        Returns
        -------
        np.array
            S1 wave at self.scattering_angles
        """
        s = self._mie.S1()

        return s[0] + 1j * s[1]

    @wrap_skif_functionfail
    def S2(self):
        """
        Returns
        -------
        np.array
            S2 wave at self.scattering_angles
        """
        s = self._mie.S2()

        return s[0] + 1j * s[1]

    @wrap_skif_functionfail
    def Pmom(self):
        """
        Returns
        -------
        np.array
            Phase moment parameters
        """
        return self._mie.PMom()

    @wrap_skif_functionfail
    def Sforward(self):
        """
        Returns
        -------
        complex
            Forward peak of S1
        """
        temp = self._mie.Sforward()
        return temp[0] + 1j * temp[1]

    @wrap_skif_functionfail
    def SBackward(self):
        """
        Returns
        -------
        complex
            Backward peak of S1
        """
        temp = self._mie.SBackward()
        return temp[0] + 1j * temp[1]

    @wrap_skif_functionfail
    def spike(self):
        """
        Returns
        -------
        float
            Mie spike amplitude
        """
        return self._mie.Spike()

    @wrap_skif_functionfail
    def TForward(self, i: int):
        """
        Parameters
        ----------
        i : int
            1 or 2

        Returns
        -------
        complex
        """
        temp = self._mie.TForward(i)
        return temp[0] + 1j * temp[1]

    @wrap_skif_functionfail
    def TBackward(self, i: int):
        """
        Parameters
        ----------
        i : int
            1 or 2

        Returns
        -------
        complex
        """
        temp = self._mie.TBackward(i)
        return temp[0] + 1j * temp[1]

    @abc.abstractmethod
    def scattering_angles(self):
        """
        Returns
        -------
        np.array
            Scattering angles in degrees
        """
        pass

    @abc.abstractmethod
    def max_legendre_moment(self):
        """
        Returns
        -------
        int
            Maximum number of legendre moments computed
        """
        pass

    @abc.abstractmethod
    def is_stable(self, wavelength, radius):
        """
        Parameters
        ----------
        wavelength : np.array like
            Wavelength in the same units as radius
        radius : np.array like
            Radius in the same units as wavelength

        Returns
        -------
        np.array
            Boolean array hinting if the Mie call is expected to be stable for these parameters. This is not a guarantee,
            only a hint.
        """
        pass


class MieWiscombe(Mie):
    def __init__(self, angles: np.array = None, max_legendre_moment: int = 64):
        """
        Implementation of the Mie scattering algorithms by Wiscombe

        Parameters
        ----------
        angles : np.array, optional
            Scattering angles in degrees.  Default is 0.1 degree resolution 5 degrees from the forward/backward
            scatter peaks, and 1 degree resolution elsewhere
        max_legendre_moment: int, optional
            Maximum legendre moment to compute. Default 64
        """
        super().__init__('WISCOMBE')

        if angles is not None:
            self._angles = angles
        else:
            self._angles = np.concatenate((np.arange(0, 5, 0.1), np.arange(5, 175, 1), np.arange(175, 180, 0.1), [180]))
        self._max_legendre_moment = max_legendre_moment

        self._mie.SetProperty('scatteringangles', self._angles)
        self._mie.SetProperty('maxlegendre', self._max_legendre_moment)

    def scattering_angles(self):
        """
        Returns
        -------
        np.array
            Scattering angles in degrees
        """
        return self._angles

    def max_legendre_moment(self):
        """
        Returns
        -------
        int
            Maximum number of legendre moments computed
        """
        return self._max_legendre_moment

    def is_stable(self, wavelength, radius):
        """
        Parameters
        ----------
        wavelength : np.array like
            Wavelength in the same units as radius
        radius : np.array like
            Radius in the same units as wavelength

        Returns
        -------
        np.array
            Boolean array hinting if the Mie call is expected to be stable for these parameters. This is not a guarantee,
            only a hint.
        """
        return (2 * np.pi * radius) / wavelength < 20000

    @wrap_skif_functionfail
    def Pmom(self):
        """
        Returns
        -------
        np.array
            Phase moment parameters
        """
        # Wiscombe code allocates extra memory for Pmom
        return self._mie.PMom()[:, :self._max_legendre_moment]
