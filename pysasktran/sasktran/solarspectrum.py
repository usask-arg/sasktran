import sasktranif.sasktranif as skif
from sasktran.exceptions import wrap_skif_functionfail
import numpy as np


class SolarSpectrum(object):
    """
    Utility class for obtaining the solar spectrum

    Parameters
    ----------
    name: str, optional
        Solar spectrum to use, valid choices are 'sao2010' and 'FONTELA_UVIS_3MICRON'.  Default is 'sao2010'

    Examples
    --------
    >>> from sasktran import SolarSpectrum
    >>> spectrum = SolarSpectrum()
    >>> spectrum.irradiance([340, 600])
    array([  1.76034000e+14,   5.38339000e+14])
    """
    @wrap_skif_functionfail
    def __init__(self, name='SAO2010'):
        self._isksolarspectrum = skif.ISKSolarSpectrum(name.upper())
        self._name = name

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_isksolarspectrum'}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._isksolarspectrum = skif.ISKSolarSpectrum(self._name.upper())

    def skif_object(self):
        """
        Get the underlying skif object

        Returns
        -------
        skif.ISKSolarSpectrum

        Examples
        --------
        >>> from sasktran import SolarSpectrum
        >>> spectrum = SolarSpectrum()
        >>> type(spectrum.skif_object())
        <class 'sasktranif.sasktranif.ISKSolarSpectrum'>
        """
        return self._isksolarspectrum

    @wrap_skif_functionfail
    def irradiance(self, wavelengths, solardistance=None, mjd=None, fwhm=None, stdev=None):
        """
        Parameters
        ----------
        wavelengths : np.array
            Wavelengths to get the irradiance at in nm
        solardistance : float, optional
            Default is None.  If not None, the solardistance is set prior to getting the irradiance.  Units of AU.
        mjd: float, optional
            Modified Julain Data. Default is None.  If not None, the solardistance is set based upon the time prior to
            getting the irradiance.
        fwhm: float, optional
            Full width at half maximum resolution in nm to convolve the solar spectrum down to.  Default is None, using the
            native resolution.  Only one of fwhm or stdev should be set.  Can either be a scalar indicating every
            wavelength has the same fwhm, or an array the same size at wavelengths.
        stdev: float, optional
            stdev in nm to convolve the solar spectrum down to.  Default is none, using the native roslution.
            Only one of fwhm or stdev should be set.  Can either be a scalar indicating every wavelength has the same
            stdev, or an array the same size at wavelengths.

        Notes
        -----
        If neither ``solardistance`` or ``mjd`` is set, then a solar distance of 1 AU is assumed.

        Returns
        -------
        np.array
            Irradiance at the requested wavelength in photons/cm2/sec/nm.  Same size as wavelengths

        Raises
        ------
        ValueError
            If both solardistance and mjd are simultaneously set.

        Examples
        --------
        >>> from sasktran import SolarSpectrum
        >>> spectrum = SolarSpectrum()
        >>> spectrum.irradiance([340, 600])
        array([  1.76034000e+14,   5.38339000e+14])
        >>> spectrum.irradiance([340, 600], mjd=54372)
        array([  1.75359594e+14,   5.36276563e+14])
        >>> spectrum.irradiance([340, 600], solardistance=1.5)
        array([  7.82373333e+13,   2.39261778e+14])
        >>> spectrum.irradiance([340, 600], solardistance=1.5, fwhm=1)
        array([  8.02407128e+13,   2.36113610e+14])
        >>> spectrum.irradiance([340, 600], solardistance=1.5, mjd=54372)
        Traceback (most recent call last):
        ...
        ValueError: Can only set one of solardistance or mjd
        """
        if solardistance is not None and mjd is not None:
            raise ValueError('Can only set one of solardistance or mjd')
        if solardistance is None and mjd is None:
            solardistance = 1

        if solardistance is not None:
            self._isksolarspectrum.SetSolarDistanceFromAU(solardistance)
        if mjd is not None:
            self._isksolarspectrum.SetSolarDistanceFromMjd(mjd)

        if fwhm is None and stdev is None:
            # Just using the native resolution, don't need to convolve down
            return self._isksolarspectrum.Irradiance(np.atleast_1d(wavelengths))[1]
        else:
            if fwhm is not None and stdev is not None:
                raise ValueError('Only one of fwhm and stdev should be set')
            if stdev is None:
                stdev = fwhm / (2 * np.sqrt(2 * np.log(2)))

            convolution_wavel = self._convolution_wavelengths(np.atleast_1d(wavelengths), stdev)
            convolution_irradiance = self._isksolarspectrum.Irradiance(convolution_wavel)[1]

            output_irradiance = np.zeros(np.shape(np.atleast_1d(wavelengths)))

            if len(np.atleast_1d(stdev)) != len(np.atleast_1d(wavelengths)):
                stdev = np.ones_like(wavelengths) * stdev

            for idx, (w, s) in enumerate(zip(np.atleast_1d(wavelengths), np.atleast_1d(stdev))):
                # Modify the standard deviation according to the underlying high resolution
                hires_stdev = self.resolution(w) / (2 * np.sqrt(2 * np.log(2)))
                mod_stdev = np.sqrt(s**2 - hires_stdev**2)

                gaussian = np.exp(-(convolution_wavel - w)**2 / (2 * mod_stdev**2))
                gaussian /= np.nansum(gaussian)

                output_irradiance[idx] = np.dot(gaussian, convolution_irradiance)

            return output_irradiance

    def _convolution_wavelengths(self, wavelengths, stdev):

        spacing = self.resolution(np.nanmean(wavelengths))

        wavels = np.arange(np.nanmin(wavelengths) - 5 * np.nanmax(stdev),
                           np.nanmax(wavelengths) + 5 * np.nanmax(stdev), spacing)

        return wavels

    @wrap_skif_functionfail
    def resolution(self, wavelengths):
        """
        Gets the resolution in fwhm for the solar spectrum at a set of specific wavelengths.
        Parameters
        ----------
        wavelengths: np.array
            Wavelengths to get the resolution at

        Returns
        -------
        np.array
            Same size as wavelengths, resolution as fwhm in nm.

        Examples
        --------
        >>> from sasktran import SolarSpectrum
        >>> spectrum = SolarSpectrum()
        >>> spectrum.resolution([340, 600])
        array([ 0.04,  0.04])
        """
        return self._isksolarspectrum.NanometerResolutionFWHM(np.atleast_1d(wavelengths))[1]

    @property
    @wrap_skif_functionfail
    def min_wavelength(self):
        """
        Returns
        -------
        float
            Minimum wavelength in nm that this specific solar spectrum is valid for

        Examples
        --------
        >>> from sasktran import SolarSpectrum
        >>> spectrum = SolarSpectrum()
        >>> spectrum.min_wavelength
        200.07
        """
        return self._isksolarspectrum.MinValidWavelength()[1]

    @property
    @wrap_skif_functionfail
    def max_wavelength(self):
        """
        Returns
        -------
        float
            Maximum wavelength in nm that this specific solar spectrum is valid for

        Examples
        --------
        >>> from sasktran import SolarSpectrum
        >>> spectrum = SolarSpectrum()
        >>> spectrum.max_wavelength
        1000.99
        """
        return self._isksolarspectrum.MaxValidWavelength()[1]
