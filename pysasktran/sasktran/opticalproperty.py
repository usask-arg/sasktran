import sasktranif.sasktranif as skif
from sasktran.climatology import Climatology, ClimatologyUserDefined
from typing import Optional
from collections import namedtuple
from copy import copy
from pathlib import Path
import numpy as np
import os
from sasktran.exceptions import wrap_skif_functionfail
from sasktran.util import to_iter
from sasktran import config


class OpticalProperty(object):
    """
    Generic optical property container that is applicable to a wide variety of optical properties.

    Parameters
    ----------
    name : str
        Name of the optical property to create

    Examples
    --------
    >>> from sasktran import OpticalProperty, Climatology
    >>> opt_prop = OpticalProperty('O3_OSIRISRES')
    >>> print(opt_prop)
    SasktranIF Optical Property: O3_OSIRISRES
    >>> atmospheric_state = Climatology('msis90')
    >>> opt_prop.calculate_cross_sections(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,\
                                          wavelengths=[350, 600])
    CrossSections(wavelengths=array([350, 600]), absorption=array([  7.09857544e-23,   5.27719024e-21]), scattering=array([ 0.,  0.]), total=array([  7.09857544e-23,   5.27719024e-21]))
    """
    _CrossSections = namedtuple('CrossSections', ['wavelengths', 'absorption', 'scattering', 'total'])

    @wrap_skif_functionfail
    def __init__(self, name: str):
        self._iskopticalproperty = skif.ISKOpticalProperty(name)
        self._name = name
        self.info = dict()

    def __getstate__(self):
        state = {k: v for k, v in self.__dict__.items() if k != '_iskopticalproperty'}
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._iskopticalproperty = skif.ISKOpticalProperty(state['_name'])

    def __repr__(self):
        return 'SasktranIF Optical Property: {}'.format(self._name)

    def skif_object(self, **kwargs):
        """

        Returns
        -------
        skif.ISKOpticalProperty
            Underlying ISKOpticalProperty object
        """
        return self._iskopticalproperty

    @wrap_skif_functionfail
    def calculate_cross_sections(self, atmospheric_state: Climatology, latitude: float, longitude: float,
                                 altitude: float, mjd: float, wavelengths: np.array):
        """
        Calculates absorption and scattering cross sections for the given optical property with a specific atmospheric
        state and location.

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

        Returns
        -------
        CrossSections
            NamedTuple with fields wavelengths, absorption, scattering, total.  wavelengths field matches the input
            wavelengths, absorption and scattering are the absorption and scattering cross sections respectively in
            :math:`\mathrm{cm^2}` and total is the sum of the two.

        Examples
        --------
        >>> from sasktran import OpticalProperty, Climatology
        >>> opt_prop = OpticalProperty('O3_OSIRISRES')
        >>> print(opt_prop)
        SasktranIF Optical Property: O3_OSIRISRES
        >>> atmospheric_state = Climatology('msis90')
        >>> opt_prop.calculate_cross_sections(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,\
                                              wavelengths=[350, 600])
        CrossSections(wavelengths=array([350, 600]), absorption=array([  7.09857544e-23,   5.27719024e-21]), scattering=array([ 0.,  0.]), total=array([  7.09857544e-23,   5.27719024e-21]))
        """
        self._iskopticalproperty.SetAtmosphericState(atmospheric_state.skif_object())
        self._iskopticalproperty.SetLocation([latitude, longitude, altitude, mjd])

        wavenumbers = np.asarray(1e7 / np.asarray(to_iter(wavelengths)))
        sort_index = np.argsort(wavenumbers)
        raw_result = self._iskopticalproperty.CalculateCrossSections(wavenumbers[sort_index])

        # applying the sorted sorting indices restores the original order
        reverse_index = np.argsort(sort_index)
        absorption = raw_result[1][reverse_index]
        total = raw_result[2][reverse_index]
        scattering = raw_result[3][reverse_index]

        return OpticalProperty._CrossSections(absorption=absorption, scattering=scattering,
                                              total=total, wavelengths=np.asarray(wavelengths))

    @wrap_skif_functionfail
    def calculate_phase_matrix(self, atmospheric_state: Climatology, latitude: float, longitude: float,
                               altitude: float, mjd: float, wavelengths: np.array, cosine_scatter_angles: np.array):
        """
        Calculates the scattering phase matrix for the given optical property with a specific atmospheric
        state and location.

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
            Wavelengths in nm. Shape (N_wavel,)
        cosine_scatter_angles : np.array
            Array of cosine of the scattering angle. Shape (N_angle,)

        Returns
        -------
        phase_matrix : np.ndarray
            Array of shape (N_wavel, N_angle, 4, 4) of phase matrices.

        Examples
        --------
        >>> import sasktran as sk
        >>> opt_prop = sk.Rayleigh()
        >>> atmospheric_state = sk.MSIS90()
        >>> opt_prop.calculate_phase_matrix(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,\
                                              wavelengths=[350], cosine_scatter_angles=[1])
        array([[[[ 1.47818063, -0.        ,  0.        ,  0.        ],
                 [-0.        ,  1.4345419 ,  0.        ,  0.        ],
                 [ 0.        ,  0.        ,  1.4345419 ,  0.        ],
                 [ 0.        ,  0.        ,  0.        ,  1.39018959]]]])

        """
        self._iskopticalproperty.SetAtmosphericState(atmospheric_state.skif_object())
        self._iskopticalproperty.SetLocation([latitude, longitude, altitude, mjd])

        wavenumbers = np.asarray(1e7 / np.asarray(to_iter(wavelengths)))
        phase_matrix = np.zeros((len(wavelengths), len(cosine_scatter_angles), 4, 4))

        for idx, w in enumerate(wavenumbers):
            for idy, cosine_scatter_angle in enumerate(cosine_scatter_angles):
                ok, raw = self._iskopticalproperty.CalculatePhaseMatrix(w, cosine_scatter_angle)

                phase_matrix[idx, idy, :, :] = raw.reshape((4, 4))

        return phase_matrix


class UserDefinedAbsorption(OpticalProperty):
    """
    User defined absorption optical property that may be temperature dependent.

    Parameters
    ----------
    wavelengths: np.ndarray
        Wavelengths the cross section is specified at. Shape (n,)
    cross_sections: np.ndarray
        Cross sections in /cm2. Shape (n,) if the cross section is not temperature dependent, and (n, m) if it is.
    temperatures: np.ndarray, optional
        Temperatures in K the cross sections are specified at, shape (m,).  May be None if the cross sections are
        not temperature dependent. Default None.

    Examples
    --------
    >>> from sasktran import UserDefinedAbsorption, MSIS90
    >>> wavelengths = [350, 800]
    >>> cross_sections = [1, 2]
    >>> opt_prop = UserDefinedAbsorption(wavelengths, cross_sections)
    >>> print(opt_prop)
    SasktranIF Optical Property: USERDEFINED_TABLES
    >>> atmospheric_state = MSIS90()
    >>> opt_prop.calculate_cross_sections(atmospheric_state, latitude=0, longitude=0, altitude=20000, mjd=54372,\
                                          wavelengths=[350, 600])
    CrossSections(wavelengths=array([350, 600]), absorption=array([ 1.        ,  1.55555556]), scattering=array([ 0.,  0.]), total=array([ 1.        ,  1.55555556]))
    """
    def __init__(self, wavelengths: np.ndarray, cross_sections: np.ndarray, temperatures: np.ndarray = None):
        super().__init__('USERDEFINED_TABLES')

        try:
            self.skif_object().SetProperty('WavelengthTruncation', 1)
        except skif._sasktranif.functionfail:  # pragma: no cover
            # This function mistakingly returns fail in older versions of SASKTRANIF so we suppress the error
            # (It still works even if it fails)
            pass

        if temperatures is None:
            # Cross section is not temperature dependent, so just cover the full temperature range
            self.skif_object().AddUserDefined(0, wavelengths, cross_sections)
            self.skif_object().AddUserDefined(1000, wavelengths, cross_sections)
        else:
            for idx, T in enumerate(temperatures):
                self.skif_object().AddUserDefined(T, wavelengths, cross_sections[:, idx])


class UserDefinedAbsorptionPressure(OpticalProperty):
    """
    A user defined absorption optical property that is a function of temperature, pressure, wavelength, and VMR
    of a secondary species that broadens the profile.

    Parameters
    ----------
    broadener_vmr: np.array
        VMR of the broadening species.  Must be in ascending order. Shape (Nv)
    pressures: np.array
        Pressure in Pa, must be in ascending order. Shape (Np)
    temperatures: np.array
        Temperatures in K at each pressure level.  Must be in ascending order. Shape (Np, Nt)
    wavel_nm: np.array
        Wavelengths in nm.  Must be in ascending order. Shape (Nw)
    xs: np.ndarray
        Cross sections in cm^2.  Shape (Nv, Np, Nt, Nw)
    broadener_clim: Optional[Climatology]
        A climatology for the broadening species, can be None if there is no broadener. Default None.
    broadener_handle: Optional[str]
        The climatology handle for the broadener. Can be None if there is no broadener. Default None.
    """
    def __init__(self,
                 broadener_vmr: np.array,
                 pressures: np.array,
                 temperatures: np.ndarray,
                 wavel_nm: np.array,
                 xs: np.ndarray,
                 broadener_clim: Optional[Climatology] = None,
                 broadener_handle: Optional[str] = None):
        super().__init__('USERDEFINED_PRESSURE')

        for cross, vmr in zip(xs, broadener_vmr):
            self.skif_object().AddUserDefinedPressure(pressures.flatten(), temperatures.flatten(), wavel_nm.flatten(), cross.flatten(), float(vmr))

        if broadener_clim is not None:
            self.skif_object().SetProperty('broadenerclimatology', broadener_clim.skif_object())
            self.skif_object().SetProperty('broadenerhandle', broadener_handle)


class UserDefinedScatterConstantHeight(OpticalProperty):
    r"""
    A user defined optical property for a scattering species that has constant optical properties in location space.

    Currently this optical property only supports specifying the phase function in terms of the Legendre expansion,
    this means that it is only useful for the SASKTRAN-DO engine.

    The Legendre expansion is defined in terms of "greek" coefficients, where the phase matrix has the standard form

    .. math::

        P(\Theta) = \begin{pmatrix}
                      a1(\Theta) & b1(\Theta) & 0 & 0 \\
                      b1(\Theta) & a2(\Theta) & 0 & 0 \\
                      0 & 0 & a3(\Theta) & b2(\Theta) \\
                      0 & 0 & -b2(\Theta) & a4(\Theta)
                    \end{pmatrix}

    The phase function is then expanded in terms of Wigner functions (generalized spherical functions related to
    Legendre polynomials/associated Legendre functions)

    .. math::

        a_1(\Theta) = \sum_l \alpha_{1,l} \, d^{l}_{00}(\Theta)

        a_2(\Theta) + a_3(\Theta) = \sum_l (\alpha_{2,l} + \alpha_{3,l}) \, d^{l}_{2,2}(\Theta)

        a_2(\Theta) - a_3(\Theta) = \sum_l (\alpha_{2,l} - \alpha_{3,l}) \, d^{l}_{2,-2}(\Theta)

        a_4(\Theta) = \sum_l \alpha_{4,l} \, d^{l}_{00}(\Theta)

        b_1(\Theta) = \sum_l \beta_{1,l} \, d^{l}_{02}(\Theta)

        b_2(\Theta) = -\sum_l \beta_{2,l} \, d^{l}_{02}(\Theta)

    Parameters
    ----------
    wavelengths : np.array
        Array of wavelengths in nm.  Shape (n,)
    xs_scat : np.array
        Scattering cross section in cm2 at wavelengths. Shape (n,)
    xs_abs : np.array
        Absorption cross section in cm2 at wavelengths. Shape (n,)
    lm_a1: np.ndarray
        Legendre coefficients :math:`\alpha_{1,l}`.  Shape (n, m).
    lm_a2: np.ndarray
        Legendre coefficients :math:`\alpha_{2,l}`.  Shape (n, m).
    lm_a3: np.ndarray
        Legendre coefficients :math:`\alpha_{3,l}`.  Shape (n, m).
    lm_a4: np.ndarray
        Legendre coefficients :math:`\alpha_{4,l}`.  Shape (n, m).
    lm_b1: np.ndarray
        Legendre coefficients :math:`\beta_{1,l}`.  Shape (n, m).
    lm_b2: np.ndarray
        Legendre coefficients :math:`\beta_{2,l}`.  Shape (n, m).
    """
    def __init__(self,
                 wavelengths: np.array,
                 xs_scat: np.array,
                 xs_abs: np.array,
                 legendre_moments: np.ndarray = None,
                 lm_a1: np.ndarray = None,
                 lm_a2: np.ndarray = None,
                 lm_a3: np.ndarray = None,
                 lm_a4: np.ndarray = None,
                 lm_b1: np.ndarray = None,
                 lm_b2: np.ndarray = None,
                 ):
        super().__init__('USERDEFINED_SCATTERCONSTANTHEIGHT')

        self.skif_object().SetProperty('wavelengths', wavelengths)
        self.skif_object().SetProperty('xs_scat', xs_scat)
        self.skif_object().SetProperty('xs_abs', xs_abs)
        if legendre_moments is not None:
            self.skif_object().SetProperty('legendremoments', legendre_moments)

        if lm_a1 is not None:
            self.skif_object().SetProperty('legendremomentsa1', lm_a1)

        if lm_a2 is not None:
            self.skif_object().SetProperty('legendremomentsa2', lm_a2)

        if lm_a3 is not None:
            self.skif_object().SetProperty('legendremomentsa3', lm_a3)

        if lm_a4 is not None:
            self.skif_object().SetProperty('legendremomentsa4', lm_a4)

        if lm_b1 is not None:
            self.skif_object().SetProperty('legendremomentsb1', lm_b1)

        if lm_b2 is not None:
            self.skif_object().SetProperty('legendremomentsb2', lm_b2)


class OpticalPropertyConvolved(OpticalProperty):
    """
    Generic optical property container that is applicable to a wide variety of optical properties.

    Parameters
    ----------
    name : str
        Name of the optical property to create
    psf_wavelength : np.ndarray
        array of wavelengths in nanometers where the point spread function is defined.
    psf : np.ndarray
        array of point spread function values in nanometers. Same length as psf_wavelength
    wavel_spacing : scalar, optional
        Spacing to calculate the high resolution wavelengths on in nanometers, default is 0.01nm
    num_stdev : scalar, optional
        Number of standard deviations to use in the gaussian for convolution, default is 4
    output_spacing : scalar, string optional
        Spacing to calculate the convolved cross section on. Either a scalar to use a fixed spacing, or None to use
        the wavelengths of the point spread function (wavel parameter). Default is none.

    Examples
    --------
    >>> from sasktran import OpticalPropertyConvolved, MSIS90, O3DBM
    >>> import numpy as np
    >>> convolved_dbm = OpticalPropertyConvolved(O3DBM(), psf_wavelength=np.linspace(250, 1000, 50),\
                                                 psf=np.linspace(0.5, 10, 50))
    >>> print(convolved_dbm)
    SasktranIF Optical Property: O3_DBM_Convolved
    >>> convolved_dbm.calculate_cross_sections(MSIS90(), latitude=0, longitude=0, altitude=20000, mjd=54372,\
                                               wavelengths=[350, 600])
    CrossSections(wavelengths=array([350, 600]), absorption=array([  2.83084701e-22,   5.01481379e-21]), scattering=array([ 0.,  0.]), total=array([  2.83084701e-22,   5.01481379e-21]))
    """

    def __init__(self, optical_property: OpticalProperty, psf_wavelength: np.ndarray, psf: np.ndarray,
                 wavel_spacing: float=0.01, num_stdev: float=4, output_spacing: np.ndarray=None):
        super().__init__(optical_property._name)

        self._unconv_optical_prop = optical_property
        self._psf_wavelength = psf_wavelength
        self._psf = psf
        self._wavel_spacing = wavel_spacing
        self._num_stdev = num_stdev
        self._output_spacing = output_spacing

        self._iskopticalproperty = self.convolved_optical_property(self._unconv_optical_prop,
                                                                   self._psf_wavelength,
                                                                   self._psf,
                                                                   self._wavel_spacing,
                                                                   self._num_stdev,
                                                                   self._output_spacing)

    def __setstate__(self, state):
        super().__setstate__(state)
        self._iskopticalproperty = self.convolved_optical_property(self._unconv_optical_prop,
                                                                   self._psf_wavelength,
                                                                   self._psf,
                                                                   self._wavel_spacing,
                                                                   self._num_stdev,
                                                                   self._output_spacing)

    def __repr__(self):
        return 'SasktranIF Optical Property: {}_Convolved'.format(self._name)

    @staticmethod
    @wrap_skif_functionfail
    def convolved_optical_property(hires_optprop: OpticalProperty, wavel: np.ndarray, psf: np.ndarray,
                                   wavel_spacing: float=0.01, num_stdev: float=4, output_spacing: np.ndarray=None):
        """
        Convolves down the hi-resolution cross sections and creates a user defined optical property. Convolution assumes
        the hi-resolution version has infinite resolution

        Parameters
        ----------
        hires_optprop : sasktran.OpticalProperty
            sasktran OpticalProperty that will be convolved down to desired resolution.
        wavel : numpy array
            Wavelengths to calculate the new cross section for in nm
        psf : numpy array
            Standard deviations of the gaussian to convolve by for each wavelength
            Can be the same size as wavel, or size 1 in which case the psf is assumed
            to be the same for all wavelengths
        wavel_spacing : scalar, optional
            Spacing to calculate the high resolution wavelengths on in nanometers, default is 0.01nm
        num_stdev : scalar, optional
            Number of standard deviations to use in the gaussian for convolution, default is 4
        output_spacing : scalar,string optional
            Spacing to calculate the convolved cross section on. Either a scalar to use a fixed spacing, or None to use
            the wavelengths of the point spread function (wavel parameter). Default is none.

        Returns
        -------
        optprop : ISKOpticalProperty
            Convolved optical property
        """

        optprop = skif.ISKOpticalProperty('USERDEFINED_TABLES')
        try:
            optprop.SetProperty('WavelengthTruncation', 1)
        except skif._sasktranif.functionfail:  # pragma: no cover
            # This function mistakingly returns fail in older versions of SASKTRANIF so we supress the error
            # (It still works even if it fails)
            pass

        # If we only have one psf value, copy it for every wavelength
        if len(np.atleast_1d(psf)) == 1:
            psf = np.repeat(psf, len(wavel))

        # Temperatures of the hi-resolution xsect
        temp = hires_optprop.info['temperatures']
        # Ranges in which each temperature has data
        temp_range = hires_optprop.info['wavelength_range']

        for T, ra in zip(temp, temp_range):
            # Get the high resolution DBM cross section
            clim_temp = ClimatologyUserDefined([0, 10000], {'SKCLIMATOLOGY_TEMPERATURE_K': [T, T]}, interp='linear')
            wavel_dbm = np.linspace(ra[0], ra[1], int(np.floor((ra[1] - ra[0]) / wavel_spacing)))
            xsect = hires_optprop.calculate_cross_sections(clim_temp, 0, 0, 1000, 54372, wavel_dbm)
            xsect = xsect.absorption

            # Only do the convolution for wavelengths in the temperature range
            if output_spacing is not None:
                wavel_to_calc_for = np.linspace(ra[0], ra[1], int(np.floor((ra[1] - ra[0]) / output_spacing)))
                psf_to_calc_for = np.interp(wavel_to_calc_for, wavel, psf)
            else:
                good = (wavel > ra[0]) & (wavel < ra[1])
                wavel_to_calc_for = wavel[good]
                psf_to_calc_for = psf[good]

                # Add a little bit to the edges to avoid nan's on the boundary
                wavel_to_calc_for = np.concatenate(([wavel_to_calc_for[0] - wavel_spacing], wavel_to_calc_for, [wavel_to_calc_for[-1] + wavel_spacing]))
                psf_to_calc_for = np.concatenate(([psf_to_calc_for[0]], psf_to_calc_for, [psf_to_calc_for[-1]]))

            wavelconv = np.zeros(np.shape(wavel_to_calc_for))
            xsectconv = np.zeros(np.shape(wavel_to_calc_for))

            for widx, (w, p) in enumerate(zip(wavel_to_calc_for, psf_to_calc_for)):
                if (w > (wavel_dbm[0] + 1 * num_stdev * p)) | (w < wavel_dbm[-1] - 1 * num_stdev * p):

                    good = (wavel_dbm > w - num_stdev * p) & (wavel_dbm < w + num_stdev * p)
                    if not np.any(good):
                        wavelconv[widx] = w
                        xsectconv[widx] = np.interp(w, wavel_dbm, xsect)
                    else:
                        wavel_filter = wavel_dbm[good] - w
                        f = np.exp(-0.5 * (wavel_filter / p)**2)
                        f /= np.sum(f)

                        wavelconv[widx] = np.dot(wavel_dbm[good], f)
                        xsectconv[widx] = np.dot(xsect[good], f)

            good = xsectconv > 0
            full_xsect = np.interp(wavel_to_calc_for, wavelconv[good], xsectconv[good], left=np.nan, right=np.nan)

            if len(wavel_to_calc_for[~np.isnan(full_xsect)]) == 0:
                continue
            optprop.AddUserDefined(T, wavel_to_calc_for[~np.isnan(full_xsect)], full_xsect[~np.isnan(full_xsect)])

        return optprop


class MieAerosol(OpticalProperty):
    """
    Specialized OpticalProperty which supports Mie Aerosol calculations.

    Parameters
    ----------
    particlesize_climatology : sasktran.Climatology

    species : str
        Molecule to use, one of ['H2SO4', 'ICE', 'WATER']
    """
    def __init__(self, particlesize_climatology: Climatology, species: str):
        super().__init__('MIEAEROSOL_' + species.upper())

        self._species = species
        self._particlesize_climatology = particlesize_climatology
        self._update_opticalproperty()

    def __setstate__(self, state):
        super().__setstate__(state)
        self._update_opticalproperty()

    @wrap_skif_functionfail
    def _update_opticalproperty(self):
        self._iskopticalproperty.SetProperty('SetParticleSizeClimatology', self._particlesize_climatology.skif_object())

    @property
    def particlesize_climatology(self):
        return self._particlesize_climatology

    @particlesize_climatology.setter
    def particlesize_climatology(self, value):
        self._particlesize_climatology = value
        self._update_opticalproperty()


class HITRANChemical(OpticalProperty):
    """
    Calculates the optical absorption and extinction of various atmospheric molecules using the Voigt line-shape and
    the HITRAN spectral line database. The object supports all of the HITRAN species specified in the HITRAN database
    file molparam.txt.  This optical property requires additional configuration, see :ref:`configuration` for more information.

    Parameters
    ----------
    chemical_name : str
        Chemical abbreviation of the molecule of interest.
    isotope_filter : int, optional
        Allows the HITRAN object to load in just one isotope of the requested molecule. The value set must match one of
        the isotope labels used for the given molecule in the HITRAN database file, molparam.txt. Note that the code
        does not adjust the line strength but uses the line strength value as written in the HITRAN database. This
        means you may have to account for and/or remove the abundance automatically built into the HITRAN database line
        strength values.
    line_tolerance : float, optional
        Allows the user to set the tolerance used to reject weak lines from the current micro-window as part of a speed optimization strategy. The default value
        is 0.0 which disables the optimization. Larger values speed up calculation of spectra but may result in choppy spectra at the smaller
        intensities, especially in extinction/absorption spectra which typically follow the log of the cross-section. A
        smaller value will reduce choppiness but increase computational speed. Only values greater
        than or equal to zero are acceptable. A value of 1.0E-09 is a good starting value for users wishing to investigate
        the effectiveness of the speed optimization.
    max_line_strength : float, optional
        Allows the user to manually set the maximum line strength within a micro-window. By default the object will take this value
        from the strongest line in the micro-window. The value is used with the line tolerance to reject weak lines from
        spectral calculations. A value of zero will disable the manual setting and reinstate usage of the default.
    use_cache: bool, optional
        If true then cross sections will be cached when used internally inside a radiative transfer calculation.
        The default is True, this should be left as True unless trying to calculate millions of wavelengths where
        memory starts to become an issue.
    """
    @wrap_skif_functionfail
    def __init__(self, chemical_name: str, line_tolerance=None, max_line_strength=None, isotope_filter=None,
                 use_cache=True):
        self._validate_registry()
        optical_property_name = "HITRANCHEMICAL_" + chemical_name
        super().__init__(optical_property_name)
        self._optical_property_name = optical_property_name

        self._line_tolerance = line_tolerance
        self._max_line_strength = max_line_strength
        self._isotope_filter = isotope_filter

        self._use_cache = use_cache
        self._cached_wavelengths = None
        self._wavenumber_range = None

        self._update_opticalproperty()

    def __setstate__(self, state):
        super().__setstate__(state)
        self._update_opticalproperty()

    def _update_opticalproperty(self):
        if self._max_line_strength is not None:
            self._iskopticalproperty.SetProperty('setmaxlinestrength', self._max_line_strength)

        if self._line_tolerance is not None:
            self._iskopticalproperty.SetProperty('setlinetolerance', self._line_tolerance)

        if self._isotope_filter is not None:
            self._iskopticalproperty.SetProperty('setisotopefilter', self._isotope_filter)

        if self._wavenumber_range is not None:
            self._iskopticalproperty.SetProperty('setwavenumberrange', self._wavenumber_range)

        if self._cached_wavelengths is not None:
            self._iskopticalproperty.SetProperty('enablecachedcrosssections', 1e7 / self._cached_wavelengths)

    @staticmethod
    def _validate_registry():
        registry_config = config.load_sasktran_registry()

        hitran_folder = registry_config['software']['usask-arg']['skopticalproperties']['hitran']['basedirectory']

        if hitran_folder.lower() == 'undefined':
            user_config = config.user_config_file_location()
            raise IOError('HITRAN Folder is not set, please add hitran_directory: xxxx to the file {}'.format(user_config))

    def _reset_skif_object(self):
        self._iskopticalproperty = skif.ISKOpticalProperty(self._optical_property_name)
        self._update_opticalproperty()

    @wrap_skif_functionfail
    def skif_object(self, **kwargs):
        # Check to see if we are being called by an engine prior
        if 'engine' in kwargs:
            wavelengths = np.asarray(kwargs['engine'].wavelengths)

            min_wavenumber = np.nanmin(1e7 / wavelengths)
            max_wavenumber = np.nanmax(1e7 / wavelengths)

            new_wavenumber_range = np.asarray([min_wavenumber, max_wavenumber])

            # If we are not using the cache, just check if the wavenumber range has changed and if so reset the
            # object
            if not self._use_cache:
                if self._wavenumber_range is None or self._wavenumber_range[0] != min_wavenumber or self._wavenumber_range[1] != max_wavenumber:
                    self._wavenumber_range = new_wavenumber_range
                    self._reset_skif_object()
            else:
                # We are using the cache, so we need to check if either the wavenumber range or the cache has changed
                if self._wavenumber_range is None or self._cached_wavelengths is None or self._wavenumber_range[0] != min_wavenumber or\
                        self._wavenumber_range[1] != max_wavenumber or not np.array_equal(wavelengths, self._cached_wavelengths):
                    self._wavenumber_range = copy(new_wavenumber_range)
                    self._cached_wavelengths = copy(wavelengths)
                    self._reset_skif_object()

        return self._iskopticalproperty


class O3DBM(OpticalProperty):
    """
    Tabulated high resolution cross-sections of O3 measured by Daumont, Brion and Malicet in the early 1990’s [1].
    The wavelength range slightly varies with temperature but covers the entire UV to NIR region, from 194.50 nm to
    830.00 nm at 0.01 to 0.02 nm resolution. The cross-section data were collected at 0.01-0.02 nm resolution and
    each wavelength/cross-section table varies in size from 22,052 to 63,501 entries. The data consists of 5 tables
    of wavelength versus cross-section for 5 temperatures.

    Notes
    -----
    Temerature Range
        Measurements are provided at 5 temperatures covering typical stratospheric and
        tropospheric conditions:

            | 218 K
            | 228 K
            | 243 K
            | 273 K
            | 295 K

    Wavelength Range
        The wavelength range of each temperature table is slightly different and is given below.
        Note that most of the temperature variation occurs in the Huggins band
        between 315 and 360 nm:

            | 218K -> 194.50nm to 650.01nm
            | 228K -> 194.50nm to 520.01nm
            | 243K -> 194.50nm to 519.01nm
            | 273K -> 299.50nm to 520.01nm
            | 295K -> 195.00nm to 830.00nm

        We looked into temperature interpolation and while DBM suggest that a quadratic interpolation scheme [3] they do
        not indicate an explicit technique. We tested several quadratic fitting routines and found that a truncated linear
        fit in temperature was visually more appealing than any of the quadratic fits and had none of the undesirable
        artifacts (excessive curvature etc.) that naturally arise with quadratic curve fitting. Consequently this object
        uses a truncated linear fit in temperature.

    Data Source
        These data are an exact replication of the data files:

            | O3_CRS_BDM_218K.dat
            | O3_CRS_BDM_228K.dat
            | O3_CRS_BDM_243K.dat
            | O3_CRS_BDM_273K.dat
            | O3_CRS_BDM_295K.dat

        Data is from the IGACO site, http://igaco-o3.fmi.fi/ACSO/files/cross_sections. The files were copied on
        July 16-July 25 2012.

    References
    ----------
    .. [1] Daumont, D., et al. "Ozone UV spectroscopy I: Absorption cross-sections at room temperature."
             Journal of Atmospheric Chemistry 15.2 (1992): 145-155.
    .. [2] Brion, J., et al. "High-resolution laboratory absorption cross section of O3. Temperature effect."
             Chemical physics letters 213.5 - 6(1993): 610-612.
    .. [3] Malicet, J., et al. "Ozone UV spectroscopy. II. Absorption cross-sections and temperature dependence."
             Journal of Atmospheric Chemistry 21.3 (1995): 263-273.
    .. [4] Brion, J., et al. "Absorption spectra measurements for the ozone molecule in the 350–830 nm region."
            Journal of Atmospheric Chemistry 30.2 (1998): 291-299.

    """
    def __init__(self):
        super().__init__('O3_DBM')
        self.info['temperatures'] = [218, 228, 243, 273, 295]
        self.info['wavelength_range'] = [[194.50, 650.01], [194.50, 520.01], [194.50, 519.01], [299.50, 520.01],
                                         [195.00, 830.00]]
        self.info['spectral sampling'] = 0.01
        self.info['spectral resolution'] = 0.02


class O3OSIRISRes(OpticalProperty):
    """"
    These are the cross-sections used in the OSIRIS level 2 analysis for Saskmart V5.07. This table is based upon
    the cross-sections of Bogumil Orphal and Burrows and has been reduced to the nominal resolution of the OSIRIS
    instrument (approx. 1.0 nm). The cross-sections range from 270nm to 820 nm in 0.1 nm steps.
    """
    def __init__(self):
        super().__init__('O3_OSIRISRES')
        self.info['temperatures'] = [203, 223, 243, 273, 293]
        self.info['wavelength_range'] = [[270.0, 820.0], [270.0, 820.0], [270.0, 820.0], [270.0, 820.0], [270.0, 820.0]]
        self.info['spectral sampling'] = 0.1
        self.info['spectral resolution'] = 1.0


class NO2Vandaele1998(OpticalProperty):
    """
    Calculates the absorption cross section of NO2 molecules from 230 nm to 1000 nm at 220 K to 294 K following [1]

    .. [1] Vandaele, Ann Carine, et al. "Measurements of the NO2 absorption cross-section from 42 000 cm− 1 to
           10 000 cm− 1 (238–1000 nm) at 220 K and 294 K." Journal of Quantitative Spectroscopy and Radiative Transfer
           59.3-5 (1998): 171-184.
    """
    def __init__(self):
        super().__init__('NO2_VANDAELE1998')
        self.info['temperatures'] = [220, 294]
        self.info['wavelength_range'] = [[238.0, 1000.0], [238.0, 1000.0]]
        self.info['spectral sampling'] = np.nan
        w = np.arange(238, 1000, 1.0)
        self.info['spectral resolution'] = 2 * w**2 / (1e7 - 2 * w)
        self.info['spectral resolution wavelengths'] = w


class NO2OSIRISRes(OpticalProperty):
    """
    Calculates the absorption cross section of NO2 molecules from 230 nm to 795 nm and 221K to 293K. The cross-sections
    have been reduced to the resolution of OSIRIS and these cross-sections have been used in the OSIRIS level 2 MART
    retrievals.
    """
    def __init__(self):
        super().__init__('NO2_OSIRISRES')
        self.info['temperatures'] = [221, 293]
        self.info['wavelength_range'] = [[230.0, 795.0], [230.0, 795.0]]
        self.info['spectral sampling'] = 0.1
        self.info['spectral resolution'] = 1.0


class Rayleigh(OpticalProperty):
    """
    An optical property object that provides Rayleigh molecular scattering in dry-air. The code closely follows the
    algorithm published by Bates 1984 and exactly replicates his cross-section calculations to the 4 significant digits
    in his Table 1. The cross-section is weighted to account for the different gas ratios in standard atmospheric
    composition. No attempt is made to track changes in CO2 composition. Note that water vapour effects are implicitly
    ignored as it only considers dry air. The only difference with Bates is that this object takes into account the
    tiny fraction of gas that is not N2, O2, Argon or CO2. Bates ignores this component while this object assumes it
    the residual gas with properties similar to Argon.
    """
    def __init__(self):
        super().__init__('RAYLEIGH')


class SimpleRayleigh(OpticalProperty):
    """
    An optical property object that provides Rayleigh molecular scattering without any corrections.
    This is very similar to the :py:class:`sasktran.Rayleigh` optical property except that the Bates correction factor
    is not included. For most calculations it is preferred to use the :py:class:`sasktran.Rayleigh` optical property
    instead.
    """
    def __init__(self):
        super().__init__('SIMPLERAYLEIGH')


class InelasticRayleigh(OpticalProperty):
    """
    An optical property class that provides Rayleigh molecular scattering in dry-air, with the elastic Cabannes line
    separated from the inelastic Raman wings. Extinction is exactly the same as the Rayleigh class, but the scattering
    cross section now represents the Cabannes line only, which is a few percent smaller than the full Rayleigh scatter.
    The Raman lines are currently inaccessible through the python interface.
    """
    def __init__(self):
        super().__init__('INELASTICRAYLEIGH')


class BaumIceCrystal(OpticalProperty):
    """
    Scattering non-spherical ice crystals based upon the database from Baum.

    Parameters
    ----------
    effective_size_microns : float
        Size of the particles, typically from 10 microns to 50 microns
    use_delta_eddington : bool, optional
        True if the delta eddington approximation is to be used.  Should be True for use
        in either the HR or DO engines. Default: True.

    """
    def __init__(self, effective_size_microns, use_delta_eddington=True):
        super().__init__('BAUM_ICECRYSTALS')
        self._validate_registry()

        self._effective_size_microns = effective_size_microns
        self._use_delta_eddington = use_delta_eddington

        self._update_opticalproperty()

    def __setstate__(self, state):
        super().__setstate__(state)
        self._update_opticalproperty()

    @wrap_skif_functionfail
    def _update_opticalproperty(self):
        if self._effective_size_microns is not None:
            self._iskopticalproperty.SetProperty('effectivesizemicrons', self._effective_size_microns)

        if self._use_delta_eddington is not None:
            self._iskopticalproperty.SetProperty('usedeltaeddington', int(self._use_delta_eddington))

    @staticmethod
    def _validate_registry():
        registry_config = config.load_sasktran_registry()

        baum_folder = registry_config['software']['usask-arg']['skopticalproperties']['baum_icecrystals']['storage']['2014database']

        if not os.path.isfile(os.path.join(baum_folder, 'GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc')):
            user_config = config.user_config_file_location()
            raise IOError(
                'Baum Folder is invalid, please add baum_directory: xxxx to the file {}'.format(user_config))


class SO2Vandaele2009(OpticalProperty):
    """
    Calculates the absorption cross section of SO2 molecules from 227 nm to 420 nm at 298 K, 318 K, 338 K, and 358 K
    following [1] and [2].

    .. [1] C. Hermans, A.C. Vandaele, and S. Fally. "Fourier transform measurements of SO2 absorption cross sections:
           I. Temperature dependence in the 24000-29000 cm-1 (345-420 nm) region," J. Quant. Spectrosc. Radiat.
           Transfer 110, 756-765 (2009); DOI: 10.1016/j.jqsrt.2009.01.031
    .. [2] A.C. Vandaele, C. Hermans, and S. Fally, "Fourier transform measurements of SO2 absorption cross sections:
           II. Temperature dependence in the 29000-44000 cm-1 (227-345 nm) region," J. Quant. Spectrosc. Radiat.
           Transfer 110, 2115-2126 (2009); DOI: 10.1016/j.jqsrt.2009.05.006
    """
    def __init__(self):
        super().__init__('SO2_VANDAELE2009')
        self.info['temperatures'] = [298, 318, 338, 358]
        self.info['wavelength_range'] = [[227, 420]] * 4
        self.info['spectral sampling'] = np.nan
        w = np.arange(227, 420, 1.0)
        self.info['spectral resolution'] = 2 * w**2 / (1e7 - 2 * w)
        self.info['spectral resolution wavelengths'] = w


class O2O2HITRAN2016(OpticalProperty):
    """
    The O2-O2 collision induced absorption cross-sections distributed in Hitran 2016 as described by Karman et al. 2019.
    It is composed of 8 spectral regions measured by several researchers that extend from 335nm to over 8 microns.
    Some of the regions are temperature dependent and others are not
    """
    def __init__(self):
        super().__init__('O2_O2_HITRAN2016')


class O2O2Thalman2013(OpticalProperty):
    """
    The collision induced absorption cross-sections measured by Thalman et al. 2013.
    The cross-sections were measured at room temperature with a Fourier Transform spectrometer in the 15000 to 29800 cm−1 region (335-667 nm) at a maximal optical path difference of 0.45 cm (resolution 2 cm-1).
    Note that an appropriate climatology is the square of the O2 number density.
    This is provided by climatologies that support SKCLIMATOLOGY_O2_O2_CM6. For example, see MSIS90 and ECMWF:
    """
    def __init__(self):
        super().__init__('O2_O2_THALMAN2013')


class O2O2Fally2000(OpticalProperty):
    """
    The collision induced absorption cross-sections measured by Fally et al. 2000.
    The cross-sections were measured at room temperature with a Fourier Transform spectrometer in the 15000 to 29800 cm−1 region (335-667 nm) at a maximal optical path difference of 0.45 cm (resolution 2 cm-1).
    Note that an appropriate climatology is the square of the O2 number density.
    This is provided by climatologies that support SKCLIMATOLOGY_O2_O2_CM6.
    """
    def __init__(self):
        super().__init__('O2_O2_FALLY2000')