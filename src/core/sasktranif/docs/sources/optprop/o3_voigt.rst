
.. _optical_o3voigt:

O3_VOIGT
========
Tabulated high resolution cross-sections of O3 measured by Voigt et al. 2001. The resolution is a 
constant 5 cm-1 (0.027 nm at 230 nm to 0.36 nm at 850 nm).

Temperatures
^^^^^^^^^^^^
The Voigt paper presents a new exponential interpolation algorithm for interpolating cross-sections in temperature. 
This class does not yet use this technique but uses truncated linear interpolation. Measurements are provided at 
the following five temperatures,

* 203 K
* 223 K
* 246 K
* 280 K
* 293 K

Example
^^^^^^^
::

   optprop = ISKOpticalProperty('O3_VOIGT')
   msis = ISKClimatology('MSIS90');
   mjd = 52393.3792987115;
   location = [0.0, 0.0, 25000.0, mjd];
   optprop.SetAtmosphericState( msis)
   optprop.SetLocation(location)
   [ok, abs, ext,sca] = optprop.CalculateCrossSections( 1.0E7/600.0, location );

Properties
^^^^^^^^^^
.. option:: SetTemperature (double n)
   
   Sets the temperature in Kelvins that will be used in the next calculation of cross-sections.

References
^^^^^^^^^^
1. Voigt S., J. Orphal, K. Bogumil, J.P. Burrows, The temperature dependence (203–293 K) of the absorption cross sections of O3 in the 230–850 nm region measured by Fourier-transform spectroscopy. Journal of Photochemistry and Photobiology A: Chemistry Vol. 143, 2001.
2.Orphal J. et al. A critical review of the absorption cross-sections of O3 and NO2 in the 240-790 nm region. ESA Technical Note MO-TN-ESA-GO-0302, 2002.

Abstract of Paper 1.
^^^^^^^^^^^^^^^^^^^^
::

   Absolute absorption cross sections of O3 were measured in the
   230–850 nm (11765–43478 cm?1) region at five different temperatures 
   (203–293 K) using a Fourier-transform spectrometer, at a spectral
   resolution of 5.0 cm-1 (corresponding to about 0.027 nm at 230 nm 
   and to about 0.36 nm at 850 nm). The spectral accuracy of the data 
   is better than 0.1 cm-1 — about 0.5 pm at 230 nm and about 7.2 pm 
   at 850 nm — validated by recording of I2 absorption spectra in the
   visible using the same experimental set-up. O3 absorption spectra 
   at different concentrations were recorded at five different sample 
   temperatures in the range 203–293 K, and at each temperature at two 
   total pressures (100 and 1000 mbar) using O2/N2 mixtures as buffer 
   gas. Within the limits of experimental uncertainties, no influence
   of total pressure on the O3 spectrum was observed in the entire region,
   as expected from the short lifetimes of the upper electronic states
   of O3. The temperature dependence of the O3 absorption cross sections
   is particularly strong in the Huggins bands between 310 and 380 nm,
   as observed in previous studies. An empirical formula is used to model
   the temperature dependence of the O3 absorption cross sections between
   236 and 362 nm, a spectral region that is particularly important for
   atmospheric remote-sensing and for photochemical modelling.
   
Header details from distributed Data Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

   ESA Study 11340/95/NL/CN UV - Cross-Sections in the UV and Visible
   O3 ABSORPTION CROSS-SECTIONS AT 203-293K
   Low pressure = 100 mbar total pressure
   High pressure = 1000 mbar total pressure
   S. Voigt, J. Orphal, and J. P. Burrows
   University of Bremen - Institute of Environmental Physics
   P. O. Box 33 04 40
   D-28334 Bremen, Germany
   Tel. + 49 (0)421 218 3526
   e-mail: Susanne.Voigt@iup.physik.uni-bremen.de
   Wavelength range: 230 - 850 nm
   Wavenumber range: 11752 - 43315 cm-1
   Spectral Resolution: 5 cm-1
   EXPERIMENTAL
   Cell: Multiple reflection quartz cell (White optics)
   optical pathlength 505 cm/120 cm (single path without White optics)

   Spectrometer:  BRUKER IFS 120HR Fourier-Transform-Spectrometer
   Wavenumber Range  Light Source      Detector
   ----------------  ------------      --------
   29500-43500 cm-1  Xe lamp     UV diode
   20000-33000 cm-1  Xe lamp     GaP diode
   12000-25000 cm-1  QTH lamp    Si diode
   Note: Spectral Regions are cutted and concatenated at 22000 cm-1, 31000 cm-1 and 33000 cm-1.
