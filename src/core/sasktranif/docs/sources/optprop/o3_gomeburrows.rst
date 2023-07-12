
.. _optical_o3gomeburrows:


O3_GOMEBURROWS
==============
Table of O3 cross-sections measured with the GOME flight instrument at medium wavelength resolution from 231 nm to 794 nm.

Spectral resolution
^^^^^^^^^^^^^^^^^^^
* 231-307 nm, 0.20 nm
* 307-316 nm, 0.20 nm
* 311-405 nm, 0.17 nm
* 405-611 nm, 0.29 nm
* 595-794 nm, 0.33 nm

Temperature Range
^^^^^^^^^^^^^^^^^
The cross-sections were measured at 5 temperatures covering normal stratospheric and
tropospheric ranges. The paper does not provide any advice on how to interpolate in temperature. 
We use standard linear interpolation. The spectra were measured at the following five temperatures.

* 202 K,
* 221 K,
* 241 K,
* 273 K
* 293 K

Example
^^^^^^^
::

   optprop = ISKOpticalProperty('O3_GOMEBURROWS')
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
Burrows J. P., A. Richter, A. Dehn, B. Deters, S. Himmelmann, S. Voigt, and J. Orphal: "Atmospheric Remote-Sensing Reference Data from GOME: 2. Temperature-Dependent Absorption Cross Sections of O3 in the 231-794 nm Range", Journal of Quantitative Spectroscopy and Radiative Transfer 61, 509-517, 1999.
