
.. _optical_o3sciabogumilv3:

O3_SCIABOGUMILV3
================
The version 3.0 cross-sections of O3 measured by Bogumil et al. with the Sciamachy flight instrument
before launch. Note that a later version, :ref:`optical_o3sciabogumilv4` is also available. The
cross-sections are measured from 230 nm to 1070 nm at five separate temperatures

* 203 K,
* 223 K,
* 243 K,
* 273 K
* 293 K

All measurements are at the intrinsic FWHM spectral resolution of Sciamachy,

* 0.32 nm below 311.7 nm
* 0.21 nm between 310.8 nm and 402 nm
* 0.52 nm between 402 nm and 597.8 nm
* 0.47 nm between 597.8 nm and 781.6 nm
* 0.62 nm between 781.6 nm and 1052.5 nm
* 1.45 nm above 1052.5 nm

Note: Spectral Regions are cutted and concatenated at 287.7 nm, 302.2 nm, 310.8 nm, 316.3 nm, 334.2 nm, 345.2 nm, 402 nm, 489.5 nm, 597.8 nm, 718.9 nm, 781.6 nm and 1052.5 nm.

Example
^^^^^^^
::

   optprop = ISKOpticalProperty('O3_SCIABOGUMILV3')
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
Bogumil K., J. Orphal, J.P. Burrows, J.M. Flaud: Vibrational progressions in the visible and near ultra-violet absorption spectrum of ozone. Chem Phys Lett, 349, pages 241-248, 2001.

