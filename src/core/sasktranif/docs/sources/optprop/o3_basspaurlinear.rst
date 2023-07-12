.. _optical_basspaurlinear:


O3_BASSPAURLINEAR
=================
Tables of the O3 cross-section at 6 different temperatures made by Bass-Paur.
A complementary class, :ref:`optical_basspaurquadratic`, provides the same data with
superior quadratic temperature interpolation. The data are a world standard but 
show signs of age and are of limited wavelength extent with no measurements at 
wavelengths longer than 343 nm. They are at high wavelength sampling and high 
wavelength resolution. 

Wavelength Range
^^^^^^^^^^^^^^^^
The table covers the wavelength range 245.0180 to 342.7800 nm in 1956 samples, of approximately 0.05 nm per step. 
The spectral resolution is better than 0.025 nm. 

Temperature Range
^^^^^^^^^^^^^^^^^
The Bass-Paur tables provide measurements at 6 temperatures:

1. 203 K
2. 223 K
3. 246 K
4. 273 K
5. 276 K 
6. 280 K

This object piecewise linearly interpolates cross-sections in temperature.

Example
^^^^^^^
::

   optprop = ISKOpticalProperty('O3_BASSPAURLINEAR')
   msis = ISKClimatology('MSIS90');
   mjd = 52393.3792987115;
   location = [0.0, 0.0, 25000.0, mjd];
   optprop.SetAtmosphericState( msis)
   optprop.SetLocation(location)
   [ok, abs, ext,sca] = optprop.CalculateCrossSections( 1.0E7/600.0, location );

Data Source
^^^^^^^^^^^
These data are an exact replication of the data files,

* bp_203clc.dat
* bp_223clc.dat
* bp_246clc.dat
* bp_273clc.dat
* bp_276clc.dat
* bp_280clc.dat

on the IGACO site, http://igaco-o3.fmi.fi/ACSO/files/cross_sections.  The files were copied on July 25 2012.
Note that the entries for 282.470 nm and 282.460 nm were swapped around as they are in the wrong order in 
the original files.


Properties
^^^^^^^^^^
.. option:: SetTemperature (double n)
   
   Sets the temperature in Kelvins that will be used in the next calculation of cross-sections.


References
^^^^^^^^^^
Bass, A. M. and Paur, R. J., The ultraviolet cross-sections of ozone, I. Measurements.  Proc. Quadrennial Ozone Symp. Halkidiki, Greece, Reidel, Dordrecht, pp. 606–610, 1984

Bass, A. M. and Paur, R. J., The ultraviolet cross-sections of ozone, II Results and temperature dependence. Proc. Quadrennial Ozone Symp. Halkidiki, Greece, Reidel, Dordrecht, pp. 611–616, 1984
