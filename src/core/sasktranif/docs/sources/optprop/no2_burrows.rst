.. _optical_no2burrows:

NO2_BURROWS
===========
Tabulates the absorption cross section of NO2 molecules at medium resolution (0.2-0.33 nm)
from 230 nm to 795 nm at 4 temperatures. The cross-sections were measured by Burrows et al. 
with the GOME instrument before flight. 

Spectral Resolution
^^^^^^^^^^^^^^^^^^^

* 231-307 nm, 0.20 nm 
* 307-316 nm, 0.20 nm
* 311-405 nm, 0.17 nm
* 405-611 nm, 0.29 nm
* 595-794 nm, 0.33 nm

Temperature Range
^^^^^^^^^^^^^^^^^
The cross-sections were measured at 4 temperatures covering normal stratospheric and
tropospheric ranges. The paper does not   provide any advice on how to interpolate in temperature. 
We use a standard (linear interpolation) technique provided by the base class. The spectra were 
measured at the following four temperatures.

1. 221K
2. 241K
3. 273K
4. 293K

Example
^^^^^^^
::

   optprop = ISKOpticalProperty('NO2_BURROWS')
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
 J. P. Burrows, A. Dehn, B. Deters, S. Himmelmann, A. Richter, S. Voigt, and J. Orphal: 
 "Atmospheric Remote-Sensing Reference Data from GOME: 1. Temperature-Dependent Absorption
 Cross Sections of NO2 in the 231-794 nm Range", Journal of Quantitative Spectroscopy
 and Radiative Transfer 60, 1025-1031, 1998.
 
Paper Abstract
^^^^^^^^^^^^^^
Absorption cross-sections of NO2 between 231-794 nm have been measured in the 221-293K temperature 
range, using the global ozone monitoring experiment (GOME) flight model (FM) satellite spectrometer.
The spectra have a resolution of about 0.2 nm below 400 nm and of about 0.3 nm above 400 nm. These
are the first reference spectra of NO2 covering at the same time the entire UV-visible-NIR 
spectral range and a broad range of relevant atmospheric temperatures. The new absorption
cross-sections are important as accurate reference data for atmospheric remote-sensing of NO2 and 
other minor trace gases.
