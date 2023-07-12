
.. _optical_no2osirisres:

NO2_OSIRISRES
=============
Calculates the absorption cross section of NO2 molecules from 230 nm to 795 nm and 221K to 293K.
The cross-sections have been reduced to the resolution of OSIRIS and these cross-sections have 
been used in the OSIRIS level 2 MART retrievals.


Example
^^^^^^^
::

   optprop = ISKOpticalProperty('NO2_OSIRISRES')
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
