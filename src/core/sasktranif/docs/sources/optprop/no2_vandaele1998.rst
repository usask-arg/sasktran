
.. _optical_no2vandaele:

NO2_VANDAELE1998
================
Calculates the absorption cross section of NO2 molecules from 230 nm to 1000 nm at 220 K to 294 K following Vandaele et al. 1998. 


Example
^^^^^^^
::

   optprop = ISKOpticalProperty('NO2_VANDAELE1998')
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
Vandaele A.C., C. Hermans, P.C. Simon, M. Carleer, R. Colin, S. Fally, M.F. Mérienne, A. Jenouvrier, 
and B. Coquart, Measurements of the NO2 absorption cross-section from 42000 cm-1 to 10000 cm-1 
(238-1000 nm) at 220 K and 294 K, JSQRT, 59, 171-184 (1998)

