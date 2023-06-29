
.. _optical_o3sciabogumilv4:

O3_SCIABOGUMILV4
================
Version 4 is a reanalysis of the actual Bogumil data used to make :ref:`optical_o3sciabogumilv3`.
A manuscript was being prepared for publication (as of July 2012) on the revised cross-sections 
(private communication with Mark Weber). The Sciamachy people had some problems in the 
Huggins band with the Bogumil data as total ozone retrievals were biased 5%.

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
