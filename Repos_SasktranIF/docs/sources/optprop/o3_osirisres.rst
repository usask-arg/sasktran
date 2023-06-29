.. _optical_o3osirisres:

O3_OSIRISRES
============
These are the cross-sections used in the OSIRIS level 2 analysis for Saskmart V5.07. 
This table is based upon the cross-sections of Bogumil Orphal and Burrows and has been 
reduced to the nominal resolution of the OSIRIS instrument (approx. 1.0 nm). 
The cross-sections range from 270nm to 820 nm in 0.1 nm steps.

Example
-------
::

   optprop = ISKOpticalProperty('O3_OSIRISRES')
   msis = ISKClimatology('MSIS90');
   mjd = 52393.3792987115;
   location = [0.0, 0.0, 25000.0, mjd];
   optprop.SetAtmosphericState( msis)
   optprop.SetLocation(location)
   [ok, abs, ext,sca] = optprop.CalculateCrossSections( 1.0E7/600.0, location );

