.. _clim_osirisl2_o3rtmodel_v507:

OSIRISL2_O3RTMODEL_V507
=======================
A climatology that provides access to individual scans from the Odin-OSIRIS ozone V507 data product.
It is only useful if you need the ozone profile for a specific scan of the Odin-OSIRIS instrument,
implying you have access to techniques to convert scan number to location and modified juliand date.

The object only uses the modified julian date to index the profile from the OSIRIS database
and only works for periods coincident with the Odin-OSIRIS scans. It cannot 
be used as a general purpose ozone climatology.

This climatology will only work if your local system has access to the 
Odin-OSIRIS Level 2, version 5.07 data products

::

   mjd = 52393.3792987115;
   climate = ISKClimatology('OSIRISL2_O3RTMODEL_V507')
   climate.UpdateCache( [0,0,0,mjd])
   [ok, value] = climate.GetParameter( SKCLIMATOLOGY_O3_CM3(), [0.0, 0.0, 25000, mjd]);

Supported Species
-----------------

* SKCLIMATOLOGY_O3_CM3

