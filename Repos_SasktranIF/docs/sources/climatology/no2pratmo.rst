.. _clim_no2pratmo:

NO2PRATMO
=========
A climatology for the NO2 number density based upon the Prather photo-chemical box model. This climatological model 
is based upon monthly mean solutions to the box model for a range of latitudes, local solar times and altitudes::

    import sasktranif as skif

    location = [52.0, -102.0, 0.0, 57005.0]
    climate = skif.ISKClimatology('NO2PRATMO')
    climate.UpdateCache(location)
    ok, NO2 = climate.GetParameter( 'SKCLIMATOLOGY_NO2_CM3', location );

Supported Species
-----------------
* SKCLIMATOLOGY_NO2_CM3


Cache Snapshot
--------------
The NO2PRATMO climatology caches the entire NO2 database climatology which provides NO2 as a function of geolocation
and day of the year.

Python extension
----------------
The NO2PRATMO climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The NO2 Pratmo climatology requires the user to download the NO2 pratmo database file. Methods to do this are provided
and described :ref:`sasktran_core`

Properties
----------
The NO2PRATMO cclimatology has no properties


