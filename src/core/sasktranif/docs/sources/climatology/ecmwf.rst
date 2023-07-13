.. _clim_ecmwf:

ECMWF
=====
A climatology that interpolates ECMWF ERA-interim data files. This class is setup to read ECMWF era-interim grib files. The files
follow a defined hierarchy that sorts files into directories baszed upon year and month. The files are of the format,
`<basedirectory>/yyyy/mm/ERA-Int_ml_yyyymmdd.grib` where `yyyy` is the 4 digit year (eg 1970, 2019), `mm` is the 2 digit
month (01-12) and `dd` is the two digit day (01-31). The climatology has been extended to include support for O2 collisional
induced absorption (CIA). ::

   import sasktranif as skif

   lat =   50.0
   lng = -102.0
   h   = 25000.0
   mjd = 53000.1234567
   location =  [lat,lng, h, mjd]
   climate = ISKClimatology('ECMWF')
   climate.UpdateCache(location)
   ok, P = climate.GetParameter( 'SKCLIMATOLOGY_PRESSURE_PA', location);

Supported Species
-----------------

==================================  ================================== ==============   ================
Sasktran Handle                     Description                        Units            Availability
==================================  ================================== ==============   ================
SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3  Air Number density                 molecules/cm3    Always
SKCLIMATOLOGY_PRESSURE_PA           Pressure in                        Pascals          Always
SKCLIMATOLOGY_TEMPERATURE_K         Temperature                        Kelvin           Always
SKCLIMATOLOGY_O2_CM3                Molecular oxygen number density    molecules/cm3    Always
SKCLIMATOLOGY_O2_O2_CM6             O2-O2 square density for CIA       molecule^2/cm6   Always
SKCLIMATOLOGY_N2_CM3                Molecular nitrogen number density  molecules/cm3    Always
==================================  ================================== ==============   ================

Cache Snapshot
--------------
The ECMWF cache internally stores two snapshots of the entire globe that straddle the time of interest. Consequently, there is
no need to invoke :meth:`ISKClimatology.UpdateCache` when position changes.

Python extension
----------------
The ECMWF climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The ECMWF climatology needs ERA-Interim files downladed from ECMWF.  preparation. The ECMWF base directory is set within the active registry. Python users
can change the base directory in the registry entry using the configuration software inside the :ref:`sasktran_core`  module.
The era interim files are assumed to be available every day at 00 UTC, 06 UTC, 12 UTC and 18 UTC. The code interpolates in latitude, longitude, height and time.

Properties
----------
The ECMWF climatology has no proeprties
.. py:module:: ECMWF

