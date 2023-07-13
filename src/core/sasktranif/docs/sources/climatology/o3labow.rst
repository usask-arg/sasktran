.. _clim_o3labow:

O3LABOW
=======
The Labow volume mixing ratio climatology of ozone. This is climatology of the volume mixing ratio
from 0 km to 60 km in 1 km increments for 18 latitudes steps of 10 degrees from -85 to +85 
for each month of the year. The climatology seems to be unpublished as all references 
refer to an unpublished piece of work by McPeters and Labow in 2002 or 2003.  
We have extended the climatology above 60 km using the differences in scale height 
between the neutrals (7 km scale height ) and ozone (4.5 km scale height) to calculate 
a new scale height to extrapolate the signal at 60 km. We have also extended the volume 
mixing ratio so it can be converted to an ozone density.  This requires an atmospheric 
number density climatology for which the MSIS90 climatology is used::

   import sasktranif as skif

   location = [50.0, -102.0, 12345.0, 53000.327869823]
   climate = skif.ISKClimatology('O3LABOW')
   climate.UpdateCache(location)
   ok, O3 = climate.GetParameter( 'SKCLIMATOLOGY_O3_CM3', location);

Supported Species
-----------------

* SKCLIMATOLOGY_O3_VMR
* SKCLIMATOLOGY_O3_CM3

Cache Snapshot
--------------
The O3LABOW climatology caches the entire LABOW ozone database which provides O3 as a function of geolocation and day of the year

Python extension
----------------
The O3LABOW climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The O3LABOW climatology needs no external preparation.

Properties
----------
The O3LABOW climatology has no special properties.

