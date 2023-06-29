.. _brdf_lambertian:


Lambertian 
==========
The Lambertian BRDF object scatters/reflects downward flux equally into all outward directions. The Lambertian object uses its albedo property to 
ratio the upward flux t the downward flux and converts the upward flux to the outbound radiance through a factor of pi. The object can be used as a 
replacement for the scalar *albedo* implement in radiative transfer engines although you must be careful to eliminate any extra divisions by pi.

Example
-------
::

   import sasktranif.sasktranif as skif
   import math
   
   brdf = skif.ISKBrdf('Lambertian')
   brdf.SetProperty('Albedo', 0.3);
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   mu_in          = math.cos( math.radians( 30.0 ) )
   mu_out         = math.cos( math.radians( 60.0 ) )
   cosdphi        = math.cos( mat.radians( 150.0 ) )
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)


Properties
----------

..  py:module:: Lambertian

Albedo
^^^^^^
.. py:function:: Albedo(double n)
   
   Sets the albedo of the surface. This should be a value between 0 and 1. The value is applied to all wavelengths

