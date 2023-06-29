.. _brdf_roujean_kernel:

Roujean Kernel
==============
This kernel is derived considering a random arrangement of rectangular blocks on a flat surface.


Example
-------
::

   import sasktranif.sasktranif as skif
   import math

   brdf = skif.ISKBrdf('ROUJEAN_KERNEL')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)


Properties
----------
The Roujean Kernel has no properties


References
^^^^^^^^^^
**Roujean, J.L.**, M. Leroy, and P.Y. Deschamps,  "A bidirectional reflectance model of the Earth's surface for the correction of remote sensing data", *J. Geophys. Res.*, **97**, D18, 20455-20468, (1992), `doi:10.1029/92JD01411 <https://doi.org/10.1029/92JD01411>`_.
