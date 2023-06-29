.. _brdf_hapke:

Hapke
=====
The Hapke BRDF depends on three parameters and is usually used
independently, not in a linear combination with other kernels. It is
similar to the Rahman kernel but the hotspot, a peak in
reflectivity back in the direction of the incident light due to the
absense of shadowing effects, is treated differently.

This is one of eight BRDF kernels listed by Spurr in
his 2002 paper.

Equation A.18 from Spurr is implemented with the following changes:

    * cos(phi) is replaced with -cos(phi)
    * the reflectance is divided by pi

Example
-------
::

   import sasktranif.sasktranif as skif
   import math

   brdf = skif.ISKBrdf('HAPKE')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)


References
-----------
**Hapke B.**, "Theory of reflectance and emittance spectroscopy". Cambridge University Press, 1993.
**Robert J. D. Spurr**, "A new approach to the retrieval of surface properties 	from earthshine measurements", *Journal of Quantitative Spectroscopy & Radiative Transfer*,
vol. 83, pp. 15-46, Oct. 9, 2002. doi:10.1016/S0022-4073(02)00283-2 	url: `http://www.sciencedirect.com/science/article/pii/S0022407302002832 <http://www.sciencedirect.com/science/article/pii/S002240730200283>`_

