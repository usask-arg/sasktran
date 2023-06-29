.. _brdf_ross_thin_kernel:

Ross Thin Kernel
================
The Ross-thin kernel is a volume-scattering kernels that
models the reflection of surface facets (leaves). It is purely geometric
and depends on no additional parameters.

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

   brdf = skif.ISKBrdf('ROSS_THIN_KERNEL')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)



References
-----------
**W. Wanner** and A. H. Strahler, "On the derivation of kernels for kernel-driven models of bidirectional reflectance," *Journal of Geophysical Research*, vol. 100, no. d10, pp. 21077-21089, Oct. 20, 1995.
doi: 10.1029/95JD02371, url: `http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract <http://onlinelibrary.wiley.com/doi/10.1029/95JD02371/abstract>`_

**Robert J. D. Spurr**, "A new approach to the retrieval of surface properties from earthshine measurements", *Journal of Quantitative Spectroscopy & Radiative Transfer*,
vol. 83, pp. 15-46, Oct. 9, 2002. doi:10.1016/S0022-4073(02)00283-2 	url: `http://www.sciencedirect.com/science/article/pii/S0022407302002832 <http://www.sciencedirect.com/science/article/pii/S002240730200283>`_
