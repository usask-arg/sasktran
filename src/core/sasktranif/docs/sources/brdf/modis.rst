.. _brdf_modis:

Modis
=====
This is the BRDF model used by MODIS for their BRDF/albedo product. It
consists of a linear combination of an isotropic (lambertian) term, a
volume scattering term (Ross-thick), and a geometric scattering term
(Li-Sparse-Reciprocal).

Example
-------
::

   import sasktranif.sasktranif as skif
   import math

   brdf = skif.ISKBrdf('MODIS')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)


Properties
-----------
..  py:module:: BRDF_MODIS

BRDFParameters
^^^^^^^^^^^^^^
..  py:function:: BRDFParameters( array param)

    A 3 element array containing the 3 parameters of the MODIS distribution.

    ================ ========================================================
      n              Setting
    ================ ========================================================
      parameters[0]  isotropic (lambertian) component
      parameters[1]  linear weight of volume scattering kernel (Ross-Thick)
      parameters[2]  linear weight of geometric scattering kernel (Li-Sparse)
    ================ ========================================================


References
-----------
**Robert J. D. Spurr**, "A new approach to the retrieval of surface properties 	from earthshine measurements", *Journal of Quantitative Spectroscopy & Radiative Transfer*,
vol. 83, pp. 15-46, Oct. 9, 2002. doi:10.1016/S0022-4073(02)00283-2 	url: `http://www.sciencedirect.com/science/article/pii/S0022407302002832 <http://www.sciencedirect.com/science/article/pii/S002240730200283>`_

