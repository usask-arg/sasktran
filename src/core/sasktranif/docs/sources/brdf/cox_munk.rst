.. _brdf_cox_munk:

Cox Munk
========
Implements the Cox Munk BRDF.	This BRDF is based on the work of Cox and Munk who
analyzed the surface of water as a function of wind using aerial
photographs in the 1950s. Their work includes dependence on wind
direction and a Gram-Charlier correction to the Gaussian distribution
of water facet orientations, but this BRDF neglects these corrections
and assumes a symmetric distribution. This BRDF depends on two
parameters, wind speed and the index of refraction of water.

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

   brdf = skif.ISKBrdf('COX_MUNK')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)


Properties
-----------
..  py:module:: COX_MUNK

BRDFParameters
^^^^^^^^^^^^^^
..  py:function:: BRDFParameters( array param)

    A 2 element array containing the wind speed and water refractive index required for the Cox Munk distribution

    ================ ===========
      n              Setting
    ================ ===========
      parameters[0]  wind speed
      parameters[1]  water refractive index
    ================ ===========

References
-----------
**Cox C**, Munk W. "The measurement of the roughness of the sea surface from photographs of the Sun glitter". *J Opt Soc Ann* 1954; 44:838–50.

**Robert J. D. Spurr**, "A new approach to the retrieval of surface properties 	from earthshine measurements", *Journal of Quantitative Spectroscopy & Radiative Transfer*,
vol. 83, pp. 15-46, Oct. 9, 2002. doi:10.1016/S0022-4073(02)00283-2 	url: `http://www.sciencedirect.com/science/article/pii/S0022407302002832 <http://www.sciencedirect.com/science/article/pii/S002240730200283>`_

