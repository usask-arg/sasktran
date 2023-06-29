.. _brdf_rahman:

Rahman
======
This is a semi-empirical model with 3 parameters, designed to model the reflectance of arbitrary natural surfaces in the visible
and near-infrared bands.

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

   brdf = skif.ISKBrdf('RAHMAN')
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)

Properties
-----------
..  py:module:: BRDF_RAHMAN

BRDFParameters
^^^^^^^^^^^^^^
..  py:function:: BRDFParameters( array param)

    A 3 element array containing the 3 parameters of the Rahman distribution.

    ================ ========================================================
      n              Setting
    ================ ========================================================
      parameters[0]  :math:`\rho`
      parameters[1]  :math:`\theta`
      parameters[2]  :math:`k`
    ================ ========================================================

References
-----------
**Rahman H**, Pinty B, Verstraete M. "Coupled surface atmosphere reflectance (CSAR) model: 2. Semi-empirical surface model usable with NOAA AVHRR data.", *J Geophys Res* 1993; 98:20,791–801.

**Robert J. D. Spurr**, "A new approach to the retrieval of surface properties 	from earthshine measurements", *Journal of Quantitative Spectroscopy & Radiative Transfer*,
vol. 83, pp. 15-46, Oct. 9, 2002. doi:10.1016/S0022-4073(02)00283-2 	url: `http://www.sciencedirect.com/science/article/pii/S0022407302002832 <http://www.sciencedirect.com/science/article/pii/S002240730200283>`_

