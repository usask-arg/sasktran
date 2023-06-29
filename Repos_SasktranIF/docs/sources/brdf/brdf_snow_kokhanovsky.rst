.. _brdf_kokhanovsky:

Kokhanovsky Snow
=================
The BRDF object that implements the algorithm for the BRDF of snow developed by Kokhanovsky and Breon in 2012. The algorithm is suitable for surface reflectance from snow covered ground.
The model has two free parameters, **L** and **M** which are either provided a-priori by the user or derived from measurements. The reader is referred to their paper for complete details. They
present results taken over Greenland in May 2006 (see their Figure 1) and use values based upon relecatnce at 1020 nm of L = 3.6 mm (3.6E06 nm) and M = 5.5E-08

Example
-------
::

   import sasktranif.sasktranif as skif
   import math

   brdf = skif.ISKBrdf('SNOW_KOKHANOVSKY2012')
   brdf.SetProperty('BRDFParameters', [3.6E06, 5.5E-08] );
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   mu_in          = math.cos( math.radians( 30.0 ) )
   mu_out         = math.cos( math.radians( 60.0 ) )
   cosdphi        = math.cos( mat.radians( 150.0 ) )
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, mu_in, mu_out, cosdphi)


Properties
----------

..  py:module:: SNOW_KOKHANOVSKY2012

BRDFParameters
^^^^^^^^^^^^^^
.. py:function:: BRDFParameters( array param )

    A 2 element array containing the 2 free parameters of the Kokhanovsky model. The Kokhanovsky 2012 paper referes to these parameters as L and M.

    ================ ====================
     Parameter         Description
    ================ ====================
     param[0]          L in nanonmeters
     param[1]          M (dimensionless)

    ================ ====================


References
^^^^^^^^^^
 - A. A. Kokhanovsky and F. M. Breon,
   "Validation of an Analytical Snow BRDF Model Using PARASOL Multi-Angular and Multispectral Observations,"
   IEEE Geoscience and Remote Sensing Letters, vol. 9, no. 5, pp. 928-932, Sept. 2012.
   doi: 10.1109/LGRS.2012.2185775

 - A. Kokhanovsky, V. V. Rozanov, T. Aoki, D. Odermatt, C. Brockmann, O. Kriger, M. Bouvet, M. Drusch, and M. Hori
   "Sizing snow grains using backscattered solar light"
   International Journal Of Remote Sensing Vol. 32 , Iss. 22,2011

