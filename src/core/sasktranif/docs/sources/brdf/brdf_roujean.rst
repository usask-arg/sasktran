.. _brdf_roujean:

Roujean
=======
The BRDF object that implements the algorithm developed by Roujean et al. 1992. The algorithm is suitable for surface reflectance from a variety of terrains and is based
upon three parameters, k0, k1 and k2 which the user must know a-priori or derive from measurements. The reader is referred to their paper for complete details.

Example
-------
::

   import sasktranif.sasktranif as skif
   import math
   
   brdf = skif.ISKBrdf('Roujean')
   brdf.SetProperty('SetPredefinedParameters', 8);
   mjd            = 52393.3792987115;
   location       = [0.0, 0.0, 25000.0, mjd];
   mu_in          = math.cos( math.radians( 30.0 ) )
   mu_out         = math.cos( math.radians( 60.0 ) )
   cosdphi        = math.cos( mat.radians( 150.0 ) )
   [ok,brdfvalue] = brdf.BRDF( 600.0, location, 0.6, 0.7, -0.8)


Properties
----------
..  py:module:: Roujean

BRDFParameters
^^^^^^^^^^^^^^
.. py:function:: BRDFParameters( array param)

    A 3 element array containing the 3 parameters of the Roujean distribution. For examples see Table 1 of Roujean et al. 1992.

    ================ ===========
      n              Setting
    ================ ===========
      parameters[0]  k0
      parameters[1]  k1
      parameters[2]  k2
    ================ ===========

SetPredefinedParameters
^^^^^^^^^^^^^^^^^^^^^^^
.. py:function:: SetPredefinedParameters( code: int )

   An integer value that loads the Roujean BRDF parameters from an internal table of values. The table replicates the values published in Table 1 of Roujean et al. 1992.
   The table provides a visible setting (580 nm to 680 nm) and a near infra-red (NIR) setting (730 nm to 1100 nm). Pass in the appropriate code
   from the table below to select the appropriate BRDF setting

    ================= =============  ===========
    Terrain            Visible Code    NIR Code
    ================= =============  ===========
    PLOWED FIELD         1             21
    ANNUAL GRASS         2             22
    HARD WHEAT           3             23
    STEPPE               4             24
    CORN                 5             25
    ORCHARD GRASS        6             26
    IRRIGATED WHEAT      7             27
    PINEFOREST           8             28
    DECIDUOUS FOREST     9             29
    SOYBEAN             10             30 
    GRASS LAWN          11             31 
    ================= =============  ===========

References
^^^^^^^^^^
**Roujean, J.L.**, M. Leroy, and P.Y. Deschamps,  "A bidirectional reflectance model of the Earth's surface for the correction of remote sensing data", *J. Geophys. Res.*, **97**, D18, 20455-20468, (1992), `doi:10.1029/92JD01411 <https://doi.org/10.1029/92JD01411>`_.
