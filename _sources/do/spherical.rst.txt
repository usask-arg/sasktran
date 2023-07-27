.. _spherical:

Spherical Approximations
************************
While SK-DO is based on the discrete ordinates technique which is designed for a plane-parallel atmosphere, SK-DO
contains several approximations to handle various levels of sphericity.

Pseudo-Spherical Approximation
------------------------------
The pseudo-spherical approximation is a common method where the incoming solar beam is treated in a spherical
atmosphere.  These transmittances are then used to initialize the discrete ordinates technique.  The pseudo-spherical
approximation is enabled by default, and there is usually no reason to disable it as it does not have any computational
overhead involved in using it.  However, it can be disabled with

.. code-block:: python

    engine.options['usepseudospherical'] = 0

if desired.

Line of Sight Spherical Corrections
-----------------------------------
SK-DO also contains the capability to include spherical effects along the line of sight.  In this mode the
line of sight ray is traced in a fully spherical atmosphere.  Single scattering source terms are calculated exactly
inside the spherical atmosphere.  Multiple scatter source terms are calculated using the discrete ordinates method
at multiple solar zenith angles along the line of sight.  These sources are then integrated along the line of sight
to get the total radiance.  This method is generally very accurate for almost every nadir viewing application, and
can be accurate for some limb-viewing applications at low tangent altitudes.

The technique can be enabled with

.. code-block:: python

    engine.viewing_mode = 'spherical'

and the main parameter controlling the accuracy of the calculation is


.. code-block:: python

    engine.num_spherical_sza = 2 # default

which controls the number of plane-parallel solutions calculated along the line of sight.
