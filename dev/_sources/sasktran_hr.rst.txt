.. _sasktran_hr:

SASKTRAN-HR
***********
Sasktran HR is the recommended default radiative transfer model that is primarily designed for limb viewing applications.


Polarization in SASKTRAN-HR
---------------------------
SASKTRAN-HR is capable of calculation the full Stokes vector.  This option is controlled by setting
two different things.  The internal model calculation mode::

    import sasktran as sk

    engine = sk.EngineHR()

    engine.polarization = 'vector'

In this mode all calculations will be done considering polarization, however the output will still be a single intensity element.
To obtain the full Stokes vector you must pass the `full_stokes_vector=True` when calculating the radiance::

    stokes = engine.calculate_radiance(full_stokes_vector=True)

See :ref:`hr_examples` for examples of usage.

Commonly Used Options
---------------------

.. autoattribute:: sasktran.EngineHR.num_orders_of_scatter
.. autoattribute:: sasktran.EngineHR.grid_spacing
.. autoattribute:: sasktran.EngineHR.surface_elevation
.. autoattribute:: sasktran.EngineHR.top_of_atmosphere_altitude
.. autoattribute:: sasktran.EngineHR.num_threads
.. autoattribute:: sasktran.EngineHR.polarization


Advanced Configuration
----------------------

.. autoattribute:: sasktran.EngineHR.include_emissions
.. autoattribute:: sasktran.EngineHR.atmosphere_dimensions
.. autoattribute:: sasktran.EngineHR.refraction


Advanced Accuracy/Performance
-----------------------------
.. autoattribute:: sasktran.EngineHR.num_diffuse_profiles

Experimental Settings
---------------------
.. autoattribute:: sasktran.EngineHR.disable_los_single_scattering
