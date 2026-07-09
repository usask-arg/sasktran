.. _sasktran_tir:

**************
The TIR Engine
**************
The ``TIR`` engine solves the radiative transfer equation in the thermal infra-region on a true spherical Earth.

.. note::

  The ``TIR`` model should be used with caution for:

  * Extremely optical thick conditions like clouds. The code does not handle the large optical depths very well.
  * Nadir geometries. They are supported but the surface emission is a simple black body emitter with scalar emissivity.

SASKTRAN-TIR is a radiative transfer model designed for limb viewing applications in the thermal infrared spectral region.
It is recommended to use this model for applications where the radiative source at the target wavelengths is dominated by thermal emissions, meaning that the contribution from solar scattering can be ignored.
Typically, this occurs for wavelengths larger than 5000 nm.
Nadir viewing is also supported where the surface is modeled as a black body with a constant user-specified emissivity.

Features
    - Scalar radiative transfer computations for Limb and Nadir lines of sight
    - Scalar weighting functions computations for species' number density
    - Supports refractive ray tracing
    - Supports varying temperature and atmospheric properties along the line of sight

Limitations
    - Basic black body surface model used in Nadir geometries
    - Currently this is no support for non-LTE effects, continuum absorption, and line-mixing

Weighting Functions
===================

SASKTRAN-TIR supports the same weighting function interface as SASKTRAN-HR, and supports multiple-species weighting
functions with a similar interface to SASKTRAN-Disco.

.. code-block:: python

    # Standard atmosphere configuration
    atmosphere = sk.Atmosphere()
    atmosphere['o3'] = sk.Species(sk.HITRANChemical('O3'), sk.Labow())
    atmosphere['no2'] = sk.Species(sk.HITRANChemical('NO2'), sk.Pratmo())

.. code-block:: python

    atmosphere.wf_species = 'o3'            # single species syntax

.. code-block:: python

    atmosphere.wf_species = ['no2', 'o3']   # multiple species syntax

.. code-block:: python

    # calculate radiance and weighting functions
    rad, wf = engine.calculate_radiance()

To compute temperature weighting functions:

.. code-block:: python

    engine.do_temperature_wf = True

Shape of ``wf``
---------------

Length definitions for this section
    ``lwav = len(engine.wavelengths)``
    ``llos = len(geometry.lines_of_sight)``
    ``lwfs = len(atmosphere.wf_species)``
    ``lwfa = len(weighting function altitudes)``

The shape of ``wf`` depends on the syntax you used to specify the weighting functions. If you used the single species
syntax to specify ``atmosphere.wf_species`` then the shape of ``wf`` is ``(lwav, llos, lwfa)``. If you used the
multiple species syntax then the shape is ``(lwav, llos, lwfs, lwfa)``.

Engine Options
==============



Weighting Functions
-------------------

.. autoattribute:: sasktran.tir.engine.EngineTIR.wf_heights
.. autoattribute:: sasktran.tir.engine.EngineTIR.wf_widths
.. autoattribute:: sasktran.tir.engine.EngineTIR.do_vmr_wf
.. autoattribute:: sasktran.tir.engine.EngineTIR.do_temperature_wf


Integration
-----------

.. autoattribute:: sasktran.tir.engine.EngineTIR.adaptive_integration
.. autoattribute:: sasktran.tir.engine.EngineTIR.max_optical_depth_of_cell
.. autoattribute:: sasktran.tir.engine.EngineTIR.min_extinction_ratio_of_cell
.. autoattribute:: sasktran.tir.engine.EngineTIR.linear_extinction


Geometry and Atmosphere
-----------------------

.. autoattribute:: sasktran.tir.engine.EngineTIR.ground_emissivity
.. autoattribute:: sasktran.tir.engine.EngineTIR.refraction
.. autoattribute:: sasktran.tir.engine.EngineTIR.surface_elevation
.. autoattribute:: sasktran.tir.engine.EngineTIR.top_of_atmosphere_altitude
.. autoattribute:: sasktran.tir.engine.EngineTIR.grid_spacing
.. autoattribute:: sasktran.tir.engine.EngineTIR.atmosphere_dimensions


Performance
-----------

.. autoattribute:: sasktran.tir.engine.EngineTIR.num_threads
.. autoattribute:: sasktran.tir.engine.EngineTIR.use_cached_cross_sections
