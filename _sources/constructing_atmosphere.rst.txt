.. _constructing_atmosphere:

Specifying the Atmosphere
=========================
One of the most important and complicated steps in constructing a radiative transfer calculation is
building a representation of the atmosphere.   To perform a radiative transfer calculation, a few things are needed:

* A set of constituent molecules which make up the atmosphere.  For each of these molecules we need to know the amount
  of the molecule, specified with :py:class:`sasktran.Climatology`, and the scattering/absorption cross sections,
  specified as a :py:class:`sasktran.OpticalProperty`.  The combination of these two objects is a
  :py:class:`sasktran.Species`.

* A :py:class:`sasktran.Climatology` that provides the pressure and temperature profiles of the atmosphere, this is
  often referred to as the atmospheric state.  These
  profiles are used when calculating cross sections of the atmospheric constitutents.  By default this is
  :py:class:`sasktran.MSIS90`

* A :py:class:`sasktran.BRDF` object which describes the reflectance of the surface of the Earth.  This could be as
  simple as a constant number between 0--1 indicating a Lambertian surface.

All of these things are combined together with the :py:class:`sasktran.Atmosphere` class which is the input
to the radiative transfer model.

Climatologies
-------------

A :py:class:`sasktran.Climatology` object is essentially a lookup table of profiles
where the tables keys are quantities name. A grid (altitude and potentially geographic)
definition is also associated.

Users can defined their own user defined climatologies using the
:py:class:`sasktran.ClimatologyUserDefined` class. Some climatologies that we use a lot are built
into SASKTRAN (see `Useful Climatology Shorthands`).

Available Climatologies
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::

    sasktran.ClimatologyUserDefined
    sasktran.ClimatologyUserDefined2D
    sasktran.ClimatologyUserDefined3D
    sasktran.Labow
    sasktran.Pratmo
    sasktran.MSIS90
    sasktran.GloSSAC


Optical Properties
------------------


Available Optical Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::

    sasktran.Rayleigh
    sasktran.O3DBM
    sasktran.O3OSIRISRes
    sasktran.NO2Vandaele1998
    sasktran.NO2OSIRISRes
    sasktran.HITRANChemical
    sasktran.MieAerosol
    sasktran.OpticalPropertyConvolved
    sasktran.UserDefinedAbsorptionPressure
    sasktran.BaumIceCrystal
    sasktran.SO2Vandaele2009


The Species Object
------------------
The :py:class:`Species <sasktran.Species>` class is a combination of a :py:class:`Climatology <sasktran.Climatology>`
and :py:class:`OpticalProperty <sasktran.OpticalProperty>`,
linking them together.  In most cases creating the Species is as simple as::

    import sasktran as sk

    ozone_species = sk.Species(climatology=sk.Labow(), optical_property=sk.O3DBM())

The Species object can then be added to the Atmosphere object through::

    atmosphere = sk.Atmosphere()

    atmosphere['ozone'] = ozone_species

However, there can be issues with this simple approach if the :py:class:`Climatology <sasktran.Climatology>` used
supports multiple species.  For example, trying::

    msis = sk.MSIS90(True)

    air = sk.Species(climatology=msis, optical_property=sk.Rayleigh())

will raise the error `ValueError: Could not automatically infer the species key, use species= one of ['SKCLIMATOLOGY_PRESSURE_PA', 'SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3', 'SKCLIMATOLOGY_TEMPERATURE_K', 'SKCLIMATOLOGY_O2_CM3']`
this is because the `msis` object supports multiple atmospheric constituents, and which one was intended to be used
could not automatically be determined.  Here we need to do::

    air = sk.Species(climatology=msis, optical_property=sk.Rayleigh(), species='SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3')

this species can then be added to the Atmosphere object in the same way::

    atmosphere['air'] = air


The Atmospheric State
---------------------
The atmospheric state is a special :py:class:`Climatology <sasktran.Climatology>` which supports both pressure
and temperature.  These pressure and temperature profiles are used to calculate the cross sections of the
atmospheric constituents, they are not used to define Rayleigh scattering within the atmosphere, that is done with
a separate Climatology and OpticalProperty.

By default the atmospheric state is set to the :py:class:`MSIS90 <sasktran.MSIS90>` climatology, however it can
be changed through::

    import sasktran as sk

    atmosphere = sk.Atmosphere()

    atmosphere.atmospheric_state = sk.MSIS90() # Default

The atmospheric state climatology must support the two keys:

* 'SKCLIMATOLOGY_PRESSURE_PA'
* 'SKCLIMATOLOGY_TEMPERATURE_K'

Thus if you are trying to set the atmospheric state as a :py:class:`UserDefinedClimatology <sasktran.UserDefinedClimatology>`
it is necessary that both of these keys are supported.

BRDF
----
The :py:class:`BRDF <sasktran.BRDF>` object controls the surface reflectance inside the radiative transfer model.
This could be as simple as a scalar which represents an isotropically reflecting Lambertian surface::

    import sasktran as sk

    atmosphere = sk.Atmosphere()

    atmosphere.brdf = 0.3

There are also many available special BRDF's::

    atmosphere.brdf = sk.Kokhanovsky()


Available BRDF's
^^^^^^^^^^^^^^^^

.. autosummary::

    sasktran.Lambertian
    sasktran.Kokhanovsky
    sasktran.Roujean
    sasktran.CoxMunk
    sasktran.Rahman
    sasktran.Hapke
    sasktran.MODIS

Kernel Based BRDF's
^^^^^^^^^^^^^^^^^^^
For advanced users it is possible to combine multiple BRDF Kernels together to create a hybrid BRDF object.
The main class which enables this is the :py:class:`LinearCombination <sasktran.LinearCombination>` object.::

    import sasktran as sk

    brdf = sk.LinearCombination()

    brdf.kernels = [sk.Lambertian(0.3), sk.RoujeanKernel()]


Available Kernels
"""""""""""""""""

.. autosummary::

    sasktran.RoujeanKernel
    sasktran.LiKernel
    sasktran.LiSparseKernel
    sasktran.LiSparseReciprocalKernel
    sasktran.LiDenseKernel
    sasktran.RossThinKernel
    sasktran.RossThickKernel
