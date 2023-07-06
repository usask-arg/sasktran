.. _emission:

Emission
********

Module: ``sasktran.emission``

Classes that describe emissions occuring within the atmosphere and are typically thermal or photochemical. The overall purpose
of the Emissions classes is to provide the volume emission rate as a function of wavelength at any location in the atmosphere.

.. autosummary::

    sasktran.EmissionTable
    sasktran.EmissionThermal
    sasktran.HITRANPhotoChemical
    sasktran.HITRANPhotoChemical_O2_ABand
    sasktran.HITRANPhotoChemical_O2_SingletDelta
    sasktran.ABandEmission

Emission Base Class
--------------------

.. autoclass:: sasktran.Emission
    :exclude-members: skif_object

EmissionTable
--------------
.. autoclass:: sasktran.EmissionTable

EmissionThermal
----------------
.. autoclass:: sasktran.EmissionThermal

HITRANPhotoChemical
--------------------
.. autoclass:: sasktran.HITRANPhotoChemical

HITRANPhotoChemical_O2_ABand
----------------------------
.. autoclass:: sasktran.HITRANPhotoChemical_O2_ABand

HITRANPhotoChemical_O2_SingletDelta
-------------------------------------
.. autoclass:: sasktran.HITRANPhotoChemical_O2_SingletDelta

ABandEmission
-------------

.. autoclass:: sasktran.ABandEmission

See Also
--------
.. include:: descriptions/atmosphere.desc
