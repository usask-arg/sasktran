.. _hitran:

HITRAN Optical Properties
*************************

SASKTRAN-TIR supports HITRAN optical properties.
However, if the :class:`sasktran.HITRANChemical` class is used, it is necessary to set the ``use_cache`` option to be false.
For example,

::

    import sasktran as sk
    hitran_o3 = sk.HITRANChemical('O3', use_cache=False)

If this option is not set and an engine object is used to calculate radiance more than once, all radiances except the first will be incorrect.
SASKTRAN-TIR includes a :py:class:`sasktran.tir.opticalproperty.HITRANChemicalTIR` class which automatically sets the use_cache parameter to be false.
This class also makes it possible to set the micro window margin during initialization.

.. autosummary::

    sasktran.tir.opticalproperty.HITRANChemicalTIR

.. autoclass:: sasktran.tir.opticalproperty.HITRANChemicalTIR
