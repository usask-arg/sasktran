.. _optical_property:

Optical Properties
******************

Module: ``sasktran.opticalproperty``

Classes that describe the optical properties of some atmospheric :ref:`species`. These objects are used to 
create :ref:`species` objects.

See our `Useful Shorthands`_ for optical properties that we use a lot.

Base Classes
------------

.. autoclass:: sasktran.OpticalProperty
    :exclude-members: skif_object

.. autoclass:: sasktran.OpticalPropertyConvolved

Useful Shorthands
-----------------

.. autoclass:: sasktran.Rayleigh

.. autoclass:: sasktran.O3DBM

.. autoclass:: sasktran.O3OSIRISRes

.. autoclass:: sasktran.NO2Vandaele1998

.. autoclass:: sasktran.NO2OSIRISRes

.. autoclass:: sasktran.HITRANChemical

.. autoclass:: sasktran.MieAerosol

.. autoclass:: sasktran.SimpleRayleigh

.. autoclass:: sasktran.BaumIceCrystal

.. autoclass:: sasktran.SO2Vandaele2009

.. autoclass:: sasktran.O2O2Fally2000

.. autoclass:: sasktran.O2O2HITRAN2016

.. autoclass:: sasktran.O2O2Thalman2013

.. autoclass:: sasktran.UserDefinedAbsorptionPressure

.. autoclass:: sasktran.UserDefinedScatterConstantHeight

See Also
--------
.. include:: descriptions/species.desc
.. include:: descriptions/climatology.desc
.. include:: descriptions/atmosphere.desc
