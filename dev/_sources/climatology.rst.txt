.. _climatology:

Climatology
***********

Module: ``sasktran.climatology``

A :py:class:`sasktran.Climatology` object is essentially a lookup table of profiles
where the tables keys are quantities name. A grid (altitude and potentially geographic) 
definition is also associated.

Users can defined their own user defined climatologies using the 
:py:class:`sasktran.ClimatologyUserDefined` class. Some climatologies that we use a lot are built 
into SASKTRAN (see `Useful Climatology Shorthands`_).

.. autosummary::

    sasktran.ClimatologyUserDefined
    sasktran.ClimatologyUserDefined2D
    sasktran.ClimatologyUserDefined3D
    sasktran.Labow
    sasktran.Pratmo
    sasktran.MSIS90
    sasktran.ECMWF
    sasktran.GloSSAC


.. autoclass:: sasktran.Climatology
    :exclude-members: skif_object


User Defined Climatologies
--------------------------

.. autoclass:: sasktran.ClimatologyUserDefined

.. autoclass:: sasktran.ClimatologyUserDefined2D

.. autoclass:: sasktran.ClimatologyUserDefined3D

.. note::

    For all climatologies the altitude grid must be given in meters.

.. note::

    All species climatologies must be in units of :math:`\mathrm{\frac{molecules}{cm^3}}`.


Useful Climatology Shorthands
-----------------------------

.. autoclass:: sasktran.Labow

.. autoclass:: sasktran.Pratmo

.. autoclass:: sasktran.MSIS90

.. autoclass:: sasktran.GloSSAC

See Also
--------
.. include:: descriptions/atmosphere.desc
.. include:: descriptions/species.desc
.. include:: descriptions/optical_property.desc
