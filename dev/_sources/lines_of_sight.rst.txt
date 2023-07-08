.. _lines_of_sight:

Lines of Sight
==============

Module: ``sasktran.lineofsight``

This is an object that describes a line of sight. A line of sight is made up of three components:

* A look vector (a unit vector in :ref:`geodetic` coordinates)
* An observer location (a position vector in :ref:`geodetic` coordinates)
* An observation time in MJD

A list of these objects are usually added to a :ref:`geometry` object to describe lines of sight 
to calculate radiances for.

.. autoclass:: sasktran.LineOfSight

See Also
--------
.. include:: descriptions/geometry.desc
.. include:: descriptions/geodetic.desc