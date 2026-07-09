.. _geometry:

Geometry
********

Module: ``sasktran.geometry``

Objects used to specify the :ref:`lines_of_sight` and other geometric specifications to your 
:ref:`engine` such as the sun's position and where you would like the engines reference point to
be. Note that at a minimum, a :py:class:`sasktran.Geometry` object must specify a list of 
:py:class:`sasktran.LineOfSight` objects (so that the engine knows what look vectors to 
calculate radiances for).

.. autoclass:: sasktran.Geometry

.. autoclass:: sasktran.VerticalImage

.. autoclass:: sasktran.NadirGeometry

See Also
--------
.. include:: descriptions/lines_of_sight.desc
.. include:: descriptions/geodetic.desc