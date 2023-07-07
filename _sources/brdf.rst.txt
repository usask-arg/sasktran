.. _brdf:

BRDF
====

Module: ``sasktran.brdf``

SASKTRAN support `BRDF surfaces
<https://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_. A surface must
be added to every :ref:`atmosphere` object. Below are the BRDF classes currently supported by
SASKTRAN.

Common BRDFs
------------

Below is the documentation for the common BRDFs that SASKTRAN supports.

.. autoclass:: sasktran.BRDF

.. autoclass:: sasktran.Lambertian

.. autoclass:: sasktran.Kokhanovsky

.. autoclass:: sasktran.Roujean

.. autoclass:: sasktran.CoxMunk

.. autoclass:: sasktran.Rahman

.. autoclass:: sasktran.Hapke

.. autoclass:: sasktran.MODIS

Kernel-Based BRDFs
------------------

Below is the documentation for the kernel-based BRDFs that SASKTRAN supports.

.. autoclass:: sasktran.LinearCombination

.. autoclass:: sasktran.RoujeanKernel

.. autoclass:: sasktran.LiKernel

.. autoclass:: sasktran.LiSparseKernel

.. autoclass:: sasktran.LiSparseReciprocalKernel

.. autoclass:: sasktran.LiDenseKernel

.. autoclass:: sasktran.RossThinKernel

.. autoclass:: sasktran.RossThickKernel

Two-Dimensional BRDFs
---------------------
SASKTRAN is capable of modelling surfaces where the BRDF varies across the surface of the Earth.
These are BRDFs that have two-dimensional variation across the surface

.. autoclass:: sasktran.LatLonBRDF



See Also
--------
.. include:: descriptions/atmosphere.desc