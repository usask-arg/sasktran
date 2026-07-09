.. _layering:

Layering and Vertical Grids in SASKTRAN-DO
******************************************
The discrete ordinates method (and by extension, SASKTRAN-DO) represents the atmosphere by a set of stacked,
homogenous layers.  Each layer is characterized by the vertical optical depth, single scattering albedo, and
Legendre coefficients of the phase function expansion.  If you are familiar with "pure scattering" radiative transfer
models such as DISORT or LIDORT these are the input quantities you may be used to, however SASKTRAN-DO operates
with a higher level interface.  This page describes internally how SASKTRAN-DO
interprets the input it is provided.

The Altitude Grid
-----------------
We can think of the primary input to SASKTRAN-DO being number density of a variety of species on a fixed altitude grid.
The altitude grid object is set through the :meth:`sasktran.EngineDO.alt_grid` method of the engine,
with the default being a linearly spaced grid from 0 km to 100 km with a spacing of 0.5 km.  Any climatology
you provide to the atmospheric state will be sampled at this altitude grid.  If you are using user defined climatology
objects we recommend matching the altitude grid in those objects to the altitude grid in the engine, although
this is not strictly necessary. When constructing homogeneous layers the model will internally integrate the atmospheric
state over this altitude grid, assuming linear interpolation between grid points.  Any weighting function with
respect to an atmospheric species will be returned with respect to this altitude grid.  The altitude grid need not be
uniformly spaced, and the bottom and top altitudes can be set to any value.

Layer Construction
------------------
While inputs are specified on an altitude grid, internally the discrete ordinate technique requires homogenous layers.
As previously mentioned, the layer properties are calculated by integrating the atmospheric state over the layer boundaries,
but there are multiple ways to determine how the layer boundaries themselves are determined.  SASKTRAN-DO
provides multiple options to do this, which method is preferred can be application dependent.  The two properties that
primarily influence the layer construction are :meth:`sasktran.EngineDO.num_layers`, which
controls the total number of homogeneous layers, and :meth:`sasktran.EngineDO.layer_construction` which
controls the technique to determine the layer boundaries.

Uniform Pressure Layers
^^^^^^^^^^^^^^^^^^^^^^^
The default method to construct layers is to place the boundaries uniform in pressure.  This is a good technique
for a variety of applications as it results in layers that are more densely spaced close to the surface where scattering is
the greatest.  It also results in layer boundaries that are independent of the input concentrations.  This is the default
setting but can be manually specified with

.. code-block:: python

    engine.layer_construction = 'uniform_pressure'

Uniform Height Layers
^^^^^^^^^^^^^^^^^^^^^
Another option is to place the layers uniformly in altitude.  This can be beneficial when operating in limb mode
where altitude variation is important or in areas of large absorption.  This can be set through

.. code-block:: python

    engine.layer_construction = 'uniform_height'

A secondary option is to match the layers directly with the internal altitude grid through

.. code-block:: python

    engine.layer_construction = 'match_altitude_grid'

note that setting this will set both the number of layers and layer construction method options.

Manually Specifying Layer Boundaries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Lastly, layers can be manual specified through

.. code-block:: python

    engine.layer_construction = array_containing_boundaries_in_m

this option will simultaneously set the layer boundaries and the number of layers.  The number of layers will be
equal to the number of layer boundaries minus one.

Accuracy Involving the Number of Layers
---------------------------------------
Including more layers both improves the models ability to represent vertical inhomogeneity, and improves the calculation
of the pseudo-spherical beam transmittance.  The number of layers required to achieve a certain level of accuracy will be
highly application dependent.  Generally in areas of the spectrum with little absorption less layers are required, while large
amounts of absorption benefits from more layers.  It is recommended to play around with the number of layers until a specific
level of accuracy has been obtained.
