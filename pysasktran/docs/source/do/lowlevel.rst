.. _lowlevel:

Overview
********
In addition to being part of the standard SASKTRAN framework, SASKTRAN-DO exposes a secondary interface
that we colloquially refer to as the low-level interface.  The low-level interface has the potential to be
faster and allow for more generic calculations, as the cost of increased complexity in setting it up.  Instead
of supplying high level constructs to the model like climatologies (number densities) and optical properties (cross sections
and phase functions), the user directly supplies low-level constructs like layer optical depth, layer single scatter albedo, and layer
Legendre coefficients.  Constructing these quantities requires significantly more developer effort and is only recommended
for advanced users who are experienced in radiative transfer.

The input to the low-level interface is four structures:

- An :class:`sktran_disco.lowlevel.Atmosphere` object which defines the layer optical depths, single scatter albedos, legendre coefficients, and other aspects to do with the layering.
- An :class:`sktran_disco.lowlevel.Config` object which sets some configuration settings
- An :class:`sktran_disco.lowlevel.ViewingGeometry` object which specifies the viewing geometry
- Optionally a :class:`sktran_disco.lowlevel.WeightingFunctions` object to specify which quantities to calculate derivatives with respect to.

Radiative transfer is performed through a single function call,
:meth:`sktran_disco.lowlevel.calculate`. We have found that the notation and definitions of various
radiative transfer quantities are rarely consistent, and as such the following sections
briefly outline
the theory of the discrete ordinates technique with a focus on how to define the necessary inputs
for the model. For details on how to construct these inputs in code, we recommend reviewing the API reference for
the input and output structures, and looking at some of the examples.

Optical Layer Quantities
------------------------
DO performs radiative transfer through a medium consisting of :math:`N_{lyr}` optically homogeneous
layers.  We use the index :math:`i` to index over each layer, :math:`i=0` always corresponds to the
layer at the Top of the Atmosphere (TOA), while :math:`i=N_{lyr}-1` is the layer touching the surface
of the Earth.  Optically each layer is defined by its vertical optical depth, :math:`\Delta_i`, single scatter albedo,
:math:`\omega_i`, and scattering Legendre coefficients, :math:`\vec{\beta}_{i,l}`.  The index
:math:`l` runs from :math:`0 \ldots N_{str}-1`, resulting in a second user defined configuration option,
:math:`N_{str}` that influences the accuracy of the calculation.

The Legendre coefficient vector is defined in terms of the phase matrix.  Suppressing
the layer index for now,
in generalized form the phase
matrix takes the shape,

.. math::

    P(\Theta) = \begin{pmatrix}
                    P_{11}(\Theta) & P_{12}(\Theta) & 0 & 0 \\
                    P_{12}(\Theta) & P_{22}(\Theta) & 0 & 0 \\
                    0 & 0 & P_{33}(\Theta) & P_{34}(\Theta) \\
                    0 & 0 & -P_{34}(\Theta) & P_{44}(\Theta)
                \end{pmatrix}

The phase matrix elements are then expanded in terms of Wigner functions (generalized spherical
functions related to Legendre polynomials/associated Legendre functions)

.. math::

    P_{11}(\Theta) = \sum_l a_{1, l} \, d^l_{00}(\Theta)

    P_{22}(\Theta) + P_{33}(\Theta) = \sum_l (a_{2, l} + a_{3, l}) \, d^l_{2,2}(\Theta)

    P_{22}(\Theta) - P_{33}(\Theta) = \sum_l (a_{2, l} - a_{3, l}) \, d^l_{2,-2}(\Theta)

    P_{44}(\Theta) = \sum_l a_{4, l} \, d^l_{00}(\Theta)

    P_{12}(\Theta) = \sum_l b_{1, l} \, d^l_{02}(\Theta)

    P_{34}(\Theta) = -\sum_l b_{2, l} \, d^l_{02}(\Theta)

The set of 6 coefficients, :math:`\vec{\beta}_l = (a_{1, l}, a_{2,l}, a_{3,l}, a_{4,l}, b_{1,l}, b_{2,l})` then
are the necessary inputs to the model. Note that in the scalar approximation (:math:`N_{stokes}=1`), we
assume every coefficient except for :math:`a_{1,l}` is zero, and thus only the :math:`a_{1,l}` coefficients
are necessary for input. The `pysasktran` package includes functionality to calculate Wigner functions
and also calculate the necessary coefficients if you only have the phase function.

Delta Scaling
^^^^^^^^^^^^^
DO includes the option to include a delta-scaling factor for the phase function.


Single Scatter Correction
^^^^^^^^^^^^^^^^^^^^^^^^^
DO implements the TMS correction which replaces the single scatter solution with an exact single scatter calculation.
In this case the first row of the phase matrix :math:`\vec{P}` must be specified for each viewing geometry.  Note that
currently this option only works in scalar mode.


Geometry Layer Quantities
-------------------------
To further define the problem we need to specify a few quantities that are purely geometry dependent.
These are:

 - The Earth Radius
 - The altitude (relative to the Earth radius) of each boundary separating the layer

These quantities are used for spherical approximations within the model, and to determine the location
of the observer inside the atmosphere. If all spherical approximations are turned
off and the observer is outside the atmosphere they are not used.

Viewing Geometry
----------------
The viewing geometry for DO is specified in terms of viewing angles and solar angles.

 - The cosine of viewing zenith angle: defined such that `cos_vza=0` corresponds to true Nadir viewing
 - The cosine of solar zenith angle: defined such that `cos_sza=0` corresponds to the sun directly overhead
 - Relative solar azimuth angle: defined as the relative azimuth angle in radians between the observer and solar direction such that `saa=0` is the backwards scattering direction.
 - Observer altitude: the altitude of the observer relative to the Earth radius, can be set to -1 to indicate that the observer is outside the atmosphere.

Currently the only way to specify the viewing geometry for DO is as a list of `(cos_vza, saa, altitude)` triplets, with a single
`cos_sza`.

Weighting Function Specification
--------------------------------
A key feature of DO is that it is analytically linearized with respect to the optical input quantities, this
includes:

- Optical depth, :math:`\Delta`
- Single scatter albedo, :math:`\omega`
- Legendre phase moments, :math:`\vec{\beta}_l`
- Delta-truncation factor, :math:`f`
- Single scattering phase function, :math:`\vec{P}`
- Surface albedo, :math:`a`

Suppose we have a quantity :math:`x` that affects some or all of these optical quantities in a single layer :math:`i`, we
can write the derivative of radiance with respect to this quantity as,

.. math::

    \frac{\partial I}{\partial x} = \frac{\partial \Delta_i}{\partial x} \frac{\partial I}{\partial \Delta_i} + \frac{\partial \omega_i}{\partial x} \frac{\partial I}{\partial \omega_i} + \sum_l \frac{\partial \vec{\beta}_l}{\partial x} \frac{\partial I}{\partial \vec{\beta}_l} + \frac{\partial \vec{P}_i}{\partial x} \frac{\partial I}{\partial \vec{P}_i} + \frac{\partial a}{\partial x} \frac{\partial I}{\partial a},

The quantities on the right hand side that involve the derivative of radiance with respect to a fundamental input quantity
(e.g. :math:`\partial I / \partial \Delta_i`) are internally calculated by DO.  The other quantities,
the derivatives of the fundamental input quantities by :math:`x` (e.g. :math:`\partial \Delta_i / \partial x`) are
inputs to the model.  These input quantities must be computed by the user and are usually a straight-forward application
of differentiation.

Quite often it is desired to calculate a derivative with respect a quantity that affects multiple layers (e.g. a total column derivative).
Currently this involves some manual work by the user.  The user must calculate how this quantity affects every layer
individually, calculate :math:`N_{lyr}` derivatives, and then sum them together.

Computational Considerations
----------------------------
For most applications DO will be called many times with the same layer configuration for multiple sets of optical properties,
for example, the same calculation performed at multiple wavelengths.  For this reason, multiple sets of optical
properties can be input to the model at the same time.  We refer to this dimension as the wavelength dimension, although
it does not necessarily have to correspond to wavelength.  In other words, rather than the layer optical depths
:math:`\Delta` being a one-dimensional array of size :math:`(N_{lyr})`, it is instead a two-dimensional array
of size :math:`(N_{lyr}, N_{wavel})`.  Every optical quantity has this additional dimension, while the geometry quantities
do not.

Specifying optical properties like this instead of the user looping over
the wavelength dimension is preferable since internally the DO engine can save computational effort by not repeating
any setup that would be common between each calculation.  Furthermore it allows the DO engine to multi-thread
over the wavelength dimension if requested.

