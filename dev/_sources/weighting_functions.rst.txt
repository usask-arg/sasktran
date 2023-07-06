.. _weighting_functions:

Weighting Functions
===================

In the SASKTRAN framework, weighting functions are returned in the form:

.. math::

    \mathrm{wf[i, j, k]} = \frac{\partial I(\lambda_\mathrm{i}, \mathrm{L_j})}{\partial 
    \mathrm{S_k}}

where

    =============================   ======================================
    Symbol                          Meaning
    =============================   ======================================
    :math:`\lambda_\mathrm{i}`      Wavelength :math:`\mathrm{i}`
    :math:`\mathrm{L_\mathrm{j}}`   Line of sight :math:`\mathrm{j}`
    :math:`\mathrm{S_\mathrm{k}}`   Species at altitude :math:`\mathrm{k}`
    =============================   ======================================

Some engines might support additional forms but by default must return the form given above. To
calculate weighting functions the :py:attr:`sasktran.Atmosphere.wf_species` must be set.

Weighting Function Altitudes
----------------------------

This is fairly self explanatory; weighting functions are calculated at specific altitudes. These 
altitudes are usually set via an engine option.

Weighting Function Widths
-------------------------

In multiple :ref:`engine`'s the concept of weighting function widths is important. 

.. figure:: images/wf_equiv_fd.svg
    :figwidth: 50%
    :align: center

    Visual of what is meant by the weighting function altitude, and width for an equivalent 
    finite difference weighting function calculation.

Weighting function widths refer to the region (around a weighting function altitude) in which a
perturbation is applied. A different way of say this is that the width refers to equivalent
region a finite difference method would need to perturb the species of interest to obtain
the same result. It should be noted that the equivalent perturbation would be linearly 
interpolated from 0 at the region boundaries to the magnitude of the perturbation at the 
weighting function altitude. Below is an illustration to clarify.


.. note::

    The user should take special care when setting the weighting function widths and altitudes.
    Note that the altitudes specify the center points and the widths specify the boundaries. 
    Special care must be taken to ensure that the perturbation region does not extend above or below
    the top of the atmosphere (TOA) or surface respectively.

See Also
--------
.. include:: descriptions/engine.desc