.. _weighting_functions:

Weighting Functions
===================

SASKTRAN-TIR supports the analytic computation of weighting functions for concentrations of gas species and temperature.
It is possible to compute weighting functions for a single parameter or for multiple parameters simultaneously.
Weighting functions computed for a single parameter have the form:

.. math::

    \mathrm{wf[i, j, k]} = \frac{\partial I(\lambda_\mathrm{i}, \mathrm{L_j})}{\partial \mathrm{S_k}}

while those computed for multiple species are given by:

.. math::

    \mathrm{wf[i, j, n, k]} = \frac{\partial I(\lambda_\mathrm{i}, \mathrm{L_j})}{\partial \mathrm{S_{n, k}}}

where

    ================================    =========================================================
    Symbol                              Meaning
    ================================    =========================================================
    :math:`\lambda_\mathrm{i}`          Wavelength :math:`\mathrm{i}`
    :math:`\mathrm{L_\mathrm{j}}`       Line of sight :math:`\mathrm{j}`
    :math:`\mathrm{S_\mathrm{n, k}}`    Species :math:`\mathrm{n}` at altitude :math:`\mathrm{k}`
    ================================    =========================================================

To calculate weighting functions for a gas species the :py:attr:`sasktran.Atmosphere.wf_species` must be set.
To compute temperature weighting functions set the :py:attr:`sasktran_tir_engine.EngineTIR.do_temperature_wf` option to True.

Weighting Function Units
------------------------

By definition, SASKTRAN-TIR weighting functions have units of :math:`\mathrm{radiance / species}`.
The radiance portion has units of :math:`\mathrm{photon}/(\mathrm{s \: cm^2 \: nm \: sr})`.
The units of the weighting function species are configurable.
By default, weighting functions are given for gas number density, in :math:`\mathrm{cm}^{-3}`.
To calculate gas species weighting functions with respect to Volume Mixing Ratio (VMR) set the
:py:attr:`sasktran_tir.engine.EngineTIR.do_vmr_wf` option to True.

Multiple Species Calculations
-----------------------------

When weighting functions are computed for more than one species, an extra dimension is added to the weighting function
array. The index :math:`n` of each species are assigned in the same order as the species were listed when setting
:py:attr:`sasktran.Atmosphere.wf_species`. If temperature weighting functions are also being computed, they are
assigned the highest species index.

For example, consider a SASKTRAN-TIR engine object with the
weighting function species configured as follows:

.. code-block:: python

    engine.atmosphere.wf_species = ['CO2', 'O3']
    engine.do_temperature_wf = True
    rad, wf = engine.calculate_radiance()

The CO2 weighting functions will be given by ``wf[:, :, 0, :]``, the O3 weighting functions by ``wf[:, :, 1, :]``, and the temperature weighting functions
by ``wf[:, :, 2, :]``.
