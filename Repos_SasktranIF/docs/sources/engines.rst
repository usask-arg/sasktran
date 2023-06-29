.. _engines:

*******
Engines
*******

The Sasktranif engines are implementations of radiative tranbsfer models that calculate radiance in the atmosphere.
The following engines have been developed and released

========================  ==================== ==========================================================
Engine                    Extension            Description
========================  ==================== ==========================================================
:ref:`HR <enginehr>`      :ref:`sasktran_core` Successive orders spherical RTM
:ref:`OCC <engineocc>`    :ref:`sasktran_core` Solar occultation RTM
DO                        :ref:`sasktran_do`   Sasktran implementation of DISORT, plane parallel RTM.
`TIR`_                    `sasktran_tir`_      Infra-red thermal emissions, spherical RTM
:ref:`SO <engineso>`      :ref:`sasktran_core` Successive orders spherical RTM (deprecated).
========================  ==================== ==========================================================

..  toctree::
    :maxdepth: 2

    engines/hr
    engines/occ
    engines/so


.. _sasktran_tir: https://arg.usask.ca/docs/sasktran_tir/index.html
.. _TIR: https://arg.usask.ca/docs/sasktran_tir/sasktran_tir.html
