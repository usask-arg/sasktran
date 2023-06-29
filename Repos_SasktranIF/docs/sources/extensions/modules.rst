.. _modules:

*******************
Sasktran Extensions
*******************

The Sasktranif package has been built so it can be extended with external modules that are compatible with the
:ref:`interfaces`. The table below is a quick guide for the extension modules that allows you to choose which python
packages you need to install on your machine.


==========================  ===========================================================================
**Module**
==========================  ===========================================================================
:ref:`sasktran_core`         Default components supplied with the standard sasktran installation
:ref:`sasktran_do`           Sasktran implementation of DISORT. A plane parallel RTM.
`sasktran_tir`_              Thermal infra-red radiative transfer engine. Suitable for Thermal IR.
`sasktran_gcm`_              Climatology interfaces to outputs from global circulation models
==========================  ===========================================================================

.. _sasktran_tir: https://arg.usask.ca/docs/sasktran_tir/index.html
.. _sasktran_gcm: https://arg.usask.ca/docs/sasktran_gcm/index.html

..  toctree::
    :maxdepth: 2

    sasktran_core
    sasktran_do



