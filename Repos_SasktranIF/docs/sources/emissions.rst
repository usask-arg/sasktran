.. _emissions:

*********
Emissions
*********
Emission objects provide the radiance emitted at locations in the atmosphere. They are created using the
:class:`ISKEmission` class.

Overview
--------
Emission objects within the sasktranif framework have the following peculiarities:

    * They provide the signal as a radiance per steradian rather than volume emission rate. This introduces a factor of :math:`4\pi`.
    * They provide the signal as photons per nanometer rather than photons per wavenumber. This will introduce a factor of :math:`\frac{\nu^2}{10^7}`
      to calculation primarily based upon wavenumber.

======================================  ===================== ==========================================================
Emission                                Extension             Description
======================================  ===================== ==========================================================
:ref:`emission_thermal`                 :ref:`sasktran_core`  Simple planck black body.
:ref:`emission_userdefined`             :ref:`sasktran_core`  User defined 2D array of emission, wavelength vs. height.
:ref:`emissions_hitran_photochemical`   :ref:`sasktran_core`  Photo-chemical emissions using the Hitran database
======================================  ===================== ==========================================================

..  toctree::
    :maxdepth: 2

    emissions/thermal
    emissions/hitran_photochemical

