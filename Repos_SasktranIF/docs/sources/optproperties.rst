.. _optproperties:

******************
Optical Properties
******************
Optical property objects provide the absorption, extinction and scattering of
individual molecules/particles in the atmosphere. They are created using the
:class:`ISKOpticalProperty` class and are used by the radiative transfer
engines to calculate scattering, extinction and absorption in the atmosphere. Many optical property
classes are solely focused on the properties of a single molecule while a few classes support
several species (eg Hitran). The majority of optical properties depend upon atmospheric state,  typically
local pressure and/or temperature, and this information is provided to the optical property object as a
climatology and location (lat, lon, height, utc).

==============================================      ======================  =====================================================
Optical Property                                    Extension               Description
==============================================      ======================  =====================================================
:ref:`optical_basspaurlinear`                       :ref:`sasktran_core`    O3, Bass-Paur Linear
:ref:`optical_basspaurquadratic`                    :ref:`sasktran_core`    O3, Bass-Paur Quadratic
:ref:`optical_o3dbm`                                :ref:`sasktran_core`    O3, Daumont, Brion and Malicet 1992.
:ref:`optical_o3gomeburrows`                        :ref:`sasktran_core`    O3, Burrows with GOME instrumenmt
:ref:`optical_o3osirisres`                          :ref:`sasktran_core`    O3, Convolved cross-section used for early OSIRIS analysis
:ref:`optical_o3sciabogumilv3`                      :ref:`sasktran_core`    O3, Bugumil for SCIA. Version 3
:ref:`optical_o3sciabogumilv4`                      :ref:`sasktran_core`    O3, Bogumil for SCIA. Version 4.
:ref:`optical_o3serdyuchenkov1`                     :ref:`sasktran_core`    O3, Serdyuchenko.
:ref:`optical_o3voigt`                              :ref:`sasktran_core`    O3, Voigt.
:ref:`optical_no2burrows`                           :ref:`sasktran_core`    NO2, Burrows with GOME instrument
:ref:`optical_no2osirisres`                         :ref:`sasktran_core`    NO2, Convolved cross-section used for OSIRIS analysis.
:ref:`optical_no2vandaele`                          :ref:`sasktran_core`    NO2, Vandaele 1998.
:ref:`so2_vandaele_2009`                            :ref:`sasktran_core`    SO2, Vandaele 2009.
:ref:`so2_bogumil_2003`                             :ref:`sasktran_core`    SO2, Bogumil 2003.
:ref:`so2_rufus_2003`                               :ref:`sasktran_core`    SO2, Rufus, 2003.
:ref:`so2_freeman1984`                              :ref:`sasktran_core`    SO2, Freeman 1984.
:ref:`o2_o2_hitran2016`                             :ref:`sasktran_core`    O2-O2 CIA. Hitran database
:ref:`o2_o2_thalman2013`                            :ref:`sasktran_core`    O2-O2 CIA. Thalman 2013.
:ref:`o2_o2_fally2000`                              :ref:`sasktran_core`    O2-O2 CIA. Fally 2000.
:ref:`optical_rayleigh`                             :ref:`sasktran_core`    Rayleigh scattering, Bates.
:ref:`optical_mieaerosol`                           :ref:`sasktran_core`    Mie aerosol
:ref:`optical_hitranchemical`                       :ref:`sasktran_core`    Hitran molecules
:ref:`optical_userdefined`                          :ref:`sasktran_core`    User-defined cross-sections
:ref:`optical_convolvedxsec`                        :ref:`sasktran_core`    Convolved cross-sections
==============================================      ======================  =====================================================



..  toctree::
    :maxdepth: 5

    optprop/o3_basspaurlinear
    optprop/o3_basspaurquadratic
    optprop/o3_dbm
    optprop/o3_gomeburrows
    optprop/o3_osirisres
    optprop/o3_sciabogumilv3
    optprop/o3_sciabogumilv4
    optprop/o3_serdyuchenkov1
    optprop/o3_voigt
    optprop/no2_burrows
    optprop/no2_osirisres
    optprop/no2_vandaele1998
    optprop/so2
    optprop/rayleigh
    optprop/mieaerosol_xxx
    optprop/hitranchemical_xxx
    optprop/o4
    optprop/userdefined_tables
    optprop/convolvedxsec
