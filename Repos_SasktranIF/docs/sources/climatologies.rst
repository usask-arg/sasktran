.. _climatologies:

*************
Climatologies
*************
Climatology objects provide the value of parameters at different locations in the atmosphere. They are created using the
:class:`ISKClimatology` class and are used by radiative transfer engines to calculate
atmospheric parameters such as number density, pressure and temperature. However, climatology objects
can be built that support any atmospheric parameter, as long as it is a scalar, and can be used as stand-alone objects
in other software projects.

==============================================      ======================  =====================================================
Climatology                                         Extension               Description
==============================================      ======================  =====================================================
:ref:`clim_msis90`                                  :ref:`sasktran_core`    Fast climatological background atmosphere using MSIS and CIRA.
:ref:`clim_no2pratmo`                               :ref:`sasktran_core`    NO2 PRATMO climatological model.
:ref:`clim_o3labow`                                 :ref:`sasktran_core`    The Labow O3 ozone climatology.
:ref:`clim_userdefined_profile`                     :ref:`sasktran_core`    User defined height profiles
:ref:`clim_userdefined_profile_table`               :ref:`sasktran_core`    User defined height profile for legacy applications
:ref:`clim_userdefined_profile_plane`               :ref:`sasktran_core`    User defined height profiles that vary along a great circle /plane.
:ref:`clim_userdefined_profile_latlonheight`        :ref:`sasktran_core`    User defined height profiles distributed across the surface of the Earth
:ref:`clim_ecmwf`                                   :ref:`sasktran_core`    The ECMWF atmosphere. Needs external database files.
:ref:`clim_linearcombo`                             :ref:`sasktran_core`    A linear combination of two other climatologies
:ref:`clim_osirisl2_aerosolrtmodel_v507`            :ref:`sasktran_core`    The OSIRIS aerosol height profile generated in Version 5.07
:ref:`clim_osirisl2_aerosolmoderadius_v600`         :ref:`sasktran_core`    The OSIRIS aerosol mode radius profile generated in Version 6.00 (unofficial)
:ref:`clim_osirisl2_o3rtmodel_v507`                 :ref:`sasktran_core`    The OSIRIS O3 height profile generated in Version 5.07.
:ref:`clim_osirisl2_no2rtmodel_v507`                :ref:`sasktran_core`    The OSIRIS NO2 height profile generated in Version 5.07.
:ref:`clim_constantvalue`                           :ref:`sasktran_core`    Climatology which is a constant value everywhere.
:ref:`clim_onepressuretemp`                         :ref:`sasktran_core`    Climatology which is a fixed pressure and temperature (useful for IR cross-sections)
`GEM-MACH`_                                         `sasktran_gcm`_         GEM-MACH regional climate model
`GEOS CHEM`_                                        `sasktran_gcm`_         GEOS-CHEM global climate model
==============================================      ======================  =====================================================

.. _sasktran_gcm: https://arg.usask.ca/docs/sasktran_gcm/index.html
.. _GEM-MACH: https://arg.usask.ca/docs/sasktran_gcm/gem_mach.html
.. _GEOS CHEM: https://arg.usask.ca/docs/sasktran_gcm/geos_chem.html

Climatology parameters are identified using standard :ref:`climatology handles <climatologyhandles>`. All climatologies
implement the concept of caching where values at atmospheric locations are calculated once for the cache and then subsequent
requests to calculate atmospheric parameters are taken from the cache rather than the full actual climatological model. This
provides a significant speed advantage for the atmospheric radiative transfer models but users must be aware of the possible
subtleties.

For example, the :ref:`clim_msis90` is able to calculate atmospheric parameters at any location in the atmosphere but
its cache is written so it only extracts a single height profile. All subsequent calls to retrieve atmospheric parameters
at any location on Earth are taken from the cached vertical profile rather than executing the entire model.

..  toctree::
    :maxdepth: 2

    climatology/msis90
    climatology/no2pratmo
    climatology/o3labow
    climatology/userdefined_profile
    climatology/userdefined_plane
    climatology/userdefined_profile3d
    climatology/userdefined_profile_table
    climatology/ecmwf
    climatology/linear_combo
    climatology/osirisl2_aerosolrtmodel_v507
    climatology/osirisl2_aerosolmoderadius_v600
    climatology/osirisl2_o3rtmodel_v507
    climatology/osirisl2_no2rtmodel_v507
    climatology/constantvalue
    climatology/onepressuretemp
