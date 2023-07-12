.. _clim_osirisl2_aerosolmoderadius_v600:

OSIRISL2_AEROSOLMODERADIUS_V600
===============================
A climatology that provides the Odin-OSIRIS V.00 aerosol mode radius. The climatology works for
all latitudes and longitudes although the interpolation is assumed to be coincident with regions
actively probed by OSIRIS measurements, i.e. the sunlit hemispheres.

The climatology actually consists of two climatologies

1. one from the ascending node measurements 
2. one from the descending mode measurements.

The difference between the two climatologies results from 
systematic errors in the aerosol particle size properties (and other parameters) used in 
the retrieval.

Users must actively select which database they wish to use before calling 
:meth:`~ISKClimatology::UpdateCache` or :meth:`~ISKClimatology::GetParameter` by setting 
property SetIsAscendingNode by calling :meth:`~ISKClimatology::SetProperty`::

   mjd = 52393.3792987115;
   climate = ISKClimatology('OSIRISL2_AEROSOLMODERADIUS_V600')
   climate.SetPropertyScalar('SetIsAscendingNode', 1)
   climate.UpdateCache( [0,0,0,mjd])
   [ok, value] = climate.GetParameter( 'SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS', [0.0, 0.0, 25000, mjd]);
   [ok, value] = climate.GetParameter( 'SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH', [0.0, 0.0, 25000, mjd]);
   
Supported Species
------------------

* SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS
* SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH

This climatology will only work if your local system has access to the Odin-OSIRIS version 6.00 aerosol retrieval products.

Cache Snapshot
--------------
The cache is the values provided.

Python extension
----------------
The OSIRISL2_AEROSOLMODERADIUS_V600 climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The OSIRISL2_AEROSOLMODERADIUS_V600 requires no external configuration.

Properties
----------

..  module:: OSIRISL2_AEROSOLMODERADIUS_V600

.. py:function:: SetIsAscendingMode

    Instructs the object to perform interpolation of the logs of the data on the height-grid.
    The default value is NaN and must be set for proper use.

    ================ ===========
      n              Setting
    ================ ===========
      0              Use the descending node database
      1              Use the ascending node database .
    ================ ===========
