.. _clim_onepressuretemp:

ONE_PRESSURE_TEMP
=================
A climatology that supports a single pressure and temperature. This is useful for cross-section calculations
that need an atmospheric state and the user wants to avoid the overhead and complication of realistic 
climatologies. In this climatology the user explicitly sets the temperature and pressure. The class supports 
pressure (pascals), temperature (kelvin) and number density per cm3::

   mjd = 52393.3792987115;
   climate = ISKClimatology(‘ONE_PRESSURE_TEMP’)
   climate.UpdateCache( [0,0,0,mjd])

Supported Species
-----------------
* SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3
* SKCLIMATOLOGY_PRESSURE_PA
* SKCLIMATOLOGY_TEMPERATURE_K


Cache Snapshot
--------------
The cache is the values provided.

Python extension
----------------
The ONE_PRESSURE_TEMP climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The ONE_PRESSURE_TEMP requires no external configuration.

Properties
----------
..  module:: ONE_PRESSURE_TEMP

.. py:function:: SetTemperature( double T)

    Sets the temperature of the climatology to the value passed in. This value is in Kelvins.

.. py:function:: SetPressure(double P)

    Sets the pressure of the climatology to the value passed in. The value is in Pascals.



