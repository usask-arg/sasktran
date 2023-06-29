.. _clim_constantvalue:

CONSTANTVALUE
==============
A climatology that supports a single value at all locations and times.
This is useful to set climatologies to 1 or 0 everywhere. The class supports all species::

   import sasktranif as skif

   mjd = 52393.3792987115;
   climate = skif.ISKClimatology(‘CONSTANTVALUE’)
   climate.SetPropertyScalar(‘SetConstantValue’, 0.0);
   climate.UpdateCache( [0,0,0,mjd])
   
Supported Species
-----------------
All species are implicitly supported


Cache Snapshot
--------------
The cahed value is the constant.

Python extension
----------------
The CONSTANTVALUE climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Properties
----------

..  module:: constantvalue

SetConstantValue
^^^^^^^^^^^^^^^^
.. py:function:: SetConstantValue( double value )

    Sets the climatology to the value passed in. This value is returned in response to all calls to
    :meth:`ISKClimatology.GetParameter`


