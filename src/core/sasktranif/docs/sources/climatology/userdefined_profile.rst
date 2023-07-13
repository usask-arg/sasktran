.. _clim_userdefined_profile:

USERDEFINED_PROFILE
===================
A climatology class that allows users to define their own height profiles of scalar
parameters.  The climatology class can hold an unlimited number of different profiles,
each indexed by a unique CLIMATOLOGY_HANDLE.  Different scalar parameters can be on
different grids within the same class instance.

Each profile consists of a height-grid definition and the values on that grid.
The height grid can either be re-used or changed for subsequent families of values.
The algorithm uses spline or piecewise linear interpolation of the values on the height grid.

The height-grid is specified as an array of heights above sea-level in ascending order
expressed in meters.  Note that radiative transfer models will
normally require the grid to cover the entire range of the atmosphere from 0.0 meters
to the top of the atmosphere, typically 100,000 meters. The height grid is set by
calling SetProperty( “Heights”, h ).

The user can request that logarithmic or linear interpolation be used by calling
SetProperty( “DoLogInterpolation”, xx ). This value is in effect
for all the grids stored in the object. The default value is linear.

The user can specify the value returned for heights above or below the heights grid
by calling SetProperty( “BadValue”, xx ). This value is in effect for all
subsequent grids until it is changed. The default value is NaN.

Values are set on the current grid by calling :meth:`~ISKClimatology.SetPropertyUserDefined`.  The current grid settings are copied and
combined with the value array and a new entry made for this grid array. If a grid
for this species already exists it is replaced.

A Matlab example is shown below. Note that the class ignores latitude, longitude
and mjd values and only uses the height value passed in by the user
to function GetParameter::

    climate = ISKClimatology('USERDEFINED_PROFILE')
    h = (0:100)*1000.0;
    v = h.*h + 1.0;
    climate.SetProperty(“DoLogInterpolation”, 1);
    climate.SetProperty( “Heights”, h );
    climate.SetPropertyUserDefined('SKCLIMATOLOGY_O3_CM3', v);
    [ok, value] = climate.GetParameter( 'SKCLIMATOLOGY_O3_CM3', [50,102, 35000.0, 53000]);

Cache Snapshot
--------------
The USERDEFINED_PROFILE caches the entire profile

Python extension
----------------
The USERDEFINED_PROFILE climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The USERDEFINED_PROFILE climatology needs no external preparation.

Properties
----------

.. py:module:: USERDEFINED_PROFILE

DoLogInterpolation
^^^^^^^^^^^^^^^^^^
.. py:function:: DoLogInterpolation (int n) [Default: 0]

    Instructs the object to perform logarithmic or linear interpolation of the scalar profile on the height-grid.
    Logarithmic interpolation is implemented as linear interpolation of the logarithm of the scalar profile followed
    by an exponentiation. This value will only be applied to profiles created by subsequent calls to :meth:`SetPropertyUserDefined`. It
    is not applied to profiles that have already been created.

    ================ ===========
      n              Setting
    ================ ===========
      0              Scalar profiles are interpolated linearly in altitude.
      1              Scalar profiles are interpolated using logarithmic interpolation.
    ================ ===========

DoPiecewiseLinear
^^^^^^^^^^^^^^^^^

.. py:function:: DoPiecewiseLinear

    Instructs the object how to interpolate between values on the height grid. This value will only be applied to 
    profiles created by subsequent calls to :meth:`SetPropertyUserDefined`. It
    is not applied to profiles that have already been created.

    ================ ===========
      n              Setting
    ================ ===========
      0              Fit a bezier spline to the height-grid points
      1              Perform piecewise linear interpolation between points.
    ================ ===========

Heights
^^^^^^^
.. py:function:: Heights

    Sets the height grid that will be used in subsequent calls to :meth:`~ISKClimatology.SetPropertyUserDefined`.
    The heights must be in ascending order and specify the height of the grid point above sea level in meters.
    The default is an empty array.
    
BadValue
^^^^^^^^
.. py:function:: BadValue

    Sets the value to be returned for heights outside the range of the height grid. The default is NaN.



 
