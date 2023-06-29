.. _clim_userdefined_profile_latlonheight:

USERDEFINED_PROFILE3D_LATLONHEIGHT
==================================
A climatology that allows users to define their own profiles of scalar parameters on a 3-D grid. 
The climatology can hold an unlimited number of different profiles, each indexed by a 
unique CLIMATOLOGY_HANDLE.  Different scalar parameters can be on different grids within 
the same class instance.

Each profile consists of a grid definition, which must be set first, followed by definition 
of the values on the grid.  The grid definition remains in effect for all subsequent value 
definitions until the grid is changed.  The algorithm uses linear interpolation of the 
values on the grid.

The 3-D grid is defined by an array of heights, an array of longitudes and an array of latitudes. 
The three arrays defining the 3 axes will normally be of different sizes.

* Heights are expressed in meters above sea-level in ascending order.
  The height grid is set by calling ``ISKClimatology.SetProperty( “Heights”, h )``.
* Longitudes are expressed in degrees in ascending order in the range 0.0 to 360.0. 
  The first element of the longitude array must be 0.0.
  The last element of the longitude array must be 360.0. 
  i.e. the longitude array must cover all longitudes as we found this significantly 
  simplifies linear interpolation across the 0 degree meridian. The longitude grid is set by calling ``ISKClimatology.SetProperty(“Longitudes”, lng )``.
* Latitudes are expressed in degrees in ascending order in the range -90.0 to +90.0.  Latitudes only need to cover the region of interest 
  although users should be aware that radiative transfer models often extend the working region by substantial amounts (eg 20 degrees). 
  The latitude grid is set by calling ``ISKClimatology.SetProperty( “Latitudes”, lat)``.  

The user can specify the value returned for points outside the grid by calling ``ISKClimatology.SetProperty( “BadValue”, xx )``.
This value is in effect for all subsequent grids until it is changed. The default value is NaN.
Values are set on the current grid by calling SetPropertyUserDefined and it is only at this time that actual action is taken.  The current grid settings are copied and combined with the value array and a new entry made for this grid array. If a grid for this species already exists it is replaced.

A Matlab example is shown below. Note that the class ignores the mjd value passed in by the user to function GetParameter::


   climate = ISKClimatology('USERDEFINED_PROFILE3D_LATLONHEIGHT')
   h = (0:100)*1000.0;
   lng = 0:360;
   lats = -50:0
   v = array( 101, 361, 51);  // Do something to assign values to the grid points
   climate.SetProperty( “Heights”,    h    );
   climate.SetProperty( “Longitudes”, lng  );
   climate.SetProperty( “Latitudes”,  lats );
   climate.SetProperty( “BadValue”,   0.0);
   climate.SetPropertyUserDefined(SKCLIMATOLOGY_O3_CM3(), v);
   [ok, value] = climate.GetParameter( SKCLIMATOLOGY_O3_CM3(), [-35,102, 35000.0, 53000]);
   
SetPropertyUserDefined
^^^^^^^^^^^^^^^^^^^^^^
Calls to method :meth:`~ISKClimatology.SetPropertyUserDefined` sets the values for the grid currently in effect in this object.
The values must be a 3-D array. The array is assumed to be stored in the order, ``values( heights, lngs, lats)``  where 
*heights* varies most rapidly in memory and *lats* varies the slowest. The code will ensure that the size of array passed in by the user
is consistent with the current grid size.

Cache Snapshot
--------------
The USERDEFINED_PROFILE3D_LATLONHEIGHT caches the entire profile

Python extension
----------------
The USERDEFINED_PROFILE3D_LATLONHEIGHT climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The USERDEFINED_PROFILE3D_LATLONHEIGHT climatology needs no external preparation.

Properties
----------

.. py:module:: USERDEFINED_PROFILE3D_LATLONHEIGHT

Heights
^^^^^^^
.. py:function:: Heights

    Sets the height grid that will be used in subsequent calls to :meth:`~ISKClimatology.SetPropertyUserDefined`.
    The heights must be in ascending order and specify the height of the grid point above sea level in meters.
    The default is an empty array.

Longitudes
^^^^^^^^^^
.. py:function:: Longitudes

    Sets the longitude grid that will be used in subsequent calls to :meth:`~ISKClimatology.SetPropertyUserDefined`.
    The longitudes are specified in ascending order in degrees in the range 0 to 360. The first element
    must be 0.0 and the last element must be 360. The default is an empty array.

Latitudes
^^^^^^^^^
.. py:function:: Latitudes

    Sets the latitude grid that will be used in subsequentm calls to :meth:`~ISKClimatology.SetPropertyUserDefined`.
    The latitudes must be in ascending order in the range -90 to +90. The default is an empty array.

BadValue
^^^^^^^^
.. py:function:: BadValue

    Sets the value to be returned for locations outside the grid. The default is NaN.


