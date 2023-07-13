.. _clim_userdefined_profile_plane:

USERDEFINED_PROFILE_PLANE
=========================
A two dimensional user defined climatology which varies both in altitude
and angle (e.g. latitude)  in a specified plane.  Interpolation is done bilinearly.  For points
not in the plane, the angle is determined by projecting the point onto the plane.
The ``values`` passed into :meth:`~ISKClimatology::SetPropertyUserDefined` must be
a 2D array stored in the order, ``values( heights, lats)``.
i.e. *heights* varies most rapidly and *lats* varies the slowest::
   
   mjd = 52393.3792987115;
   climate = ISKClimatology(‘USERDEFINED_PROFILE_PLANE’)
   climate.UpdateCache( [0,0,0,mjd])


Cache Snapshot
--------------
The USERDEFINED_PROFILE_PLANE caches values across the plane

Python extension
----------------
The USERDEFINED_PROFILE_PLANE climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The USERDEFINED_PROFILE_PLANE climatology needs no external preparation.

Properties
----------

.. py:module:: USERDEFINED_PROFILE_PLANE

Heights
^^^^^^^
.. py:function:: Heights

    Sets the height grid that will be used in subsequent calls to :meth:`~ISKClimatology.SetPropertyUserDefined`.
    The heights must be in ascending order and specify the height of the grid point above sea
    level in meters.The default is an empty array                           |

Angles
^^^^^^
.. py:function:: Angles

    Sets the angular grid used for the climatology. This must be set before adding any species to the
    climatology. The angles are specified in degrees. An angle of 0 corresponds to the reference
    vector of the plane and angles increase towards the direction of normal cross reference.

NormalAndReference
^^^^^^^^^^^^^^^^^^
.. py:function:: NormalAndReference

    A 6 element array. The first 3 elements define the unit vector normal to the plane (X,Y,Z).
    The last 3 elements define a reference vector in the plane that defines the direction of
    angle zero. i.e defines the x-axis of the plane.

DoLogInterpolation
^^^^^^^^^^^^^^^^^^
.. py:function:: DoLogInterpolation

    Instructs the object to perform interpolation of the logs of the data on the height-grid.

    ================ ===========
      n              Setting
    ================ ===========
      0              Scalar profiles are interpolated linearly in altitude.
      1              Scalar profiles are interpolated using logarithmic interpolation.
    ================ ===========
