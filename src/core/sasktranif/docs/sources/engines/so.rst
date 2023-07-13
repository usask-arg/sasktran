.. _engineso:

SO
===

The ``SO`` engine solves the scalar radiative transfer equation using the successive orders technique. It was the
first engine developed for the sasktran framework. The **SO** engine is deprecated and is no longer under development.
It has been superceded by the :ref:`enginehr` engine which also uses a successive orders algorithm but has significant
accuracy improvements as well as additional important features such as weighting function generation and adaptive ray
splitting.

The theoretical basis for the **SO** technique can be found in [soBourassa2008]_. The ``SO`` engine, which was developed for
ozone retrievals from odin-osiris measurements, is optimized for computational speed at the expense of memory.
The model solves the radiative transfer over a geographic region surrounding a location known as the *Reference Point*;
this is typically the location where the solution is required (eg location of limb scan etc.) The model
approximates the Earth as a sphere with the radius matching the North-South curvature of the true Earth
at the *Reference Point*.

The model works in two stages. The first stage is a *geometry* stage where all of the rays used by the model are traced
through the atmosphere and split into manageable segments throughout the atmosphere. Many geometrical factors and calculations generated
during the ray-tracing are internally cached and tabulated for later use by the second stage. The second stage peforms all the
*wavelength specific* calculations.

The primary advantage of the ``SO`` engine is that the ray tracing is performed
once and is shared between subsequent wavelength calculations. This saving can be significant within atmospheric retrievals where
the optical properties of the atmosphere are modified multiple times as part of the fitting process.

On the other hand this strength is also a weakness as (i) the ray tracing is not adaptive and cannot split rays when the
optical depth is too large and (ii) the internal memory footprint can get large if many diffuse profiles are required across
the region of interest. The model works well for scenarios that avoids these two problem conditions, for example
uv-visible calculations with the Sun well above the horizon and small horizontal gradients in atmospheric properties. The
:ref:`enginehr` addresses these issues and is the recommended engine for successive orders calculations.

Properties
----------

..  module:: SO

NumOrdersOfScatter
^^^^^^^^^^^^^^^^^^
..  py:function:: NumOrdersOfScatter(int n)

    Sets the number of orders of scatter used in the successive orders approximation. This number should be set high
    enough to ensure the solution has converged. We have found that values as high as 50 are often required at UV-VIS
    wavelengths close to 350 nm to ensure convergence.

NumDiffuseProfiles
^^^^^^^^^^^^^^^^^^
..  py:function:: NumDiffuseProfiles(int n)

    Sets the number of diffuse profiles used in the calculation. The default is 1, implying there is no spatial variation
    variation in the diffuse signal across the calculation region. We recommend setting odd values for this number so the
    middle profile is centered on the reference point. Typical values are the odd numbers up to ~11 although smaller
    values (3 or 5) probably capture the variation of the diffuse signal across the scene for most situations.  This
    parameter is not intended to accurately capture the variation in the diffuse profile across the terminator

    :param int n:
        The number of diffuse profiles.

ConfigureUserDefinedShells
^^^^^^^^^^^^^^^^^^^^^^^^^^
..  py:function:: ConfigureUserDefinedShells( array heights )

    Configures the model to use the altitude shells defined by the user. The
    shell altitudes must be in ascending order and are specified in meters above sea level.
    By default, shell heights are evenly spaced at 1000 m intervals from 0 to 100,000 m.
    The diffuse points, by default, are placed in the middle of the shells and the optical
    properties are placed both in the middle and on the boundary of each shell.

    :param array[n] heights:
        The array of shell altitudes in meters above sea level.

ManualDiffuseHeights
^^^^^^^^^^^^^^^^^^^^
..  py:function:: ManualDiffuseHeights( array heights )

    Manually overides the default placement of diffuse points. This is an advanced option and most users do not need to
    call this property as the default values are usually sufficient.

    :param array[n] heights:
        The array of diffuse altitudes specified in meters above sea level. Diffuse points are placed at each altitude
        on each diffuse profile.

DiffuseIncomingResolution
^^^^^^^^^^^^^^^^^^^^^^^^^
..  py:function:: DiffuseIncomingResolution( array[4] parameters)

    Specifies four parameters that specify the distribution of azimuth and zenith angles across the incoming unit sphere of
    each diffuse point. This distribution of points is used to calculate the incoming signal which is then scattered to the
    outbound unit sphere of each diffuse point.The algorithm breaks the zenith angles of the incoming unit sphere of each
    diffuse point into 3 regions:

        * Upward incoming signal from the lower (ground) regions.
        * Horizontal incoming signalfrom the angles close to the limb
        * Downward signal from altitudes above the diffuse point.

    This is an alternative method to using :py:func:`~SO.ConfigureIncomingZenithAngles` and :py:func:`~SO.ConfigureIncomingAzimuthAngles`
    to set the incoming zenith and azimuth angles. Note that the default values are set in these other functions.

    :param  array[4] parameters:
        An array of 4 integers that configure the distribution of incoming points on the incoming unit sphere of each diffuse
        point.

            * [0] Ground resolution: number of zenith angles between 180 and 100 degrees on the incoming diffuse unit sphere.
            * [1] Horizon resolution: number of zenith angles between 100 and 80 degrees on the incoming diffuse unit sphere.
            * [2] Atmosphere resolution: number of zenith angles between 0 and 80 degrees on the incoming diffuse unit sphere.
            * [3] Number of azimuths: the number of horizontal azimuths around the incoming diffuse unit sphere.


ConfigureIncomingZenithAngles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  py:function:: ConfigureIncomingZenithAngles( array(n) zenithangles)

    Configures the incoming zenith angles for the incoming unit sphere of each diffuse point. This is an array of ascending
    zenith angles in degrees. Each element of the array specifies the zenith angle at which the solid angle of the
    corresponding segment finishes. I.E. The first value should not be 0.0 as that implies the solid angle segment starts
    and ends at 0.0. Likewise the last value should be 180.0 as the last segment will finish at 180.0 zenith angle.

    :param array[n] zenithangles:
        The array of zenith angles in degrees in ascending order. The default value is, [15, 30, 40, 50, 60, 70, 75, 80,
        85, 87, 89, 90, 91, 93, 95, 100, 105, 110, 120, 130, 140, 150, 165, 180].

ConfigureIncomingAzimuthAngles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  py:function:: ConfigureIncomingAzimuthAngles( array(n) azimuthangles)

    Configures the incoming azimuth angles for the incoming unit sphere of each diffuse point.

    :param array[n] azimuthangles:
        This is an array of ascending azimuth angles in degrees. The default value is [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]

SetReferencePoint_TargetAltitude
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  py:function:: SetReferencePoint_TargetAltitude( double value)

    Sets the target altitude to be used when determining the reference point. The target altitude is used with the
    target range variable in limb viewing geometries to apply an increased weighting to
    lines of sight tangential in the vicinity of the target altitude. This encourages
    the reference point to be closer to the lines of sight which are tangential in the
    region of the reference point.

SetReferencePoint_TargetRange
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  py:function:: SetReferencePoint_TargetRange(double value)

    Sets the reference point target range parameter. Specifies the altitude range above and below
    the target altitude for enhanced weighting of limb viewing lines of sight, default 15000 meters.

SetGroundAltitude
^^^^^^^^^^^^^^^^^
..  py:function:: SetGroundAltitude(double value)

    Sets the altitude in meters of the ground shell above the oblate spheroid. Default is 0.0

SetUpperBoundAltitude
^^^^^^^^^^^^^^^^^^^^^
..  py:function:: SetUpperBoundAltitude(double value)

    the maximum altitude of the atmosphere in meters used when considering the reference point

SetLowerBoundAltitude
^^^^^^^^^^^^^^^^^^^^^
..  py:function:: SetLowerBoundAltitude(double value)

    the minimum altitude of the atmosphere in meters used when considering reference point

SetReferencePoint
^^^^^^^^^^^^^^^^^
..  py:function:: SetReferencePoint ( array position )

    Manually set the reference point to this location.

    :param array[4] position:
        A three element array specifing [latitude, longitude, height, mjd]

SetSun
^^^^^^
..  py:function:: SetSun( array position )

     Manually set the unit vector from the Earth to the Sun. A three element array specifing the unit vector from the
     Earth to the Sun [x,y,z] in the global geographic coordinate system.

References
----------
.. [soBourassa2008] | A.E. Bourassa, D.A. Degenstein, E.J. Llewellyn, SASKTRAN: A spherical geometry radiative transfer code for efficient estimation of limb scattered sunlight, Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 109, Issue 1, January 2008, Pages 52-73, ISSN 0022-4073, http://dx.doi.org/10.1016/j.jqsrt.2007.07.007.

