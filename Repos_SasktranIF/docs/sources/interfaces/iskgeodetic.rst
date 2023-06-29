.. _ISKGeodetic:

***********
ISKGeodetic
***********
A helper class for SasktranIF users. The class contains an oblate spheroid model of Earth and provides
methods to convert geocentric X, Y, Z vectors to geodetic latitude, longitude and altitude coordinates.
The class also provides methods for defining lines of sight for the SasktranIF models. The coordinate
conversion methods use a 2 step process. The first step sets a location using various techniques; 
the location is internally stored in the class instance. The second step calls various get methods
to convert the internal location to the desired output coordinates.

Overview 
-------------------
This is the overview section

Geocentric Coordinates
^^^^^^^^^^^^^^^^^^^^^^
The oblate spheroid uses two primary coordinate systems. This first is the geocentric X, Y, Z system. 
This coordinate system is a 3-D Cartesian system with the origin located at the center of the Earth. 
All dimensions are specified in meters. The X axis is in the equatorial plane and points to the
Greenwich meridian at 0° longitude. The Y axis is in the equatorial plane and points to the
meridian at 90° E. The Z axis is parallel to the northward pointing spin axis of the Earth.

Geodetic Coordinates
^^^^^^^^^^^^^^^^^^^^
The second coordinate system is the geodetic latitude, longitude and altitude system. Users
should note that geodetic coordinates are not the same as geocentric latitude, longitude
and altitude. Geodetic coordinates refer to the local vertical at any position on the 
oblate spheroid and this unit vector due to the oblateness of the Earth is inclined to the 
geocentric radius.  Geodetic coordinates are the coordinate system used on most maps and 
GPS systems and are generally the coordinates normally referred to when people say 
latitude, longitude and altitude.

Geoid Model
^^^^^^^^^^^
Different models for the oblate Earth exist and the ISKGeodetic class by default uses the
IAU 1976 model for the Earth.  Other standard models are available through method 
SetProperty. Oblate spheroid models attempt to model the geopotential surface of
the Earth and are generally accurate to about 100 meters. Users requiring higher
precision models of the Earth’s geopotential must use other techniques.

Example
^^^^^^^
A simple example::

   import sasktranif.sasktranif as skif

   geoid    = skif.ISKGeodetic();
   ok       = geoid.SetLocationLatLonAlt( 52.1315, -106.6335, 530.0)
   location = geoid.GetLocationXYZ();
   lat      = geoid.GetLatitude();
   lng      = geoid.GetLongitude();
   H        = geoid.GetAlt();


.. py:class:: ISKGeodetic()
   Creates an instance of ISKGeodetic

SetLocationLatLonAlt
^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKGeodetic.SetLocationLatLonAlt( lat, lon, alt) -> ok

      Sets the internal location to the point specified by the geodetic latitude, longitude and altitude.
      The 3 coordinates can be used to cover the entire 3-D space from the center of the Earth to infinity::

         ok = geoid.SetLocationLatLonAlt( lat, lon, alt )

      :param double lat:
         The geodetic latitude in degrees (-90 to + 90).

      :param double lon:
         The geodetic longitude in degrees (-180 to 360)

      :param double alt:
         The altitude above mean sea-level in meters.
         
      :param bool ok:
         Returns true if successful
         
      :return: returns true if successful

SetLocationXYZ
^^^^^^^^^^^^^^
.. py:method:: ISKGeodetic.SetLocationXYZ( geocentric ) -> ok

      Sets the internal location to the point specified by the geocentric vector::

         ok = geoid.SetLocationXYZ( geocentric )

      :param nxVector geocentric:
         The geocentric coordinates of the location specified in meters from the centre of the Earth.
         nxVector is represented as a 3 element array (X,Y,Z)

      :param bool ok:
         Returns true if successful
         
      :return: returns true if successful

SetLocationFromTangentPoint
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKGeodetic.SetLocationFromTangentPoint( observer, look ) -> ok
      
      Sets the internal location from the tangent point implied by the straight-line
      from the given *observer* in the direction of *look* away from the *observer*.
      The tangent point is the location where a straight line passing through 
      the observer’s location is parallel to the surface of the oblate Earth.
      The algorithm is robust for any observer looking in any direction as there 
      are always two points on the surface of the Earth parallel to the ray. 
      The tangent point is chosen to be the point closest to the ray. 
      Note that the tangent point calculated by the algorithm may be behind the observer::

         ok = geoid.SetLocationFromTangentPoint( observer, look)

      :param nxVector observer:
         The geocentric location of the observer in meters. nxVector is expressed as a 3 element array (X, Y, Z).
         
      :param nxVector lookvector:
         The look vector away from the observer in geocentric coordinates. This should be a unit vector. 
         nxVector is expressed as a 3 element array (X, Y, Z).
         
      :param bool ok:
         Returns true if successful
         
      :return: returns true if successful

SetLocationFromTangentAltitude
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKGeodetic.SetLocationFromTangentAltitude( required_altitude, observer, boresight_direction)-> ok,look

      A method useful for emulating satellite limb-scanning measurements.  The code sets the location 
      from the tangent point derived from the observer looking at a limb tangent at the *required_altitude*.
      The code returns the required look vector.
      
      At any given observer location there is a circle of solutions (or horizon) that are tangent
      at the required altitude. The tangent point is selected by choosing the point that intersects
      the plane formed by the boresight_direction and the observer position. In other words the 
      code looks for a tangent point in the bore-sight direction but allows the look vector 
      to rotate vertically up or down to meet the required altitude.

      The look vector is constrained to lie in a plane defined by the local vertical and
      the *boresight_direction*. The tangent point is the location where a straight line
      that passes through the observer’s location is parallel to the surface of the Earth.
      The algorithm is robust for any observer looking in any direction as there are always
      two points on the surface of the Earth parallel to the ray. The tangent point is chosen
      to be the point closest to the ray. Note that the tangent point calculated by
      the algorithm may be behind the observer::

         ok, required_look = geoid.SetLocationFromTangentAltitude( required_altitude, observer, boresight_direction)

      :param double required_altitude:
         The altitude of the required tangent point in meters above sea level.  The observer and the 
         returned *look* vector will be looking at this altitude.
         
      :param nxVector observer:
         The geocentric location of the observer in meters. The observer should be located 
         above the required altitude. nxVector is represented as a 3 element array.
         
      :param nxVector boresight_direction:
         The boresight direction of the instrument. The look vector towards the tangent point
         will lie in the plane defined by the local vertical and this bore-sight direction. This vector
         should not be in the lcoal vertical direction.
         
      :param nxVector required_look:
         Returns the look unitvector from the observer towards the tangent point at the desired altitude. This is
         a unitvector.
         
      :param bool ok:
         Returns true if successful
         
      :return: returns two element list (ok, look)

