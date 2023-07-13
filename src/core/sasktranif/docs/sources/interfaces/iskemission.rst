.. _ISKEmission:

***********
ISKEmission
***********


.. py:class:: ISKEmission(name)
   
   A class used to represent thermal and photo-chemical emissions in the atmosphere. The ISKEmission class
   was originally created in early 2015 to model photochemical emissions in the upper mesosphere. 
   It is intended that the class could also be used to describe thermal and auroral emissions in the atmosphere. 
   The class only considers isotropic emissions::
   
      import sasktranif.sasktranif as skif
      
      o  = skif.ISKEmission('OH')
      
   :param str name:
      The name of the emission to be created. The name must correspond to 
      an installed sasktranif emission object. 

IsValidObject
^^^^^^^^^^^^^
 .. py:method:: ISKEmission.IsValidObject() -> ok

      Used to identify if the underlying C++ optical property is properly created.
      The function is primarily intended for internal usage::
      
         ok = climate.IsValidObject();
      
      :param boolean ok:
         The return value, true if the underlying C++ object is properly constructed otherwise false.

      :return: returns true if successful

UpdateLocation
^^^^^^^^^^^^^^
.. py:method:: ISKEmission.UpdateLocation( location, isgorund ) -> ok

      Sets the location that will be used to determine the isotropic radiance of
      the object at the next call to IsotropicRadiance. The interface provides the option to inform
      the emission object if this new location is a ground point. This can be used
      by the emission object to provide one isotropic radiance for the ground and
      another value for the air just above the ground. This is useful for modelling
      the ground discontinuity in thermal emissions::

         ok = emission.UpdateLocation( location, isground )

      :param GEODETIC_INSTANT location:
         The location for which the emission is required. The GEODETIC_INSTANT is a 4 element array [latitude, longitude, height_meters, mjd].
         Latitude and longitude are geodetic coordinates in degrees, height_meters is height above sea-level in meters and mjd is Modified Julian Date expressed in days.

      :param boolean isground:
         True if this location is a ground point. False if it is a point in the atmosphere.

      :param boolean ok:
         The return value, true if successful.

      :return: returns true if successful

UpdateCache
^^^^^^^^^^^
.. py:method:: ISKEmission.UpdateCache( location ) -> ok

      Provides the emission object an opportunity to update the caches on any climatologies it
      uses internally::

         ok = emission.UpdateCache( location, isground )

      :param GEODETIC_INSTANT location:
         The location for which the cache is required. The GEODETIC_INSTANT is a 4 element array [latitude, longitude, height_meters, mjd].
         Latitude and longitude are geodetic coordinates in degrees, height_meters is height above sea-level in meters and mjd is Modified Julian Date expressed in days.

      :param boolean ok:
         The return value, true if successful otherwise false.

      :return: returns true if successful

IsotropicRadiance
^^^^^^^^^^^^^^^^^
.. py:method:: ISKEmission.IsotropicRadiance( wavenumber ) -> ok,isotropicradiance

      Calculates the isotropic radiance emitted at the requested wavenumber at the location specified by the last call to
      SetAtmosphericState.  The radiance is returned as photons/sec/nm/ster/cm2. The incoming wavenumber can refer to either
      vacuum or air STP and this depends upon the context of the optical property object. UV-VIS optical property objects
      normally report cross-sections at STP while HITRAN /IR optical property object normally report cross-section in vacuum::

         ok,isotropicradiance= emission.IsotropicRadiance( wavenumber )

      :param float/array wavenumber:
         The wavenumber (cm-1) at which the cross-sections are required. If wavenumber is a
         scalar double then it will call the C++ scalar version. If wavenumber is an array
         then the code will call the C++ array version. The array version may be significantly more
         efficient for thermal emissions involving Voigt function lineshapes.
         
      :param float/array isotropicradiance:
         Returns the isotropic radiance emitted at the given location in photons/sec/nm/ster/cm2. May be NaN if there were errors.
         *isotropicradiance* type and size will match parameter wavenumber.
         
      :param boolean ok:
         The return value, true if successful otherwise false.

      :return: returns true if successful

SetProperty
^^^^^^^^^^^
.. py:method:: ISKEmission.SetProperty( propertyname, value) -> ok
   
      Set custom properties of the emission. The user must refer to 
      documentation about the specific emission object to see what properties it supports::
   
         ok = emission.SetProperty(propertyname, value)
   
      :param string propertyname:
         The name of the custom property to be modified.
   
      :param double/array/object/string value:
         The new value of the property.
         
      :param boolean ok:
         returns true if successful

      :return: returns true if successful

         
         

