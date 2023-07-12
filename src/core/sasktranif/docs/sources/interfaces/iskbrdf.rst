.. _ISKBrdf:

*******
ISKBrdf
*******

.. py:class:: ISKBrdf(name:str)
   
   The ISKBrdf interface provides implementations of Bidirectional Reflectance
   Distribution Function. The brdf is currently implemented in a scalar
   form and for the moment the polarized version should be assumed 
   to be unpolarised which is reasonable for the surface reflectances implemented to date (2016-12-20) 
   but would not be true for special, surface reflectances such as sea-glint.
   
   An example is given below::
   
      import sasktranif.sasktranif as skif
      
      sun  = skif.ISKBrdf('Lambertian')
   
      :param str name:
      The name of the BRDF object to be created. The name must correspond to 
      an installed sasktranif BRDF. 

IsValidObject
^^^^^^^^^^^^^
.. py:method:: ISKBrdf.IsValidObject() -> ok

      Used to identify if the underlying C++ ISKBRDF is properly created.
      The function is primarily intended for internal usage::
      
         ok = brdfobject.IsValidObject();
      
      :param boolean ok:
         The return value, true if the C++ object is properly constructed otehrwise false.

      :return: returns true if successful

BRDF
^^^^
.. py:method:: ISKBrdf.BRDF( wavelen_nm, pt, MU_in, MU_out, COSDPHI) -> ok, brdf
   
      Returns the (scalar) bidirectional reflectance distribution function at given wavelength and location
      for a given incoming and outgoing direction. solar irradiance at the top of the atmosphere at the current solar distance
      for the given wavelengths::
      
         ok,brdf = brdfobject.BRDF( wavelength, pt, cosincoming, cosoutgoing, cosazimuth)
   
      :param float wavelen_nm:
         The wavelength in nanonmeters at which to calculate the BRDF.

      :param GEODETIC_INSTANT location:
         The location and instant in time for which the BRDF is required. The GEODETIC_INSTANT is a
         4 element array [latitude, longitude, height_meters, mjd]. Latitude and longitude are
         geodetic coordinates in degrees, height_meters is height above sea-level in meters and mjd
         is Modified Julian Date expressed in days.

      :param float MU_in:
         The cosine of the zenith angle of the incoming ray. Note we consider the angle between
         the local vertical and the outbound direction of the ray, i.e. this will be a value between 
         0 and 1 for regular surfaces.
         
      :param float MU_out:
         The cosine of the zenith angle of the scattered ray. Note we consider the angle between
         the local vertical and the outbound direction of the ray, i.e. this will be a value between 
         0 and 1 for regular surfaces.
         
      :param float COSDPHI:
         The cosine of the azimuthal angle between the incoming and outgoing ray.  Note that the
         azimuthal angle is between the incoming ray **pointing into** the surface and the outgoing ray 
         **pointing out** of the surface.  Forward scatter occurs at an azimuth of 0 degrees (cosine = 1)
         and backscatter occurs at an azimuth of 180 degrees (cosine = -1). COSDPHI  will have a value 
         between -1 and 1. 
      
         
      :param boolean ok:
         returns true if successful

      :return: the value of the brdf for this configuration.
      :rtype: float

SetProperty
^^^^^^^^^^^
.. py:method:: ISKBrdf.SetProperty( propertyname, value) -> ok

      Set custom properties of the solar spectrum object. The user must refer to 
      documentation about the specific solar spectrum object to see what properties it supports::

         ok = sun.SetProperty(propertyname, value)

      :param string propertyname:
         The name of the custom property to be modified.

      :param double/array/object value:
         The new value of the property. The value must be a scalar double, array of doubles or a SasktranIF object

      :param boolean ok:
         returns true if successful

      :return: returns true if successful

