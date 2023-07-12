.. _ISKSolarSpectrum:

****************
ISKSolarSpectrum
****************

.. py:class:: ISKSolarSpectrum(name)
   
   The ISKClimatology interface is used to describe climatologies used in the
   radiative transfer models. The climatologies typically specify the number
   density of a given optically active species at various locations in the Earth’
   s atmosphere and when combined with the cross-section data from the
   ISKOpticalProperty interface allows the optical depth of atmospheric rays to
   be calculated.
   
   The ISKClimatology interface is much more flexible than just
   providing number densities, it can be used and extended to provide any scalar
   value at any location in the atmosphere. In addition one ISKClimatology object
   can provide more than one scalar, a frequent example is that models supporting
   atmospheric state will provide temperature, pressure and number density.
   
   The ISKClimatology objects implement the concept of caching. Almost all
   climatological models and databases incur some overhead penalty when the user
   requests information from different times and locations in the atmosphere. In
   some models, such as user supplied tables, the overhead is insignificant but
   in other models, such as raw ECMWF files, the overhead can be very large. The
   ISKClimatology interface manages this issue by requesting that each
   ISKClimatology model update its internal cache before returning any
   atmospheric parameters. The size and nature of the cache is left to specific
   ISKClimatology object to decide but it is expected that all caches will store
   at least a vertical profile at the requested location::
   
      import sasktranif.sasktranif as skif
      
      sun  = skif.ISKSolarSpectrum('SAO2010')
   
      :param str name:
      The name of the solar spectrum to be created. The name must correspond to 
      an installed sasktranif solar spectrum. 

IsValidObject
^^^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.IsValidObject() -> ok

      Used to identify if the underlying C++ solar spectrum is properly created.
      The function is primarily intended for internal usage::
      
         ok = climate.IsValidObject();
      
      :param boolean ok:
         The return value, true if the C++ object is properly constructed otehrwise false.

      :return: returns true if successful

Irradiance
^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.Irradiance( wavelengths) -> ok,irradiance
   
      Return the solar irradiance at the top of the atmosphere at the current solar distance
      for the given wavelengths. The wavelengths are specified in nanometers in vacuum.  The
      current solar distance is set by calls to :meth:`~ISKSolarSpectrum.SetSolarDistanceFromMjd`
      or :meth:`~ISKSolarSpectrum.SetSolarDistanceFromAU`::
      
         ok,irradiance = sun.Irradiance( wavelengths)
   
      :param scalar/array wavelengths:
         The values (scalar or array) of wavelengths at which to calculate the top of the atmosphere
         solar irradiance.

      :param scalar/array irradiance:
         Returns the values of irradiance at the requested wavelength in photons/cm2/sec/nm.
         
      :param boolean ok:
         returns true if successful

      :return: returns the 2 element list (ok, irradiance)

IrradianceAt1AU
^^^^^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.IrradianceAt1AU( wavelengths) -> ok,irradiance

      Return the solar irradiance at the top of the atmosphere at 1 AU for the given wavelengths.
      The wavelengths are specified in nanometers in vacuum::

         ok,irradiance = sun.IrradianceAt1AU( wavelengths)

      :param scalar/array wavelengths:
         The values (scalar or array) of wavelengths at which to calculate the 1 AU solar
         irradiance.

      :param scalar/array irradiance:
         Returns the values of irradiance at the requested wavelength in photons/cm2/sec/nm.

      :param boolean ok:
         returns true if successful

      :return: returns the 2 element list (ok, irradiance)

SetSolarDistanceFromMjd
^^^^^^^^^^^^^^^^^^^^^^^
 .. py:method:: ISKSolarSpectrum.SetSolarDistanceFromMjd( mjd) -> ok

      Set the distance of the Earth from the Sun using the modified julian date given in *mjd*.
      This distance will be used in subsequent calls to `~ISKSolarSpectrum.Irradiance` to calculate 
      the top-of-the-atmosphere irradiance. No attempt is made to correct for the distance of the
      observer/TOA from the centre of the Earth.:: 

         ok = sun.SetSolarDistanceFromMjd (mjd )

      :param double mjd:
         The modified Julian Date used to determine the position of the Earth.

      :param boolean ok:
         Returns true if successful otherwise false

      :return: Returns true if successful otherwise false

SetSolarDistanceFromAU
^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.SetSolarDistanceFromAU( au) -> ok

      Set the distance of the Earth from the Sun to the number of astronomical units specified by ``au``.
      This distance will be used in subsequent calls to `~ISKSolarSpectrum.Irradiance` to calculate 
      the top-of-the-atmospher irradiance. One astronomical unit is the mean distance of the Earth from the Sun. :: 

         ok = sun.SetSolarDistanceFromAU (au)

      :param double au:
         The distance of the Earth from the Sun in astronomical units.

      :param boolean ok:
         Returns true if successful otherwise false

      :return: Returns true if successful otherwise false

MinValidWavelength
^^^^^^^^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.MinValidWavelength() -> ok, minwavelength

      Returns the minimum wavelength supported by this solar spectrum::

         ok,minwavelen = sun.MinValidWavelength()

      :param double minwavelen:
         Returns the minimum wavelength in nanometers in vacuum supported by this
         solar spectrum object.

      :param boolean ok:
         Returns true if successful otherwise false

      :return: Returns the list (ok, minwavelen)

MaxValidWavelength
^^^^^^^^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.MaxValidWavelength() -> ok, minwavelength

      Returns the maximum wavelength supported by this solar spectrum::

         ok,maxwavelen = sun.MaxValidWavelength()

      :param double maxwavelen:
         Returns the maximum wavelength in nanometers in vacuum supported by this
         solar spectrum object.

      :param boolean ok:
         Returns true if successful otherwise false

      :return: Returns the list (ok, minwavelen)

NanometerResolutionFWHM
^^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.NanometerResolutionFWHM( wavelengths) -> ok,fwhm

      Return the spectral resolution of this solar spectrum object at the specified wavelengths::

         ok,fwhm = sun.NanometerResolutionFWHM( wavelengths)

      :param scalar/array wavelengths:
          The values (scalar or array) of wavelengths in nanometers at which to calculate the spectral resolution o fthe solar spectrum.

      :param scalar/array fwhm:
          Returns the spectral resolution of the solar spectrum at the specified wavelengths as Full-Width-Half Maximum expressed in nanometers.

      :param boolean ok:
          returns true if successful

      :return: returns the 2 element list (ok, fwhm)

SetProperty
^^^^^^^^^^^
.. py:method:: ISKSolarSpectrum.SetProperty( propertyname, value) -> ok

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

