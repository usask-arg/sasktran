.. _ISKClimatology:

***************
ISKClimatology
***************
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
   at least a vertical profile at the requested location.

.. py:class:: ISKClimatology(name)

      The ISKClimatology constructor::

      >>> import sasktranif.sasktranif as skif
      >>>
      >>> climate  = skif.ISKClimatology('MSIS90')

   :param str name:
      The name of the climatology to be created. The name must correspond to
      an installed sasktranif climatology.


IsValidObject
^^^^^^^^^^^^^
.. py:method:: ISKClimatology.IsValidObject() -> ok

      Used to identify if the underlying C++ climatology is properly created.
      The function is primarily intended for internal usage::

         ok = climate.IsValidObject();

      :param boolean ok:
         The return value, true if the C++ object is properly constructed otehrwise false.

      :return: returns true if successful

UpdateCache
^^^^^^^^^^^
.. py:method:: ISKClimatology.UpdateCache(location) -> ok
   
      Instructs the climatology to update its internal cache given a specific location and time.
      All subsequent calls to GetParameter will use the internal cache. The size and structure
      of the internal cache varies from one climatology to another. Simple user-defined height
      profiles, for example, don’t do anything in response to UpdateCache request as the 
      height profile “is” the cache. At the other end of the spectrum ECMWF climatologies load
      in two snapshots of the entire globe that straddle the time of interest at all pressure
      levels within the model. Many intermediate climatologies load in a single profile at the
      requested location and time and it is desirable that single profile is used for all
      subsequent calls to GetParameter::
      
         ok = climate.UpdateCache ( location );
   
      :param GEODETIC_INSTANT location:
         The location and instant in time for which the cache is required. The GEODETIC_INSTANT is a
         4 element array [latitude, longitude, height_meters, mjd]. Latitude and longitude are
         geodetic coordinates in degrees, height_meters is height above sea-level in meters and mjd
         is Modified Julian Date expressed in days.
         
      :param boolean ok:
         returns true if successful

      :return: returns true if successful
      
SetProperty
^^^^^^^^^^^
.. py:method:: ISKClimatology.SetProperty( propertyname, value) -> ok
   
      Set custom properties of the climatology. The user must refer to 
      documentation about the specific climatology object to see what properties it supports::
   
         ok = climate.SetProperty(propertyname, value)
   
      :param string propertyname:
         The name of the custom property to be modified.
   
      :param double/array/object value:
         The new value of the property. The value must be a scalar double, array of doubles or a SasktranIF object
         
      :param boolean ok:
         returns true if successful

      :return: returns true if successful

GetParameter
^^^^^^^^^^^^
.. py:method:: ISKClimatology.GetParameter( species, location) -> value
      
      Fetch the value of the *species* at *location*. The climatology instance must support the requested *species*
      and it will extract the value from the cached data (see UpdateCache):: 

         ok,value = climate.GetParameter (species, location)

      :param string  species:
         The identification code of the required species. The climatology must support this species.
      
      :param GEODETIC_INSTANT location:
         The location in space and time for which the value is required.  See UpdateCache for a description of
         GEODETIC_INSTANT. Note that the value will be extracted by appropriate interpolation of the existing cache
         but does not update the current cache.
         
      :param boolean ok:
         The first value of the returned list. true if successful otehrwise false

      :param double value:
         The second value of the returned list. The value of the parameter at the requested location.
         Climatologies may return this value as NaN if the interpolation request was inappropriate.
       
      :return: A two element list [*ok, value*]
      
         
GetHeightProfile
^^^^^^^^^^^^^^^^
.. py:method:: ISKClimatology.GetHeightProfile( species, location, altitudes) -> value
      
      Fetch a height profile of the *species* at *location*. The climatology instance must support the requested *species*
      and it will extract the value from the cached data (see UpdateCache):: 

         ok,profile = climate.GetParameter (species, location, altitudes)

      :param string  species:
         The identification code of the required species. The climatology must support this species.
      
      :param GEODETIC_INSTANT location:
         The location in space and time for which the value is required.  See UpdateCache for a description of
         GEODETIC_INSTANT. Note that the value will be extracted by appropriate interpolation of the existing cache
         but does not update the current cache.
         
      :param 1d-array altitudes:
         An array of altitudes expressed in meters. The code will resturn the value of the species on the given altitude grid. 
         See ::meth:`UpdateCache` for a description of GEODETIC_INSTANT. Note that the value will be extracted by appropriate 
         interpolation of the existing cache but does not update the current cache.
       
      :return: A two element list [*ok, value*]
      
          **ok** True if successful
          
          **value** The height profile of the requested species on the given altitudes. This is returned either as a numpy array if altitudes was an array
          or as a scalar float if altitudes was a scalar float
         
Create_New_ClimatologyName
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKClimatology.Create_New_ClimatologyName( speciesname ) -> ok

      Allows a user to define a new species name that can be used in subsequent calls to SetPropertyUserDefined. This
      allows users to make user-defined climatologies of species which are not already defined inside the SasktranIF libraries::

         ok = climate.Create_New_ClimatologyName( speciesname)
      
      :param string  speciesname:
               The name of the new species to be added to the internal library of supported species. Note this name is added to a global library
               inside SasktranIF and is avaialble to all SasktranIF objects. Note that new entries are only made if the species name does not 
               currently exist in the internal library. The species name is case insensitive
      
      :return: returns true if successful

      
SetPropertyUserDefined
^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKClimatology.SetPropertyUserDefined( species, location) -> value

      A special property reserved for user-defined climatoglies
      that allows users to specify their own height profiles for a given species.  This functionality
      is required for retrieval algorithms. You typically have to call SetProperty('Heights',h), see :ref:`clim_userdefined_profile`, before
      calling this function to set the required height grid.::
      
         ok = climate.SetPropertyUserDefined( speciesname, profile)
         
      :param string speciesname:
         The name of the custom profile to be modified. Note that the names must exist in the internal global table (to avoid typing errors).  New names
         can be added to the global table with method :meth Create_New_ClimatologyName,

      :param array profile:
         The array of values for the climatology of locations previously defined by appropriate calls to :meth SetProperty( 'Heights' )

      :param boolean ok:
         returns true if successful

      :return: returns true if successful

      An example::
      
         climate = skif.ISKClimatology(‘USERDEFINED_PROFILE’)
         h = (0:100)*1000.0;
         v = h.*h + 1.0;
         climate.SetPropertyScalar('DoLogInterpolation', 1);
         climate.SetPropertyArray( 'Heights', h );
         climate.SetPropertyUserDefined( skif.SKCLIMATOLOGY_O3_CM3(), v);
         ok, value = climate.GetParameter( skif.SKCLIMATOLOGY_O3_CM3(), [50,102, 35000.0, 53000]);
