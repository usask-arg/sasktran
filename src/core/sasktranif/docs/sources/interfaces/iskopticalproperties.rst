.. _ISKOpticalProperty:

******************
ISKOpticalProperty
******************
   The ISKOpticalProperty interface is used to describe the absorption,
   extinction and scattering cross-sections of individual atoms, molecules or particles and includes
   support for calculating the scattering phase matrix. The optical properties are often combined
   with number density climatologies to produce extinction, absorption and scattering coefficents within the atmossphere. Many cross-section
   depend upon the prevailing atmopsheric state and this aspect is supported by the interface

   Note,
   
   * Atmospheric thermal and airglow emissions are not supported by this interface but by interface ISKEmission. 
   * Exotic atmospheric parameters such as aerosol mode radius which typically apply to only one optical property object may be defined through the optprop.SetProperty interface o fthe specific otpical property object.

.. py:class:: ISKOpticalProperty(name)
   The ISKOpticalProperty constructor.

   :param str name:
      The name of the optical property to be created. The name must correspond to
      an installed sasktranif optical property. If no name is provided then a blank
      ISKOpticalProperty is created. This blank value can be useful for updating
      climatologies with ISKEngine method :meth:`~ISKEngine.AddSpecies`

AddUserDefined
^^^^^^^^^^^^^^
.. py:method:: ISKOpticalProperty.AddUserDefined( species, location) -> value

  A method to set user-supplied, custom cross-section data.  This method is only
  supported by optical property object USERDEFINED_TABLES. The object allows you to tabulate
  cross-section as a function of temperature and wavelength, similar to many UV/Visible cross-sections databases.
  By default::

     ok = userdefined_optprop.AddUserDefined( temperatures, wavelen_nm, crosssection) -> ok

  :param double temperature:
     The temperature in Kelvins to be associated with this cross-section data.

  :param array wavelen_nm:
     The array of wavelengths in nanometers of the cross-section data.

  :param array crosssection
     The array of cross-section data. Cross-sections must be provided as cm2 per atom/molecule/”average particle”.

  :param bool ok:
     returns true if successful

  :return: returns true if successful


CalculateCrossSections
^^^^^^^^^^^^^^^^^^^^^^

.. py:method:: ISKOpticalProperty.CalculateCrossSections(wavenumber)->(ok, absxs, extxs, scatxs)

  Calculates the absorption, extinction and scattering cross-sections of one atom, molecule
  or “average” particle of the optical property object at the wavenumber. The cross-sectional areas are returned
  as cm2/particle. Users can set the location at which the cross-section is required with a prior call to
  SetAtmosphericState.

  This method will calculates cross-sections for either a single scalar wavenumber or for an array of wavenumbers which must
  be in ascending wave-number order. The python and matlab code identifies if the incoming wavenumber is an array or not.

  The incoming wavenumber is currently ambiguous as it refers to vacuum for some optical property objects
  and to air at STP for others. We know this is an issue and we will address it in the future. As a
  rule of thumb the UV/Visible optical property objects usually report cross-sections at STP
  while HITRAN and IR optical property objects report cross-section in vacuum::

     ok,absxs, extxs, scatxs = optprop.CalculateCrossSections( wavenumber )

  :param double wavenumber:
     The scalar or ascending array of wavenumber (cm-1) at which the cross-sections are required. As discussed above the wavenumber may refer to
     air at STP or vacuum depending upon context.

  :param double absxs:
     Returns the absorption cross-section in cm2/(atom-molecule-particle) for the corresponding wavenumbers. May be NaN if there were errors.

  :param double extxs:
     Returns the extinction cross-section in cm2/(atom-molecule-particle) for the corresponding wavenumbers. May be NaN if there were errors.

  :param double scatxs:
     Returns the scattering cross-section in cm2/(atom-molecule-particle) for the corresponding wavenumbers. May be NaN if there were errors.

  :param bool ok:
     returns true if successful

  :return: returns the 4 element list [ok, absxs, extxs, scatxs]

InternalClimatology_UpdateCache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKOpticalProperty.InternalClimatology_UpdateCache( location )

    This is a speciality function that most users will not need to use. Updates the caches of any internal
    climatologies used by the optical property object. For example, ice and aerosol optical properties store
    internal climatologies of particle size distribution, typically mode radius and mode width. This method
    reloads the internal caches of the internal climatologies::

        ok = optprop.InternalClimatology_UpdateCache (location );

    :param GEODETIC_INSTANT location:
         The location of the new cache. GEODETIC_INSTANT is a 4 element array [latitude, longitude, height_meters, mjd]. Latitude
         and longitude are geodetic coordinates in degrees, height_meters is height above sea-level in meters and mjd is Modified
         Julian Date expressed in days.

    :param bool ok:
        Returns true if successful.

    :return: returns true if successful
   
IsValidObject
^^^^^^^^^^^^^
.. py:method:: ISKOpticalProperty.IsValidObject() -> ok
    :noindex:

    Used to identify if the underlying C++ optical property is properly created.
    The function is primarily intended for internal usage::

         ok = climate.IsValidObject();

    :param boolean ok:
         The return value, true if the underlying C++ object is properly constructed otherwise false.

    :return: returns true if successful

SetAtmosphericState
^^^^^^^^^^^^^^^^^^^
.. py:method:: ISKOpticalProperty.SetAtmosphericState( atmosphere ) -> ok

    Sets the ISKClimatology object that will be used as the background atmosphere in
    subsequent calls to :meth:`~ISKOpticalProperty.CalculateCrossSections`. Many optical property objects need object ``atmosphere`` to
    support atmospheric temperature and several need atmospheric pressure; in prinicple the list of required parameters is unlimited but in
    practice we have not yet encountered a need for anything other than **P** and **T**. However, in the final analysis it is the callers
    responsibility to ensure that the background atmospheric climatology supports all the parameters needed by the optical property object
    (you should see warnings if problems occur).

    It is the responsibility of the calling code to ensure the cache of the ``atmosphere`` object has been previously loaded with a call to
    UpdateCache before :meth:`~ISKOpticalProperty.SetLocation` or :meth:`~ISKOpticalProperty.CalculateCrossSections` are invoked. Note that
    :class:`ISKEngine` objects make these calls for the user. Optical property objects avoid invocations of :meth:`~ISKOpticalProperty.UpdateCache`
    for reasons of speed and efficiency and only call :meth:`~ISKOpticalProperty.GetParameter`::

            ok = optprop.SetAtmosphericState ( atmosphere )

    :param ISKClimatology atmosphere:
     The ISKClimatology object used to calculate atmospheric state parameters by the optical property object on subsequent calls to CalculateCrossSection.
     The ISKClimatology object will typically support calculation of pressure and temperature.

    :return: returns true if successful
      
   
SetLocation
^^^^^^^^^^^
.. py:method:: ISKOpticalProperty.SetLocation( location )

      Sets the location in the atmosphere for subsequent cross section calculations. This normally translates to the location used in
      the background atmosphere, see :meth:`ISKOpticalProperty.SetAtmopshericState`, to calculate **P** and **T**. Note that
      :class:`ISKEngine` objects make this call for the user::

         ok = optprop.SetLocation ( location )

      :param Tuple[lat,lng,heightm, mjd] location:
         The location and time in the atmosphere used to perform the next cross-section calculation. The location is a 4 element array [latitude, longitude, height_meters, mjd].
         Latitude and longitude are geodetic coordinates in degrees, height_meters is height above sea-level in meters and mjd is Modified Julian Date expressed in days.

      :return: returns true if successful

SetProperty
^^^^^^^^^^^
.. py:method:: ISKOpticalProperty.SetProperty( propertyname, value) -> ok
       :noindex:

      Set custom properties of the optical property object. The user must refer to
      documentation about the specific optical property object to see what properties it supports::

         ok = optprop.SetProperty(propertyname, value)

      :param string propertyname:
         The name of the custom property to be modified.

      :param double/array/object value:
         The new value of the property. The value must be a scalar double, array of doubles or a SasktranIF object

      :return: returns true if successful
