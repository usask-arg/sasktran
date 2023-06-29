.. _optical_userdefined:


USERDEFINED_TABLES
==================
A special optical property object that allows a user to define their own tables of cross-section with wavelength and temperature.
The method ISKOpticalProperty::AddUserDefined is valid for this object and can be used to define tables of absorption
(or scattering) cross-sections as a function of wavelength. Each table is specified for a given temperature. Users should ensure
the tables span the useful range of wavelengths as cross-sections outside the wavelength ranges are set to zero. Cross-sections 
requested outside the range of temperatures are truncated to the nearest temperature table.

Properties
^^^^^^^^^^
.. option:: SetTemperature (double n)
   
   Sets the temperature in Kelvins that will be used in the next calculation of cross-sections.

.. option:: SetIsScatterer (int n)

   A value of zero will set the object as a pure absorber. A non-zero value will set the object as a pure scatterer.

