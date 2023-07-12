
..  _emission_thermal:

THERMAL
=======
Implements a Planck black body spectral emission for any location in the atmopshere. It uses a climatology provided
by the user to obtain the temperature of the atmosphere.

Example
-------
::

   import sasktranif.sasktranif as skif
   sun = skif.ISKEmission('THERMAL')



Sasktran Availability
---------------------
This object is available from `sasktran` as class ``sasktran.EmissionThermal``. The object is implemented in the :ref:`sasktran_core`
extension and is part of the default sasktran installation.

Properties
----------

..  py:module::EMISSION_THERMAL

.. py:function:: emissivity( float value )

    Sets the emissivity of the Planck black body. Normally a value between 0 and 1.

.. py:function:: atmospheric_state (object climatology)

    Sets the climatology used to determine the temperature of the atmosphere.


..  _emission_userdefined:

USERDEFINED_WAVELENGTHHEIGHT
============================
Implements a user-defined emission object as a 2D array of wavelength (nm) and height (meters).

Example
-------
::

   import sasktranif.sasktranif as skif
   sun = skif.ISKEmission('USERDEFINED_WAVELENGTHHEIGHT')

Properties
-----------
..  py:module::EMISSION_USERDEFINED_WAVELENGTHHEIGHT

.. py:function:: heights( array heightm )

    Sets the values for the *height* axis of the 2-D array of emissions. The heights are given as a 1-D array in ascending
    order and are specified in meters. This property must be called before calling property `emissiontable`.

.. py:function:: wavelengths( array wavelen_nm )

    Sets the values for the *wavelength* axis of the 2-D array of emissions. The wavelengths are given as a 1D array in ascending order
    and are specified in nanometers. This property must be called before calling property `emissiontable`.

..  py:function:: emissiontable( array emissiontable )

    Sets the values for the 2-D array of emission values. This must be a contiguous array of values whose dimensions match
    sizes in previous calls to *heights* and *wavelengths*. The *wavelength* values changing the most rapidly in memory. The
    emission values are given as units of photons/cm2/sec/steradian/nm
