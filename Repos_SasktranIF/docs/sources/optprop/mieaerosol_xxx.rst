
.. _optical_mieaerosol:

MIE AEROSOL
===========
Calculate the optical properties of spherical aerosol particles. The class models the sulphates
as a log-normal distribution. The user must provide an ISKClimatology of the log-normal parameters

* SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS()
* SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH()

The log normal parameters climatology must be given to the MIEAEROSOL_XXX object using 
:meth:`SetProperty( ‘SetParticleSizeClimatology’, object) <ISKOpticalProperty.SetProperty>`

All `MIEAEROSOL` optical properties function similarly, the only difference is the species index
of refraction.

==================================  ===========
 SasktranIF Name                    Notes
==================================  ===========
MIEAEROSOL_H2SO4                                   
MIEAEROSOL_DUST                     
MIEAEROSOL_WATER
MIEAEROSOL_ICE
==================================  ===========

Local Caching and Slow Performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Calculation of optical properties of spherical particle distributions using Mie code
is numerically intensive and will be quite slow. Consequently, the MIEAEROSOL_XXX object 
cache optical properties for different aerosol parameters (mode radius, mode width and wavelength) 
on the local machine.

Note that the code will be slow each time the code encounters a new particle distribution or wavelength configuration.
The location of the cache folder can be found on Windows and Linux machines at the registry entry
HKEY_LOCAL_MACHINE\SOFTWARE\Usask-ARG\SkOpticalProperties\MieAerosol\Aerosol_Cache_Directory.


Example
^^^^^^^
::

   optprop     = ISKOpticalProperty('MIEAEROSOL_H2SO4')
   modeclimate = ISKClimatology('USERDEFINED_PROFILE');
   modeclimate.SetPropertyArray('Heights', [0, 100000] );
   modeclimate.SetPropertyUserDefined( SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS(), [0.08, 0.08]);
   modeclimate.SetPropertyUserDefined( SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH(),          [1.6,  1.6]);
   aerosol.SetProperty('SetParticleSizeClimatology', modeclimate);
   msis = ISKClimatology('MSIS90');
   mjd = 52393.3792987115;
   location = [0.0, 0.0, 25000.0, mjd];
   optprop.SetAtmosphericState( msis)
   optprop.SetLocation(location)
   [ok, abs, ext,sca] = optprop.CalculateCrossSections( 1.0E7/600.0, location );


Properties
^^^^^^^^^^
.. option:: SetParticleSizeClimatology (object modeclimate)
   
   Sets the climatology used to provide the mode radius and mode width of the aerosol log normal distribution
   at any point in the atmosphere. The user will always set this property. This property should be set as part 
   of the initialization before either :meth:`~ISKOpticalProperty.SetAtmosphericState`, 
   :meth:`~ISKOpticalProperty.SetLocation` or :meth:`~ISKOpticalProperty.CalculateCrossSections` are called. 
   The climatology must support
   
   * SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS()
   * SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH()



