.. _changelog:

Changelog
*********
1.7.1
=====
* Development of the model has moved over to github, see https://github.com/usask-arg/sasktran
* Fixed several bugs involving running sasktran on non-anaconda based python environments

1.7.0
=====
* Added support for python 3.10/3.11 builds on both Windows and Linux
* Added capability for HR to calculate polarized weighting functions

1.6.0
=====
* Added the sk.UserDefinedAbsorptionPressure optical property which is a user defined optical property as a function of pressure, temperature, wavelength, and a broadener
* Added the sk.UserDefinedScatterConstantHeight optical property which is used to set scattering properties for the SASKTRAN-DO engine.
* Added a few osculating spheroid functions and the ability for the NadirGeometry class to adjust lines of sight based on the geoid used
* Fix a bug with HRSSApprox and HITRAN
* Considerably speed up the loading of the Baum table
* Added internal functionality to allow for better interfacing between the optical properties and SASKTRAN-DO
* Added Mie interface wrapper, still to be documented.
* Add option to BaumSpecies to set the name.
* Added some internal methods in SASKTRAN-HR to better specify diffuse profile locations

1.5.0
=====
* Add the ability for HR to calculate weighting functions for multiple species simultaneously
* Added a spectral varying BRDF option
* Added the option for HR to calculate weighting functions with respect to surface reflectance when directly viewing the surface

1.4.0
=====
* Fix HR WF calculation for scattering species when operating in vector mode
* Add an option to the HR engine to calculate cell optical depths
* Fix a bug where the default settings for HR vector mode were not precise
* Multithreading improvements when running on a large number of threads
* Added EngineHRSSapprox which uses HR to calculate radiance using an SS approximation
* Added support for Python 3.9 (still in experimental phases)
* Added documentation for the occultation and monte carlo engines in EngineOCC and EngineMC specifically
* Added support for calculating box AMFs with the MC engine
* Added support for inelastic Ramen scattering with the MC engine 
* Added experimenal support for solar refraction in the HR engine
* Improved support and documentation for Baum ice crystal calculations in HR

1.3.0
=====
* Build the package for Python 3.8
* Large update to HITRAN internals to speed up calculations
* Updates to the internal MSIS90 model to extend altitude range up to 200-300 km
* Added O2O2 collisional cross sections
* Added an experimental option to disable single scattering sources along the LOS.
* Bug fixes for the Baum ice crystal database

1.2.0
=====
* Large update on the underlying package orientation, sasktranif is now part of the sasktran wheel.
* Update HITRAN to 2016
* Add options to download HITRAN and PRATMO binary files

1.1.2
=====
* Fixed support for older pyyaml versions due to a regression in 1.1.1


1.1.1
=====
* Fixed a bug where the polarization basis returned back by both MC and HR engines were not correct
* Fixed a deprecation warning in pyyaml 5.1
* Fixed a deprecation warning in numpy 1.16


1.1.0
=====
* Added BaumIceCrystal optical property and SpeciesBaumIceCloud species
* Fix a bug where output_format='xarray' was not working correctly when engine.wavelengths was a scalar
* Added A Band emission model


1.0.0
=====
* Added refraction option to HR
* SolarSpectrum class can now convolve to lower resolution
* Update required SASKTRANIF version to 4.3.1
* Added SimpleRayleigh optical property for use in comparisons with other models
* Added a function to OpticalProperty to output the phase matrix
* Added an interface for polarization within HR
* Package is now supported and tested on Python 3.7


0.2.0
=====
* Added examples showing how to obtain diagnostic information from SASKTRAN-HR
* Added support for basic thermal emissions
* Improved the xarray output radiance structure when weighting functions are calculated
* Updated required SASKTRANIF to version 4.2.6, SKTRAN_Disco to 0.1.12

0.1.3
=====
* Fixed a bug with the GloSSAC species not working correctly with SASKTRAN
* Updated SASKTRANIF to version 4.2.5