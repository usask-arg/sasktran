.. _clim_msis90:

MSIS90
======
A climatology that implements the MSIS-90 atmospheric model, (Hedin 1991). This is typically used in Sasktran applications
as a quick and robust background atmospheric state. The MSIS-90 model was built by the ionospheric/thermospheric community
using mass spectrometer and incohorent radar scatter data to provide background atmospheric state in ther thermosphere under
varying geomagntic conditions. The model was coupled to CIRA-86 (Chandra 1990 and Fleming 1990) to provide atmospheric state
for altitudes between the ground and ~120 km.

Most sasktran applications only require atmospheric state below 100 km and only
utilize the CIRA-86 part of the MSIS model. Thus we have configured the default implementation of the MSIS-90 model to only provide
the 3 basic atmospheric state parameters, pressure, temperature and number density and a fourth parameter, molecular oxygen number density,
between 0 km and 120 km.

Users may configure the MSIS-90 climatology using the objects properties outlined below to fetch six other species over
a larger height range with specific geomagnetic and solar flux conditions. Users are referered to Hedin's 1991 publication
for details on how to configure the parameters.  The current sasktran MSIS model does not provide any method to access the `TSELEC` function
describes in the fortran code.

Supported Species
-----------------
The MSIS-90 object supports up to 10 species which are listed in the table below. Four of the species are always loaded
while the remaining six are only loaded if explicitly requested by calling property `AddSpecies`. We recommend
that extra parameters are requested shortly after construction as part of the initialization process. Once a
species has been added to a specific MSIS instance it cannot be removed.

==================================  ================================== ==============   ================
Sasktran Handle                     Description                        Units            Availability
==================================  ================================== ==============   ================
SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3  Air Number density                 molecules/cm3    Always
SKCLIMATOLOGY_PRESSURE_PA           Pressure in                        Pascals          Always
SKCLIMATOLOGY_TEMPERATURE_K         Temperature                        Kelvin           Always
SKCLIMATOLOGY_O2_CM3                Molecular oxygen number density    molecules/cm3    Always
SKCLIMATOLOGY_O2_O2_CM6             O2-O2 square density for CIA       molecule^2/cm6   Always
SKCLIMATOLOGY_O_CM3                 Atomic oxygen number density       molecules/cm3    *AddSpecies*
SKCLIMATOLOGY_N2_CM3                Molecular nitrogen number density  molecules/cm3    *AddSpecies*
SKCLIMATOLOGY_N_CM3                 Atomic nitrogen number density     molecules/cm3    *AddSpecies*
SKCLIMATOLOGY_He_CM3                Helium number density              molecules/cm3    *AddSpecies*
SKCLIMATOLOGY_Ar_CM3                Argon number density               molecules/cm3    *AddSpecies*
SKCLIMATOLOGY_H_CM3                 Atomic hydrogen number density     molecules/cm3    *AddSpecies*
==================================  ================================== ==============   ================

Cache Snapshot
--------------
The MSIS model caches a single vertical profile when :meth:`ISKClimatology.UpdateCache` is called.  Consequently the MSIS
climatology object does not provide the horizontal variation of parameters within Sasktran calculations which only
extract values from the cached profile.

Python extension
----------------
The MSIS90 climatology is in the :ref:`sasktran_core` extension which is part of the default sasktran installation.

Configuration
-------------
The MSIS90 climatology needs no external preparation.

Properties
----------

.. py:module:: MSIS90

F10.7Avg
^^^^^^^^^^
.. py:function:: F10.7Avg(double flux) [Default: flux=150.0]

    Sets the 3 month average of the F10.7 flux. The default value is 150.0

F10.7
^^^^^
.. py:function:: F10.7(double flux) [Default: flux=150.0]

    Sets the daily F10.7 flux for the previous day. The default value is 150.0

Ap
^^
.. py:function:: Ap( array ap)[ Default: array[] = 4.0]

    A 7 element array that defines the prevailing geomagnetic conditions used by the MSIS model. The default value
    is a 7 element array with each element equal to 4.0

        #. Daily Ap.
        #. 3 hour Ap index for the current time.
        #. 3 hour Ap index for 3 hours before the current time.
        #. 3 hour Ap index for 6 hours before the current time.
        #. 3 hour Ap index for 9 hours before the current time.
        #. Average of eight 3 hour Ap indicies from 12 to 33 hours prior to current time.
        #. Average of eight 3 hour Ap indicies from 36 to 59 hours prior to current time.

MaxHeightKMS
^^^^^^^^^^^^
.. py:function:: MaxHeightKMS(double h)[Default h = 120.0]

    Sets the nominal maximum height that will calculated by the model. This option defaults to 120.0 which is a suitable
    value for most sasktran applications.

HeightSpacingKMS
^^^^^^^^^^^^^^^^
.. py:function:: HeightSpacingKMS(double s)[Default s = 1.0]

    Sets the spacing in kilometers between sample points internally stored by the model. THe default value is 1.0. There is
    probably very little benefit in changing this parameter from its default. It is provided for completeness. Once a
    species has been added to a specific MSIS instance it cannot be removed.

AddSpecies
^^^^^^^^^^
.. py:function:: AddSpecies( str handle )

    Requests that the MSIS climatology load the requested species. The variable `handle` must be one of the 10 supported
    variables and represented as a string, e.g. 'SKCLIMATOLOGY_O_CM3'. Note that handles are case sensitive.


Examples
--------

Using MSIS below 120 km
^^^^^^^^^^^^^^^^^^^^^^^^
In this example we show how to fetch the pressure at an altitude of 25 km at 50N, -102E on mjd 53000.45::

    import sasktranif as skif

    climate = skif.ISKClimatology('MSIS90')
    location = [50.0, -102.0, 25000.0, 53000.45]
    ok, pressure = climate.GetParameter( 'SKCLIMATOLOGY_PRESSURE_PA', location );


Using MSIS within the thermosphere
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this example we demonstrate how to use the MSIS climatology object to get height profiles of various atmospheric
components up to 300 km. The main  ::

    import sasktranif as skif
    import numpy as np
    import matplotlib.pyplot as plt

    skmsis   = skif.ISKClimatology('MSIS90')
    maxheight= 300.0
    f10p7    = 150.0
    f10p7avg = 150.0
    Ap       = [ 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0]

    skmsis.SetProperty( 'F10.7Avg',     f10p7)
    skmsis.SetProperty( 'F10.7',        f10p7avg)
    skmsis.SetProperty( 'Ap',           Ap)
    skmsis.SetProperty( 'MaxHeightKMS', maxheight)
    skmsis.SetProperty( 'AddSpecies',  'SKCLIMATOLOGY_O_CM3')       # Request that MSIS loads these species.
    skmsis.SetProperty( 'AddSpecies',  'SKCLIMATOLOGY_He_CM3')
    skmsis.SetProperty( 'AddSpecies',  'SKCLIMATOLOGY_N2_CM3')
    skmsis.SetProperty( 'AddSpecies',  'SKCLIMATOLOGY_Ar_CM3')
    skmsis.SetProperty( 'AddSpecies',  'SKCLIMATOLOGY_H_CM3')
    skmsis.SetProperty( 'AddSpecies',  'SKCLIMATOLOGY_N_CM3')


    location = [52.0, -102.0, 0.0, 57005.0]
    h       = np.arange( 300)
    skmsis.UpdateCache( location )

    ok,O  = skmsis.GetHeightProfile('SKCLIMATOLOGY_O_CM3',  location, h*1000.0 )
    ok,O2 = skmsis.GetHeightProfile('SKCLIMATOLOGY_O2_CM3', location, h*1000.0 )
    ok,He = skmsis.GetHeightProfile('SKCLIMATOLOGY_He_CM3', location, h*1000.0 )
    ok,N2 = skmsis.GetHeightProfile('SKCLIMATOLOGY_N2_CM3', location, h*1000.0 )
    ok,Ar = skmsis.GetHeightProfile('SKCLIMATOLOGY_Ar_CM3', location, h*1000.0 )
    ok,N  = skmsis.GetHeightProfile('SKCLIMATOLOGY_N_CM3',  location, h*1000.0 )
    ok,H  = skmsis.GetHeightProfile('SKCLIMATOLOGY_H_CM3',  location, h*1000.0 )
    ok,ND = skmsis.GetHeightProfile('SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3',  location, h*1000.0 )
    ok,P  = skmsis.GetHeightProfile('SKCLIMATOLOGY_PRESSURE_PA',  location, h*1000.0 )
    ok,T  = skmsis.GetHeightProfile('SKCLIMATOLOGY_TEMPERATURE_K',  location, h*1000.0 )

    plt.figure(1)
    plt.plot( np.log10(O),  h,  'k.-')
    plt.plot( np.log10(O2), h,  'r.-')
    plt.plot( np.log10(He), h,  'g.-')
    plt.plot( np.log10(N2), h,  'b.-')
    plt.plot( np.log10(Ar), h,  'c.-')
    plt.plot( np.log10(N),  h,  'm.-')
    plt.plot( np.log10(H),  h,  'y.-')
    plt.plot( np.log10(ND), h,  'y.-')
    plt.plot( np.log10(P),  h,  'r.-')
    plt.plot( np.log10(T),  h,  'g.-')


    plt.legend( ['O','O2','He','N2','Ar','N','H','Num Density','Pressure','Temperature'] )
    plt.ylabel('Height Kms')
    plt.xlabel('Log10(Number density)');
    plt.title('MSIS Example height profiles')
    plt.show()



References
----------
**Fleming, E.L.**, S. Chandra, J.J. Barnett and M. Corney (1990), Zonal mean temperature, pressure, zonal wind, and geopotential height as functions of latitude, COSPAR International Reference Atmosphere: 1986, Part II: Middle Atmosphere Models,
*Adv. Space Res.*, **10**, 12, 11-59, `doi:10.1016/0273-1177(90)90386-E <http://dx.doi.org/10.1016/0273-1177%2890%2990386-E>`_.

**Chandra, S.**, E.L. Fleming, M.R. Schoeberl, J.J. Barnett, (1990), Monthly mean global climatology of temperature, wind, geopotential height and pressure for 0–120 km,
*Advances in Space Research*, **10**, 6, 3-12, `doi.org/10.1016/0273-1177(90)90230-W <https://doi.org/10.1016/0273-1177(90)90230-W>`_.

**Hedin, A. E.** (1991), Extension of the MSIS Thermosphere Model into the middle and lower atmosphere, *J. Geophys. Res.*, **96** ( A2), 1159– 1172, `doi:10.1029/90JA02125 <https://doi.org/10.1029/90JA02125>`_

