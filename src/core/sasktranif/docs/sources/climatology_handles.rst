.. _climatologyhandles:

*************************
Climatology Handles/GUID
*************************

Climatology handles are a mechanism used in the Sasktran framework to label physical quantities such as pressure,
temperature and number density. Standard names are used so that components, such as climatology objects, can be used
interchangeably within the framework. Within the C++ Sasktran code the climatology handles
are implemented as 128 bit structures called GUIDs which are guaranteed to be `globally unique identifiers <https://www.guidgenerator.com/>`_
while SasktranIF implements the GUIDs as user-friendly strings. Keep in mind that most climatology objects only support a
small sub-set of the physical quantities listed below, see the documentation for each climatology object.

In Python the climatology handle is simply the name of the handle as a string, eg, "SKCLIMATOLOGY_PRESSURE_PA". The string
is not arbitrary and python users cannot make up their own values as the string is internally converted by the C++ code
to the GUID representation, errors will occur if the string cannot be converted.

Users can only use the GUIDs listed below, they are not yet able to define their own GUIDs within the framework.

Atmospheric State GUIDS
-----------------------
    ===========================================  =====================================  =======================
    GUID String                                  Species                                Units
    ===========================================  =====================================  =======================
    SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3           Air number density                     Molecules per cm3.
    SKCLIMATOLOGY_PRESSURE_PA                    Pressure                               Pascals
    SKCLIMATOLOGY_TEMPERATURE_K                  Temperature                            Kelvins
    SKCLIMATOLOGY_POTENTIAL_TEMPERATURE_K        Potential temperature                  Kelvins
    SKCLIMATOLOGY_ALBEDO                         Albedo                                 No units
    SKCLIMATOLOGY_EPV                            Ertels Potential Vorticity
    SKCLIMATOLOGY_GEOMETRIC_HEIGHT				 Geometric height.
    SKCLIMATOLOGY_GEOPOTENTIAL_HEIGHT		     Geopotential height.
    SKCLIMATOLOGY_SURFACE_GEOPOTENTIAL_HEIGHT    Surface Geopotential Height
    SKCLIMATOLOGY_SURFACE_PRESSURE_PA            Surface Pressure
    SKCLIMATOLOGY_CLOUD_FRACTION                 Cloud fraction
    SKCLIMATOLOGY_QI_MMR                         Mass fraction of cloud ice water
    SKCLIMATOLOGY_QL_MMR                         Mass fraction of cloud liquid water
    SKCLIMATOLOGY_QV                             Specific Humidity
    SKCLIMATOLOGY_RH                             Relative Humidity
    SKCLIMATOLOGY_AOA_DAYS                       Age of Air                             Days
    ===========================================  =====================================  =======================

Aerosol GUIDS
-------------

    ==========================================      =====================================   ===============================================
    GUID String                                     Species                                 Units
    ==========================================      =====================================   ===============================================
    SKCLIMATOLOGY_AEROSOL_CM3                       Aerosol.                                Mean particle number density per cm\ :sup:`3`\
    SKCLIMATOLOGY_AEROSOLH2SO4_CM3                  Sulphate aerosol.                       Mean particle number density per cm\ :sup:`3`\
    SKCLIMATOLOGY_AEROSOLDUST_CM3                   Dust aerosol.                           Mean particle number density per cm\ :sup:`3`\
    SKCLIMATOLOGY_AEROSOLICE_CM3                    Ice aerosol.                            Mean particle number density per cm\ :sup:`3`\
    SKCLIMATOLOGY_AEROSOLWATER_CM3                  Water aerosol                           Mean particle number density per cm\ :sup:`3`\
    SKCLIMATOLOGY_ICE_CM3                           Ice                                     Ice crystal number density   per cm\ :sup:`3`\
    SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM           Extinction                              Extinction per km
    SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS      Log-Normal Mode Radius                  Microns
    SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH               Log-Normal Mode width                   No Units
    SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS             Aerosol Effective Size                  Microns
    SKCLIMATOLOGY_AEROSOLSURFACEAREA_UM2PerCM3      Aerosol Surface Area                    um\ :sup:`2`\ cm\ :sup:`-3`\
    SKCLIMATOLOGY_DUST_0p7mu                        Dust Aerosol, R\ :sub:`eff`\ 0.7 um
    SKCLIMATOLOGY_DUST_1p4mu                        Dust Aerosol, R\ :sub:`eff`\ 1.4 um
    SKCLIMATOLOGY_DUST_2p4mu                        Dust Aerosol, R\ :sub:`eff`\ 2.4 um
    SKCLIMATOLOGY_DUST_4p5mu                        Dust Aerosol, R\ :sub:`eff`\ 4.5 um
    SKCLIMATOLOGY_BCPI                              Geos Chem Black Carbon I
    SKCLIMATOLOGY_BCPO                              Geos Chem Black Carbon O
    SKCLIMATOLOGY_SALA                              Geos Chem Sea Salt aerosol Accum
    SKCLIMATOLOGY_SALC                              Geos Chem Sea Salt aerosol Coarse
    SKCLIMATOLOGY_OCPI                              Geos Chem Organic Carbon Aerosol I
    SKCLIMATOLOGY_OCPO                              Geos Chem Organic Carbon Aerosol O
    ==========================================      =====================================   ===============================================


Molecular Number Density GUIDS
------------------------------

    ============================ ======================================================  ==============================
    GUID String                  Species                                                 Units
    ============================ ======================================================  ==============================
    SKCLIMATOLOGY_Ar_CM3         Argon                                                   Atoms per cm\ :sup:`3`\
    SKCLIMATOLOGY_BRCL_CM3       BrCl                                                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_BRNO3_CM3      BrNO\ :sub:`3`\                                         Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_BRO_CM3        BrO                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_BRX_CM3        BrX                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_BRY_CM3        BrY                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_BR_CM3         Br                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_C2H2_CM3       C\ :sub:`2`\ H\ :sub:`2`\                               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_C2H4_CM3       C\ :sub:`2`\ H\ :sub:`4`\                               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_C2H6_CM3       C\ :sub:`2`\ H\ :sub:`6`\                               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_C3H6O_CM3      C\ :sub:`3`\ H\ :sub:`6`\ O  Acetone                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_C5H8_CM3       C\ :sub:`5`\ H\ :sub:`8`\    Isoprene                   Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CCL4_CM3       CCl\ :sub:`4`\                                          Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CF2CL2_CM3     CF\ :sub:`2`\ Cl\ :sub:`2`\                             Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CF4_CM3        CF\ :sub:`4`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CFCL3_CM3      CFCl\ :sub:`3`\                                         Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH2O_CM3       CH\ :sub:`2`\ O              Formaldehyde               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH3BR_CM3      CH\ :sub:`3`\ Br                                        Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH3CL_CM3      CH\ :sub:`3`\ Cl                                        Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH3CN_CM3      CH\ :sub:`3`\ CN                                        Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH3I_CM3       CH\ :sub:`3`\ I              Methyl iodide              Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH3OH_CM3      CH\ :sub:`3`\ OH                                        Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CH4_CM3        CH\ :sub:`4`\                Methane                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CL2O2_CM3      Cl\ :sub:`2`\ O\ :sub:`2`\                              Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CL2_CM3        Cl\ :sub:`2`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CLNO3_CM3      CLNO\ :sub:`3`\                                         Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CLONO2_CM3     ClONO\ :sub:`2`\                                        Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CLO_CM3        ClO                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CLY_CM3        ClY                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CL_CM3         Cl                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CO2_CM3        CO\ :sub:`2`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_COF2_CM3       COF\ :sub:`2`\                                          Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_CO_CM3         CO                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_H2CO_CM3       H\ :sub:`2`\ CO                                         Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_H2O2_CM3       H\ :sub:`2`\ O\ :sub:`2`\                               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_H2O_CM3        H\ :sub:`2`\ O                                          Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_H2S_CM3        H\ :sub:`2`\ S                                          Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_H2_CM3         H\ :sub:`2`\                                            Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HBR_CM3        HBr                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HCL_CM3        HCl                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HCN_CM3        HCN                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HCOOH_CM3      HCOOH                                                   Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HF_CM3         HF                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HI_CM3         HI                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HNO2_CM3       HNO\ :sub:`2`\               Nitrous Acid               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HNO3_CM3       HNO\ :sub:`3`\               Nitric Acid                Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HNO4_CM3       HNO\ :sub:`4`\                                          Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HO2_CM3        HO\ :sub:`2`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HOBR_CM3       HOBr                                                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_HOCL_CM3       HOCl                                                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_H_CM3          H                                                       Atoms per cm\ :sup:`3`\
    SKCLIMATOLOGY_He_CM3         He                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_MECL_CM3       MECl                                                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_N_CM3          N                                                       Atoms per cm\ :sup:`3`\
    SKCLIMATOLOGY_N2O5_CM3       N\ :sub:`2`\ O\ :sub:`5`\                               Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_N2O_CM3        N\ :sub:`2`\ O                                          Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_N2_CM3         N\ :sub:`2`\                                            Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_NH3_CM3        NH\ :sub:`3`\               Ammonia                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_NITS           Inorganic Nitrates
    SKCLIMATOLOGY_NO2_CM3        NO\ :sub:`2`                Nitrogen Dioxide            Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_NO3_CM3        NO\ :sub:`3`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_NOPLUS_CM3     NO\ :sup:`+`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_NOY_CM3        NOY                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_NO_CM3         NO                                                      Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_O2_CM3         O\ :sub:`2`                 Molecular oxygen            Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_O3_CM3         O\ :sub:`3`                 Ozone                       Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_OCLO_CM3       OClO                                                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_OCS_CM3        OCS                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_OH_CM3         OH                          Hydroxyl                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_O_CM3          O                           Atomic oxygen               Atoms per cm\ :sup:`3`\
    SKCLIMATOLOGY_PAN_CM3        PAN                         Peroxy acetyl nitrate       Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_PH3_CM3        PH\ :sub:`3`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_ROOH_CM3       ROOH                                                    Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_ROO_CM3        ROO                                                     Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_SF6_CM3        SF\ :sub:`6`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_SO2_CM3        SO\ :sub:`2`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_SO4_CM3        SO\ :sub:`4`\                                           Molecules per cm\ :sup:`3`\
    SKCLIMATOLOGY_XXX_CM3        XXX                                                     Molecules per cm\ :sup:`3`\
    ============================ ======================================================  ==============================




Molecular Volume Mixing Ratio GUIDS
-----------------------------------

    ==========================================      ==============================  ===============================================
    GUID String                                     Species                         Units
    ==========================================      ==============================  ===============================================
    SKCLIMATOLOGY_BRX_VMR                           BrX                             Volume Mixing Ratio
    SKCLIMATOLOGY_BRY_VMR                           BrY                             Volume Mixing Ratio
    SKCLIMATOLOGY_BRO_VMR                           BrO                             Volume Mixing Ratio
    SKCLIMATOLOGY_CO2_VMR                           CO\ :sub:`2`\                   Volume Mixing Ratio
    SKCLIMATOLOGY_C3H6O_VMR                         C\ :sub:`3`\ H\ :sub:`6`\ O     Volume Mixing Ratio
    SKCLIMATOLOGY_C5H8_VMR                          C\ :sub:`5`\ H\ :sub:`8`\       Volume Mixing Ratio
    SKCLIMATOLOGY_CCL4_VMR                          CCl\ :sub:`4`\                  Volume Mixing Ratio
    SKCLIMATOLOGY_CF2CL2_VMR                        CF\ :sub:`2`\ Cl\ :sub:`2`\     Volume Mixing Ratio
    SKCLIMATOLOGY_CFCL3_VMR                         CFCl\ :sub:`3`\                 Volume Mixing Ratio
    SKCLIMATOLOGY_CH2O_VMR                          CH\ :sub:`2`\ O                 Volume Mixing Ratio
    SKCLIMATOLOGY_CH3BR_VMR                         CH\ :sub:`3`\ Br                Volume Mixing Ratio
    SKCLIMATOLOGY_CH3CL_VMR                         CH\ :sub:`3`\ Cl                Volume Mixing Ratio
    SKCLIMATOLOGY_CH3I_VMR                          CH\ :sub:`3`\ I                 Volume Mixing Ratio
    SKCLIMATOLOGY_CH4_VMR                           CH\ :sub:`4`\                   Volume Mixing Ratio
    SKCLIMATOLOGY_CLY_VMR                           ClY                             Volume Mixing Ratio
    SKCLIMATOLOGY_CO_VMR                            CO                              Volume Mixing Ratio
    SKCLIMATOLOGY_CO2_VMR                           CO\ :sub:`2`\                   Volume Mixing Ratio
    SKCLIMATOLOGY_H2_VMR                            H\ :sub:`2`\                    Volume Mixing Ratio
    SKCLIMATOLOGY_H2O_VMR                           H\ :sub:`2`\ O                  Volume Mixing Ratio
    SKCLIMATOLOGY_HNO2_VMR                          HNO\ :sub:`2`\                  Volume Mixing Ratio
    SKCLIMATOLOGY_HNO3_VMR                          HNO\ :sub:`3`\                  Volume Mixing Ratio
    SKCLIMATOLOGY_MECL_VMR                          MECl                            Volume Mixing Ratio
    SKCLIMATOLOGY_N2_VMR                            N\ :sub:`2`\                    Volume Mixing Ratio
    SKCLIMATOLOGY_N2O_VMR                           N\ :sub:`2`\ O                  Volume Mixing Ratio
    SKCLIMATOLOGY_NO2_VMR                           NO\ :sub:`2`\                   Volume Mixing Ratio
    SKCLIMATOLOGY_NO_VMR                            NO                              Volume Mixing Ratio
    SKCLIMATOLOGY_NOY_VMR                           NOY                             Volume Mixing Ratio
    SKCLIMATOLOGY_NH3_VMR                           NH\ :sub:`3`                    Volume Mixing Ratio
    SKCLIMATOLOGY_O3_VMR                            O\ :sub:`3`\                    Volume Mixing Ratio
    SKCLIMATOLOGY_O2_VMR                            O\ :sub:`2`\                    Volume Mixing Ratio
    SKCLIMATOLOGY_PAN_VMR                           PAN                             Volume Mixing Ratio
    SKCLIMATOLOGY_SO2_VMR                           SO\ :sub:`2`\                   Volume Mixing Ratio
    SKCLIMATOLOGY_SO4_VMR                           SO\ :sub:`4`\                   Volume Mixing Ratio
    SKCLIMATOLOGY_XXX_VMR                           XXX                             Volume Mixing Ratio
    ==========================================      ==============================  ===============================================


Photochemical emission  GUIDS
-----------------------------

    ==========================================      ==============================  ===============================================
    GUID String                                     Species                         Units
    ==========================================      ==============================  ===============================================
    SKEMISSION_PHOTOCHEMICAL_0                      User defined species 0
    SKEMISSION_PHOTOCHEMICAL_1                      User defined species 1
    SKEMISSION_PHOTOCHEMICAL_2                      User defined species 2
    SKEMISSION_PHOTOCHEMICAL_3                      User defined species 3
    SKEMISSION_PHOTOCHEMICAL_4                      User defined species 4
    SKEMISSION_PHOTOCHEMICAL_5                      User defined species 5
    SKEMISSION_PHOTOCHEMICAL_6                      User defined species 6
    SKEMISSION_PHOTOCHEMICAL_7                      User defined species 7
    SKEMISSION_PHOTOCHEMICAL_8                      User defined species 8
    SKEMISSION_PHOTOCHEMICAL_9                      User defined species 9
    SKEMISSION_PHOTOCHEMICAL_O2                     O2
    SKEMISSION_PHOTOCHEMICAL_OH                     OH
    SKEMISSION_PHOTOCHEMICAL_O3                     O3
    SKEMISSION_THERMAL                              Thermal emission
    ==========================================      ==============================  ===============================================

Collisional Induced Absorption (CIA)
-------------------------------------

    ==========================================      ==============================  =======================================================
    GUID String                                     Species                         Description
    ==========================================      ==============================  =======================================================
    SKCLIMATOLOGY_O2_O2_CM6                         square of O2-O2 number density  Suitable for O2/O2 collisional induced absorption (CIA)
    ==========================================      ==============================  =======================================================



Miscellaneous GUIDS
-------------------

    ==========================================      ==============================
    GUID String                                     Description
    ==========================================      ==============================
    SKCLIMATOLOGY_UNDEFINED                         Undefined quantity
    SKCLIMATOLOGY_JH2O                              Used in Pratmo
    ==========================================      ==============================

