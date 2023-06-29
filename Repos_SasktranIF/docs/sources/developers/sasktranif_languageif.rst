.. _sasktranis_languageif:

*****************************************
SasktranIF, C++ Interface
*****************************************

The ``SasktranIF`` module provides a standard interface so it can be (ii) integrated into different high level programming
languages such as Python, Matlab, IDL etc and (ii) load sasktran objects from  multiple sasktran C++ extensions. The interface
is still a work in progress but has been used by the ARG group for several years withouit issue. Ideally it woulkd be fully ABI
compliant and not not pass class definitions across the shareable object/DLL boundary. In reality we pass simple C++ classes
across the boundary,  which are simple, almost purely virtual. They are similar to the COM interfaces from software programming
at the turn of the century. We have not experienced any problems on Windows or Linux.

Each extension module exposes a limited set of C/C++ functions to the SasktranIF module. Some of the functions allow SasktranIF
to create new sasktran objects, these are the *factory* objects. Other functions allow SasktranIF to configure and query the
extension during its loading and initialization. Finally there is a third component which is the global registry. Each C/C++
extension must upon first-time initialization create entries in the registry to indicate which sasktran objects it creates.
The registry entries are the only way that SasktranIF knows what objects are created by which extension.


Factory Functions
------------------
The following factory functions create sasktran objects for use by the high-level language. One function is provided
for each type of sasktran object. All functions are optional. If one of the functions is not implemented then the
extension does not support creation of that type of sasktran object.  Note that these functions must be exported by the
dynamic link library of shareable object.

======================================= ================================================================================
C++ Factory Functions                   Description
======================================= ================================================================================
:ref:`SKTRANIF_CreateEngine2`           Optional. Factory function for ISKEngine objects within C++ extension
:ref:`SKTRANIF_CreateClimatology2`      Optional. Factory function for ISKClimatology objects within C++ extension
:ref:`SKTRANIF_CreateOpticalProperty2`  Optional. Factory function for ISKOpticalProperty objects within C++ extension
:ref:`SKTRANIF_CreateEmission2`         Optional. Factory function for ISKEmission objects within C++ extension
:ref:`SKTRANIF_CreateBRDF2`             Optional. Factory function for ISKBrdf objects within C++ extension
:ref:`SKTRANIF_CreateSolarSpectrum2`    Optional. Factory function for ISKSolarSpectrum objects within C++ extension
:ref:`SKTRANIF_CreateGeodetic2`         Optional. Factory function for ISKGeodetic objects within C++ extension
======================================= ================================================================================

Logistical Functions
--------------------
The following functions provide the logistics required to integrate a given sasktran extension into the sasktran framework.

=================================================   =============================================================================
Logistical Functions                                Description
=================================================   =============================================================================
:ref:`SKTRAN_IFSetRegistryDirectoryInChildDLL`      Notify extension of folder location of the global YAML based registry
:ref:`SKTRAN_IFCreateRegistryEntriesForChildDLL`    Request extension to write entries into the global registry
:ref:`SKTRAN_IFInitializeLogger`                    Notify extension of the object used to log messages back to high level language
:ref:`SKTRAN_IFGlobalHandleTable`                   Optional. Request extension to provide tables of climatology handles created by this extension
=================================================   =============================================================================

Functions
---------

..  _SKTRANIF_CreateEngine2:

SKTRANIF_CreateEngine2
^^^^^^^^^^^^^^^^^^^^^^

Factory function to create the requested sasktran engine object. A Sasktran C++ module
only implements and exports this function if it provides radiative transfer engines. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateEngine2( const char* userenginename, ISKEngine_Stub**  engine)

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateEngine2( const char* userenginename, ISKEngine_Stub**  engine)
    {
        nxString    enginename(userenginename);

        enginename.MakeUpper();
        if      (enginename == "SO")  *engine = new ISKEngine_Stub_SO;
        else if (enginename == "HR")  *engine = new ISKEngine_Stub_HR;
        else if (enginename == "MC")  *engine = new ISKEngine_Stub_MC;
        else if (enginename == "OCC") *engine = new ISKEngine_Stub_OCC;
        else
        {
            *engine = nullptr;
            nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateEngine, engine [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) userenginename);

        }
        return (*engine != nullptr);
    }

The C++ extension must also create appropriate registry entries within `software/usask-arg/sasktranif/engines` during its
one-time initialization in function :ref:`SKTRAN_IFCreateRegistryEntriesForChildDLL`, for example::

     software:
      usask-arg:
        sasktranif:
         engines:
            hr:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            mc:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            so:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll


..  _SKTRANIF_CreateClimatology2:

SKTRANIF_CreateClimatology2
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Requests that a binary module create a new instance of a climatology object for the user.  A Sasktran C++ module
implements and exports this function only if it provides climatologies. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateClimatology2( const char* userclimatename , ISKClimatology_Stub**  climatology )

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateClimatology2( const char* userclimatename , ISKClimatology_Stub**  climatology )
    {
        bool        ok = false;
        nxString    name(userclimatename);

        name.MakeUpper();
        if      (name == "MSIS90"       )                       *climatology = new ISKClimatology_Stub_MSIS                             ( new skClimatology_MSIS90);
        else if (name == "ECMWF"        )                       *climatology = new ISKClimatology_Stub_Base                             ( new skClimatology_Ecmwf );
        else if (name == "O3LABOW"      )                       *climatology = new ISKClimatology_Stub_Base                             ( new skClimatology_LabowOzoneVMR );
        else if (name == "NO2PRATMO"    )                       *climatology = new ISKClimatology_Stub_Base                             ( new skClimatology_Pratmo );
        else if (name == "OSIRISL2_O3RTMODEL_V507" )            *climatology = new ISKClimatology_Stub_Base                             ( new skClimatology_OsirisL2_O3RTModel_V507 );
        else if (name == "OSIRISL2_NO2RTMODEL_V507")            *climatology = new ISKClimatology_Stub_Base                             ( new skClimatology_OsirisL2_NO2RTModel_V507 );
        else if (name == "OSIRISL2_AEROSOLRTMODEL_V507")        *climatology = new ISKClimatology_Stub_Base                             ( new skClimatology_OsirisL2_AerosolRTModel_V507 );
        else if (name == "OSIRISL2_AEROSOLMODERADIUS_V600")     *climatology = new ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600  ( new skClimatology_OsirisAerosolModeRadiusV600 );
        else if (name == "USERDEFINED_PROFILE" )                *climatology = new ISKClimatology_Stub_UserDefined                      ( new skClimatology_UserTableSpline );
        else if (name == "USERDEFINED_PROFILE3D_LATLONHEIGHT" ) *climatology = new ISKClimatology_Stub_UserDefined3D                    ( new skClimatology_UserDefined3D_LatLonHeight );
        else if (name == "USERDEFINED_PROFILE_TABLE" )          *climatology = new ISKClimatology_Stub_UserDefinedTable                 ( new skClimatology_UserDefinedTable );
        else if (name == "USERDEFINED_PROFILE_PLANE" )          *climatology = new ISKClimatology_Stub_UserDefinedPlane                 ( new skClimatology_UserDefinedPlane );
        else if (name == "ONE_PRESSURE_TEMP" )                  *climatology = new ISKClimatology_Stub_OnePressureTemp                  ( new skClimatology_OneTemperatureAndPressure );
        else if (name == "CONSTANTVALUE" )                      *climatology = new ISKClimatology_Stub_Constant                         ( new skClimatology_Constant );
        else if (name == "LINEARCOMBO" )                        *climatology = new ISKClimatology_Stub_LinearCombination                ( new skClimatologyLinearCombination );
        else
        {
            *climatology = nullptr;
             nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateClimatology, climatology [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) userclimatename);
        }
        return (*climatology != nullptr);
    }

The C++ extension must also create appropriate registry entries within `software/usask-arg/sasktranif/climatology` during its
one-time initialization in function :ref:`SKTRAN_IFCreateRegistryEntriesForChildDLL`, for example::

    software:
      usask-arg:
        sasktranif:
          climatology:
            constantvalue:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            ecmwf:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            linearcombo:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            msis90:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            no2pratmo:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3labow:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            one_pressure_temp:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            osirisl2_aerosolmoderadius_v600:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            osirisl2_aerosolrtmodel_v507:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            osirisl2_no2rtmodel_v507:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            osirisl2_o3rtmodel_v507:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_profile:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_profile3d_latlonheight:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_profile_plane:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_profile_table:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll


..  _SKTRANIF_CreateOpticalProperty2:

SKTRANIF_CreateOpticalProperty2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Requests that a binary module create a new instance of an optical property object for the user.  A Sasktran C++ module
implements and exports this function only if it provides optical properties. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateOpticalProperty2( const char* useroptpropname, ISKOpticalProperty_Stub**  optprop )

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateOpticalProperty2( const char* useroptpropname, ISKOpticalProperty_Stub**  optprop )
    {
        nxString    name(useroptpropname);
        bool        ok = true;

        name.MakeUpper();
        if      (name == "RAYLEIGH")             *optprop = new ISKOpticalProperty_Stub_Base        ( new skOpticalProperties_RayleighDryAir);
        else if (name == "SIMPLERAYLEIGH")       *optprop = new ISKOpticalProperty_Stub_Base        ( new skOpticalProperties_SimpleRayleigh);
        else if (name == "NO2_BURROWS")          *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_NO2_Burrows98);
        else if (name == "NO2_OSIRISRES")        *optprop = new ISKOpticalProperty_Stub_Base        ( new skOpticalProperties_NO2_OSIRISRes);
        else if (name == "NO2_VANDAELE1998")     *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_NO2_Vandaele1998);
        else if (name == "SO2_VANDAELE2009")     *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_SO2_Vandaele2009);
        else if (name == "SO2_FREEMAN1984")      *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_SO2_Freeman1984);
        else if (name == "SO2_RUFUS2003")        *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_SO2_Rufus2003);
        else if (name == "SO2_BOGUMIL2003")      *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_SO2_Bogumil2003);
        else if (name == "O3_BASSPAURLINEAR")    *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_BassPaur);
        else if (name == "O3_BASSPAURQUADRATIC") *optprop = new ISKOpticalProperty_Stub_Base        ( new skOpticalProperties_O3_BassPaurQuadratic);
        else if (name == "O3_DBM")               *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_DaumontBrionMalicet);
        else if (name == "O3_VOIGT")             *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_FTSVoigt);
        else if (name == "O3_GOMEBURROWS")       *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_GomeBurrows);
        else if (name == "O3_OSIRISRES")         *optprop = new ISKOpticalProperty_Stub_Base        ( new skOpticalProperties_O3_OSIRISRes);
        else if (name == "O3_SCIABOGUMILV3")     *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_SciaBogumilV3);
        else if (name == "O3_SCIABOGUMILV4")     *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_SciaBogumilV4);
        else if (name == "O3_SERDYUCHENKOV1")    *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O3_SerdyuchenkoV1);
        else if (name == "O2_O2_FALLY2000")      *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O4_Fally2000);
        else if (name == "O2_O2_THALMAN2013")    *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_O4_Thalman2013);
        else if (name == "O2_O2_HITRAN2016" )    *optprop = new ISKOpticalProperty_Stub_Base        ( new skOpticalProperties_O4_Hitran2016);
        else if (name == "MIEAEROSOL_H2SO4")     *optprop = new ISKOpticalProperty_Stub_Aerosol     ( new skOpticalProperties_AerosolProfileH2SO4 );
        else if (name == "MIEAEROSOL_DUST")      *optprop = new ISKOpticalProperty_Stub_Aerosol     ( new skOpticalProperties_AerosolProfileDust );
        else if (name == "MIEAEROSOL_WATER")     *optprop = new ISKOpticalProperty_Stub_Aerosol     ( new skOpticalProperties_AerosolProfileWater);
        else if (name == "MIEAEROSOL_ICE")       *optprop = new ISKOpticalProperty_Stub_Aerosol     ( new skOpticalProperties_AerosolProfileIce_Mie);
        else if (name == "TMATRIXAEROSOL_ICE")   *optprop = new ISKOpticalProperty_Stub_Aerosol     ( new skOpticalProperties_AerosolProfileIce );
        else if (name == "BAUM_ICECRYSTALS")     *optprop = new ISKOpticalProperty_Stub_Baum        ( new skOpticalProperties_BaumIceCrystals2014 );
        else if (name == "USERDEFINED_TABLES")   *optprop = new ISKOpticalProperty_Stub_UserDefined ( new skOpticalProperties_UserDefinedAbsorption );
        else if (name == "CONVOLVED_CROSSSECTION")*optprop = new ISKOpticalProperty_Stub_ConvolvedFixedFWHM ( new skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM );
        else if (name == "MART_HYBRIDPROFILE")   *optprop  = new ISKOpticalProperty_Stub_MartHybridProfile( new skOpticalProperties_MartHybridProfile );
        else if (name.Find("HITRANCHEMICAL_") == 0)
        {
            nxString    chemicalname;
            chemicalname = name.Right( name.GetLength()-15 );
            *optprop = new ISKOpticalProperty_Stub_Hitran( new skOpticalProperties_HitranChemical(chemicalname, 0, 500000.0 ) );
        }
        else
        {
            ok = false;
            *optprop = nullptr;
            nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateOpticalProperty, opticalproperty [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) useroptpropname);

        }
        ok = ok && nullptr!=*optprop;
        return ok;
    }

The C++ extension must also create appropriate registry entries within `software/usask-arg/sasktranif/opticalproperty` during its
one-time initialization in function :ref:`SKTRAN_IFCreateRegistryEntriesForChildDLL`, for example::

    software:
      usask-arg:
        sasktranif:
          opticalproperty:
            baum_icecrystals:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            convolved_crosssection:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_c2h2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_c2h4:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_c2h6:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_c2n2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_c4h2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_cf4:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ch3br:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ch3cl:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ch3cn:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ch3oh:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ch4:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_clo:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_clono2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_co:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_co2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_cocl2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_cof2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_cs:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_h2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_h2co:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_h2o:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_h2o2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_h2s:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hbr:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hc3n:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hcl:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hcn:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hcooh:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hf:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hi:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hno3:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ho2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hobr:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_hocl:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_n2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_n2o:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_nh3:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_no:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_no2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_noplus:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_o:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_o2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_o3:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ocs:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_oh:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_ph3:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_sf6:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_so2:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hitranchemical_so3:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            mieaerosol_dust:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            mieaerosol_h2so4:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            mieaerosol_ice:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            mieaerosol_water:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            no2_burrows:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            no2_osirisres:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            no2_vandaele1998:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o2_o2_fally2000:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o2_o2_hitran2016:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o2_o2_thalman2013:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_basspaurlinear:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_basspaurquadratic:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_dbm:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_gomeburrows:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_osirisres:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_sciabogumilv3:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_sciabogumilv4:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_serdyuchenkov1:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            o3_voigt:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            osiris_convolveddiscretewavelen:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            rayleigh:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            simplerayleigh:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            so2_bogumil2003:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            so2_freeman1984:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            so2_rufus2003:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            so2_vandaele2009:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            tmatrixaerosol_ice:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_tables:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll


..  _SKTRANIF_CreateEmission2:

SKTRANIF_CreateEmission2
^^^^^^^^^^^^^^^^^^^^^^^^

Requests that a binary module create a new instance of an emission object for the user.  A Sasktran C++ module
implements and exports this function only if it provides emission properties. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateEmission2( const char* useroptpropname, ISKEmission_Stub**  emission )

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateEmission2( const char* useroptpropname, ISKEmission_Stub**  emission )
    {
        nxString    name(useroptpropname);
        bool        ok = false;

        name.MakeUpper();
        if      (name == "USERDEFINED_WAVELENGTHHEIGHT") *emission = new ISKEmission_Stub_Tabulated_HeightWavelength( new skEmission_Tabulated_HeightWavelength);
        else if (name == "THERMAL")                      *emission = new ISKEmission_Stub_Thermal                   ( new skEmission_Thermal);
        else if (name == "HITRAN_PHOTOCHEMICAL")         *emission = new ISKEmission_Stub_HitranChemical            ( new skEmission_HitranChemical);
        else
        {
            *emission = nullptr;
            nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateEmission, emission [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) useroptpropname);

        }
        return (*emission != nullptr);
    }

The C++ extension must also create appropriate registry entries within `software/usask-arg/sasktranif/emission` during its
one-time initialization in function :ref:`SKTRAN_IFCreateRegistryEntriesForChildDLL`, for example::

    software:
      usask-arg:
        sasktranif:
          emission:
            hitran_photochemical:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            thermal:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_wavelengthheight:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll

..  _SKTRANIF_CreateBRDF2:

SKTRANIF_CreateBRDF2
^^^^^^^^^^^^^^^^^^^^

Requests that a binary module create a new instance of an BRDF object for the user.  A Sasktran C++ module
implements and exports this function only if it provides BRDF object. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateBRDF2( const char* brdfname, ISKBrdf_Stub**  brdf )

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateBRDF2( const char* brdfname, ISKBrdf_Stub**  brdf )
    {
        nxString    name(brdfname);
        bool        ok = false;

        *brdf = nullptr;
        name.MakeUpper();
        if      (name == "LAMBERTIAN")                  *brdf = new ISKBrdf_Stub_LambertianAlbedo                   ( new SKTRAN_BRDF_Lambertian );
        else if (name == "SNOW_KOKHANOVSKY2012")        *brdf = new ISKBrdf_Stub_Snow_Kokhanovsky2012               ( new SKTRAN_BRDF_Snow_Kokhanovsky2012);
        else if (name == "ROUJEAN")                     *brdf = new ISKBrdf_Stub_Roujean                            ( new SKTRAN_BRDF_Roujean);
        else if (name == "ROUJEAN_KERNEL")              *brdf = new ISKBrdf_Stub_Roujean_Kernel                     ( new SKTRAN_BRDF_Roujean_Kernel);
        else if (name == "LI_SPARSE_KERNEL")            *brdf = new ISKBrdf_Stub_Li_Kernel                          ( new SKTRAN_BRDF_LiSparse_Kernel);
        else if (name == "LI_DENSE_KERNEL")             *brdf = new ISKBrdf_Stub_Li_Kernel                          ( new SKTRAN_BRDF_LiDense_Kernel);
        else if (name == "LI_SPARSE_RECIPROCAL_KERNEL") *brdf = new ISKBrdf_Stub_Li_Kernel                          ( new SKTRAN_BRDF_LiSparseReciprocal_Kernel);
        else if (name == "ROSS_THIN_KERNEL")            *brdf = new ISKBrdf_Stub_Ross_Kernel                        ( new SKTRAN_BRDF_RossThin_Kernel);
        else if (name == "ROSS_THICK_KERNEL")           *brdf = new ISKBrdf_Stub_Ross_Kernel                        ( new SKTRAN_BRDF_RossThick_Kernel);
        else if (name == "COX_MUNK")                    *brdf = new ISKBrdf_Stub_Cox_Munk                           ( new SKTRAN_BRDF_CoxMunk);
        else if (name == "RAHMAN")                      *brdf = new ISKBrdf_Stub_Rahman                             ( new SKTRAN_BRDF_Rahman);
        else if (name == "HAPKE")                       *brdf = new ISKBrdf_Stub_Hapke                              ( new SKTRAN_BRDF_Hapke);
        else if (name == "LINEAR_COMBINATION")          *brdf = new ISKBrdf_Stub_LinearCombination                  ( new SKTRAN_BRDF_LinearCombination);
        else if (name == "MODIS")                       *brdf = new ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal  ( new SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal);
        else if (name == "USERDEFINED_LATLON")          *brdf = new ISKBrdf_Stub_UserDefinedLatLon                  ( new SKTRAN_BRDF_UserDefinedLatLon );
        else
        {
            *brdf = nullptr;
            nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateBRDF2, BRDF [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) name);
        }
        return (*brdf != nullptr);
    }

The C++ extension must also create appropriate registry entries within `software/usask-arg/sasktranif/brdf` during its
one-time initialization in function :ref:`SKTRAN_IFCreateRegistryEntriesForChildDLL`, for example::

    software:
      usask-arg:
        sasktranif:
          brdf:
            cox_munk:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            hapke:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            lambertian:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            li_dense_kernel:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            li_sparse_kernel:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            li_sparse_reciprocal_kernel:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            linear_combination:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            modis:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            rahman:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            ross_thick_kernel:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            ross_thin_kernel:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            roujean:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            roujean_kernel:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            snow_kokhanovsky2012:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll
            userdefined_latlon:
              dllname: C:/Users/nickl/anaconda3/envs/sasktran/lib/site-packages/sasktran_core/_sasktran_core_internals.dll


..  _SKTRANIF_CreateSolarSpectrum2:

SKTRANIF_CreateSolarSpectrum2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Requests that a binary module create a new instance of an SolarSpectrum object for the user.  A Sasktran C++ module
implements and exports this function only if it provides SolarSpectrum objects. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateSolarSpectrum2( const char* solarname, ISKSolarSpectrum_Stub**  solar )

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateSolarSpectrum2( const char* solarname, ISKSolarSpectrum_Stub**  solar )
    {
        nxString    name(solarname);
        bool        ok = false;

        *solar = nullptr;
        name.MakeUpper();
        if      (name == "SAO2010")                 *solar = new ISKSolarSpectrum_Stub_Base( new skSolarSpectrum_SAO2010);
        else if (name == "FONTELA_UVIS_3MICRON")    *solar = new ISKSolarSpectrum_Stub_Base( new skSolarSpectrum_FontelaUVIS3Micron);
        else if (name == "FONTELA_UVIS_100MICRON")  *solar = new ISKSolarSpectrum_Stub_Base( new skSolarSpectrum_FontelaUVIS100Micron);
        else
        {
            *solar = nullptr;
            nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateSolarSpectrum, solarspectrum [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) name);
        }
        return (*solar != nullptr);
    }



..  _SKTRANIF_CreateGeodetic2:

SKTRANIF_CreateGeodetic2
^^^^^^^^^^^^^^^^^^^^^^^^

Requests that a binary module create a new instance of an geodetic object for the user.  A Sasktran C++ module
implements and exports this function only if it provides a geodetic object. Note that teh core components distribute
working version of a geodetic object. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateGeodetic2( const char* geoidname, ISKGeodetic_Stub**  geoid)

An example piece of code taken from the sasktran_core implementation::

    extern "C" DLL_PUBLIC bool SKTRANIF_CreateGeodetic2( const char* geoidname, ISKGeodetic_Stub**  geoid)
    {
        nxString    enginename(geoidname);
        bool        ok = false;

        enginename.MakeUpper();
        if (enginename == "STANDARD")
        {
            *geoid = new ISKGeodetic_Stub_std;
            ok = (geoid != nullptr);
        }
        else
        {
            *geoid = nullptr;
            nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateGeodetic, geoid [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) enginename);

        }
        return ok;
    }



..  _SKTRAN_IFSetRegistryDirectoryInChildDLL:

SKTRAN_IFSetRegistryDirectoryInChildDLL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This function must be implemented by the C++ binary module. It informs the C++ module of the location of a YAML based registry scheme.
This scheme works well with the ``nxRegistry`` classes used in the sasktran core components.
We strongly recommend that new C++ modules manage their registry settings with the same ``nxRegistry`` classes. The format of
the function is::

    extern "C" DLL_PUBLIC bool SKTRAN_IFSetRegistryDirectoryInChildDLL( const char* registrydirname)

..  _SKTRAN_IFCreateRegistryEntriesForChildDLL:

SKTRAN_IFCreateRegistryEntriesForChildDLL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This function must be implemented by the C++ binary module. It instructs the module to register its sasktran objects in
the current registry. The high level languages configure the SasktranIF with registry information and SasktranIF passes
this information onto the modules. Its a bit complicated but it works. We strongly recommend that new C++ modules manage their
registry settings with the same ``nxRegistry`` classes. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRAN_IFCreateRegistryEntriesForChildDLL( const char* paramstr )

..  _SKTRAN_IFInitializeLogger:

SKTRAN_IFInitializeLogger
^^^^^^^^^^^^^^^^^^^^^^^^^

This function should  be implemented by the C++ binary module. It allows the  module to send its logging information,
via ``nxLog::Record``, to the high level language logging system. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRAN_IFInitializeLogger( InxLog* logger)

..  _SKTRAN_IFGlobalHandleTable:

SKTRAN_IFGlobalHandleTable
^^^^^^^^^^^^^^^^^^^^^^^^^^

This function is optional for the C++ extension. It allows the  extension to notify SasktranIF of :ref:`climatology handles <climatologyhandles>`
created in this extension. Python users cannot use :ref:`climatology handles <climatologyhandles>` defined in extensions unless they are registered
with SasktranIF. The format of the function is::

    extern "C" DLL_PUBLIC bool SKTRAN_IFGlobalHandleTable( GlobalClimatologyHandleTable** entry, int* numpoints )
