// SasktranImpl.cpp : Defines the exported functions for the DLL application.
//

#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"


/*-----------------------------------------------------------------------------
 *					SetModule_SandboxRegistryDirectory		 2016- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SetModule_SandboxRegistryDirectory( const char* registrydir)
{
	nxString	regdir(registrydir);

	regdir.MakeDirectorySeparatorsOSConsistent();
	nxDirectory::CreateADirectory( regdir);
	nxRegistryKey::RegistryLocationVar()->Set_BaseDirectory(regdir);
//	printf("SetModule_SandboxRegistryDirectory, The directory registry is <%s>", (const char*)registrydir);
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateClimatology		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateClimatology2( const char* userclimatename , ISKClimatology_Stub**  climatology )
{
	bool		ok = false;
	nxString	name(userclimatename);

	name.MakeUpper();
	if      (name == "MSIS90"       )						*climatology = new ISKClimatology_Stub_MSIS								( new skClimatology_MSIS90);
	else if (name == "O3LABOW"      )						*climatology = new ISKClimatology_Stub_Base								( new skClimatology_LabowOzoneVMR );
	else if (name == "NO2PRATMO"    )						*climatology = new ISKClimatology_Stub_Base								( new skClimatology_Pratmo );
	else if (name == "USERDEFINED_PROFILE" )				*climatology = new ISKClimatology_Stub_UserDefined						( new skClimatology_UserTableSpline );
	else if (name == "USERDEFINED_PROFILE3D_LATLONHEIGHT" )	*climatology = new ISKClimatology_Stub_UserDefined3D					( new skClimatology_UserDefined3D_LatLonHeight );
	else if (name == "USERDEFINED_PROFILE_TABLE" )			*climatology = new ISKClimatology_Stub_UserDefinedTable					( new skClimatology_UserDefinedTable );
	else if (name == "USERDEFINED_PROFILE_PLANE" )			*climatology = new ISKClimatology_Stub_UserDefinedPlane					( new skClimatology_UserDefinedPlane );
	else if (name == "ONE_PRESSURE_TEMP" )					*climatology = new ISKClimatology_Stub_OnePressureTemp					( new skClimatology_OneTemperatureAndPressure );
	else if (name == "CONSTANTVALUE" )						*climatology = new ISKClimatology_Stub_Constant							( new skClimatology_Constant );
	else if (name == "LINEARCOMBO" )						*climatology = new ISKClimatology_Stub_LinearCombination				( new skClimatologyLinearCombination );
	else
	{
		*climatology = nullptr;
		 nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateClimatology, climatology [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) userclimatename);
	}
	return (*climatology != nullptr);
}

/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateOpticalProperty		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateOpticalProperty2( const char* useroptpropname, ISKOpticalProperty_Stub**  optprop )
{
	nxString	name(useroptpropname);
	bool		ok = true;

	name.MakeUpper();
	if (name == "RAYLEIGH")             *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_RayleighDryAir);
	else if (name == "SIMPLERAYLEIGH")       *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_SimpleRayleigh);
	else if (name == "INELASTICRAYLEIGH")	 *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_RayleighDryAir_Inelastic);
	else if (name == "NO2_BURROWS")          *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_NO2_Burrows98);
	else if (name == "NO2_OSIRISRES")        *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_NO2_OSIRISRes);
	else if (name == "NO2_VANDAELE1998")     *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_NO2_Vandaele1998);
	else if (name == "SO2_VANDAELE2009")     *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_SO2_Vandaele2009);
	else if (name == "SO2_FREEMAN1984")      *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_SO2_Freeman1984);
	else if (name == "SO2_RUFUS2003")        *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_SO2_Rufus2003);
	else if (name == "SO2_BOGUMIL2003")      *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_SO2_Bogumil2003);
	else if (name == "O3_BASSPAURLINEAR")    *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_BassPaur);
	else if (name == "O3_BASSPAURQUADRATIC") *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_O3_BassPaurQuadratic);
	else if (name == "O3_DBM")		         *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_DaumontBrionMalicet);
	else if (name == "O3_VOIGT")		     *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_FTSVoigt);
	else if (name == "O3_GOMEBURROWS")       *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_GomeBurrows);
	else if (name == "O3_OSIRISRES")         *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_O3_OSIRISRes);
	else if (name == "O3_SCIABOGUMILV3")     *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_SciaBogumilV3);
	else if (name == "O3_SCIABOGUMILV4")     *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_SciaBogumilV4);
	else if (name == "O3_SERDYUCHENKOV1")    *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O3_SerdyuchenkoV1);
	else if (name == "O2_O2_FALLY2000")      *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O4_Fally2000);
	else if (name == "O2_O2_THALMAN2013")    *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_O4_Thalman2013);
	else if (name == "O2_O2_HITRAN2016")    *optprop = new ISKOpticalProperty_Stub_Base(new skOpticalProperties_O4_Hitran2016);
	else if (name == "MIEAEROSOL_H2SO4")     *optprop = new ISKOpticalProperty_Stub_Aerosol(new skOpticalProperties_AerosolProfileH2SO4);
	else if (name == "MIEAEROSOL_DUST")      *optprop = new ISKOpticalProperty_Stub_Aerosol(new skOpticalProperties_AerosolProfileDust);
	else if (name == "MIEAEROSOL_WATER")	 *optprop = new ISKOpticalProperty_Stub_Aerosol(new skOpticalProperties_AerosolProfileWater);
	else if (name == "MIEAEROSOL_ICE")       *optprop = new ISKOpticalProperty_Stub_Aerosol(new skOpticalProperties_AerosolProfileIce_Mie);
	else if (name == "TMATRIXAEROSOL_ICE")   *optprop = new ISKOpticalProperty_Stub_Aerosol(new skOpticalProperties_AerosolProfileIce);
	else if (name == "BAUM_ICECRYSTALS")     *optprop = new ISKOpticalProperty_Stub_Baum(new skOpticalProperties_BaumIceCrystals2014);
	else if (name == "USERDEFINED_TABLES")   *optprop = new ISKOpticalProperty_Stub_UserDefined(new skOpticalProperties_UserDefinedAbsorption);
	else if (name == "USERDEFINED_PRESSURE") *optprop = new ISKOpticalProperty_Stub_UserDefinedPressure(new skOpticalProperties_UserDefinedAbsorptionPressure);
	else if (name == "USERDEFINED_SCATTERCONSTANTHEIGHT") *optprop = new ISKOpticalProperty_Stub_UserDefinedScatterConstantHeight(new skOpticalProperties_UserDefinedScatterConstantHeight);
	else if (name == "CONVOLVED_CROSSSECTION")*optprop = new ISKOpticalProperty_Stub_ConvolvedFixedFWHM ( new skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM );
	else if (name.Find("HITRANCHEMICAL_") == 0)
	{
		nxString	chemicalname;
		chemicalname = name.Right( name.GetLength()-15 );
		*optprop = new ISKOpticalProperty_Stub_Hitran( new skOpticalProperties_HitranChemical(chemicalname	) );
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

/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateOpticalProperty		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateEmission2( const char* useroptpropname, ISKEmission_Stub**  emission )
{
	nxString	name(useroptpropname);
	bool		ok = false;

	name.MakeUpper();
	if      (name == "USERDEFINED_WAVELENGTHHEIGHT") *emission = new ISKEmission_Stub_Tabulated_HeightWavelength( new skEmission_Tabulated_HeightWavelength);
	else if (name == "THERMAL")						 *emission = new ISKEmission_Stub_Thermal					( new skEmission_Thermal);
	else if (name == "HITRAN_PHOTOCHEMICAL")		 *emission = new ISKEmission_Stub_HitranChemical			( new skEmission_HitranChemical);
	//else if (name == "XXXXXX")                     *emission = new ISKOpticalProperty_Stub_UserDefined	    ( new skOpticalProperties_NO2_Burrows98);
	else
	{
		*emission = nullptr;
		nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateEmission, emission [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) useroptpropname);

	}
	return (*emission != nullptr);
}


/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateSolarSpectrum		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateSolarSpectrum2( const char* solarname, ISKSolarSpectrum_Stub**  solar )
{
	nxString	name(solarname);
	bool		ok = false;

	*solar = nullptr;
	name.MakeUpper();
	if      (name == "SAO2010")					*solar = new ISKSolarSpectrum_Stub_Base( new skSolarSpectrum_SAO2010);
	else if (name == "FONTELA_UVIS_3MICRON")	*solar = new ISKSolarSpectrum_Stub_Base( new skSolarSpectrum_FontelaUVIS3Micron);
	else if (name == "FONTELA_UVIS_100MICRON")	*solar = new ISKSolarSpectrum_Stub_Base( new skSolarSpectrum_FontelaUVIS100Micron);
	else
	{
		*solar = nullptr;
		nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateSolarSpectrum, solarspectrum [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) name);
	}
	return (*solar != nullptr);
}


/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateSolarSpectrum		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateBRDF2( const char* brdfname, ISKBrdf_Stub**  brdf )
{
	nxString	name(brdfname);
	bool		ok = false;

	*brdf = nullptr;
	name.MakeUpper();
	if      (name == "LAMBERTIAN")					*brdf = new ISKBrdf_Stub_LambertianAlbedo					( new SKTRAN_BRDF_Lambertian );
	else if (name == "SNOW_KOKHANOVSKY2012")		*brdf = new ISKBrdf_Stub_Snow_Kokhanovsky2012				( new SKTRAN_BRDF_Snow_Kokhanovsky2012);
	else if (name == "ROUJEAN")						*brdf = new ISKBrdf_Stub_Roujean							( new SKTRAN_BRDF_Roujean);
	else if (name == "ROUJEAN_KERNEL")				*brdf = new ISKBrdf_Stub_Roujean_Kernel						( new SKTRAN_BRDF_Roujean_Kernel);
	else if (name == "LI_SPARSE_KERNEL")			*brdf = new ISKBrdf_Stub_Li_Kernel							( new SKTRAN_BRDF_LiSparse_Kernel);
	else if (name == "LI_DENSE_KERNEL")				*brdf = new ISKBrdf_Stub_Li_Kernel							( new SKTRAN_BRDF_LiDense_Kernel);
	else if (name == "LI_SPARSE_RECIPROCAL_KERNEL") *brdf = new ISKBrdf_Stub_Li_Kernel							( new SKTRAN_BRDF_LiSparseReciprocal_Kernel);
	else if (name == "ROSS_THIN_KERNEL")			*brdf = new ISKBrdf_Stub_Ross_Kernel						( new SKTRAN_BRDF_RossThin_Kernel);
	else if (name == "ROSS_THICK_KERNEL")			*brdf = new ISKBrdf_Stub_Ross_Kernel						( new SKTRAN_BRDF_RossThick_Kernel);
	else if (name == "COX_MUNK")					*brdf = new ISKBrdf_Stub_Cox_Munk							( new SKTRAN_BRDF_CoxMunk);
	else if (name == "RAHMAN")						*brdf = new ISKBrdf_Stub_Rahman								( new SKTRAN_BRDF_Rahman);
	else if (name == "HAPKE")						*brdf = new ISKBrdf_Stub_Hapke								( new SKTRAN_BRDF_Hapke);
	else if (name == "LINEAR_COMBINATION")			*brdf = new ISKBrdf_Stub_LinearCombination					( new SKTRAN_BRDF_LinearCombination);
	else if (name == "MODIS")						*brdf = new ISKBrdf_Stub_MODIS_RossThickLiSparseReciprocal	( new SKTRAN_BRDF_MODIS_RossThickLiSparseReciprocal);
	else if (name == "USERDEFINED_LATLON")          *brdf = new ISKBrdf_Stub_UserDefinedLatLon                  ( new SKTRAN_BRDF_UserDefinedLatLon );
	else if (name == "SPECTRAL_VARYING")            *brdf = new ISKBrdf_Stub_SpectralVarying                    ( new SKTRAN_BRDF_SpectralVarying );
	else
	{
		*brdf = nullptr;
		nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateBRDF2, BRDF [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) name);
	}
	return (*brdf != nullptr);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_IFSetRegistryDirectory		 2016- 11- 14*/
/** This will force the DLL to use a text based YAML file for the registry 
 *	if registrydirname is not NULL.
**/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRAN_IFSetRegistryDirectoryInChildDLL( const char* registrydirname)
{
	bool	ok = true;
	
	if (registrydirname == nullptr) 
	{
		ok = nxRegistryKey::RegistryLocation().UseNativeRegistry();
		if (!ok)
		{
			printf("SasktranV3 Internal Registry Initialization::SKTRAN_IFSetRegistryDirectoryInChildDLL, the caller has requested using the native registry but that is not available on this build");
		}
	}
	else
	{

		nxString   regdir(registrydirname);
		ok =       regdir.MakeDirectorySeparatorsOSConsistent();
		ok = nxRegistryKey::RegistryLocationVar()->Set_BaseDirectory(registrydirname);
		if (!ok)
		{
			printf("SasktranV3 Internal Registry Initialization::SKTRAN_IFSetRegistryDirectoryInChildDLL, there were errors setting the registry to create and use base directory <%s>", (const char*)registrydirname);
		}
	}
	return ok;
}


extern bool SKTRAN_IFCreateRegistryEntriesForChildDLLImplementation( const char* paramstr, const char* modulename );

/*-----------------------------------------------------------------------------
 *					SKTRAN_IFCreateRegistryEntriesForChildDLL		 2016- 11- 14*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRAN_IFCreateRegistryEntriesForChildDLL( const char* paramstr )
{
	bool ok;

	ok = SKTRAN_IFCreateRegistryEntriesForChildDLLImplementation( paramstr,"sasktran_core" );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_IFInitializeLogger		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRAN_IFInitializeLogger( InxLog* logger)
{
	if  (logger != nullptr) nxLog::SetAsDefaultLogger(logger);
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_IFGlobalHandleTable		 2015- 11- 19*/
/** This is a hook where DLL's can add extra CLIMATOLOGY_HANDLES to the Sasktran IF.
	The handles will be avilable to Matlab/Python/IDL as strings.
 **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRAN_IFGlobalHandleTable( GlobalClimatologyHandleTable** entry, int* numpoints )
{
	*entry = NULL;
	*numpoints = 0;
	return true;
}

extern "C" DLL_PUBLIC bool SKTRAN_IFSetParentHandleTable(std::map<nxString, CLIMATOLOGY_HANDLE>* parenttable)
{
	SetParentHandleTable(parenttable);
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateEngine		2014-2-8*/
/** DLL Exorted function that is used by the Sasktran User Interface to
 *	export engines to users.
 **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateEngine2( const char* userenginename, ISKEngine_Stub**  engine)
{
	nxString	enginename(userenginename);
	bool		ok = false;

	enginename.MakeUpper();
	if (enginename == "SO")
	{
		*engine = new ISKEngine_Stub_SO;
		ok = (engine != nullptr);
	}
	else if (enginename == "HR")
	{
		*engine = new ISKEngine_Stub_HR;
		ok = (engine != nullptr);
	}
	else if (enginename == "MC")
	{
		*engine = new ISKEngine_Stub_MC;
		ok = (engine != nullptr);
	}
	else if (enginename == "OCC")
	{
		*engine = new ISKEngine_Stub_OCC;
		ok = (engine != nullptr);
	}
	else
	{
		*engine = nullptr;
		nxLog::Record(NXLOG_WARNING,"SKTRANIF_CreateEngine, engine [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*) userenginename);

	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANIF_CreateGeodetic		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" DLL_PUBLIC bool SKTRANIF_CreateGeodetic2( const char* geoidname, ISKGeodetic_Stub**  geoid)
{
	nxString	enginename(geoidname);
	bool		ok = false;

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

extern "C" DLL_PUBLIC bool SKTRANIF_CreateMie2(const char* geoidname, ISKMie_Stub * *geoid)
{
	nxString	enginename(geoidname);
	bool		ok = false;

	enginename.MakeUpper();
	if (enginename == "WISCOMBE")
	{
		*geoid = new ISKMie_Stub_Wiscombe;
		ok = (geoid != nullptr);
	}
	else
	{
		*geoid = nullptr;
		nxLog::Record(NXLOG_WARNING, "SKTRANIF_CreateMie, mie [%s] is not available in this DLL/shareable object. This may mean your registry settings are damaged", (const char*)enginename);

	}
	return ok;
}
