#include "../dllimplementation/stdafx.h"

#include "yaml-cpp/yaml.h"

/*-----------------------------------------------------------------------------
 *					ConfigureSasktranIFEntry		 2015- 8- 11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ConfigureSasktranIFEntry( const char* sasktranifentry, const char* dllname, bool createentries)
{
	bool	ok;
	nxRegistryConfiguration		register_entry( "USask-ARG", sasktranifentry,       nxRegistryConfiguration::GLOBAL_INI, false);

	if (createentries)
	{
		ok  = register_entry.SetPath( "DLLName", dllname);
		if (!ok)
		{
			printf("ERROR, There was an error configuring the SasktranIF Registry settings. That is a problem. Check pathnames and priveleges"); 
		}
	}
	else
	{
		nxString	value;
		ok  = register_entry.GetPath( "DLLName", &value);
		printf("%s:DLLName = %s\n", (const char*)sasktranifentry, (const char*)value   ); 
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ConfigureSasktranIFEntries		 2015- 8- 11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ConfigureSasktranIFEntries( const char* dllname, bool createentries)
{
	bool	ok1;
	bool	ok = true;

	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Engines/SO",										dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Engines/HR",										dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Engines/MC",										dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Engines/OCC",										dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Geodetic/STANDARD",								dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/SolarSpectrum/SAO2010",							dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/SolarSpectrum/FONTELA_UVIS_3MICRON",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/SolarSpectrum/FONTELA_UVIS_100MICRON",			dllname, createentries); ok = ok && ok1;

	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/MSIS90",								dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/ECMWF",								dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/O3LABOW",								dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/NO2PRATMO",							dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/OSIRISL2_O3RTMODEL_V507",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/OSIRISL2_NO2RTMODEL_V507",			dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/OSIRISL2_AEROSOLRTMODEL_V507",		dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/OSIRISL2_AEROSOLMODERADIUS_V600",		dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/USERDEFINED_PROFILE",					dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/USERDEFINED_PROFILE3D_LATLONHEIGHT",	dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/USERDEFINED_PROFILE_TABLE",			dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/USERDEFINED_PROFILE_PLANE",			dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/ONE_PRESSURE_TEMP",					dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/CONSTANTVALUE",						dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Climatology/LINEARCOMBO",							dllname, createentries); ok = ok && ok1;
     
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/RAYLEIGH",						dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/SIMPLERAYLEIGH",                   dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/INELASTICRAYLEIGH",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/NO2_BURROWS",						dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/NO2_OSIRISRES",					dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/NO2_VANDAELE1998",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_BASSPAURLINEAR",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_BASSPAURQUADRATIC",			dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_DBM",							dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_VOIGT",						dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_GOMEBURROWS",					dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_OSIRISRES",					dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_SCIABOGUMILV3",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_SCIABOGUMILV4",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O3_SERDYUCHENKOV1",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O2_O2_FALLY2000",					dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O2_O2_THALMAN2013",			    dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/O2_O2_HITRAN2016",			    dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/SO2_FREEMAN1984",					dllname, createentries); ok = ok && ok1; 
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/SO2_RUFUS2003",					dllname, createentries); ok = ok && ok1;    
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/SO2_BOGUMIL2003",					dllname, createentries); ok = ok && ok1;    
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/SO2_VANDAELE2009",				dllname, createentries); ok = ok && ok1;

 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/MIEAEROSOL_H2SO4",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/MIEAEROSOL_WATER",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/MIEAEROSOL_ICE",					dllname, createentries); ok = ok && ok1;

 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/MIEAEROSOL_DUST",					dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/TMATRIXAEROSOL_ICE",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/BAUM_ICECRYSTALS",				dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/USERDEFINED_TABLES",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/USERDEFINED_PRESSURE",            dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/USERDEFINED_SCATTERCONSTANTHEIGHT", dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/CONVOLVED_CROSSSECTION",			dllname, createentries); ok = ok && ok1;
 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/OSIRIS_CONVOLVEDDISCRETEWAVELEN",	dllname, createentries); ok = ok && ok1;

 	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_H2O",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CO2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_O3",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_N2O",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CO",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CH4",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_O2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_NO",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_SO2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_NO2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_NH3",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HNO3",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_OH",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HF",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HCL",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HBR",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HI",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CLO",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_OCS",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_H2CO",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HOCL",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_N2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HCN",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CH3CL",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_H2O2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_C2H2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_C2H6",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_PH3",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_COF2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_SF6",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_H2S",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HCOOH",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HO2",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_O",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CLONO2",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_NOPLUS",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HOBR",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_C2H4",				dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CH3OH",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CH3BR",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CH3CN",			dllname, createentries); ok = ok && ok1;
  	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CF4",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_C4H2",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_HC3N",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_H2",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_CS",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_SO3",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_C2N2",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/OpticalProperty/HITRANCHEMICAL_COCL2",			dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Emission/USERDEFINED_WAVELENGTHHEIGHT",			dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Emission/THERMAL",								dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Emission/HITRAN_PHOTOCHEMICAL",					dllname, createentries); ok = ok && ok1;

	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/LAMBERTIAN",									dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/SNOW_KOKHANOVSKY2012",						dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/ROUJEAN",									dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/ROUJEAN_KERNEL",								dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/LI_SPARSE_KERNEL",							dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/LI_DENSE_KERNEL",							dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/LI_SPARSE_RECIPROCAL_KERNEL",				dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/ROSS_THIN_KERNEL",							dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/ROSS_THICK_KERNEL",							dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/COX_MUNK",									dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/RAHMAN",										dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/HAPKE",										dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/LINEAR_COMBINATION",							dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/MODIS",										dllname, createentries); ok = ok && ok1;	
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/USERDEFINED_LATLON",                         dllname, createentries); ok = ok && ok1;
	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/BRDF/SPECTRAL_VARYING",                           dllname, createentries); ok = ok && ok1;

	ok1 = ConfigureSasktranIFEntry ( "/SasktranIF/Mie/WISCOMBE",                                    dllname, createentries); ok = ok && ok1;


	return ok;
}

/*-----------------------------------------------------------------------------
 *					ConfigureSasktranRegistry		2009-12-15*/
/** THIS IS NOW DEPRECATED. IT IS IMPLEMENTED INSIDE THE PYTHON CODE**/
/*---------------------------------------------------------------------------*/

/*
bool ConfigureSasktranRegistry( const char* aerosolcachedir, 
								const char* icecrystalcachedir, 
								const char* baumdir,
								const char* hitrandir,
								const char* pratmodir, 
								const char*	erainterim_ecmwfdir,
								const char*	gribdefs_folder,
								const char*	gribsamps_folder)

{
	bool		ok1;
	bool		ok2;
	bool		ok3;
	bool		ok4;
	bool		ok5;
	bool		ok8;
	bool		ok9;
	bool		ok10;
	bool		ok;

	nxRegistryConfiguration		register_mieaerosol   ( "USask-ARG", "skOpticalProperties/MieAerosol/",				nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_mieiceaerosol( "USask-ARG", "skOpticalProperties/Mie_IceCrystals",			nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_baumice      ( "USask-ARG", "skOpticalProperties/Baum_IceCrystals/Storage",nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_hitran       ( "USask-ARG", "skOpticalProperties/Hitran",					nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_pratmo       ( "USask-ARG", "/Climatology/Pratmo/",						nxRegistryConfiguration::GLOBAL_INI, false);
    nxRegistryConfiguration     register_ecmwf        ( "USask-ARG", "/Climatology/ECMWF/Storage",					nxRegistryConfiguration::GLOBAL_INI, false);
    nxRegistryConfiguration     register_grib         ( "USask-ARG", "/Grib",										nxRegistryConfiguration::GLOBAL_INI, false);
	
	ok1  = register_mieaerosol.   SetPath ( "Aerosol_Cache_Directory",     aerosolcachedir); 
    ok2  = register_mieiceaerosol.SetPath ( "Ice_Crystal_Cache_Directory", icecrystalcachedir); 
    ok3  = register_baumice.      SetPath ( "2014Database",                baumdir); 
    ok4  = register_hitran.       SetPath ( "BaseDirectory",               hitrandir);
	ok5  = register_pratmo.       SetPath ( "Storage",                     pratmodir    );
 	ok8  = register_ecmwf.        SetPath ( "EraInterimDir",               erainterim_ecmwfdir );

#if defined(NX_WINDOWS)
 	ok9  = register_grib.         SetPath ( "definitions_folder",          gribdefs_folder );
 	ok10 = register_grib.         SetPath ( "samples_folder",              gribsamps_folder );
#else
	ok9  = true;
	ok10 = true;
#endif

	ok =    ok1 && ok2 && ok3 && ok4 && ok5  && ok8 && ok9 && ok10;
	if (!ok)
	{
		printf("ERROR, There was an error configuring the Sasktran Core Components Registry settings. That is a problem. Check pathnames and file I/O priveleges");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					ConfigureSasktranRegistry		2009-12-15*/
/** **/
/*---------------------------------------------------------------------------*/

/* bool PrintSasktranRegistry( )
{
	bool		ok1, ok2, ok3, ok4, ok5, ok8, ok9, ok10;
	bool		ok;
	nxString	aerosolcachedir;
	nxString	icecrystalcachedir;
	nxString	baumdir;
	nxString	hitrandir;
	nxString	ames_ecmwfdir;
	nxString	odinnwp_ecmwfdir;
	nxString	erainterim_ecmwfdir;
	nxString	gribdefs_folder;
	nxString	gribsamps_folder;
	nxString	pratmodir;

	nxRegistryConfiguration		register_mieaerosol   ( "USask-ARG", "skOpticalProperties/MieAerosol/",				nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_mieiceaerosol( "USask-ARG", "skOpticalProperties/Mie_IceCrystals",			nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_baumice      ( "USask-ARG", "skOpticalProperties/Baum_IceCrystals/Storage",nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_hitran       ( "USask-ARG", "skOpticalProperties/Hitran",					nxRegistryConfiguration::GLOBAL_INI, false);
	nxRegistryConfiguration		register_pratmo       ( "USask-ARG", "Climatology/Pratmo/",						    nxRegistryConfiguration::GLOBAL_INI, false);
    nxRegistryConfiguration     register_ecmwf        ( "USask-ARG", "Climatology/ECMWF/Storage",					nxRegistryConfiguration::GLOBAL_INI, false);
    nxRegistryConfiguration     register_grib         ( "USask-ARG", "Grib",										nxRegistryConfiguration::GLOBAL_INI, false);
	
	ok1  = register_mieaerosol.   GetPath ( "Aerosol_Cache_Directory",     &aerosolcachedir); 
    ok2  = register_mieiceaerosol.GetPath ( "Ice_Crystal_Cache_Directory", &icecrystalcachedir); 
    ok3  = register_baumice.      GetPath ( "2014Database",                &baumdir); 
    ok4  = register_hitran.       GetPath ( "BaseDirectory",               &hitrandir);
	ok5  = register_pratmo.       GetPath ( "Storage",                     &pratmodir    );
 	ok8  = register_ecmwf.        GetPath ( "EraInterimDir",               &erainterim_ecmwfdir );

#if defined (NX_WINDOWS)
 	ok9  = register_grib.         GetPath ( "definitions_folder",          &gribdefs_folder );
 	ok10 = register_grib.         GetPath ( "samples_folder",              &gribsamps_folder );
#else
	ok9 = true;
	ok10 = true;
#endif

	printf("Mie Aerosol cache directory          = %s\n", (const char*) aerosolcachedir); 
	printf("Ice crystal cache directory          = %s\n", (const char*) icecrystalcachedir);
	printf("Baum 2014 Crystal Database directory = %s\n", (const char*) baumdir);
	printf("Hitran Folder location               = %s\n", (const char*) hitrandir);
	printf("Pratmo Mean NO2 climatology file     = %s\n", (const char*) pratmodir);
	printf("ECMWF Era-Interim Grib directory     = %s\n", (const char*) erainterim_ecmwfdir );
#if defined (NX_WINDOWS)
	printf("Grib definitions folder              = %s\n", (const char*) gribdefs_folder );
	printf("Grib samples folder                  = %s\n", (const char*) gribsamps_folder );
#endif


	ok = ok1 && ok2 && ok3 && ok4 && ok5 && ok8 && ok9 && ok10;
	if (!ok)
	{
		printf("ERROR, There was errors reading the Sasktran Core Components Registry settings. That is a problem. Check pathnames and file I/O privileges\n");
	}
	return ok;
}

*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_IFCreateRegistryEntriesForChildDLLImplementation		 2016- 11- 15*/
/** Creates the initial registry settings for the Sasktran Core Components typically
 *  for a Python distribution.
 *
 *	Parameters:
 *	\param paramstr
 *		A list of parameters separated by ";" on Windows or ":" on Linux. 
 *
 *		param[0] =  full path name of this DLL.
 **/	
 /*---------------------------------------------------------------------------*/


bool SKTRAN_IFCreateRegistryEntriesForChildDLLImplementation( const char* paramstr, const char* modulename )
{
	nxStringArray	argv;
	int				argc;
	nxString		dllname;
	bool			ok2, ok3, ok5;
	bool			ok;

	
#if defined(NX_WINDOWS)
	argv.Strtok(paramstr, ";");
#else
	argv.Strtok(paramstr, ":");
#endif

	argc = argv.GetSize();
	ok = (argc >= 1);
	if (ok)
	{
		dllname             = argv.GetAt(0);
		ok2 = ConfigureSasktranIFEntries( (const char*)dllname, true );
		ok3 = nxRegistryConfiguration::FlushRegistry(nxRegistryConfiguration::GLOBAL_INI );
		ok =  ok2 && ok3;
	}
	if (!ok)
	{
		printf("\nERROR, There was an error configuring the Sasktran Registry settings for %s\n",(const char*)modulename);
		printf("\n\n------- registry details for the module --------------\n\n");
		ok5 = ConfigureSasktranIFEntries("", false);
		printf("\n\n------- specific component setting --------------\n\n");
		//ok4 = PrintSasktranRegistry();
	}
	return ok;
}

