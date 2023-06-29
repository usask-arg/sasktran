#include <sasktran.h>
//#include <sasktranv21.h>
#include "include/compare_sktran_to_hr.h"
#include <boost/timer/timer.hpp>

bool LoadOzoneTextCLimatology(const char* foldername, std::vector<double>& alt, std::vector<double>& lat, std::vector<double>& lon, nx3dArray<double>& profile)
{
	bool ok = true;
	char filename[100];
	double temp;

	sprintf(filename, "%s/alts.txt", foldername);
	std::ifstream infile;
	infile.open( filename );
	while(infile.good())
	{
		infile >> temp;
		alt.push_back(temp);
	}
	alt.pop_back();
	infile.close();
	infile.clear();
	sprintf(filename, "%s/lats.txt", foldername);
	infile.open(filename);
	do
	{
		infile >> temp;
		lat.push_back(temp);
	} while(infile.good());
	lat.pop_back();
	infile.close();
	infile.clear();
	sprintf(filename, "%s/lons.txt", foldername);
	infile.open(filename);
	do
	{
		infile >> temp;
		lon.push_back(temp + 180);
	} while(infile.good());
	lon.pop_back();
	infile.close();
	infile.clear();
	profile.SetSize( lat.size(), lon.size(), alt.size() );
	sprintf(filename, "%s/profile.txt", foldername);
	infile.open(filename);
	size_t count = 0;
	size_t altidx, lonidx, latidx;
	do
	{
		infile >> temp;
		latidx = count % lat.size();
		altidx = (count / lat.size()) % alt.size();
		lonidx = count / ( lat.size() * alt.size() );
		if( count < alt.size() * lon.size() * lat.size() )
			profile.At( latidx, lonidx, altidx ) = temp;
		++count;
	} while(infile.good());
	infile.close();


	return ok;
}

bool ConfigureStratAerosolOpticalPropsSK( skOpticalProperties_AerosolProfileH2SO4* optaerosol, const char* climatologyDatabaseDir  )
{
    nxString	        aersizeparamfilename;           // file where the aerosol parameters are stroed
    nx2dArray<double>       aerosolparams;                  // the actual array of aerosol parameters
    nx1dArray<double>       buffer;                         // temporary buffer
    nx1dArray<double>       ah;                             // altitudes
    nx1dArray<double>       moderad;                        // mode radius
    nx1dArray<double>       modewidth;                      // mode width
    bool                    ok;

    // get the actual file name + path of the aerosol parameters
    aersizeparamfilename.sprintf( "%s/StratAerosolSizeParameters.txt", climatologyDatabaseDir );

    // load in the data
    ok = aerosolparams.InputColumnMajorText( aersizeparamfilename, 3, 0 );
    if (ok)
    {
        // read in the data to the specified variables
        ah.SetSize       ( aerosolparams.XSize( ) );
        moderad.SetSize	 ( aerosolparams.XSize( ) );
        modewidth.SetSize( aerosolparams.XSize( ) );
        for (size_t altctr=0; altctr<aerosolparams.XSize(); altctr++ )
        {
            ah.At       (altctr) = 1000.0*aerosolparams.At(altctr,0);
            moderad.At  (altctr) = aerosolparams.At(altctr,1);
            modewidth.At(altctr) = aerosolparams.At(altctr,2);
        }
        // set the log normal climatology
        ok = optaerosol->SetLogNormalProfileClimatology( ah.UnsafeArrayBasePtr(), moderad.UnsafeArrayBasePtr(), modewidth.UnsafeArrayBasePtr(), ah.size() );
    }

    if (!ok)
    {
        nxLog::Record(NXLOG_WARNING,"ConfigureAerosolOpticalProps, Error configuring aerosol mie scattering height profile from file <%s>", (const char*)aersizeparamfilename );
    }
    return ok;
}
bool MakeSingleScatterComparisonLinesOfSight( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* sun = NULL )
{
	bool 	ok;
	double	mjd 	= 54832.0;
	double  lat 	= 0;
	double	lng 	= 70;
	double 	sza 	= 60;
	double  saa 	= 157.5;
	double  rayazi 	= 90;
	//double	tanheights_meters[] = {4000.0,8000.0,10000.0,12000.0,14000.0,16000.0,18000.0,20000.0,24000.0,28000.0,32000.0,36000.0,40000.0};
	double tanheights_meters[36];
	//tanheights_meters[0] = 10000;
	for (int i=0; i < 36; i++) tanheights_meters[i] = i*2000.0 + 4000;	
	ok = linesofsight.SetRaysFromTangentHeightArray( mjd, lat, lng, sza, saa, rayazi, tanheights_meters, 36, 600000.0, sun );
	
	
	//for( int i = 0; i < 13; i++ )
	//{
	//	ok = linesofsight.AddEquatorialLineOfSight( sza, saa, tanheights_meters[i], 600000, mjd );
	//}
	
	
	
	//ok = linesofsight.AddEquatorialLineOfSight( sza, saa, 16000, 600000, mjd );
	//ok = linesofsight.AddEquatorialLineOfSight

	return ok;
}

/*-----------------------------------------------------------------------------
 *					MakeSpeciesList		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool CreateOpticalState( SKTRAN_AtmosphericOpticalState_V21* opticalstate )
{
/*	std::vector<double> alt;
	std::vector<double> lat;
	std::vector<double> lon;
	nx3dArray<double>	profile;

	skOpticalProperties_RayleighDryAir*			rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;	
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skClimatology_UserDefinedTable*			aerosol;
	skClimatology_UserDefined_LatLon_Table* o3three;
	skOpticalProperties_AerosolProfileH2SO4* aero_optprop;
	skOpticalProperties_O3_OSIRISRes*			o3_opticalprops;			// optical properties of one O3 molecule
//	skClimatology_Ecmwf*					ecmwf;	
//	skRTExtinction_MieAerosol*				mieAerosol;					// optical properties of one mie aerosol molecule
//	skRTExtinction_AerosolProfile*			aerosolProfile;	
	bool									ok;

	rayleigh         = new skOpticalProperties_RayleighDryAir;
	msis90           = new skClimatology_MSIS90;
	o3numberdensity  = new skClimatology_LabowOzoneVMR;
	o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;
	aerosol			 = new skClimatology_UserDefinedTable(SKCLIMATOLOGY_AEROSOLH2SO4_CM3, "C:/development/Repos_SasktranV3/input/Aerosol_6432012.txt");
	aero_optprop     = new skOpticalProperties_AerosolProfileH2SO4;
	o3three          = new skClimatology_UserDefined_LatLon_Table;
//	mieAerosol		 = new skRTExtinction_MieAerosol;
//	aerosolProfile   = new skRTExtinction_AerosolProfile;
//	ecmwf			 = new skClimatology_Ecmwf;

	
	LoadOzoneTextCLimatology( "C:/development/Repos_SasktranV3/input/", alt, lat, lon, profile );
	ConfigureStratAerosolOpticalPropsSK( aero_optprop, "C:/development/Repos_SasktranV3/input" );
	o3three->LoadProfileFromData( profile, alt, lat, lon, SKCLIMATOLOGY_O3_CM3 );

	opticalstate->erase();
	ok =       (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90,          rayleigh );
	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_PRESSURE_PA, msis90, rayleigh );
	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_TEMPERATURE_K, msis90, rayleigh );
	//ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_AEROSOLH2SO4_CM3, aerosol, aero_optprop );
	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_O3_CM3,               o3three, o3_opticalprops);
	ok = ok && opticalstate->SetAlbedo( 1.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list");
	}
	return ok;
	*/
	return false;
}

void CompareSingleScatter( )
{
	SKTRAN_LineOfSightArray_V21			linesofsight;
	nx2dArray<double>					hr_radiance;
	nx2dArray<double>					sk_radiance;
	nx2dArray<double>					mc_radiance;
	std::vector<double>					wavelen;
	SKTRAN_AtmosphericOpticalState_V21	opticalstate;
	double								saa, sza;
	saa = 157.5;
	sza = 86;
	//nxVector								sun ( 4.999999998815835e-001,   0.000000000000000e+000,   8.660254038528065e-001 );		// The unit vector from Earth to sun 
	//nxVector							sun ( 1.0, 0, 0 );
	//nxVector							sun ( nxmath::sind(sza) * nxmath::cosd(saa), nxmath::sind(sza) * nxmath::sind(saa), nxmath::cosd(sza) );
	nxVector sun;

	CreateOpticalState( &opticalstate );
	//make wavelenegths
	/*
	for(int i = 280; i < 800; i+=5)
	{
		wavelen.push_back(i);
	}
	*/
	
	wavelen.push_back(280);
	wavelen.push_back(300);
	wavelen.push_back(320);
	wavelen.push_back(340);
	wavelen.push_back(360);
	wavelen.push_back(380);
	wavelen.push_back(400);
	wavelen.push_back(450);
	wavelen.push_back(500);
	wavelen.push_back(550);
	wavelen.push_back(600);
	wavelen.push_back(650);
	wavelen.push_back(700);
	wavelen.push_back(750);
	wavelen.push_back(800);
	//wavelen.push_back( 280 );
	MakeSingleScatterComparisonLinesOfSight( linesofsight, &sun );
	CreateOpticalState( &opticalstate );

	//SK_SingleScatter( wavelen, linesofsight, opticalstate, sun, sk_radiance );
	//HR_SingleScatter( wavelen, linesofsight, opticalstate, sun, hr_radiance, true,"C:/development/private_repos/sktran_hr_comparisons/hr_curve.txt" );
	HR_SingleScatter( wavelen, linesofsight, opticalstate, sun, hr_radiance, false,"C:/development/private_repos/sktran_hr_comparisons/hr_3d_equator_ss_ray90_MT_TableFill.txt" );
	//MC_SingleScatter( wavelen, linesofsight, opticalstate, sun, mc_radiance, 0.001 );

	double maxpercdiff = 0;
	double temp;
	for( size_t i = 0; i < sk_radiance.XSize(); i++ )
	{
		for( size_t j = 0; j < sk_radiance.YSize(); j++ )
		{
			temp = fabs(sk_radiance.At(i,j) - hr_radiance.At(i,j)) / sk_radiance.At(i,j) * 100;
			if( temp > maxpercdiff )
			{
				maxpercdiff = temp;
			}
		}
	}
	printf("Maximum difference between HR and SK single scatter is %f%%\n",maxpercdiff);


}

void SK_SingleScatter( const std::vector<double>& wavelen, const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance )
{
#if SKTRAN_HR_VERBOSE_TIMING
	boost::timer::auto_cpu_timer t;
#endif
	SKTRANSO_Engine					sk_engine;
	std::vector<SKTRAN_StokesScalar>    sk_radiance_temp;
	SKTRANSO_SpecificationsUser_Legacy	specs;

	specs.RayTracingRegionManagerVar()->SetSun( sun );
	specs.DiffuseSpecificationsVar()->ConfigureNumberDiffuseProfiles( 1 );
	sk_engine.ConfigureModel( specs, linesofsight, 4 );
	sk_radiance_temp.resize( linesofsight.NumRays() );
	radiance.SetSize( linesofsight.NumRays(), wavelen.size() );
	for( size_t i = 0; i < wavelen.size(); i++ )
	{
		sk_engine.CalculateRadiance(&sk_radiance_temp, wavelen[i], 50, &opticalstate );
		for(size_t j = 0; j < sk_radiance_temp.size(); j++ )
		{
			radiance.At(j,i) = sk_radiance_temp[j];
		}
	}
	
	radiance.WriteToTextFile( "C:/development/private_repos/sktran_hr_comparisons/sk_scatter.txt", true , "%e" );

}

void HR_SingleScatter( const std::vector<double>& wavelen, const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance, bool usecurve, std::string filename )
{
#if SKTRAN_HR_VERBOSE_TIMING
	boost::timer::auto_cpu_timer t;
#endif

	SKTRAN_HR_Engine					hr_engine;
	std::vector<SKTRAN_StokesScalar>    hr_radiance_temp;
	SKTRAN_HR_Specs_User				hr_specs;

	//hr_specs.RayTracingSpecs().SetLinesOfSightType( SKTRAN_HR_RayTracer_Curved );
	//hr_specs.RayTracingSpecs().SetDiffuseType( SKTRAN_HR_RayTracer_Curved );
	hr_specs.IntegratorSpecs().UseLegacySasktran21Technique( false );
	hr_specs.IntegratorSpecs().SetMaxOpticalDepth( 10000 );
	hr_specs.RayTracingSpecs().SetCurvedSeparation( 1000 );
	hr_specs.DiffuseSpecs().SetNumProfiles( 1 );
	hr_specs.OpticalPropertiesSpecs().SetNumProfiles( 100 );
	//hr_specs.OpticalPropertiesSpecs().SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere );
	//hr_specs.DiffuseSpecs().SetNumBeforeHoriz( 12 );
	//hr_specs.DiffuseSpecs().SetNumHoriz( 20 );
	//hr_specs.DiffuseSpecs().SetNumAfterHoriz( 16 );
	hr_specs.RayTracingSpecs().SetShellSpacing( 500 );
	// hr_specs.DiffuseSpecs().SetNumoffLook( 3 );
	// hr_specs.DiffuseSpecs().SetNumOutgoing( 600 );

	hr_engine.SetSun( sun );
	hr_engine.ConfigureModel( hr_specs, linesofsight, 4 );
	hr_radiance_temp.resize( linesofsight.NumRays() );
	radiance.SetSize( linesofsight.NumRays(), wavelen.size() );
	for( size_t i = 0; i < wavelen.size(); i++ )
	{
		hr_engine.CalculateRadiance(&hr_radiance_temp, wavelen[i], 1, &opticalstate );
		//hr_engine.CalculateSecondOrderRadiance(&hr_radiance_temp, wavelen[i],2,&opticalstate );
		for(size_t j = 0; j < hr_radiance_temp.size(); j++ )
		{
			radiance.At(j,i) = hr_radiance_temp[j];
		}
	}
	
	radiance.WriteToTextFile( filename.c_str(), true, "%e" );
	//radiance.WriteToTextFile( "C:/development/private_repos/sktran_hr_comparisons/hr_scatter.txt", true , "%e" );

}

/*
void MC_SingleScatter( const std::vector<double>& wavelen, const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance, const double& mc_precision )
{
#if SKTRAN_HR_VERBOSE_TIMING
	std::timer::auto_cpu_timer			t;
#endif

	SKTRAN_Engine_MC_V21					mc_engine;
	std::vector<SKTRAN_StokesScalar>		mc_radiance_temp;
	SKTRAN_Specifications_MC				mc_specs;
	
	mc_specs.SetSun( sun );
	mc_engine.ConfigureModel( mc_specs, linesofsight, 0 );
	mc_engine.SetNumPhotonsPerLOS((size_t)(1/(mc_precision*mc_precision)));
	mc_engine.SetMinimumRelativeWeight(mc_precision/3.0);

	mc_radiance_temp.resize( linesofsight.NumRays() );
	radiance.SetSize( linesofsight.NumRays(), wavelen.size() );
	for( int i = 0; i < wavelen.size(); i++ )
	{
		mc_engine.CalculateRadiance(&mc_radiance_temp, wavelen[i], 1, &opticalstate );
		for(int j = 0; j < mc_radiance_temp.size(); j++ )
		{
			radiance.At(j,i) = mc_radiance_temp[j];
		}
	}

	radiance.WriteToTextFile( "C:/development/private_repos/sktran_hr_comparisons/mc_singlescatter_600.txt", true , "%e" );
	
}
*/