#include <sasktran.h>
#include "include/test_hr.h"

SKTRAN_HR_Test::SKTRAN_HR_Test( bool verbose, double tolerance )
{
	m_hardcodefolder = "../../../../Application/TestInstallation/hardcodetestvalues/";
	m_verbose = verbose;
	m_tol = tolerance;
}

void SKTRAN_HR_Test::MakeAtmosphericState()
{
	skOpticalProperties_RayleighDryAir*		rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;	
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*		o3_opticalprops;			// optical properties of one O3 molecule
	bool									ok = true;

	m_opticalstate = std::unique_ptr<SKTRAN_AtmosphericOpticalState_V21> ( new SKTRAN_AtmosphericOpticalState_V21 );

	rayleigh         = new skOpticalProperties_RayleighDryAir;
	msis90           = new skClimatology_MSIS90;
	o3numberdensity  = new skClimatology_LabowOzoneVMR;
	o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;

	ok =       (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && m_opticalstate->AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90,          rayleigh );
	ok = ok && m_opticalstate->AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
	ok = ok && m_opticalstate->SetAlbedo( 1.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list");
	}
}

void SKTRAN_HR_Test::MakeLinesOfSight()
{
	double	mjd 	= 54832.0;
	double  lat 	= -57.5;
	double	lng 	= 70;
	double 	sza 	= 60;
	double  saa 	= 157.5;
	double  rayazi 	= 0;
	double tanheights_meters[4];
	const int numHeights = 4;

	double deltaHeight  = 2000.0;
	double lowestHeight = 10000.0;

	m_linesofsight = std::unique_ptr<SKTRAN_LineOfSightArray_V21>( new SKTRAN_LineOfSightArray_V21);

	// make the lines of sight
	for (int i=0; i < numHeights; i++)
	{
		tanheights_meters[i] = i*deltaHeight + lowestHeight;	
	}
	m_linesofsight->SetRaysFromTangentHeightArray( mjd, lat, lng, sza, saa, rayazi, tanheights_meters, numHeights, 600000.0, &m_sun );
}

bool SKTRAN_HR_Test::RunCurrentTest()
{
	nx2dArray<double> radiance;
	std::vector<SKTRAN_StokesScalar>	radiance_temp;
	radiance.SetSize( m_linesofsight->NumRays(), m_wavel.size() );

	m_engine->SetSun( m_sun );
	m_engine->ConfigureModel( *m_specs, *m_linesofsight, 0 );
	std::cout << "Running test " << m_testname << ":";
	for( size_t i = 0; i < m_wavel.size(); i++ )
	{
		m_engine->CalculateRadiance(&radiance_temp, m_wavel[i], m_scatterorder, m_opticalstate.get() );
		for( size_t j = 0; j < radiance_temp.size(); j++ )
		{
			radiance.At(j,i) = radiance_temp[j];
		}
	}
	// uncomment to remake the hardcoded values
	//radiance.WriteToTextFile( (m_hardcodefolder + m_testname).c_str(), true, "%.15e");
	for( size_t i = 0; i < m_wavel.size(); i++ )
	{
		for( size_t j = 0; j < radiance_temp.size(); j++ )
		{
			double percdiff;
			percdiff = (m_hardcodevalues.At(j,i) - radiance.At(j,i))/ m_hardcodevalues.At(j,i) * 100;
			if( abs(percdiff) > m_tol )
			{
				std::cout << " Failed" << std::endl;
				DisplayFullTestInfo( radiance );
				return false;
			}
		}
	}
	std::cout << " Passed" << std::endl;
	if( m_verbose )
		DisplayFullTestInfo( radiance );
	return true;
}

void SKTRAN_HR_Test::SetupShortTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine> ( new SKTRAN_HR_Engine );
	m_specs  = std::unique_ptr<SKTRAN_HR_Specs_User> ( new SKTRAN_HR_Specs_User );
	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 50;
	m_testname = "hr_default";

	m_hardcodevalues.SetSize( 4,2 );
	m_hardcodevalues.InputColumnMajorText( (m_hardcodefolder + m_testname).c_str() );
	m_wavel.clear();
	m_wavel.push_back( 340 );
	m_wavel.push_back( 600 );
	
}

void SKTRAN_HR_Test::SetupReplicateSOTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine> ( new SKTRAN_HR_Engine );
	m_specs  = std::unique_ptr<SKTRAN_HR_Specs_User> ( new SKTRAN_HR_Specs_User );
	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 50;
	m_testname = "hr_so-settings";

	m_hardcodevalues.SetSize( 4,2 );
	m_hardcodevalues.InputColumnMajorText( (m_hardcodefolder + m_testname).c_str() );
	m_wavel.clear();
	m_wavel.push_back( 340 );
	m_wavel.push_back( 600 );


	m_specs->IntegratorSpecs().UseLegacySasktran21Technique( true );
	m_specs->DiffuseSpecs().SetForceV21IncomingSphere( true );
	m_specs->RayTracingSpecs().SetShellSpacing( 1000 );
	m_specs->RayTracingSpecs().SetSolarShellSpacing( 1000 );
	m_specs->RayTracingSpecs().SetDiffuseType( SKTRAN_HR_RayTracer_Shells );
	m_specs->RayTracingSpecs().SetLinesOfSightType( SKTRAN_HR_RayTracer_Shells );
	m_specs->RayTracingSpecs().SetSolarType( SKTRAN_HR_RayTracer_Shells );

}

void SKTRAN_HR_Test::SetupSingleScatterTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine> ( new SKTRAN_HR_Engine );
	m_specs  = std::unique_ptr<SKTRAN_HR_Specs_User> ( new SKTRAN_HR_Specs_User );
	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 1;
	m_testname = "hr_singlescatter";

	m_wavel.clear();
	for( double wav = 340; wav < 800; wav += 5 )
		m_wavel.push_back( wav );
	m_hardcodevalues.SetSize( 4, m_wavel.size() );
	m_hardcodevalues.InputColumnMajorText( (m_hardcodefolder + m_testname).c_str() );
}

void SKTRAN_HR_Test::SetupMultipleDiffuseProfileTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine> ( new SKTRAN_HR_Engine );
	m_specs  = std::unique_ptr<SKTRAN_HR_Specs_User> ( new SKTRAN_HR_Specs_User );
	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 50;
	m_testname = "hr_multiplediffuseprofile";

	m_hardcodevalues.SetSize( 4,1 );
	m_hardcodevalues.InputColumnMajorText( (m_hardcodefolder + m_testname).c_str() );
	m_wavel.clear();
	m_wavel.push_back( 340 );

	m_specs->DiffuseSpecs().SetNumProfiles( 5 );
}

void SKTRAN_HR_Test::SetupTwoDimTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine> ( new SKTRAN_HR_Engine );
	m_specs  = std::unique_ptr<SKTRAN_HR_Specs_User> ( new SKTRAN_HR_Specs_User );
	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 50;
	m_testname = "hr_twodim";

	m_hardcodevalues.SetSize( 4,1 );
	m_hardcodevalues.InputColumnMajorText( (m_hardcodefolder + m_testname).c_str() );
	m_wavel.clear();
	m_wavel.push_back( 340 );

	m_specs->OpticalPropertiesSpecs().SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType_LOSPlane );
	m_specs->OpticalPropertiesSpecs().SetNumProfiles( 625 );
}

void SKTRAN_HR_Test::DisplayFullTestInfo( const nx2dArray<double>& radiance )
{
	double percdiff;
	if( m_verbose )
		printf("----------------------- Values -------------------------\n");
	else
		printf("------------------- Failed Values ----------------------\n");
	for( size_t waveidx = 0; waveidx < m_wavel.size(); waveidx++ )
	{
		for( size_t altidx = 0; altidx < m_linesofsight->NumRays(); altidx++ )
		{
			percdiff = (radiance.At(altidx, waveidx) - m_hardcodevalues.At(altidx, waveidx)) / m_hardcodevalues.At(altidx, waveidx) * 100;
			if( abs(percdiff) > m_tol || m_verbose )
				printf("Wavelength %f nm, losidx %i: %19.15f %%\n", m_wavel[waveidx], (unsigned int)altidx, percdiff);
		}
	}
	printf("--------------------------------------------------------\n");
}

void SKTRAN_HR_Test::RunShortTests()
{
	SetupShortTest();
	RunCurrentTest();
}

void SKTRAN_HR_Test::RunAllTests()
{
	SetupShortTest();
	RunCurrentTest();

	SetupReplicateSOTest();
	RunCurrentTest();

	SetupSingleScatterTest();
	RunCurrentTest();

	SetupMultipleDiffuseProfileTest();
	RunCurrentTest();

	SetupTwoDimTest();
	RunCurrentTest();

}