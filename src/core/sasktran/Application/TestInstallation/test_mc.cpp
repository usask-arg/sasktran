#include <sasktran.h>
#include "include/test_mc.h"
#include <float.h>


/*---------------------------------------------------------------------------
 *                 SKTRAN_MC_Test::SKTRAN_MC_Test                 2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_MC_Test::SKTRAN_MC_Test( bool verbose, double tolerance )
{
	m_hardcodefolder = "../../../../Application/TestInstallation/hardcodetestvalues/";
	m_verbose = verbose;
	m_tol = tolerance;
}


/*---------------------------------------------------------------------------
 *              SKTRAN_MC_Test::MakeAtmosphericState              2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_MC_Test::MakeAtmosphericState()
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
	//ok = ok && m_opticalstate->AddSpecies( SKCLIMATOLOGY_PRESSURE_PA, msis90, rayleigh );
	//ok = ok && m_opticalstate->AddSpecies( SKCLIMATOLOGY_TEMPERATURE_K, msis90, rayleigh );
	ok = ok && m_opticalstate->AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
	ok = ok && m_opticalstate->SetAlbedo( 1.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list");
	}
}


/*---------------------------------------------------------------------------
 *                SKTRAN_MC_Test::MakeLinesOfSight                2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_MC_Test::MakeLinesOfSight()
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

bool SKTRAN_MC_Test::RunCurrentTest()
{
	nx2dArray<double> radiance;
	nx2dArray<double> variance;
	nx2dArray<double> angles;
	nx1dArray<skRTStokesVector> stokes;
	std::vector<SKTRAN_StokesScalar>	radiance_temp;
	radiance.SetSize( m_linesofsight->NumRays(), m_wavel.size() );
	variance.SetSize( m_linesofsight->NumRays(), m_wavel.size() );
	angles  .SetSize( m_linesofsight->NumRays(), m_wavel.size() );

	m_specs->SetSunGeographicPosition( m_sun );
	m_engine->ConfigureModel( *m_specs, *m_linesofsight, 1 );
	std::cout << "Running test " << m_testname << ":";
	for( size_t i = 0; i < m_wavel.size(); i++ )
	{
		m_engine->CalculateRadiance(&radiance_temp, m_wavel[i], m_scatterorder, m_opticalstate.get() );
		m_engine->GetStokesVectors( stokes );
		for( size_t j = 0; j < radiance_temp.size(); j++ )
		{
			radiance.At(j,i) = radiance_temp[j];
			angles  .At(j,i) = std::atan2(stokes.At(j).At(2),stokes.At(j).At(1));
			m_engine->GetMeasurementVariance(j,variance.At(j,i));
		}
	}
	bool writeVals = false;
	if( writeVals ){
		radiance.WriteToTextFile( (m_hardcodefolder + m_testname         ).c_str(), true, "%1.16e");
		variance.WriteToTextFile( (m_hardcodefolder + m_testname + "_var").c_str(), true, "%1.16e");
		angles  .WriteToTextFile( (m_hardcodefolder + m_testname + "_ang").c_str(), true, "%1.16e");
	}
	for( size_t i = 0; i < m_wavel.size(); i++ )
	{
		for( size_t j = 0; j < radiance_temp.size(); j++ )
		{
			double raddiff, vardiff, angdiff;
			raddiff = (m_hardcoderadiance.At(j,i) - radiance.At(j,i))/ m_hardcoderadiance.At(j,i) * 100.0;
			vardiff = (m_hardcodevariance.At(j,i) - variance.At(j,i))/ m_hardcodevariance.At(j,i) * 100.0;
			angdiff = std::fmod( m_hardcodeangles.At(j,i) - angles.At(j,i), nxmath::TwoPi ) / nxmath::TwoPi * 100.0;
			if( !(raddiff==raddiff&&vardiff==vardiff&&angdiff==angdiff) || abs(raddiff) > m_tol ||abs(vardiff) > m_tol || abs(angdiff)>m_tol )
			{
				std::cout << " Failed" << std::endl;
				DisplayFullTestInfo( radiance, variance, angles );
				return false;
			}
		}
	}
	std::cout << " Passed" << std::endl;
	if( m_verbose )
		DisplayFullTestInfo( radiance, variance, angles );
	return true;
}



void SKTRAN_MC_Test::SetUpDefault ( )
{
	m_engine = std::unique_ptr<SKTRAN_Engine_MC_V21> ( new SKTRAN_Engine_MC_V21 );
	m_specs  = std::unique_ptr<SKTRAN_Specifications_MC> ( new SKTRAN_Specifications_MC );
	
	size_t numPhotonsPerProc = 500;
	size_t numPhotons        = numPhotonsPerProc * 1;// omp_get_max_threads(); // Need to compute 500 photons per processor to get hard coded value

	MakeAtmosphericState();
	MakeLinesOfSight();
	
	m_specs->ConfigureDefaults();
	m_specs->SetSunGeographicPosition ( m_sun );
	m_specs->SetNumPhotonsPerLOS      ( numPhotons );
	m_specs->SetPrecisionMC           ( 0.0 );
	m_specs->SetMinimumRelPathWeight  ( 0.0 );
	m_specs->SetMinFractionHigherOrder( 1.0 );
	m_specs->SetNumRayTracingAlts     ( 100 + 1 );
	m_specs->SetLOSRayTracerType      ( SKTRAN_Specifications_MC::RayTracerType::shell );
	m_specs->SetMSRayTracerType       ( SKTRAN_Specifications_MC::RayTracerType::shell );
	m_specs->SetSolarRayTracerType    ( SKTRAN_Specifications_MC::RayTracerType::shell );
	m_specs->SetSunType               ( SKTRAN_Specifications_MC::SunType::point );
	m_specs->SetSolarTableType        ( SKTRAN_Specifications_MC::SolarTableType::noTable );
	m_specs->SetOptPropIntType        ( SKTRAN_Specifications_MC::OptPropIntType::straight );
	m_specs->SetSunType               ( SKTRAN_Specifications_MC::SunType::point );
	m_specs->SetOptTableType          ( SKTRAN_Specifications_MC::OptTableType::dim1 );
	m_specs->SetPolType               ( SKTRAN_Specifications_MC::PolType::none );

	m_specs->SetAllowDynamicThreads   ( false );
	//m_specs->SetChunkSize             ( 0 );
	m_specs->SetRngSeed               ( 1234 );
	m_specs->SetScatterPositionRes    ( 50.0 );

	m_scatterorder = 50;
	
	m_wavel.clear();
	m_wavel.push_back( 340 );
	m_wavel.push_back( 600 ); 
}

void SKTRAN_MC_Test::SetUpShortTest ( )
{
	m_testname = "MC_ShortTest";

	SetUpDefault ( );
	m_specs->FinalizeSpecs ( );
	
	LoadHardcodedValues( 4, 2 );
}


void SKTRAN_MC_Test::SetUpSingleScattTest ( )
{
	m_testname = "MC_SingleScatter";

	SetUpDefault ( );
	m_scatterorder = 1;
	m_specs->FinalizeSpecs ( );
	
	LoadHardcodedValues( 4, 2 );
}

void SKTRAN_MC_Test::SetUpSecondOrderTest ( )
{
	m_testname = "MC_SecondOrder";

	SetUpDefault ( );
	m_scatterorder = 2;
	m_specs->FinalizeSpecs ( );
	
	LoadHardcodedValues( 4, 2 );
}

void SKTRAN_MC_Test::SetUpPseudoPolTest ( )
{
	m_testname = "MC_PseudoPol";

	SetUpDefault ( );
	m_specs->SetPolType ( SKTRAN_Specifications_MC::PolType::pv1 );
	m_specs->FinalizeSpecs ( );

	LoadHardcodedValues( 4, 2 );
}

void SKTRAN_MC_Test::SetUpPolarizationTest ( )
{
	m_testname = "MC_Polarized";

	SetUpDefault ( );
	m_specs->SetPolType ( SKTRAN_Specifications_MC::PolType::pol );
	m_specs->FinalizeSpecs ( );

	LoadHardcodedValues( 4, 2 );
}

void SKTRAN_MC_Test::SetUpHighPrecisionTest ( )
{

	m_testname = "MC_HighPrecision";

	SetUpDefault ( );
	
	size_t numPhotonsPerProc = 2000;
	size_t numPhotons        = numPhotonsPerProc * 1;// omp_get_max_threads(); // Need to compute 500 photons per processor to get hard coded value
	m_specs->SetNumPhotonsPerLOS      ( numPhotons );
	m_specs->SetPolType ( SKTRAN_Specifications_MC::PolType::pol );
	m_specs->FinalizeSpecs ( );

	LoadHardcodedValues( 4, 2 );
}

void SKTRAN_MC_Test::LoadHardcodedValues( int numLOS, int numWavs )
{
	bool ok;

	m_hardcoderadiance .SetSize( numLOS, numWavs );
	ok = m_hardcoderadiance .InputColumnMajorText( (m_hardcodefolder + m_testname         ).c_str() );
	if(!ok){ m_hardcoderadiance.SetSize ( numLOS, numWavs ); m_hardcoderadiance.SetTo(std::numeric_limits<double>::infinity());}
	m_hardcodevariance .SetSize ( numLOS, numWavs );
	ok = m_hardcodevariance .InputColumnMajorText( (m_hardcodefolder + m_testname + "_var").c_str() );
	if(!ok){ m_hardcodevariance.SetSize ( numLOS, numWavs ); m_hardcodevariance.SetTo(std::numeric_limits<double>::infinity());}
	m_hardcodeangles   .SetSize ( numLOS, numWavs );
	ok = m_hardcodeangles   .InputColumnMajorText( (m_hardcodefolder + m_testname + "_ang").c_str() );
	if(!ok){ m_hardcodeangles  .SetSize ( numLOS, numWavs ); m_hardcodeangles  .SetTo(std::numeric_limits<double>::infinity());}
}


void SKTRAN_MC_Test::DisplayFullTestInfo( const nx2dArray<double>& radiance, const nx2dArray<double>& variance, const nx2dArray<double>& angles )
{
	double raddiff, vardiff, angdiff;

	if( m_verbose )
		printf("------------------------------- Values ---------------------------------\n");
	else
		printf("--------------------------- Failed Values ------------------------------\n");
	for( size_t waveidx = 0; waveidx < m_wavel.size(); waveidx++ )
	{
		for( size_t altidx = 0; altidx < m_linesofsight->NumRays(); altidx++ )
		{
			raddiff = (radiance.At(altidx, waveidx) - m_hardcoderadiance.At(altidx, waveidx)) / m_hardcoderadiance.At(altidx, waveidx) * 100;
			vardiff = (variance.At(altidx, waveidx) - m_hardcodevariance.At(altidx, waveidx)) / m_hardcodevariance.At(altidx, waveidx) * 100;
			angdiff = std::fmod( m_hardcodeangles.At(altidx,waveidx) - angles.At(altidx,waveidx), nxmath::TwoPi ) / nxmath::TwoPi * 100.0;
			if( !(raddiff==raddiff&&vardiff==vardiff&&angdiff==angdiff) || abs(raddiff) > m_tol ||abs(vardiff) > m_tol || abs(angdiff)>m_tol || m_verbose )
				printf("Wavelength %f nm, losidx %i:\n   r:%1.15e v:%1.15e a:%1.15f (%%)\n) ", m_wavel[waveidx], (unsigned int)altidx, raddiff, vardiff, angdiff);
		}
	}
	printf("------------------------------------------------------------------------\n");
}

void SKTRAN_MC_Test::RunShortTests()
{
	SetUpShortTest();
	RunCurrentTest();
}

void SKTRAN_MC_Test::RunAllTests()
{
	SetUpShortTest();
	RunCurrentTest();
	
	SetUpSingleScattTest();
	RunCurrentTest();
	
	SetUpSecondOrderTest();
	RunCurrentTest();

	SetUpPseudoPolTest();
	RunCurrentTest();

	SetUpPolarizationTest();
	RunCurrentTest();
	
	SetUpHighPrecisionTest();
	RunCurrentTest();


}