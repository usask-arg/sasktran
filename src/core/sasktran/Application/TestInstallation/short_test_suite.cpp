#include <sasktran.h>
#include "include/short_test_suite.h"
#include <boost/timer/timer.hpp>

/*---------------------------------------------------------------------------
 *               SKTRAN_Short_Test_Base::Initialize               2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::Initialize( GeometryType geometrytype )
{
	bool ok = true;

	m_geometrytype = geometrytype;
	ok = ok && MakeOpticalState();
	switch (m_geometrytype)
	{
		case GeometryType::limb:
		{
			ok = ok && MakeLimbBasedLinesOfSight();
			break;
		}
		case GeometryType::ground:
		{
			ok = ok && MakeGroundBasedLinesOfSight();
			break;
		}
		case GeometryType::nadir:
		{
			ok = ok && MakeNadirBasedLinesOfSight();
			break;
		}
		default:
		{
			ok = false;
		}
	}

	ok = ok && MakeWavelengths();
	ok = ok && MakeSpecs();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::MakeOpticalState		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::MakeOpticalState( )
{
	skOpticalProperties_RayleighDryAir*		rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;	
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*		o3_opticalprops;			// optical properties of one O3 molecule
	bool									ok = true;

	m_opticalstate.erase();
	rayleigh         = new skOpticalProperties_RayleighDryAir;
	msis90           = new skClimatology_MSIS90;
	o3numberdensity  = new skClimatology_LabowOzoneVMR;
	o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;

	ok =       (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90,     rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity,			o3_opticalprops);
	ok = ok && m_opticalstate.SetAlbedo( 1.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list");
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *       SKTRAN_Short_Test_Base::MakeLimbBasedLinesOfSight        2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::MakeLimbBasedLinesOfSight()
{
	bool 	ok;
	double	mjd 	= 54832.0;
	double  lat 	= -57.5;
	double	lng 	= 70;
	double 	sza 	= 60;
	double  saa 	= 157.5;
	double  rayazi 	= 0;
	double tanheights_meters[4];
	const SKTRAN_LineOfSightEntry_V2* entry;
	const int numHeights = 4;

	double deltaHeight  = 2000.0;
	double lowestHeight = 10000.0;

	// make the lines of sight
	for (int i=0; i < numHeights; i++) tanheights_meters[i] = i*deltaHeight + lowestHeight;	
	ok = m_linesofsight.SetRaysFromTangentHeightArray( mjd, lat, lng, sza, saa, rayazi, tanheights_meters, numHeights, 600000.0, &m_sun );
	
	//printf("Test Installation Sun location   = [ %25.18e, %25.18e, %25.18e ]\n", (double)m_sun.X(), (double)m_sun.Y(), (double)m_sun.Z() );
	//for (int i=0; i < numHeights; i++)
	//{
	//	m_linesofsight.GetRay(i,&entry);
	//	printf("Test Installation, Line of sight = %25.18e [ %25.18e, %25.18e, %25.18e ] [ %25.18e, %25.18e, %25.18e ] \n",(double)entry->Mjd(), (double)entry->Observer().X(), (double)entry->Observer().Y(), (double)entry->Observer().Z(), (double)entry->Look().X(), (double)entry->Look().Y(), (double)entry->Look().Z() );
	//}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::MakeGroundBasedLinesOfSight		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::MakeGroundBasedLinesOfSight()
{
	double			lat			= 48.4667;
	double			lng			= 81.3333;
	double			height		= 2000.0;
	double			solarazi    = 180;				// Sun is in the South 
	double			solarzen    = 60.0;				//     at 60 degrees zenitha angle
	double			obsazi      = 90.0;				// Observer is look East
	double			obszen;
	double			mjd			= 54832.0;
	double			azi;
	nxGeodetic		geoid;
	nxVector		observer;
	nxVector		look;
	nxVector		west, south, up;
	double			Re;

	geoid.FromGeodetic( lat,lng,height );
	observer = geoid.Location();
	Re       = observer.Magnitude() - height;
	geoid.GetGeodeticWestSouthUp( &west, &south, &up );
	azi = 270.0-solarazi;										// Convert compass azimuth to west, south, up coordinate system
	m_sun =  (nxmath::sind(solarzen)*nxmath::cosd(azi))*west  + (nxmath::sind(solarzen)*nxmath::sind(azi))*south + nxmath::cosd(solarzen)*up;

	azi = 270.0 - obsazi;
	for(int i = 0; i <= 30; i++ )
	{
		obszen = 3*i;
		look   =  (nxmath::sind(obszen)*nxmath::cosd(azi))*west  + (nxmath::sind(obszen)*nxmath::sind(azi))*south + nxmath::cosd(obszen)*up;
		m_linesofsight.AddLineOfSight( observer, look, mjd );
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::MakeGroundBasedLinesOfSight		2014-4-7*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::MakeNadirBasedLinesOfSight()
{
	double			lat = 48.4667;
	double			lng = 81.3333;
	double			solarazi = 180;				// Sun is in the South 
	double			solarzen = 60.0;				//     at 60 degrees zenitha angle
	double			obsazi = 90.0;				// Observer is look East
	double			mjd = 54832.0;
	double			azi;
	nxGeodetic		geoid;
	nxVector		observer;
	nxVector		ground;
	nxVector		look;
	nxVector		west, south, up;
	double			Re;

	geoid.FromGeodetic(lat, lng, 0.0);
	ground = geoid.Location();
	Re = ground.Magnitude();
	geoid.GetGeodeticWestSouthUp(&west, &south, &up);
	azi = 270.0 - solarazi;										// Convert compass azimuth to west, south, up coordinate system
	m_sun = (nxmath::sind(solarzen) * nxmath::cosd(azi)) * west + (nxmath::sind(solarzen) * nxmath::sind(azi)) * south + nxmath::cosd(solarzen) * up;

	azi = 270.0 - obsazi;
	for (auto&& obszen : std::vector<double>({ 89.0 }))//, 85.0, 80.0, 72.0, 60.0, 30.0}))
	{
		look = -(nxmath::sind(obszen) * nxmath::cosd(azi)) * west - (nxmath::sind(obszen) * nxmath::sind(azi)) * south - nxmath::cosd(obszen) * up;
		observer = ground - 2e6 * look;
		m_linesofsight.AddLineOfSight(observer, look, mjd);
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test_Base::MakeWavelengths		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::MakeWavelengths()
{
	m_wavelen.clear();

	m_wavelen.push_back( 600 );
	m_wavelen.push_back( 340 );

	return true;
}



/*---------------------------------------------------------------------------
 *                  SKTRAN_HR_Test::MakeHRSpecs                   2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_HR::MakeSpecs()
{
	bool ok = true;

	bool debugmode = false;
	if(debugmode)
	{
		if (m_hrspecs == nullptr) m_hrspecs = new SKTRAN_HR_Specs_User;
		m_hrspecs->IntegratorSpecs().SetMaxOpticalDepth( 10000 );
		m_hrspecs->IntegratorSpecs().UseLegacySasktran21Technique( false );
		m_hrspecs->RayTracingSpecs().SetShellSpacing( 5000.0 );
		m_hrspecs->RayTracingSpecs().SetSolarShellSpacing( 500 );
		m_hrspecs->DiffuseSpecs().SetHeightRes(5000.0);
		m_hrspecs->DiffuseSpecs().SetNumAzi(5);
		m_hrspecs->DiffuseSpecs().SetNumHoriz(5);
		m_hrspecs->DiffuseSpecs().SetNumAfterHoriz(5);
		m_hrspecs->DiffuseSpecs().SetNumBeforeHoriz(5);
	} else
	{
		if (m_hrspecs == nullptr) m_hrspecs = new SKTRAN_HR_Specs_User;
		m_hrspecs->IntegratorSpecs().SetMaxOpticalDepth( 10000 );
		m_hrspecs->IntegratorSpecs().UseLegacySasktran21Technique( false );
		m_hrspecs->RayTracingSpecs().SetShellSpacing( 1000 );
		m_hrspecs->RayTracingSpecs().SetSolarShellSpacing( 500 );
	}
	m_hrspecs->OpticalPropertiesSpecs().SetMaxPolarizationOrder( 0 );
	m_hrspecs->OpticalPropertiesSpecs().SetPolarizationHigherOrderBehaviour( SKTRAN_HR_PolHOType::unpolarized );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::HR_Scatter		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_HR::Scatter( nx2dArray<double>& radiance, size_t scattorder)
{
	bool								ok = true;
	bool								ok1;
	SKTRAN_HR_Engine					hr_engine;
	std::vector<SKTRAN_StokesScalar>	hr_radiance_temp;

	ok = radiance.SetSize( m_linesofsight.NumRays(), m_wavelen.size() );
	hr_radiance_temp.resize( m_linesofsight.NumRays() );
	ok =       hr_engine.SetSun( m_sun );
	ok = ok && hr_engine.ConfigureModel( *m_hrspecs, m_linesofsight, 0 );
	if (ok)
	{
		for( int i = 0; i < (int)m_wavelen.size(); i++ )
		{
			ok1 = hr_engine.CalculateRadiance(&hr_radiance_temp, m_wavelen[i], scattorder, &m_opticalstate );
			ok  = ok && ok1;
			if (ok1)
			{
				for(int j = 0; j < (int)hr_radiance_temp.size(); j++ )
				{
					radiance.At(j,i) = hr_radiance_temp[j];
				}
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Short_Test::HR_Scatter, There were errors executing the HR model. This should be investigated in more detail.");
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *              SKTRAN_HR_Test::LoadHardCodedValues               2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_HR::LoadHardCodedValues( nx2dArray<double>& hardcode )
{
	bool ok = true;

//  double hr[] = { 5.5766149230950877e-002, 4.5255524380402937e-002, 3.6227758391763336e-002, 2.8765117789293491e-002, 
//                  1.5636573909540583e-001, 1.5547102189648537e-001, 1.5436932643191056e-001, 1.5266761467383869e-001 
//                };
//	
//	double hr[] = { 5.5766147577812335e-002, 4.5255523495638238e-002, 3.6227756700478939e-002, 2.8765116972648257e-002, 
//					1.5636572012370423e-001, 1.5547100354938806e-001, 1.5436930906847379e-001, 1.5266759717037634e-001 
//	};

	double hr[] = { 5.5732634669672404e-02, 4.5223961154029879e-02, 3.6201236594918210e-02, 2.8744647800144005e-02, 
                    1.5626052284053002e-01, 1.5536235331861920e-01, 1.5425566694176107e-01, 1.5255063507284350e-01
                  };


	hardcode    = nx2dArray<double>(4,2,hr);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::MakeMCSpecs		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_MC::MakeSpecs()
{
	bool ok = true;


	size_t numPhotonsPerProc = 500;
	size_t numPhotons        = numPhotonsPerProc * 1; // Need to compute 500 photons per processor to get hard coded value

	if (m_mcspecs == nullptr) m_mcspecs = new SKTRAN_Specifications_MC;
	ok = ok && m_mcspecs->ConfigureDefaults();
	ok = ok && m_mcspecs->SetSunGeographicPosition ( m_sun );
	ok = ok && m_mcspecs->SetNumPhotonsPerLOS      ( numPhotons );
	ok = ok && m_mcspecs->SetPrecisionMC           ( 0.0 );
	ok = ok && m_mcspecs->SetMinimumRelPathWeight  ( 0.0 );
	ok = ok && m_mcspecs->SetNumRayTracingAlts     ( 100 + 1 );
	ok = ok && m_mcspecs->SetLOSRayTracerType      ( SKTRAN_Specifications_MC::RayTracerType::shell );
	ok = ok && m_mcspecs->SetMSRayTracerType       ( SKTRAN_Specifications_MC::RayTracerType::shell );
	ok = ok && m_mcspecs->SetSolarRayTracerType    ( SKTRAN_Specifications_MC::RayTracerType::shell );
	ok = ok && m_mcspecs->SetSunType               ( SKTRAN_Specifications_MC::SunType::point );
	ok = ok && m_mcspecs->SetSolarTableType        ( SKTRAN_Specifications_MC::SolarTableType::noTable );
	ok = ok && m_mcspecs->SetOptPropIntType        ( SKTRAN_Specifications_MC::OptPropIntType::straight );
	ok = ok && m_mcspecs->SetOptTableType          ( SKTRAN_Specifications_MC::OptTableType::dim1 );
	ok = ok && m_mcspecs->SetAllowDynamicThreads   ( false );
	//ok = ok && m_mcspecs->SetChunkSize             ( 0 );
	ok = ok && m_mcspecs->SetRngSeed               ( 1234 );
	ok = ok && m_mcspecs->SetScatterPositionRes    ( 50.0 );
	ok = ok && m_mcspecs->SetMinFractionHigherOrder( 1.0 );

	ok = ok && m_mcspecs->SetPolType( SKTRAN_Specifications_MC::PolType::none );
	ok = ok && m_mcspecs->FinalizeSpecs();

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::MC_Scatter		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_MC::Scatter( nx2dArray<double>& radiance, size_t scattorder )
{
	bool								ok = true;
		bool								ok1;
		SKTRAN_Engine_MC_V21				mc_engine;
		std::vector<SKTRAN_StokesScalar>	mc_radiance_temp;
		size_t numThreads = 1;

		std::vector<skRTStokesVector> svecs;
		ok = radiance.SetSize( m_linesofsight.NumRays(), m_wavelen.size() );
		mc_radiance_temp.resize( m_linesofsight.NumRays() );
		ok =       m_mcspecs->SetSunGeographicPosition( m_sun );
		ok = ok && m_mcspecs->FinalizeSpecs();
		ok = ok && mc_engine.ConfigureModel( *m_mcspecs, m_linesofsight, numThreads ); 
		ok = ok && radiance.SetSize( m_linesofsight.NumRays(), m_wavelen.size() );
		if (ok)
		{
			for( size_t i = 0; i < m_wavelen.size(); i++ )
			{
				ok1 = mc_engine.CalculateRadiance(&mc_radiance_temp, m_wavelen[i], scattorder, &m_opticalstate, &svecs );
				ok = ok && ok1;
				if (ok1)
				{
					for( size_t j = 0; j < mc_radiance_temp.size(); j++ )
					{
						radiance.At(j,i) = mc_radiance_temp[j];
					}
				}
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_Short_Test::MC_Scatter, There were errors executing the MC model. This should be investigated in more detail.");
		}

		return ok;
	}



	/*---------------------------------------------------------------------------
	 *              SKTRAN_MC_Test::LoadHardCodedValues               2020-01-31 */
	/** **/
	/*---------------------------------------------------------------------------*/

	bool SKTRAN_Short_Test_MC::LoadHardCodedValues( nx2dArray<double>& hardcode )
	{
		bool ok = true;


		//double mc[] = { 5.5724447594630168e-002, 4.7369386898081893e-002, 3.6310971905386479e-002, 2.8946355947833770e-002, 
		//               1.6274666281905237e-001, 1.6661422591543568e-001, 1.5977541781659052e-001, 1.5845057932349851e-001
		//};
		double mc[] = { 5.5447117636213766e-02, 4.4223540943734159e-02, 3.7519314370382577e-02, 3.0172066272158241e-02,
						1.6468969860965310e-01, 1.6016495844526352e-01, 1.6321000676050257e-01, 1.5193692132980371e-01
		};

		hardcode = nx2dArray<double>(4,2,mc);

		return ok;
	}

/*---------------------------------------------------------------------------
 *          SKTRAN_Short_Test_Base::RunGroundBasedTests           2020-01-31 */
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_Short_Test_Base::RunGroundBasedTests()
{
	bool					ok;
	nx2dArray<double>		sk_radiance;

	ok =	   Initialize( GeometryType::ground );
	ok = ok && sk_radiance.SetSize( m_linesofsight.NumRays(), m_wavelen.size() );
	sk_radiance.SetTo(0.0);

	printf("---- Starting ground based Test ----\n");
	ok = ok && Scatter( sk_radiance, 10 );
	

	printf("------------------------------------------------\n");
	for( size_t waveidx = 0; waveidx < m_wavelen.size(); waveidx++ )
	{
		printf("Wavelength %f nm\n", m_wavelen[waveidx]);
		printf("------ SK ---------------- HR ---------------- MC ------  |\n");
		for( size_t altidx = 0; altidx < m_linesofsight.NumRays(); altidx++ )
		{
			printf("%.10e \n", (double)sk_radiance.At(altidx, waveidx) );
		}
		printf("----------------------------------------------------------\n");
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_Short_Test::RunGroundBasedTests(), There were errors executing the ground based tests. This should be investigated");
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *          SKTRAN_Short_Test_Base::RunNadirBasedTests           2022-06-30 */
 /** **/
 /*---------------------------------------------------------------------------*/
bool SKTRAN_Short_Test_Base::RunNadirBasedTests()
{
	bool					ok;
	nx2dArray<double>		sk_radiance;

	ok = Initialize( GeometryType::nadir );
	ok = ok && sk_radiance.SetSize(m_linesofsight.NumRays(), m_wavelen.size());
	sk_radiance.SetTo(0.0);

	printf("---- Starting ground based Test ----\n");
	ok = ok && Scatter(sk_radiance, 10);


	printf("------------------------------------------------\n");
	for (size_t waveidx = 0; waveidx < m_wavelen.size(); waveidx++)
	{
		printf("Wavelength %f nm\n", m_wavelen[waveidx]);
		printf("------ SK ---------------- HR ---------------- MC ------  |\n");
		for (size_t altidx = 0; altidx < m_linesofsight.NumRays(); altidx++)
		{
			printf("%.10e \n", (double)sk_radiance.At(altidx, waveidx));
		}
		printf("----------------------------------------------------------\n");
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Short_Test::RunNadirBasedTests(), There were errors executing the nadir based tests. This should be investigated");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::TestHR		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Base::RunStandardTest( )
{
	bool					ok = true;
	nx2dArray<double>		sk_radiance;
	nx2dArray<double>		sk_hardcode;
	double					skpercdiff;

	ok = ok && Initialize( GeometryType::limb );
	ok = ok && LoadHardCodedValues( sk_hardcode );
	ok = ok && sk_radiance.SetSize(sk_hardcode.XSize(),sk_hardcode.YSize());
	sk_radiance.SetTo(0.0);

	printf("%s Short Limb Test: tolerance = %7.3e percent\n", (const char*)Name(), (double)m_errortolerance );
	ok = ok && Scatter( sk_radiance, 50 );
    if(!ok) printf("\n\nCOULD NOT RUN TESTS.\nAll tests should run cleanly for release code.\n\n");

	for( size_t waveidx = 0; waveidx < m_wavelen.size(); waveidx++ )
	{
		printf("  Wavelength %f nm\n", m_wavelen[waveidx]);
 		for( size_t altidx = 0; altidx < m_linesofsight.NumRays(); altidx++ )
		{
			skpercdiff = (sk_radiance.At(altidx, waveidx) - sk_hardcode.At(altidx, waveidx)) / sk_hardcode.At(altidx, waveidx) * 100;
			printf("   [%2d] this = %19.14e hardcode= %19.14e percent diff = %19.15f\n", (int)altidx, (double)sk_radiance.At(altidx, waveidx), (double)sk_hardcode.At(altidx, waveidx),  (double)skpercdiff);
			ok =    fabs(skpercdiff)  < m_errortolerance;
		}
		printf( ok ? "   PASS\n" : "   FAIL\n");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::SKTRAN_Short_Test		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_Short_Test::SKTRAN_Short_Test ()
{
	m_dolimb = true;
	m_doHR   = true;
	m_doSO   = true;
	m_doMC   = true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::~SKTRAN_Short_Test		2014-4-9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_Short_Test::~SKTRAN_Short_Test()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test::RunTests		2014-4-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test::RunTests( bool doso, bool dohr, bool domc, double errortolerance)
{
	SKTRAN_Short_Test_Base*	shorttest = nullptr;
	bool					ok = true;

	for (int i = 0; i < 3; i++)
	{
		switch  (i)
		{
		case 1 : shorttest = dohr ? new SKTRAN_Short_Test_HR(errortolerance) : nullptr; break;
		case 2 : shorttest = domc ? new SKTRAN_Short_Test_MC(errortolerance) : nullptr; break;
		};
		if (shorttest != nullptr)
		{
			ok = ok && shorttest->RunStandardTest();
			delete shorttest;
			shorttest = nullptr;
		}
	}
	return ok;
}

bool SKTRAN_Short_Test::RunGroundBasedTests	(bool doso, bool dohr, bool domc)
{
		return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test_Inelastic_MC::MakeMCSpecs		2020-04-03*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Inelastic_MC::MakeSpecs()
{
	bool ok = true;


	size_t numPhotonsPerProc = 500;
	size_t numPhotons = numPhotonsPerProc * 1; // Need to compute 500 photons per processor to get hard coded value

	if (m_mcspecs == nullptr) m_mcspecs = new SKTRAN_Specifications_MC;
	ok = ok && m_mcspecs->ConfigureDefaults();
	ok = ok && m_mcspecs->SetSunGeographicPosition(m_sun);
	ok = ok && m_mcspecs->SetNumPhotonsPerLOS(numPhotons);
	ok = ok && m_mcspecs->SetPrecisionMC(0.0);
	ok = ok && m_mcspecs->SetMinimumRelPathWeight(0.0);
	ok = ok && m_mcspecs->SetNumRayTracingAlts(100 + 1);
	ok = ok && m_mcspecs->SetLOSRayTracerType(SKTRAN_Specifications_MC::RayTracerType::shell);
	ok = ok && m_mcspecs->SetMSRayTracerType(SKTRAN_Specifications_MC::RayTracerType::shell);
	ok = ok && m_mcspecs->SetSolarRayTracerType(SKTRAN_Specifications_MC::RayTracerType::shell);
	ok = ok && m_mcspecs->SetSunType(SKTRAN_Specifications_MC::SunType::point);
	ok = ok && m_mcspecs->SetSolarTableType(SKTRAN_Specifications_MC::SolarTableType::noTable);
	ok = ok && m_mcspecs->SetOptPropIntType(SKTRAN_Specifications_MC::OptPropIntType::straight);
	ok = ok && m_mcspecs->SetOptTableType(SKTRAN_Specifications_MC::OptTableType::dim1);
	ok = ok && m_mcspecs->SetAllowDynamicThreads(false);
	//ok = ok && m_mcspecs->SetChunkSize             ( 0 );
	ok = ok && m_mcspecs->SetRngSeed(1234);
	ok = ok && m_mcspecs->SetScatterPositionRes(50.0);
	ok = ok && m_mcspecs->SetMinFractionHigherOrder(1.0);

	ok = ok && m_mcspecs->SetPolType(SKTRAN_Specifications_MC::PolType::none);

	if (m_simultaneous)
	{
		m_mcspecs->SetWavelengthType(SKTRAN_Specifications_MC::WavelengthType::simultaneous);
		m_mcspecs->SetRadianceWavelengths(m_wavelen);
		m_mcspecs->SetPrimaryWavelength(m_wavelen[0]);
	}

	if (!m_optimized)
	{
		if (!m_ring) 
			m_mcspecs->SetScatterType(SKTRAN_Specifications_MC::ScatterType::inelastic);
		else 
			m_mcspecs->SetScatterType(SKTRAN_Specifications_MC::ScatterType::both);
	}
	else
	{
		if (!m_ring) 
			m_mcspecs->SetScatterType(SKTRAN_Specifications_MC::ScatterType::manualInelastic);
		else 
			m_mcspecs->SetScatterType(SKTRAN_Specifications_MC::ScatterType::manualBoth);

		std::vector<size_t> maxRamanOrders = { 3, 2 };
		ok = ok && m_mcspecs->SetMaxRamanOrders(maxRamanOrders);
		std::vector<double> minFrac = { 0.05, 0.05, 0.05 };
		ok = ok && m_mcspecs->SetMinFractionHigherOrder(minFrac);
	}

	if (m_ring)
		ok = ok && m_mcspecs->SetSecondaryOutput(SKTRAN_Specifications_MC::SecondaryOutput::ringSpectrum);

	std::vector<double> optwavelengths = { 330, 340, 350, 590, 600, 610 };
	ok = ok && m_mcspecs->SetOpticalPropertiesWavelengths(optwavelengths);

	// set solar spectrum
	std::vector<double> solarSpectrum = { 1, 1, 1, 1, 1, 1 };
	ok = ok && m_mcspecs->SetSolarSpectrum(solarSpectrum);

	ok = ok && m_mcspecs->FinalizeSpecs();

	// replace the air species
	skOpticalProperties_RayleighDryAir_Inelastic*	rayleigh;
	skClimatology_MSIS90*							msis90;

	rayleigh = new skOpticalProperties_RayleighDryAir_Inelastic;
	msis90 = new skClimatology_MSIS90;

	ok = (rayleigh != NULL) && (msis90 != NULL);

	ok = ok && m_opticalstate.RemoveSpecies(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3);
	ok = ok && m_opticalstate.AddSpecies(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90, rayleigh);

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Short_Test_Inelastic_MC::MC_Scatter		2020-04-03*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Inelastic_MC::Scatter(nx2dArray<double>& radiance, size_t scattorder)
{
	bool								ok = true;
	bool								ok1;
	SKTRAN_Engine_MC_V21				mc_engine;
	std::vector<SKTRAN_StokesScalar>	mc_radiance_temp;
	size_t numThreads = 1;

	std::vector<skRTStokesVector> svecs;
	ok = radiance.SetSize(m_linesofsight.NumRays(), m_wavelen.size());
	mc_radiance_temp.resize(m_linesofsight.NumRays());
	ok = m_mcspecs->SetSunGeographicPosition(m_sun);
	ok = ok && m_mcspecs->FinalizeSpecs();
	ok = ok && mc_engine.ConfigureModel(*m_mcspecs, m_linesofsight, numThreads);
	ok = ok && radiance.SetSize(m_linesofsight.NumRays(), m_wavelen.size());
	if (ok)
	{
		for (size_t i = 0; i < m_wavelen.size(); i++)
		{
			ok1 = mc_engine.CalculateRadiance(&mc_radiance_temp, m_wavelen[i], scattorder, &m_opticalstate, &svecs);
			ok = ok && ok1;
			if (ok1)
			{
				for (size_t j = 0; j < mc_radiance_temp.size(); j++)
				{
					if (!m_ring) radiance.At(j, i) = mc_radiance_temp[j];
					else ok = ok && mc_engine.GetSecondaryMeasurement(j, radiance.At(j, i));
				}
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Short_Test::MC_Scatter, There were errors executing the MC model. This should be investigated in more detail.");
	}

	return ok;
}



/*---------------------------------------------------------------------------
 *         SKTRAN_Short_Test_Inelastic_MC::LoadHardCodedValues     2020-04-03 */
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_Short_Test_Inelastic_MC::LoadHardCodedValues(nx2dArray<double>& hardcode)
{
	bool ok = true;
	int testcase = (m_ring ? 4 : 0) + (m_optimized ? 2 : 0) + (m_simultaneous ? 1 : 0);
	std::vector<double> mc;
	switch (testcase)
	{
		case 0: // inelastic
			mc  = { 5.83989324623629e-02, 4.12098259733346e-02, 3.64124212670311e-02, 3.00017491987263e-02,
					1.69330026855362e-01, 1.72766794921033e-01, 1.64013299129608e-01, 1.60269179309271e-01 };
			break;
		case 1: // simultaneous
			mc  = { 5.83989324623629e-02, 4.12098259733346e-02, 3.64124212670311e-02, 3.00017491987263e-02,
					1.58678857101968e-01, 1.59190394844242e-01, 1.39370388355339e-01, 1.69689307323443e-01 };
			break;
		case 2:  // optimized
			mc  = { 5.68244276676006e-02, 5.13637379438321e-02, 3.19911268134996e-02, 2.95360491504307e-02,
					1.59156007299124e-01, 1.63058630311870e-01, 1.61819376865440e-01, 1.50074846614705e-01 };
			break;
		case 3: // optimized simultaneous
			mc  = { 5.68244276676006e-02, 5.13637379438321e-02, 3.19911268134996e-02, 2.95360491504307e-02,
					1.61535158610863e-01, 1.64237525384473e-01, 1.38946123097123e-01, 1.35347155516999e-01 };
			break;
		case 4: // ring
			mc  = { 7.32114553518655e-05, 6.57406917755582e-05, 8.72182577022348e-05, 3.74015499189707e-05,
					3.80490720726684e-05, 4.46691118247620e-04, - 4.23168093617869e-07, 1.54040173070888e-04 };
			break;
		case 5: // simultaneous ring
			mc  = { 7.32114553518655e-05, 6.57406917755582e-05, 8.72182577022348e-05, 3.74015499189707e-05,
					- 1.49882548270708e-04, - 1.79780855648580e-06, - 1.19333875387083e-04, - 8.45186108756882e-04 };
			break;
		case 6: // optimized ring
			mc  = { 6.40423107195561e-05, 9.00503515726087e-05, 7.17288311565054e-05, 6.05103163474209e-05,
					- 4.36030617852743e-05, - 3.36399761179754e-05, - 3.85896774051736e-05, - 3.44004270794588e-05 };
			break;
		case 7:// optimized simultaneous ring
			mc  = { 6.40423107195561e-05, 9.00503515726087e-05, 7.17288311565054e-05, 6.05103163474209e-05,
					- 3.56401502763240e-05, - 9.05373136546410e-05, - 3.17115580209637e-05, - 5.91555494448382e-05 };
			break;
	}
	   	 
	hardcode = nx2dArray<double>(4, 2, &mc[0]);

	return ok;
}

const char * SKTRAN_Short_Test_Inelastic_MC::Name() const
{
	int testcase = (m_ring ? 4 : 0) + (m_optimized ? 2 : 0) + (m_simultaneous ? 1 : 0);
	switch (testcase)
	{
	case 0: // inelastic
		return "Inelastic MC";
	case 1: // simultaneous
		return "Simultaneous MC";
	case 2:  // optimized
		return "Optimized MC";
	case 3: // optimized simultaneous
		return "Optimized Simultaneous MC";
	case 4: // ring
		return "Ring MC";
	case 5: // simultaneous ring
		return "Simultaneous Ring MC";
	case 6: // optimized ring
		return "Optimized Ring MC";
	case 7:// optimized simultaneous ring
		return "Optimized Simultaneous Ring MC";
	}
}