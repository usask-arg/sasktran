#include "include/sktran_montecarlo_internals.h"

#include <vector>
#include <limits>
#include <boost/timer/timer.hpp>

/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::SKTRAN_Engine_MC_V21		2011-07-12
 */
/*---------------------------------------------------------------------------*/
SKTRAN_Engine_MC_V21::SKTRAN_Engine_MC_V21()
{ 
	InitializeEngine( );
}

void SKTRAN_Engine_MC_V21::InitializeEngine( ){
	curveAllRays = true;
	curveLOSRays = true;

	m_prevWavelen				= 0.0;
	m_opticalpropertiestable	= NULL;
	m_amfopticalpropertiestable = NULL;
	m_mcconfig					= nullptr;
	m_opticalpropsintegrator	= NULL;
	m_amfopticalpropsintegrator = NULL;
    m_mcOptTable                = NULL;

	m_rayFactory_los			= NULL;
	m_rayFactory_secondary		= NULL;
	m_rayFactory_solar			= NULL;
	m_scatterop                 = NULL;

	m_aveKernel                 = NULL;

    m_emissionTable             = nullptr;
	m_solarTransmissionTable	= NULL;

	m_numPhotons.resize(1);
	m_numPhotons[0]             = 10000;	// Default 1% precision
	m_targetStd                 = 0.01;
	m_minimumRelativeWeight	    = 0.0;
	m_scatterPosResolution		= 50.0;

	m_minFractionHigherOrder = std::vector<double>(1, 0.1);

}


/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::~SKTRAN_Engine_MC_V21		2011-07-12
 */
/*---------------------------------------------------------------------------*/
SKTRAN_Engine_MC_V21::~SKTRAN_Engine_MC_V21()
{
	ReleaseResources();
}

GEODETIC_INSTANT SKTRAN_Engine_MC_V21::ReferencePoint() const
{
	GEODETIC_INSTANT point;

	point.latitude  = m_mcconfig->CoordinateSystemPtr()->ReferencePtLatitude();
	point.longitude = m_mcconfig->CoordinateSystemPtr()->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd		= m_linesofsight.MeanMJD();

	return point;
}


/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::ReleaseResources		2011-07-12
 */
/*---------------------------------------------------------------------------*/
void SKTRAN_Engine_MC_V21::ReleaseResources()
{
    // Release everything else
	if(NULL!=m_opticalpropsintegrator)	m_opticalpropsintegrator	->Release();	m_opticalpropsintegrator	= NULL;
	if(NULL!=m_solarTransmissionTable)	m_solarTransmissionTable	->Release();	m_solarTransmissionTable	= NULL;
    if(nullptr!=m_emissionTable)        m_emissionTable             ->Release();    m_emissionTable             = nullptr;
	if (NULL != m_amfopticalpropertiestable)	m_amfopticalpropertiestable->Release(); m_amfopticalpropertiestable = NULL;
	if (NULL != m_amfopticalpropsintegrator)	m_amfopticalpropsintegrator->Release(); m_amfopticalpropsintegrator = NULL;


    m_scatterop = nullptr; 

	m_mcconfig = nullptr;
	
	// Release internally-stored LOS
	for(size_t losIdx=0; losIdx<m_firstOrderPhotons.size(); losIdx++)
	{
		std::unique_ptr<SKTRAN_RayOptical_Base>	nullpointer;
		m_firstOrderPhotons[losIdx]->SetOpticalRay(std::move(nullpointer));
	}

	m_aveKernel = NULL;

    // delete m_mcOptTable
    m_mcOptTable = nullptr;
    if(nullptr!=m_opticalpropertiestable){
        for(;m_opticalpropertiestable->Release() > 0;);
        m_opticalpropertiestable = nullptr;
    }

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Engine_MC_V21::SetRayFactory_LOS		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Engine_MC_V21::SetRayFactory_LOS( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_los)
{
	m_rayFactory_los = rayFactory_los;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Engine_MC_V21::SetRayFactory_SOLAR		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Engine_MC_V21::SetRayFactory_SOLAR( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_solar)
{
	m_rayFactory_solar = rayFactory_solar;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_Engine_MC_V21::SetRayFactory_SECONDARY		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Engine_MC_V21::SetRayFactory_SECONDARY( std::shared_ptr<SKTRAN_RayFactory_Base>	rayFactory_secondary)
{
	m_rayFactory_secondary = rayFactory_secondary;
	return true;
}


/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::ConfigureModel		2011-07-12*/
/**
 *	Create threads, create random number generators, configure geometry 
 *	for tables used in radiance calculations. 
**/
/*---------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::ConfigureMinorResolutionUpdates( SKTRAN_SpecsUser_Base& userspecs_base )
{
	bool ok = true;
	const SKTRAN_Specifications_MC* mcspecs = dynamic_cast<const SKTRAN_Specifications_MC*>(&userspecs_base);
    mcspecs->GetNumPhotonsPerLOS(m_numPhotons);
	m_targetStd = mcspecs->GetPrecisionMC();

	m_minimumRelativeWeight = mcspecs->GetMinimumRelPathWeight();
	m_scatterPosResolution = mcspecs->GetScatterPositionRes();
	return ok;
}


bool SKTRAN_Engine_MC_V21::ConfigureModel( SKTRAN_SpecsUser_Base& userspecs_base, const SKTRAN_LineOfSightArray_V21& linesofsight, size_t numthreads )
{
	bool		ok = true;
	nxVector    sunDir;

	ReleaseResources();

	SKTRAN_Specifications_MC* mcspecs = dynamic_cast<SKTRAN_Specifications_MC*>(&userspecs_base);
	ok = ok && NULL!=mcspecs;
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, MC engine must be configured using MC specs.");
	if(!ok) nxLog::Record(NXLOG_ERROR," SKTRAN_Engine_MC_V21::ConfigureModel, user specs must be configured to cache solar transmissions for the entire atmosphere -- solar transmission values will be wrong! Use SKTRAN_GridSpecificationsLegacy_MC_V21 or call userspecs.ConfigureForMonteCarlo(true) BEFORE specs are configured.");
	ok = ok && mcspecs->HasBeenFinalized();
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Specs have not been finalized!");

	ok = ok && ConfigureModel_SetThreads( mcspecs, numthreads );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set up threads.");


	// Get ray tracing shells from user
	ok = ok && mcspecs->CreateConfigurationManager( m_mcconfig );
	ok = ok && mcspecs->GetSun(&sunDir);
	ok = ok && m_mcconfig->ConfigureCoordinateTransform( sunDir, linesofsight, mcspecs->GetGroundShiftAlt(), mcspecs->GetTOAHeight(), mcspecs->GetReferencePoint(), mcspecs->GetNadirReferencePointOnGround() );

	ok = ok && ConfigureMinorResolutionUpdates( userspecs_base );

	ok = ok && InitializeRandomNumberGenerators ( mcspecs->GetRngSeed() );
	if(!ok)  nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not initialize random number generators.");
	ok = ok && CreateOpticalPropertyTables      ( mcspecs );
	if(ok) m_mcOptTable->MakeThreadsafeFor( m_internalNumThreads );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create optical properties table.");
	ok = ok && mcspecs->CreateScatterOperator  ( m_scatterop );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create MC scatter operator.");
	ok = ok && m_scatterop->SetOpticalProperties(m_opticalpropertiestable);
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Scatter operator would not accept optical properties table.");
	ok = ok && m_linesofsight.DeepCopy( linesofsight );
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not configure lines of sight, that's a problem!");

    std::unique_ptr< SKTRAN_Sun_Base > sun(nullptr); // Change references to m_sun to sun, give ownership to solarTransmissionTable 
	ok = ok && mcspecs->CreateAveragingKernel(&m_aveKernel);
	if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create averaging kernel.");
	ok = ok && mcspecs->CreateSun( sun, m_randGens, m_internalNumThreads );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create sun object.");
	ok = ok && mcspecs->CreateOpticalPropsIntegrator( &m_opticalpropsintegrator );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not configure optical properties integrator.");
	ok = ok && mcspecs->CreateRayTracers( this );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create ray tracers.");
	ok = ok && mcspecs->CreateSolarTransmissionTable( m_mcconfig->CoordinateSystemObjectVar(), &m_solarTransmissionTable, m_internalNumThreads );
	ok = ok && m_solarTransmissionTable->SetSun(sun);
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not configure solar transmission table.");
	ok = ok && mcspecs->CreateEmissionTable( m_mcconfig->CoordinateSystemObjectVar(), &m_emissionTable );
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not configure emission table.");
	ok = ok && mcspecs->SetRayTracers(this, m_mcconfig->CoordinateSystemObjectVar());
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set ray tracers.");
	
	// Should really get rid of these awful flags
	curveLOSRays = SKTRAN_Specifications_MC::RayTracerType::curved==mcspecs->GetLOSRayTracerType();
	curveAllRays = SKTRAN_Specifications_MC::RayTracerType::curved==mcspecs->GetMSRayTracerType(); 

	ok = ok && m_opticalpropsintegrator->SetOpticalProps(m_opticalpropertiestable);
	if(!ok)	nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21:ConfigureModel, Could not set optical properties table in integrator.");
	ok = ok && m_solarTransmissionTable->ConfigureOptical(m_rayFactory_solar, m_opticalpropsintegrator);
	ok = ok && m_solarTransmissionTable->MakeThreadSafeFor(m_internalNumThreads);
	if(!ok)	nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21:ConfigureModel, Could not configure ray tracer and integrator in solar transmission table.");

	m_minFractionHigherOrder = mcspecs->GetMinFractionHigherOrder ( );

	ok = ok && mcspecs->CreateAirMassFactorCalculator(m_amfcalculator);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create the air mass factor calculator.");
	ok = ok && mcspecs->CreateAirMassFactorOpticalPropertiesTable(m_opticalpropertiestable, &m_amfopticalpropertiestable);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create the air mass factor optical property table.");
	if (m_amfopticalpropertiestable != nullptr)	ok = ok && m_amfopticalpropertiestable->SetCoords(m_mcconfig->CoordinateSystemObjectVar());
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set coordinates for the air mass factor optical property table.");
	ok = ok && mcspecs->CreateAirMassFactorOpticalPropsIntegrator(m_opticalpropsintegrator, m_amfopticalpropertiestable, &m_amfopticalpropsintegrator);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create the air mass factor optical property integrator.");
	ok = ok && mcspecs->SetAirMassFactorRayTracers(m_amfcalculator, m_mcconfig->CoordinateSystemObjectVar());
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set the ray tracers for the air mass factor calculator.");
	ok = ok && m_amfcalculator->SetOpticalPropertiesTable(m_amfopticalpropertiestable);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set the optical property table for the air mass factor calculator.");
	ok = ok && m_amfcalculator->SetOpticalPropertiesIntegrator(m_amfopticalpropsintegrator);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set the optical property integrator for the air mass factor calculator.");
	ok = ok && m_amfcalculator->SetRayTracingShells(mcspecs->GetAMFRayTracingShells()); 
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set the ray tracing shells for the air mass factor calculator.");
	ok = ok && m_amfcalculator->SetCoords(m_mcconfig->CoordinateSystemObjectVar());
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not set coordinates for the air mass factor calculator.");

	ok = ok && mcspecs->ConfigureSimultaneousWavelengths(m_simwl);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not configure the simultaneous wavelength manager.");
	ok = ok && mcspecs->CreateOptimalScatterSequenceManager(m_seqManager);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::ConfigureModel, Could not create the optimal scatter sequence manager.");

	ok = ok && mcspecs->CreatePhotonTemplate(m_photonTemplate);
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Enging_MC_V21::ConfigureModel, Could not create the MCPhoton template.");

	m_chunkSize = mcspecs->GetChunkSize();

	omp_set_num_threads(m_externalNumThreads);		// This is not a safe way to handle this problem

	return ok;
}

bool SKTRAN_Engine_MC_V21::ConfigureModel_SetThreads( const SKTRAN_Specifications_MC* mcspecs, size_t numthreads )
{
	bool ok = true;

	ok = ok && 0==omp_in_parallel();
	if(ok)
	{
		m_externalNumThreads = omp_get_max_threads();

		if(0 == numthreads)
		{
//#if defined(NXDEBUG)
//			omp_set_num_threads(1);
//#else
			omp_set_num_threads(nxmax(omp_get_num_procs(),omp_get_max_threads()));
//#endif
		} else{
				omp_set_num_threads((int)numthreads);
			}
			m_internalNumThreads = omp_get_max_threads();
			if(omp_get_max_threads() < omp_get_num_procs()){
				NXTRACE_ONCEONLY(firsttime,("SKTRAN_Engine_MC_V21::ConfigureModel, user-defined number of threads is less than the number of processors on this machine. Performance may be increased by increasing number of threads (OpenMP detects %d processors).", omp_get_num_procs() ));
			}
			// TODO: not really a windows thing, more of an omp version thing
			#ifdef NX_WINDOWS
				omp_set_nested(0);
			#else 
				omp_set_max_active_levels(1);										     // Can only have one level of parallelism -- some arrays are accessed by thread number
			#endif
			omp_set_dynamic(mcspecs->GetAllowDynamicThreads() ? 1 : 0);	 // True if RTE can adjust number of threads to optimize performance
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_Engine_MC_V21::ConfigureModel, Specs may be called in parallelized region -- thread settings ignored!");
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *		SKTRAN_Engine_MC_V21::CalculateOpticalPropertiesTable	2011-07-12*/
/**
 *	Fill optical property tables. 
**/
/*---------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::CalculateOpticalPropertiesTable( double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology )
{
	bool                                  ok = true;

	GEODETIC_INSTANT                      point;
    const SKTRAN_CoordinateTransform_V2*  coords;

	// Make sure optical property table creates cumulative distribution function
	ok = NULL!=m_mcOptTable && NULL!=m_opticalpropertiestable;
	if( !ok ){
		nxLog::Record(NXLOG_ERROR,"SKTRANSO_Engine::CalculateOpticalPropertiesTable, cannot fill optical property tables, make sure table creates scatter cdf" );
		return false;
	}

	// Get reference point for filling optical property tables
    coords          = m_mcconfig->CoordinateSystemPtr();
	point.latitude  = coords->ReferencePtLatitude();
	point.longitude = coords->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd       = coords->ReferencePointMJD();
	NXASSERT((point.mjd > 10000.0));
	if (point.mjd < 10000.0)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateOpticalPropertiesTable the mjd being used for the climatologies is probably out of range. Its value is %e", (double)point.mjd );
	}

	// Fill optical property tables
	ok  = ok && opticalstate->SetTimeAndLocation( point, userupdateclimatology );	
	ok  = ok && m_opticalpropertiestable->ConfigureOptical( wavelen, *opticalstate );
	//ok = ok && m_inelasticOptTable->ConfigureOptical( wavelen, *opticalstate );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateRadiance, Error calculating the optical properties table. Thats not good");
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *		SKTRAN_Engine_MC_V21::CalculateAMFOpticalPropertiesTable	2018-09-20 */
 /**
  *	Fill amf optical property tables.
 **/
 /*---------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::CalculateAMFOpticalPropertiesTable(double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology)
{
	bool                                  ok = true;

	GEODETIC_INSTANT                      point;
	const SKTRAN_CoordinateTransform_V2*  coords;

	// Get reference point for filling optical property tables
	coords = m_mcconfig->CoordinateSystemPtr();
	point.latitude = coords->ReferencePtLatitude();
	point.longitude = coords->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd = coords->ReferencePointMJD();
	NXASSERT((point.mjd > 10000.0));
	if (point.mjd < 10000.0)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_Engine::CalculateAMFOpticalPropertiesTable the mjd being used for the climatologies is probably out of range. Its value is %e", (double)point.mjd);
	}

	// Fill optical property tables
	ok = ok && opticalstate->SetTimeAndLocation(point, userupdateclimatology);
	ok = ok && m_amfopticalpropertiestable->ConfigureOptical(wavelen, *opticalstate);

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_Engine::CalculateRadiance, Error calculating the AMF optical properties table. Thats not good");
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *		SKTRAN_Engine_MC_V21::TraceMotherRays	2012-07-12*/
/**
 *	Create rays for lines of sight representing observer lines of sight. 
 *	Ray tracing for all photons starts from these rays; their "children"
 *	represent the multiple-scattered photons. Multiple-scatter rays are 
 *	created and destroyed during the calculation; the first-order scatter
 *	rays are cached since #m_numPhotons are scattered along these LOSs.	
**/
/*---------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::TraceMotherRays( ){

	bool								ok = true;
	bool								ok1;

    const SKTRAN_CoordinateTransform_V2* coords = m_mcconfig->CoordinateSystemPtr();
    const SKTRAN_LineOfSightEntry_V2*	 entry;

	std::vector<double> sigmak;
	std::vector<double> sigmaf;

//	if(curveLOSRays) m_raytracer_shells_curved->SetUseCurve(curveLOSRays);
	
	m_firstOrderPhotons.resize(m_linesofsight.NumRays());

	for(size_t losIdx=0; losIdx<m_linesofsight.NumRays(); losIdx++)
	{
		std::unique_ptr<SKTRAN_RayOptical_Base> tempRay;
		std::unique_ptr<SKTRAN_MCPhoton_Base> tempPhoton(m_photonTemplate->Clone());

		//coords = m_mcconfig->CoordinateSystemPtr();

		// Create mother ray
		m_firstOrderPhotons[losIdx] = std::move(tempPhoton);
		entry =       m_linesofsight.Entry(losIdx);
		ok1    =       (nullptr !=entry);
		ok1    = ok1 && m_rayFactory_los->CreateRayObject(&tempRay);
		//ok1    = ok1 && temp->SetWavelength(m_prevWavelen);

		ok1 = ok1 && m_firstOrderPhotons[losIdx]->SetOpticalRay( std::move(tempRay) );
		//ok1 = ok1 && m_firstOrderPhotons[losIdx]->Configure( m_photonTemplate );
		//ok1 = ok1 && m_firstOrderPhotons[losIdx]->SetWavelengths(m_simwl.GetWavelengths());
		ok1 = ok1 && m_firstOrderPhotons[losIdx]->SetCurrentWavelength(m_simwl.GetPrimaryWavelength());
		ok1 = ok1 && m_firstOrderPhotons[losIdx]->photonOptical()->MoveObserver( coords->GeographicToHelio(coords->TranslateGeoidToOsculatingSphere(entry->Observer())), coords->GeographicToHelio(entry->Look()).UnitVector() );
		ok1 = ok1 && m_firstOrderPhotons[losIdx]->DefineRayBasis();
		ok1 = ok1 && m_firstOrderPhotons[losIdx]->TraceRays(m_opticalpropsintegrator, curveLOSRays);
////		ok1 = ok1 && m_rayFactory_los->TraceRay( m_firstOrderPhotons[losIdx].photonOptical() );
//		ok1 = ok1 && m_firstOrderPhotons[losIdx].photonOptical()->TraceRay_NewMethod() ;
//		
//		if(true==curveLOSRays)
//		{
//			sigmak.resize(m_firstOrderPhotons[losIdx].photonOptical()->GetNumQuadraturePoints());
//			sigmaf.resize(m_firstOrderPhotons[losIdx].photonOptical()->GetNumQuadraturePoints());
//			ok1 = ok1 && m_opticalpropsintegrator->CalculateRayScalarTransmissionVector(m_firstOrderPhotons[losIdx].photonOptical(),NULL, false, true, &sigmak, &sigmaf);
//		} 
//		else
//		{
//			ok1 = ok1 && m_opticalpropsintegrator->CalculateRayScalarTransmission_withMinContainer(m_firstOrderPhotons[losIdx].photonOptical(), NULL, false, true);
//		}
		if (!ok1)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_Engine_MC_V21::TraceMotherRays, There were errors tracing motherray  (%d)", (int) losIdx);
		}
		ok = ok && ok1;
	}

	if(!ok) nxLog::Record(NXLOG_ERROR,"SKTRAN_Engine_MC_V21::TraceMotherRays, error creating mother rays");

	return ok;
}

void SKTRAN_Engine_MC_V21::DeleteMotherRays()
{
	for(size_t midx=0; midx<m_firstOrderPhotons.size(); midx++)
	{
		std::unique_ptr<SKTRAN_RayOptical_Base>	nullpointer;
		m_firstOrderPhotons[midx]->SetOpticalRay( std::move(nullpointer));
	}
}


/*-----------------------------------------------------------------------------
 *		SKTRAN_Engine_MC_V21::InitializeRandomNumberGenerators	2011-07-12*/
/**
 *	Random number generators are created and initialized. Since the monte
 *	carlo code is multithreaded and mulitple threads may require random
 *	numbers simultaneously, it's easier to create several random number
 *	generators than to provide mutex protection for a single generator. 
**/
/*---------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::InitializeRandomNumberGenerators( const std::uint32_t seed )
{

	bool ok = true;

	if (!ok){
		nxLog::Record(NXLOG_WARNING," SKTRAN_Engine_MC_V21::InitializerandomNumberGenerators, m_randGens must be cleared before initialization! ");
	}

	// Create array of random number generators
	m_randGens.resize( m_internalNumThreads );
	
	// Initialize one random number generator for each thread
	if(ok){
		for(int i = 0; i < m_internalNumThreads; i++){
			if(0==seed){
				m_randGens[i].SetSeed((const std::uint32_t)(std::time(NULL)+i) );	// Each RNG needs its own engine and needs to be initialized to a different number
			} else{
				m_randGens[i].SetSeed( seed );                                      // Initialize each RNG to the same value, probably for debugging/testing
			}
		}
	}


	if (!ok){
		nxLog::Record(NXLOG_ERROR," SKTRAN_Engine_MC_V21::InitializerandomNumberGenerators, Error creating random number generators. ");
	}

	return ok;
}



/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::CalculateRadiance		2011-07-12*/
/** 
 *	Calculate optical property tables for the specified wavelength/atmospheric
 *	state and perform multiple-scatter monte carlo radiative transfer simulation. 
**/
/*---------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::CalculateRadiance( std::vector<SKTRAN_StokesScalar>  *losradiance, 
											 double                              wavelen,
											 size_t                              numordersofscatter, 
											 SKTRAN_AtmosphericOpticalState_V21* opticalstate,
											 std::vector<skRTStokesVector>*      losvector,
											 bool                                updateclimatology, 
											 SKTRAN_DiagnosticInterface         *diag )
{
	bool					ok = true;
	GEODETIC_INSTANT		pt;

	m_simwl.SetCurrentWavelength(wavelen);

	//SKTRAN_AtmosphericOpticalState_V21 amfopticalstate;
	//skClimatology* amfclimatology;
	//skOpticalProperties* amfopticalproperties;

	if (!m_simwl.Active() || !m_simwl.Complete())
	{
		omp_set_num_threads(m_internalNumThreads);
		// Update maximum number of scatters per photon
		if (0 == numordersofscatter) {
			m_maxNumScatters = 1000;
		}
		else {
			m_maxNumScatters = numordersofscatter;
		}
		m_seqManager->SetMaxOrder(numordersofscatter);

		// Update optical properties caches if optical state has changed since last calculation
	//	ok = ok && (wavelen==m_prevWavelen && !m_optStateAltered);	// Detection of change of atmospheric state should be automatic
		ok = ok && CalculateOpticalPropertiesTable(wavelen, opticalstate, true);	// Optical properties need to be cached for new wavelength
		if (ok) {
			m_prevWavelen = wavelen;
		}
		else {
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21:CalculateRadiance, Could not calculate optical properties table.");
		}

		// calculate the AMF optical properties (contains only the target AMF species)
		ok = ok && m_amfcalculator->CalculateOpticalPropertiesTable(wavelen, opticalstate, true);

		// The solar transmission table is configured, now fill it
		ok = ok && m_solarTransmissionTable->FillTable();
		m_scatterop->AddSourceTerm(m_solarTransmissionTable);
		ok = ok && m_scatterop->SetCoordinateSystem(m_mcconfig->CoordinateSystemObjectVar());
		ok = ok && m_emissionTable->ConfigureOptical(wavelen, opticalstate->EmissionObjectVar(), m_opticalpropertiestable);
		if (ok) {
			m_scatterop->AddSourceTerm(m_emissionTable);
		}
		else {
			nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21:CalculateRadiance, Could not configure solar and/or emission tables.");
		}

		//if(NULL!=m_ThermalEmission){
		//	pt = opticalstate->GetTimeAndLocation();
		//	ok = ok && m_ThermalEmission->ConfigureOptical(wavelen, 4.8e14, opticalstate, pt, m_mcconfig->CoordinateSystemObjectVar() );
		//	m_scatterop->AddSourceTerm( m_ThermalEmission );
		//}

		// Trace the rays representing the satellite's lines of sight
		ok = ok && TraceMotherRays();
		if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::CalculateRadiance, Could not calculate transmission along satellite lines of sight.");

		// Create internal array to hold radiances output by MC simulation
		m_radiances.SetSize(m_simwl.GetNumWavelengths(), m_firstOrderPhotons.size());
		m_radiances.SetTo(0.0);
		m_svecs.SetSize(m_simwl.GetNumWavelengths(), m_firstOrderPhotons.size());
		m_variances.SetSize(m_simwl.GetNumWavelengths(), m_firstOrderPhotons.size());
		m_variances.SetTo(0.0);

		if (m_seqManager->SecondaryMeasurement())
		{
			m_secondaryMeasurements.SetSize(m_simwl.GetNumWavelengths(), m_firstOrderPhotons.size());
			m_secondaryMeasurements.SetTo(0.0);
			m_secondaryMeasurementVariances.SetSize(m_simwl.GetNumWavelengths(), m_firstOrderPhotons.size());
			m_secondaryMeasurementVariances.SetTo(0.0);
		}

		m_airMassFactors.SetSize(m_firstOrderPhotons.size(), m_amfcalculator->NumAMFCells());
		m_airMassFactorVariances.SetSize(m_firstOrderPhotons.size(), m_amfcalculator->NumAMFCells());
		m_airMassFactors.SetTo(0.0);
		m_airMassFactorVariances.SetTo(0.0);

		m_aveKernel->ConfigureKernel(0.0, 200, 500.0, -1.0, 201, 0.01, m_internalNumThreads);

		m_lastRunTiming.SetSize(m_firstOrderPhotons.size(), 3);

		// Perform MC simulation on each LOS
		for (size_t losIdx = 0; losIdx < m_firstOrderPhotons.size(); losIdx++) {
			boost::timer::cpu_timer t;
			t.start();
			ok = ok && MonteCarloMultiScatter(losIdx);
			auto tdata = t.elapsed();
			m_lastRunTiming.At(losIdx, 0) = (double)(tdata.system);
			m_lastRunTiming.At(losIdx, 1) = (double)(tdata.user);
			m_lastRunTiming.At(losIdx, 2) = (double)(tdata.wall);
		}



		m_scatterop->ClearSourceTerms(); // Really need a better time to do this

		DeleteMotherRays();

		// Return to external number of threads in case another sasktran engine wants to run using omp multithreading 
		omp_set_num_threads(m_externalNumThreads);

		m_simwl.SetComplete(true);
	}

	// Pass simulation results to user. 
	losradiance->resize(m_firstOrderPhotons.size());
	size_t wlIdx = m_simwl.GetCurrentWavelengthIndex();
	if (nullptr != losvector)	losvector->resize(m_svecs.size());
	for (size_t losIdx = 0; losIdx < (size_t)m_firstOrderPhotons.size(); losIdx++) {
		losradiance->at(losIdx) = m_radiances.At(wlIdx, losIdx);
		if (nullptr != losvector) losvector->at(losIdx) = m_svecs.At(wlIdx, losIdx);
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::MonteCarloMultiScatter		2011-07-12*/
/** 
 *	Environment in which multiple-scatter simulation is performed. For the 
 *	specified line of sight, #m_numPhotons photon paths are simulated. For
 *	each path, first a scatter location p0 is chosen along the observer
 *	line of sight for a once-scattered photon. The radiance due to that once-
 *	scattered photon is calculated. It is then assumed that p0 is the second
 *	scatter location for a twice-scattered photon; an incoming direction for
 *	the twice-scattered photon is calculated and its first scatter
 *	location p1 is chosen along that direction. The radiance at the observer
 *	for a photon that scatters at p1 and then at p0 is then calculated and
 *	added to the total radiance due to photons travelling along "this" path.
 *	This continues until
 *		- there are #m_maxNumScatters scatter locations, or
 *		- the contribution to total radiance along this path by each higher
 *		  order of scattering is guaranteed to be less than #m_minimumRelativeWeight
 *	The total radiance along the observer line of sight is taken to be the average
 *	radiance along all #m_numPhotons paths.
 *	
 *	Monte Carlo integration applies statistics to evaluate an integral discretely
 *	over some space; it is useful when the space is especially large (so
 *	that computation of the integrand is computationally intensive, as in 
 *	the general case of evaluation of the RTE), or when construction of an a priori 
 *	set of discrete points at which to evaluate the integrand is difficult (as in
 *	cases where the gradients of atmospheric optical properties are strong functions
 *	of space or wavelength). This radiative transfer model uses statistics to determine
 *	the points p{n} which are then used to determine limb radiance. That is,
 *		- The integrand is the radiative transfer equation
 *		- The space is the set of all possible photon paths
 *		- The integral gives limb radiance for some atmospheric state, observer, and look direction
 *
**/
/*-----------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::MonteCarloMultiScatter(size_t losIdx){

	bool					ok = true;
	const SKTRAN_MCPhoton_Base * const motherPhoton = m_firstOrderPhotons[losIdx].get();
	const SKTRAN_RayOptical_Base * const mother = motherPhoton->photonOptical();
	
	// Each thread saves its results into a temporary array
	const size_t numthreads = omp_get_max_threads();

	m_higherOrderPhotons.resize(numthreads);
	for(size_t iray=0; iray < m_higherOrderPhotons.size(); iray++)
	{
		std::unique_ptr<SKTRAN_MCPhoton_Base> tempPhoton(m_photonTemplate->Clone());
		std::unique_ptr<SKTRAN_RayOptical_Base> tempRay;
		m_higherOrderPhotons[iray] = std::move(tempPhoton);
		ok = ok && m_rayFactory_secondary->CreateRayObject( &tempRay );
		ok = ok && m_higherOrderPhotons[iray]->SetOpticalRay( std::move(tempRay) );
		//ok = ok && m_higherOrderPhotons[iray]->SetWavelengths(m_simwl.GetWavelengths());
		ok = ok && m_higherOrderPhotons[iray]->SetCurrentWavelength(m_simwl.GetPrimaryWavelength());
	}

	if(curveAllRays){
		m_sigmaks.resize( numthreads );
		m_sigmafs.resize( numthreads );
	}

	m_aveKernel->WipeKernel();
	m_aveKernel->ConfigureObserverGeometry( *mother, mother->Coordinates()->AltitudeToRadius(0.0) );
	
	// Simulate user-specified number of photons 
	const int numPhotonsInt = (int)m_numPhotons[std::min(losIdx,m_numPhotons.size()-1)];	    // omp needs this number to be an integer
	//const int chunkSize( (int) ComputeChunkSize(losIdx));  // omp needs this to be constant
    //const int chunkSize = 500;
	const int minNumRaysPerThread = min( 300, (int)(numPhotonsInt/m_internalNumThreads ) );

	SKTRAN_MCScatterOperatorContainer scatterOpContainer(m_scatterop);

	SKTRAN_MCThreadRadianceLogger     radianceLogger; 
	radianceLogger.ConfigureWavelengths(m_simwl.GetWavelengths(), m_simwl.GetCurrentWavelength(), m_simwl.GetPrimaryWavelengthIndex());
	radianceLogger.ConfigureOptimalScatterSequence(m_seqManager.get());
	radianceLogger.SetChunkSize(m_chunkSize);
	SKTRAN_MCVarianceLogger           varianceLogger; varianceLogger.SetNumThreads(m_internalNumThreads);

	size_t minNumRays = min( 300, numPhotonsInt/m_internalNumThreads );

	double debugDistance = 0.0;
	//debugDistance = -1.0 * mother->GeometryRay()->GetObserver().Magnitude() * ( mother->GeometryRay()->GetObserver().UnitVector() & mother->GeometryRay()->LookVector() ); // tobs
	//debugDistance = 429756.879357503; // (One) quadrature point distance where HR and MC disagree for Elash tests
	debugDistance = 2808220.9111485; // distance for shortTestSuite debugging
	srdLog.SetDebugDistance( false ? debugDistance : -1.0 );

	radianceLogger.SetMinFractionHigherOrder(m_minFractionHigherOrder);
	radianceLogger.SetMaxOrder(m_maxNumScatters);

	// configure AMF resources
	SKTRAN_MCAirMassFactorLogger amfLogger;
	m_amfcalculator->InitializeLogger(&amfLogger);
	m_amfcalculator->AllocatePhotons(m_higherOrderPhotons);
	m_amfcalculator->AllocateRayOptical(m_internalNumThreads);
	m_amfcalculator->TraceMotherRay(mother);

    //// Don't uncomment this one //#pragma omp parallel for schedule(guided, chunkSize) firstprivate(scatterOpContainer, radianceLogger) reduction(&&:ok) // VS15 doesn't like #chunkSize when it comes from call to ComputeChunkSize
    #pragma omp parallel for firstprivate(scatterOpContainer, radianceLogger, amfLogger) reduction(&&:ok)
	for(int photonCount=0; photonCount<numPhotonsInt; photonCount++){
		double measEst = radianceLogger.TargetMeasurement();
		double userStd = measEst * m_targetStd;
		double stdEstimate = std::sqrt( varianceLogger.GetTotalVariance() );

        // Rarely we can get nans in the std estimate, in these cases set it high so we don't stop
        if(stdEstimate != stdEstimate) {
            stdEstimate = 1e99;
        }

		if ( !radianceLogger.RunningSums().minSamplesComplete || userStd < stdEstimate ){
			bool ok1 = true;

			// Initialize photon path parameters for scatter simulation
			const size_t       threadid      = omp_get_thread_num();
			int                order         = 1;
			SKTRAN_MCPhoton_Base*   mcphoton = m_higherOrderPhotons[threadid].get();
			SKTRAN_MCScatterOperator_Base* scatterop = scatterOpContainer.Op();

			scatterop->DeclareNewRay();

			mcphoton->m_isGroundScatter = false;
			//mcphoton->m_scatterWeight   = 1.0;
			std::fill(mcphoton->ScatterWeights().begin(), mcphoton->ScatterWeights().end(), 1.0);

		    mcphoton->photonOptical()->MoveObserver(mother->GetObserver(), mother->LookVector());
			//mcphoton->photonOptical()->SetWavelength(m_prevWavelen);
		    mcphoton->ResetRadiance();
		    mcphoton->DefineRayBasis();

			//std::fill(mcphoton->m_solarSlantColumns.begin(), mcphoton->m_solarSlantColumns.end(), 0.0);
			//std::fill(mcphoton->m_scatterSlantColumns.begin(), mcphoton->m_scatterSlantColumns.end(), 0.0);
			m_amfcalculator->ClearPhoton(mcphoton);

			double r = m_randGens[threadid]();
		    //const size_t maxOrderScatter = radianceLogger.OptimalMaxOrderScatter( r ); //, m_maxNumScatters );
			size_t scatterSequenceIndex = 0, maxOrder = 0;
			//ok = ok && m_seqManager->OptimalScatterSequenceIndex( radianceLogger.RunningSums(), r, scatterSequenceIndex, radianceLogger.RunningSums().minSamplesComplete);
			//ok = ok && m_seqManager->Order(scatterSequenceIndex, maxOrder);
			ok = ok && radianceLogger.OptimalScatterSequenceIndex(r, scatterSequenceIndex, maxOrder);
			//printf("optimal scatter index: %d\n", (int)maxOrderScatter);
	//		while(order <= m_maxNumScatters ){
			//while (ok && order <= (int)maxOrderScatter && (mcphoton->photonRadiance().GetScalar()*m_minimumRela tiveWeight) < mcphoton->m_scatterWeight) {
			while(ok && ok1 && order <= (int)maxOrder && (mcphoton->photonRadiance().GetScalar()*m_minimumRelativeWeight)< mcphoton->ScatterWeight() ){
				ok = ok && scatterop->InputScatterOrder(order);
				ok1 = ok1 && MultipleScatterContribution(scatterop, motherPhoton, mcphoton, mcphoton->m_scatterVector, scatterSequenceIndex, order, threadid);
				//ok = ok && m_aveKernel->AddToKernel(mcphoton->m_scatterVector, mcphoton->photonRadiance().GetRecentContribVec(), order-1, mcphoton->m_isGroundScatter, threadid);
				if( ok && ok1 ) ok = ok && m_aveKernel->AddToKernel( mcphoton, order-1, threadid );
				if( ok && ok1 ) radianceLogger.Submit( scatterSequenceIndex, order-1, mcphoton );
				if( ok && ok1 ) amfLogger.Submit( order-1, mcphoton );
			}

			if ( ok1 ) 
			{
				ok = ok && radianceLogger.DeclareRayDone( scatterSequenceIndex, varianceLogger, threadid );
				if ( ok ) amfLogger.DeclareRayDone( );
			}
			else // if something went wrong, discard the samples instead of crashing the whole calculation (sometimes some funny business at the top of atmosphere causes the program to crash, and it's hard to reproduce)
			{
				ok = ok && radianceLogger.DiscardRay( );
				if ( ok ) amfLogger.DiscardRay( );
				nxLog::Record(NXLOG_WARNING, "SKTRAN_Engine_MC_V21::MonteCarloMultiScatter, photon number %d was discarded due to an error at scatter order %d. If this happens many times there could be a problem.", photonCount, order);
			}
			//varianceLogger.UpdateThreadVariance(radianceLogger.Variance(), radianceLogger.n(1), threadid);
	
		}
	}

	radianceLogger.ExportStatistics(losIdx, m_simwl.GetCurrentWavelength());
	
	size_t wlIdx = m_simwl.GetCurrentWavelengthIndex();
	for (size_t wlIdx = 0; wlIdx < m_simwl.GetNumWavelengths(); wlIdx++)
	{
		m_variances.At(wlIdx, losIdx) = radianceLogger.Variance(wlIdx);
		m_radiances.At(wlIdx, losIdx) = radianceLogger.TotalMeasurement(wlIdx).I();
		SKTRAN_Stokes_NC svec_temp(radianceLogger.TotalMeasurement(wlIdx));
		m_svecs.At(wlIdx, losIdx).SetTo(svec_temp.I(), svec_temp.Q(), svec_temp.U(), svec_temp.V());
	}

	if (radianceLogger.SecondaryMeasurementExists())
	{
		for (size_t wlIdx = 0; wlIdx < m_simwl.GetNumWavelengths(); wlIdx++)
		{
			m_secondaryMeasurements.At(wlIdx, losIdx) = radianceLogger.SecondaryMeasurement(wlIdx);
			m_secondaryMeasurementVariances.At(wlIdx, losIdx) = radianceLogger.SecondaryVariance(wlIdx);
		}
	}
	
	for (size_t amfidx = 0; amfidx < m_airMassFactors.YSize(); amfidx++)
	{
		m_airMassFactors.At(losIdx, amfidx) = amfLogger.AirMassFactor(amfidx);
		m_airMassFactorVariances.At(losIdx, amfidx) = amfLogger.AirMassFactorVariance(amfidx);
	}
	//if (m_calcAMF) {
	//	std::vector<double> amf = amfLogger.AirMassFactor();
	//	std::vector<double> amfvar = amfLogger.AirMassFactorVariance();

	//	for (size_t amfidx = 0; amfidx < m_shellGrid_amf->NumCells(); ++amfidx) {
	//		m_airMassFactors.At(losIdx, amfidx) = amf[amfidx];
	//		m_airMassFactorVariances.At(losIdx, amfidx) = amfvar[amfidx];
	//	}
	//}

	//printf("\n\n");
	//printf("mc(%2u,:,:) = [...\n", losIdx+1 );
	//for(size_t oidx=1; oidx<=radianceLogger.GetNumDistinctOrders(); ++oidx){
	//	printf("% 8.5e   % 8.5e   % 8.5e \n", radianceLogger.Measurement(oidx).At(1), radianceLogger.Measurement(oidx).At(2), radianceLogger.Measurement(oidx).At(3) );
	//}
	//printf("];");
	//printf("\n");


	for(size_t iray=0; iray<m_higherOrderPhotons.size(); iray++)
	{
		std::unique_ptr<SKTRAN_RayOptical_Base>	nullpointer;
		m_higherOrderPhotons[iray]->SetOpticalRay( std::move(nullpointer) );
	}

	std::stringstream kernelStream;
	kernelStream << "C:/Users/srd740/Documents/MATLAB/elashFix/alongLosMC/photons/sts_sp_los";
	if(0== ((int)(losIdx/100))) kernelStream << "0";
	if(0== ((int)(losIdx/ 10))) kernelStream << "0";
	kernelStream << losIdx;
	//<< std::sqrt(std::pow(mother->GeometryRay()->GetObserver().Magnitude(),2) - std::pow(mother->GeometryRay()->GetObserver() & mother->GeometryRay()->LookVector(), 2) );
	ok = ok && m_aveKernel->PrintKernel( kernelStream.str() );

	return ok;
}



/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::MultipleScatterContribution		2012-07-12*/
/** 
 *	For a photon travelling backwards along #scatterRay and undergoing
 *	#orderOfScatter-order scattering at #scatterVector, return the radiance 
 *	at some observer due to the next order of scattering taking into account 
 *	that the photon is attenuated by factor #scatterWeight on subsequent scatters
 *	on its path to the observer. Update input parameters to reflect the 
 *	geometry of this next order of scatter. 
 *	Input #groundScatter is true iff the next order of scatter occurs at the ground
 *	(i.e. iff input #scatterVector is a ground point).
 *
**/
/*-----------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::MultipleScatterContribution( SKTRAN_MCScatterOperator_Base* scatterop, SKTRAN_MCPhoton_Base const* const motherPhoton, SKTRAN_MCPhoton_Base* mcphoton, HELIODETIC_VECTOR& prevScatterVector, size_t scatterSequenceIndex, int& orderOfScatter, size_t threadid)
{

	bool							ok = true;
	double							r = 0.0;

	HELIODETIC_POINT				scatterPoint;
	HELIODETIC_POINT				prevScatterPoint;
	HELIODETIC_UNITVECTOR			sun;
	std::vector<double>				solarSlantColumns;


	const SKTRAN_CoordinateTransform_V2*	coords		= motherPhoton->photonOptical()->Coordinates();//m_secondaryRays[threadid]->GeometryRay()->Coordinates();
	//SKTRAN_RayOptical_Base const*			scatterRay	= NULL;
	SKTRAN_MCPhoton_Base const*					incomingPhoton = NULL;



	// Set #scatterRay as the LOS ray for this order of scatter
	if(1==orderOfScatter){
		incomingPhoton = motherPhoton;
	} else{	
		mcphoton->ResetFactors();
		ok = ok && coords->HelioVectorToHelioPoint(prevScatterVector, &prevScatterPoint);
		ok = ok && ChangeScatterDirection( scatterop, prevScatterPoint, mcphoton, threadid, orderOfScatter );
		incomingPhoton = mcphoton;
		
	}
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::MultipleScatterContribution, Could not create scatter ray.");

	ok = ok && m_seqManager->ElasticScatter(scatterSequenceIndex, orderOfScatter, mcphoton->m_elasticScatter); // choose elastic/inelastic for the source term AND the next scatter (optimized mode only)
	ok = ok && ChooseScatterPoint( incomingPhoton, mcphoton, scatterPoint, orderOfScatter );

	ok = ok && mcphoton->PreGenerateRandNum(m_randGens[threadid]); // pre-generate the random number for wavelength selection for the source term AND the next scatter
	if (ok) mcphoton->ManualScatterFactor() = 1.0; // reset the manual selection factor (it should be set during the source term collection, if needed)

	// Gather contributions from all source terms
	//if( 0.0 < mcphoton->m_scatterWeight ){
	if( 0.0 < mcphoton->ScatterWeight() ){

		ok = ok && scatterop->CollectSourceTerms(scatterPoint, mcphoton, *coords, threadid);
		if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::MultipleScatterContribution, Could not evaluate source term.");
	} else{
        SKTRAN_Stokes_NC zeroVector; 
		zeroVector.SetTo(0.0);
        for (auto&& pr : mcphoton->photonRadiances()) pr.AddToVector( zeroVector ); // No contribution from this scatter
	} 

	// calculate slant columns along the new scatter segment and along the corresponding solar ray
	ok = ok && m_amfcalculator->CalculateSlantContribution(orderOfScatter, incomingPhoton->photonOptical(), scatterPoint, mcphoton, threadid);

	ok = ok && scatterop->AcceptScatterPoint( mcphoton );

	// update order of scatter
	++orderOfScatter;

	return ok;

}


/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::ChangeScatterDirection		2014-04-24*/
/** 
 *  Scatter
 *  SHOULD BE MADE CONST -- remove m_sigmaks stuff, etc.
**/
/*-----------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::ChangeScatterDirection( const SKTRAN_MCScatterOperator_Base* scatterop,  const HELIODETIC_POINT& prevScatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid, size_t orderOfScatter ) 
{
	bool ok = true;

	ok = ok && scatterop->RandomScatter( prevScatterPoint, mcphoton, m_randGens[threadid], (int)orderOfScatter );
	
//	ok = ok && m_rayFactory_secondary->TraceRay( mcphoton->photonOptical() );
	//ok = ok && mcphoton->photonOptical()->TraceRay_NewMethod( );
	ok = ok && mcphoton->TraceRays(m_opticalpropsintegrator, curveAllRays);
	
	//if(curveAllRays){
	//	if(m_sigmaks[threadid].size() < mcphoton->photonOptical()->GetNumQuadraturePoints()){
	//		m_sigmaks[threadid].resize( mcphoton->photonOptical()->GetNumQuadraturePoints()) ;
	//		m_sigmafs[threadid].resize( mcphoton->photonOptical()->GetNumQuadraturePoints()) ;
	//	}
	//	ok = ok && m_opticalpropsintegrator->CalculateRayScalarTransmissionVector( mcphoton->photonOptical(),NULL,false,true, &m_sigmaks[threadid], &m_sigmafs[threadid] );
	//} else{
	//	ok = ok && m_opticalpropsintegrator->CalculateRayScalarTransmission_withMinContainer( mcphoton->photonOptical(), NULL, false, true );
	//}
	return ok;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::ChooseScatterPoint		2014-04-24*/
/** 
 *  Choose a point along #scatterRay. the point is #scatterPoint; the effect
 *  of the scatter is multiplied onto #mcphoton, which is updated to reflect
 *  that the scatter happens there. 
**/
/*-----------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::ChooseScatterPoint( SKTRAN_MCPhoton_Base const* incomingPhoton, SKTRAN_MCPhoton_Base* mcphoton, HELIODETIC_POINT& scatterPoint, size_t order ) const
{
	bool ok = true;
#define SKTRAN_ERRORVAL 9999.0;
	double r, forcedScatterCorrFactor;
	bool forcedScatter;
	std::vector<double> slantLengths, slantColumns;
	HELIODETIC_POINT obs;
	const SKTRAN_RayOptical_Base* scatterRay = incomingPhoton->photonOptical();


	// Correct for forced scatter only if it is possible for photon to travel entire LOS and exit atmosphere OR if we are forcing an inelastic scatter
	forcedScatter = !scatterRay->Storage()->GroundIsHit() || (mcphoton->m_manualScatter && !mcphoton->m_elasticScatter);
	if (forcedScatter)
	{
		forcedScatterCorrFactor = 1.0 - exp(-scatterRay->TotalOpticalDepth());	// Correction factor is required; note \tau is in meters
	}
	else
	{
		forcedScatterCorrFactor = 1.0;		// Scatter _will_ happen, either in atmo or from ground; correction factor built into algorithm below
	}

	// Find the next scattering point for the photon being traced
	r = m_randGens[omp_get_thread_num()]();
	if(	!forcedScatter && ( -std::log(r) > scatterRay->TotalOpticalDepth())) {
		// scatter is from ground
		scatterRay->Storage()->LocationOfPoint( scatterRay->Storage()->NumCells(), &scatterPoint);
		mcphoton->m_scatterVector.SetCoords( scatterPoint.UnitVector(), scatterPoint.Radius() );
		//ok = ok && m_opticalpropertiestable->Get_AlbedoForDeprecatedLegacyCode(scatterPoint, &mcphoton->m_albedo);
		
		ok = ok && mcphoton->CalculateTransmissionsGroundScatter(incomingPhoton, scatterPoint);

		auto wavelength = mcphoton->CurrentWavelengths().cbegin();
		auto eAlbedo = mcphoton->Albedos(true).begin();
		for (auto&& albedo : mcphoton->Albedos(false))
		{
			m_opticalpropertiestable->GetBRDF(*(wavelength++), scatterPoint, 0., 0., 0., &albedo);
			albedo *= nxmath::Pi; // assumes Lambertian
			*(eAlbedo++) = albedo;
		}
		
		mcphoton->m_isGroundScatter = true;
		//mcphoton->m_distanceProb = exp(-scatterRay->TotalOpticalDepth());
		mcphoton->m_distanceProb = mcphoton->Transmission();
		mcphoton->m_targetTau = SKTRAN_ERRORVAL;
	} else if(1<scatterRay->GetNumQuadraturePoints() && 1e-6<scatterRay->Storage()->DistanceOfPointFromOrigin( scatterRay->GetNumQuadraturePoints()-1 ) ){
		// scatter is from atmosphere
		r = m_randGens[omp_get_thread_num()]();
		ok = ok && m_opticalpropsintegrator->FindNextScatterPosition(r, m_scatterPosResolution, m_opticalpropertiestable, scatterRay, mcphoton->m_scatterVector, mcphoton->m_distanceProb, mcphoton->m_targetTau, 1==order ? srdLog.GetDebugDistance() : -1.0 );		// Find how far along scatter was
		ok = ok && scatterRay->Coordinates()->HelioVectorToHelioPoint(mcphoton->m_scatterVector, &scatterPoint);
		ok = ok && mcphoton->CalculateTransmissionsAtmoScatter(m_opticalpropertiestable, incomingPhoton, scatterPoint);
		ok = ok && mcphoton->CalculateAlbedo(m_opticalpropertiestable, m_mcOptTable, scatterPoint);
	} else{
		// Condition where scatter is at exactly the edge of the atmosphere looking out
		//mcphoton->m_scatterWeight = 0.0;	
		//mcphoton->m_albedo = 0.0;
		std::fill(mcphoton->ScatterWeights().begin(), mcphoton->ScatterWeights().end(), 0.0);
		std::fill(mcphoton->Albedos().begin(), mcphoton->Albedos().end(), 0.0);

		scatterRay->Coordinates()->HelioVectorToHelioPoint( scatterRay->GetObserver(), &scatterPoint );
		mcphoton->m_scatterVector = scatterPoint.Vector();
		mcphoton->m_isGroundScatter = false;
		mcphoton->m_distanceProb = 0.0;
		mcphoton->m_targetTau = 0.0;
	}

	// Update scattered photon weight
	ok = ok && mcphoton->UpdateScatterWeight(forcedScatterCorrFactor);

	//mcphoton->m_weightFactor = forcedScatterCorrFactor * mcphoton->m_albedo;
	//mcphoton->m_scatterWeight *= mcphoton->m_weightFactor;

	//obs.FromVector(scatterRay->GetObserver(), Coordinates());
	//if (m_calcAMF) {
	//	//ok = ok && AMFCellSlantColumnOneSegment(obs, scatterRay->LookVector(), (scatterRay->GetObserver() - scatterPoint.Vector()).Magnitude(), slantColumns);
	//	ok = ok && m_amfintegrator->AMFCellSlantColumnOneSegment(obs, scatterRay->LookVector(), (scatterRay->GetObserver() - scatterPoint.Vector()).Magnitude(), slantColumns);
	//	
	//	for (size_t amfidx = 0; amfidx < m_shellGrid_amf->NumCells(); amfidx++) {
	//		mcphoton->m_slantColumns[amfidx] += slantColumns[amfidx];
	//	}
	//}

	return ok;
}

/*-----------------------------------------------------------------------------
 *			SKTRAN_Engine_MC_V21::CreateOpticalPropertyTables	2011-07-12*/
 /**
  *	Create tables to cache extinction values, phase functions, and scatter CDFs.
 **/
 /*-----------------------------------------------------------------------------*/
bool SKTRAN_Engine_MC_V21::CreateOpticalPropertyTables(SKTRAN_Specifications_MC* mcspecs)
{
	bool ok = true;

	//if(nullptr!=m_opticalpropertiestable) m_opticalpropertiestable->Release();
	ok = ok && mcspecs->CreateOpticalPropertyTables(&m_opticalpropertiestable, &m_mcOptTable, *m_mcconfig->CoordinateSystemPtr());
	ok = ok &m_opticalpropertiestable->SetCoords(m_mcconfig->CoordinateSystemObjectVar());

	return ok;
}


bool SKTRAN_Engine_MC_V21::GetStokesVectors(nx1dArray<skRTStokesVector>& dest) const
{
	bool ok = true;

	size_t wlidx = m_simwl.GetCurrentWavelengthIndex();
	dest.SetSize(m_svecs.YSize());
	for (size_t svidx = 0; svidx < dest.size(); ++svidx) dest.At(svidx) = m_svecs.At(wlidx, svidx);

	return ok;
}


bool SKTRAN_Engine_MC_V21::GetMeasurementVariance(size_t losidx, double& variance) const
{
	bool ok = true;

	ok = ok && losidx < m_variances.size();
	size_t wlidx = m_simwl.GetCurrentWavelengthIndex();
	if (ok) {
		variance = m_variances.At(wlidx, losidx);
	}
	else {
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetMeasurementVariance, losidx=%i is out of range.", losidx);
	}

	return ok;
}

bool SKTRAN_Engine_MC_V21::GetAirMassFactors(size_t losidx, std::vector<double>& amf) const
{
	bool ok = true;

	ok = ok && losidx < m_airMassFactors.XSize();
	amf.resize(m_amfcalculator->NumAMFCells(), 0.0);
	
	if (ok) {
		for (size_t amfidx = 0; amfidx < m_amfcalculator->NumAMFCells(); ++amfidx) {
			ok = ok && amfidx < m_airMassFactors.YSize();
			if (ok) amf[amfidx] = m_airMassFactors.At(losidx, amfidx);
			else nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetAirMassFActors, amfidx=%i is out of range.", amfidx);
		}
	}
	else {
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetAirMassFactors, losidx=%i is out of range.", losidx);
	}

	return ok;
}

bool SKTRAN_Engine_MC_V21::GetAirMassFactorVariance(size_t losidx, std::vector<double>& amfvar) const
{
	bool ok = true;

	ok = ok && losidx < m_airMassFactorVariances.XSize();
	amfvar.resize(m_amfcalculator->NumAMFCells(), 0.0);

	if (ok) {
		for (size_t amfidx = 0; amfidx < m_amfcalculator->NumAMFCells(); ++amfidx) {
			ok = ok && amfidx < m_airMassFactorVariances.YSize();
			if (ok) amfvar[amfidx] = m_airMassFactorVariances.At(losidx, amfidx);
			else nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetAirMassFactorVariances, amfidx=%i is out of range.", amfidx);
		}
	}
	else {
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetAirMassFactorVariances, losidx=%i is out of range.", losidx);
	}

	return ok;
}

bool SKTRAN_Engine_MC_V21::GetSecondaryMeasurement(size_t losidx, double & measurement) const
{
	bool ok = true;

	ok = ok && losidx < m_secondaryMeasurements.size();
	if (ok) {
		size_t wlidx = m_simwl.GetCurrentWavelengthIndex();
		measurement = m_secondaryMeasurements.At(wlidx, losidx);
	}
	else {
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetSecondaryMeasurement, losidx=%i is out of range.", losidx);
	}

	return ok;
}

bool SKTRAN_Engine_MC_V21::GetSecondaryVariance(size_t losidx, double & variance) const
{
	bool ok = true;

	ok = ok && losidx < m_secondaryMeasurementVariances.size();
	if (ok) {
		size_t wlidx = m_simwl.GetCurrentWavelengthIndex();
		variance = m_secondaryMeasurementVariances.At(wlidx, losidx);
	}
	else {
		nxLog::Record(NXLOG_ERROR, "SKTRAN_Engine_MC_V21::GetSecondaryVariance, losidx=%i is out of range.", losidx);
	}

	return ok;
}


bool SKTRAN_Engine_MC_V21::GetTimingData(nx2dArray<double>& t) const
{
	bool ok = true;

	t = m_lastRunTiming;

	return ok;
}