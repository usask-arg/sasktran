#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables::SKTRAN_EngineDiffuseTables		2010-4-5*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_EngineDiffuseTables::SKTRAN_EngineDiffuseTables()
{
	m_diffusetable			  = nullptr;
	m_solartransmissiontable  = nullptr;
	m_diffusegroundpointtable = nullptr;
	m_emissiontable           = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables::~SKTRAN_EngineDiffuseTables		2010-4-5*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_EngineDiffuseTables::~SKTRAN_EngineDiffuseTables()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables::ReleaseResources		2010-4-5*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_EngineDiffuseTables::ReleaseResources()
{
	if ( m_diffusetable             != nullptr ) m_diffusetable->Release();
	if ( m_solartransmissiontable   != nullptr ) m_solartransmissiontable->Release();
	if ( m_diffusegroundpointtable  != nullptr ) m_diffusegroundpointtable->Release();
	if (m_emissiontable             != nullptr ) m_emissiontable->Release();
	m_diffusetable			  = nullptr;
	m_solartransmissiontable  = nullptr;
	m_diffusegroundpointtable = nullptr;
	m_emissiontable           = nullptr;

}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::DiffuseGeometryTablesAlreadyExist		2010-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_EngineDiffuseTables::DiffuseTablesAlreadyExist()
{
	return     ( m_diffusetable             != nullptr )			// and we dopnt need to create 
		    && ( m_solartransmissiontable   != nullptr )
			&& ( m_diffusegroundpointtable	!= nullptr )
			&& ( m_emissiontable            != nullptr);						//NOte the emission table is optional, it may be null
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables::ConfigureOpticalTables		2010-4-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_EngineDiffuseTables::ConfigureOpticalTables(  bool										singlescatter,
														  const SKTRAN_TableOpticalProperties_V21*	opticalproptable,
														  SKTRAN_ThreadManager*						threadmanager )
{
	bool	ok;
	bool	ok1;
	bool	ok2;

	NXASSERT(( DiffuseTablesAlreadyExist() ));

	ok1 = singlescatter && m_solartransmissiontable->NotRequiredForSingleScatter();
	if (!ok1)
	{
		ok1 =  m_solartransmissiontable ->ConfigureOptical( singlescatter, opticalproptable, threadmanager );
	}

	ok2 = singlescatter;
	if (!ok2)
	{

		ok2 =        m_diffusetable           ->ConfigureOptical( singlescatter, opticalproptable, threadmanager );
		ok2 = ok2 && m_diffusegroundpointtable->ConfigureOptical( singlescatter, opticalproptable, threadmanager );
	}

	ok  = ok1 && ok2;
		
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables::ConfigureOpticalTables		2010-4-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_EngineDiffuseTables::ConfigureOpticalEmissionTables(  double wavelen, const SKTRAN_CoordinateTransform_V2* coords, SKTRAN_AtmosphericEmission* atmosphericemissions)
{
	bool	ok;
	ok = m_emissiontable->ConfigureOptical( wavelen, coords, atmosphericemissions);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineDiffuseTables::CreateEmptyDiffuseTables		2010-5-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_EngineDiffuseTables::CreateEmptyDiffuseTables(  const SKTRAN_SpecsInternal_V21*		modelspecifications )
{
	bool									ok;

	ReleaseResources();
	ok = modelspecifications->DiffuseSpecs()->CreateEmptyDiffuseTables( &m_diffusetable, &m_solartransmissiontable,&m_diffusegroundpointtable, &m_emissiontable);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_EngineDiffuseTables::CreateEmptyDiffuseTables, Error creating empty diffuse table");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CreateDiffuseTables		2010-3-30*/
/** Creates teh diffuse tables. The tables will use the thread manager to
 *	split the work up between multiple threads
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_EngineDiffuseTables::ConfigureDiffuseGeometryTables( bool									singlescatter, 
													             SKTRAN_SpecsInternal_V21*		modelspecifications,
													             SKTRAN_ThreadManager*					threadmanager )
{
	bool ok;
	bool ok0;
	bool ok1;
	bool ok2;

	ok = DiffuseTablesAlreadyExist();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_EngineDiffuseTables::CreateDiffuseGeometryTables, There was a problem creating the empty diffuse geometry tables");
	}
	else
	{
		ok0 = m_emissiontable->IsEmpty();
		if (!ok0)
		{
			ok0 = m_emissiontable->ConfigureGeometry(modelspecifications);
		}

		ok1 = m_solartransmissiontable->IsDefined();																	// Is the solar transmission table already initialized
		if (!ok1)																										// IF it is not
		{																												// then 
			ok1 = ( singlescatter && (modelspecifications->DiffuseSpecs()->LOSSingleScatterTableFactory() != nullptr));	// we dont need the solar transmission table if single scatter and each line of sight ray has its own single scatter table
			if (!ok1)
			{
				ok1 =         m_solartransmissiontable  ->ConfigureGeometry ( modelspecifications, threadmanager);
				ok1 = ok1 &&  m_solartransmissiontable  ->AttachOpticalToGeometry( );
			}
		}

		ok2 = m_diffusetable->IsDefined() || singlescatter;															// We dont need to configure the diffuse tables if they are
		if (!ok2)																									// already defined or we are just doing single scatter calcs
		{
			ok2 =        modelspecifications->DiffuseSpecsVar()->FirstTimeInitializeDiffuseObjects();
			ok2 = ok2 && m_diffusetable             ->ConfigureGeometry ( modelspecifications, threadmanager);
			ok2 = ok2 && m_diffusegroundpointtable  ->ConfigureGeometry ( modelspecifications, threadmanager);
			ok2 = ok2 && m_diffusetable             ->CreateJIndexTables_RayIntegral       ( threadmanager );
		    ok2 = ok2 && m_diffusegroundpointtable  ->CreateJIndexTables_HemisphereIntegral( threadmanager );
			ok2 = ok2 &&  m_diffusetable            ->AttachOpticalToGeometry( );
			ok2 = ok2 &&  m_diffusegroundpointtable ->AttachOpticalToGeometry( );
		}

		ok = ok0 && ok1 && ok2;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_EngineDiffuseTables::CreateDiffuseGeometryTables, There was a problem initializing the diffuse geometry tables");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::SKTRANSO_Engine		 2015- 3- 6*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_Engine::SKTRANSO_Engine()
{
	m_lastopticalstateobject        = nullptr;					// Add a few variables
	m_lastatmosphericemissionobject = nullptr;					// To make sure uses call ConfigureModel and CalculateRadiance within the rules;
	m_isfirsttimeafterinit          = false;

	m_opticalpropertiestable = nullptr;
	m_updateclimatologycache = true;
//	m_threadmanager = SKTRAN_ThreadManagerBoost_CreateNewInstance();			// Create a boost implementation of the thread manager by default
	m_threadmanager = SKTRAN_ThreadManagerOpenMP_CreateNewInstance();
	m_threadmanager->AddRef();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::~SKTRANSO_Engine		 2015- 3- 6*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_Engine::~SKTRANSO_Engine()
{
	if (m_threadmanager != nullptr)
	{
		m_threadmanager->Release();
		m_threadmanager = nullptr;
	}
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::ReleaseResources		 2015- 3- 6*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_Engine::ReleaseResources()
{
	m_tables.ReleaseResources();
	m_linesofsight.ReleaseResources();
	if (m_opticalpropertiestable != nullptr) m_opticalpropertiestable->Release();
	m_opticalpropertiestable = nullptr;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::ConfigureOpticalTables		2010-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::ConfigureOpticalTables(bool singlescatter, const SKTRAN_TableOpticalProperties_V21*  opticalproperties )
{
	bool	ok;

	ok =       m_threadmanager->SetOpticalProps( opticalproperties );
	ok = ok && m_tables.ConfigureOpticalTables ( singlescatter, opticalproperties, m_threadmanager);
	ok = ok && m_linesofsight.ConfigureOptical ( singlescatter, opticalproperties, m_threadmanager);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," SKTRANSO_Engine::ConfigureOpticalTables, There was an error configuring the optical tables");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::ConfigureModel		2007-12-12*/
/** Configures the radiative transfer model. Each call to this function
 *	resets all of the internal structures and tables.  This will force the
 *	all internal tables to recreate themselves upon the next call to
 *	CalculateRadiance.
 *
 *	The code also sets a request for the climatologies
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::ConfigureModel( SKTRAN_SpecsUser_Base& auserspecifications, const SKTRAN_LineOfSightArray_V21& linesofsight, size_t numthreads )
{
	bool	ok;
	SKTRANSO_SpecificationsUser*	userspecifications;
	
	ReleaseResources();

	userspecifications = dynamic_cast<SKTRANSO_SpecificationsUser*>(&auserspecifications);
	
	ok =       userspecifications->UpdateUndefinedParametersFromLinesOfSight ( linesofsight );		// Set the unit vector from Earth to Sun (Note: only one sun is used for all lines of sight)
	ok = ok && m_modelspecifications.Initialize ( userspecifications );								// setup the new model specifications
	ok = ok && m_tables.CreateEmptyDiffuseTables( &m_modelspecifications );
	ok = ok && m_threadmanager->CreateThreads   ( numthreads, &m_modelspecifications, &m_tables );
	ok = ok && CreateOpticalPropertyTables      ( );
	ok = ok && m_linesofsight.SetLinesOfSight   ( linesofsight, m_modelspecifications.CoordinateSystemPtr(), userspecifications->RayTracingRegionManagerVar() );							// setup the new lines of sight
	if (!ok)																			// see if it worked
	{																					// if it did not then log an error
		nxLog::Record(NXLOG_WARNING," SKTRANSO_Engine::ConfigureModel, there was an error configuring the lines of sight, thats a problem");
	}
	m_updateclimatologycache = true;
	m_lastopticalstateobject        = nullptr;					// Add a few variables
	m_lastatmosphericemissionobject = nullptr;					// To make sure uses call ConfigureModel and CalculateRadiance within the rules;
	m_isfirsttimeafterinit          = true;

	return ok;																			// and return the status;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CheckGeometry		2010-4-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::CheckGeometryTables( bool singlescatter )
{
	bool	ok;

	ok =       m_tables.ConfigureDiffuseGeometryTables  ( singlescatter, &m_modelspecifications, m_threadmanager );
	ok = ok && m_linesofsight.ConfigureGeometry         ( singlescatter, &m_modelspecifications, m_threadmanager );
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CreateOpticalPropertyTables		2010-5-25*/
/** **/
/*---------------------------------------------------------------------------*/
/*
bool SKTRANSO_Engine::CreateOpticalPropertyTables( )
{
	bool	ok;

	if (m_opticalpropertiestable != nullptr) m_opticalpropertiestable->Release();
	ok =       m_modelspecifications.OpticalTableSpecs()->CreateEmptyOpticalPropertiesTable( &m_opticalpropertiestable );
	ok = ok && m_opticalpropertiestable->ConfigureGeometry( &m_modelspecifications );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CreateOpticalPropertyTables, Error creating Optical properties table");
	}
	return ok;
}
*/



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CalculateOpticalPropertiesTable_SingleThreaded		2010-6-1*/
/** Calculate the optical properties table used by the model. This calculation
	is executed in a single thread as it calls all of the climatologies and
	optical properties. We founbd it much simpler to let all of this run
	in one thread rather than try to ensure that all climatologies and optical
	properties object are thread safe.
**/
/*---------------------------------------------------------------------------*/
/**
bool SKTRANSO_Engine::CalculateOpticalPropertiesTable_SingleThreaded( double wavelen, SKTRAN_AtmosphericOpticalState_V21*	 opticalstate, bool userupdateclimatology )
{
	GEODETIC_INSTANT							point;
	const SKTRAN_CoordinateTransform_V2*		coords;
	bool										ok;

	coords          = m_modelspecifications.CoordinateSystem();
	point.latitude  = coords->ReferencePtLatitude();
	point.longitude = coords->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd       = m_linesofsight.LinesOfSight()->MeanMJDWithCache();
	NXASSERT((point.mjd > 10000.0));
	if (point.mjd < 10000.0)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateOpticalPropertiesTable_SingleThreaded, the mjd being used for the climatologies is probably out of range. Its value is %e", (double)point.mjd );
	}
	ok  =        opticalstate->SetTimeAndLocation( point, m_updateclimatologycache || userupdateclimatology );	
	ok  = ok &&  m_opticalpropertiestable->ConfigureOptical( wavelen, *opticalstate );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateRadiance, Error calculating the optical properties table. Thats not good");
	}
	m_updateclimatologycache = !ok;
	return ok;
}
*/



/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineOptical_V2::CalculateIncomingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::CalculateIncomingRadianceAndScatter(bool ignorehighaltdiffuse)
{
	bool	ok;

	ok =       m_threadmanager->DiffusePointsTable_EvaluateIncomingRays   ( m_tables.DiffusePointsTableVar()->OpticalTableVar(), ignorehighaltdiffuse);
	ok = ok && m_threadmanager->DiffusePointsTable_ScatterIncomingRadiance( m_tables.DiffusePointsTableVar()->OpticalTableVar(), ignorehighaltdiffuse);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineOptical_V2::ScatterToOutgoingRadiation		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/


bool SKTRANSO_Engine::ScatterIncomingRadianceAtGround(bool ignorehighaltdiffuse)
{
	bool	ok;

//	ok =       m_threadmanager->DiffusePointsTable_ScatterRays( m_tables.DiffusePointsTableVar()->OpticalTableVar(), ignorehighaltdiffuse);
	NXTRACE_ONCEONLY(firsttime, ("*** SKTRANSO_Engine::ScatterIncomingRadiance need to multithread ground point scattering"));
	ok = m_tables.DiffuseGroundPointTableVar()->ScatterIncomingRadiance( m_threadmanager );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_EngineOptical_V2::InitializeFirstOrderIncomingRadiances		2008-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::InitializeFirstOrderIncomingRadiancesAndScatter()
{
	bool	ok;

	ok =       m_threadmanager->DiffusePointsTable_FirstOrderInitialize   ( m_tables.DiffusePointsTableVar()->OpticalTableVar() );
	ok = ok && m_threadmanager->DiffusePointsTable_ScatterIncomingRadiance( m_tables.DiffusePointsTableVar()->OpticalTableVar(), true);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::GetRayTracingGrid		2011-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::GetRayTracingGrid( std::vector<double>* shells ) const
{
	bool	ok;
	const	SKTRAN_SpecsInternal_RayTracing_V21*	raytracingspecs;
	size_t	numheights;
	size_t	idx;

	raytracingspecs = m_modelspecifications.RayTracingSpecs();
	ok = (raytracingspecs != nullptr);
	if (ok)
	{
		numheights = raytracingspecs->RayTracingShells()->NumGridPoints();
		shells->resize( numheights );
		for (idx = 0; idx < numheights; idx++)
		{
			shells->at(idx) = raytracingspecs->RayTracingShells()->At(idx);
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::ReferencePoint		2011-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

GEODETIC_INSTANT SKTRANSO_Engine::ReferencePoint() const
{
	GEODETIC_INSTANT point;

	point.latitude  = m_modelspecifications.CoordinateSystemPtr()->ReferencePtLatitude();
	point.longitude = m_modelspecifications.CoordinateSystemPtr()->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd       = m_modelspecifications.CoordinateSystemPtr()->ReferencePointMJD();
	return point;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CalculateOpticalPropertiesTable		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::CalculateOpticalPropertiesTable(double								 wavelen,
														SKTRAN_AtmosphericOpticalState_V21*  opticalstate,
														bool								 userupdateclimatology )
{

	return CalculateOpticalPropertiesTable_SingleThreaded(wavelen, opticalstate, userupdateclimatology);

}

/*-----------------------------------------------------------------------------
 *		SKTRANSO_Engine::CalculateOpticalPropertiesTable_SingleThreaded		2010-6-1*/
/** Calculate the optical properties table used by the model. This calculation
	is executed in a single thread as it calls all of the climatologies and
	optical properties. We founbd it much simpler to let all of this run
	in one thread rather than try to ensure that all climatologies and optical
	properties object are thread safe.
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::CalculateOpticalPropertiesTable_SingleThreaded( double wavelen, SKTRAN_AtmosphericOpticalState_V21*	 opticalstate, bool userupdateclimatology )
{
	GEODETIC_INSTANT							point;
	const SKTRAN_CoordinateTransform_V2*		coords;
	bool										ok;

	coords          = m_modelspecifications.CoordinateSystemPtr();
	point.latitude  = coords->ReferencePtLatitude();
	point.longitude = coords->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd       = coords->ReferencePointMJD();			// Reference point MJD is now a weighted mean, weighted towards target altitude, it makes a tiny difference due to floating roundoff
//	point.mjd       = m_linesofsight.LinesOfSight()->MeanMJD();

//	nxLog::Record(NXLOG_INFO, "mjdrefpoint = %25.18f, lines mean mjd = %25.18f", (double)coords->ReferencePointMJD(), (double)m_linesofsight.LinesOfSight()->MeanMJD());

	NXASSERT((point.mjd > 10000.0));
	if (point.mjd < 10000.0)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateOpticalPropertiesTable_SingleThreaded, the mjd being used for the climatologies is probably out of range. Its value is %e", (double)point.mjd );
	}
	ok  =        opticalstate->SetTimeAndLocation( point, m_updateclimatologycache || userupdateclimatology );	
	ok  = ok &&  m_opticalpropertiestable->ConfigureOptical( wavelen, *opticalstate );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateRadiance, Error calculating the optical properties table. Thats not good");
	}
	m_updateclimatologycache = !ok;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::ConfigureOpticalEmissionTables		 2015- 3- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::ConfigureOpticalEmissionTables(  double wavelen, SKTRAN_AtmosphericEmission* atmosphericemissions)
{
	bool										ok;

	ok = m_tables.EmissionsTableVar()->ConfigureOptical( wavelen, m_modelspecifications.CoordinateSystemPtr() , atmosphericemissions);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CreateOpticalPropertyTables		2010-5-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::CreateOpticalPropertyTables( )
{
	bool	ok;

	if (m_opticalpropertiestable != nullptr) m_opticalpropertiestable->Release();
	ok =       m_modelspecifications.OpticalTableSpecs()->CreateEmptyOpticalPropertiesTable( &m_opticalpropertiestable );
	ok = ok && m_opticalpropertiestable->ConfigureGeometry( &m_modelspecifications );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CreateOpticalPropertyTables, Error creating Optical properties table");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::CalculateRadiance		2010-4-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_Engine::CalculateRadiance( std::vector<SKTRAN_StokesScalar>*		  losradiance,
										   double								  wavelen,
										   size_t								  numordersofscatter,
										   SKTRAN_AtmosphericOpticalState_V21*	  opticalstate,
										   std::vector<skRTStokesVector>*         losvector,
										   bool									  userupdateclimatology,
										   SKTRAN_DiagnosticInterface*		      diag)


{
	SKTRAN_AtmosphericEmission*		atmosphericemissions;
	bool							ok= true;
	size_t							idx;
	size_t							num1;
//	SKTRAN_CodeTimer				s1;
	bool							singlescatter;
	size_t							orderofscatter;
	double*							buffer;
	bool							ignorehighalt;

	atmosphericemissions = opticalstate->EmissionObjectVar();
	if ( m_isfirsttimeafterinit	 )											// Is this the first call after ConfigureModel
	{																		// So save pointers 
		m_lastopticalstateobject        = opticalstate;						// So we can see if the user changes objects between calls to CalculateRadiance
		m_lastatmosphericemissionobject = atmosphericemissions;				// So we can checck that calls to ConfigureModel and CalculateRadiance within the rules;
		m_isfirsttimeafterinit          = false;
	}
	else																	// This is a subsequent call to CalculateRadiance
	{																		// Lets make sure the user is playing by the rules
		ok = (m_lastopticalstateobject == opticalstate) && (opticalstate != nullptr) && ( m_lastatmosphericemissionobject == atmosphericemissions);
		if (!ok)															// Technically we are okay not checking opticalstate but changes to atmosphericemissions can cause big problems
		{
			nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateRadiance, The opticalstate or atmospheric emissions object has changed since the last cal or opticalstate is null. This can create internal problems. You should call ConfigureModel first ");
		}
	}

	ok = m_tables.DiffuseTablesAlreadyExist();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_Engine::CalculateRadiance, You must successfully call ConfigureModel before calling CalculateRdaiance");
		losradiance->resize(0);
	}
	else
	{

		m_tables.EmissionsTableVar()->SetEmpty( atmosphericemissions == nullptr );					// Flag that The emissions table is empty to avoid making unecessary Jindex/JValue values
		ok =  CalculateOpticalPropertiesTable( wavelen, opticalstate, userupdateclimatology );

		if (ok)
		{
			singlescatter = (numordersofscatter < 2);
			ok = CheckGeometryTables( singlescatter );
			orderofscatter = 0;																// reset the number of orders of scatter currently executed.
			ok =       ConfigureOpticalEmissionTables(  wavelen, atmosphericemissions);
			ok = ok && ConfigureOpticalTables        ( singlescatter, OpticalProperties() );

			if (diag != nullptr) diag->Diagnose(SKTRAN_DIAGNOSE_CONFIGUREOPTICAL, orderofscatter, ok );			// Let external user diagnose this stage of processing

			if ( !singlescatter && ok )
			{
				ok =       InitializeFirstOrderIncomingRadiancesAndScatter();									// Initialize tables with first order incoming radiance
				ok = ok && ScatterIncomingRadianceAtGround ( true );											// Scatter the signal at ground, using the incoming rays from the last iteration
				orderofscatter++;																				// Increment to first order of scatter
				if (diag != nullptr) diag->Diagnose(SKTRAN_DIAGNOSE_FIRSTORDERINCOMING, orderofscatter, ok );		// Let external user diagnose this stage of processing

				if (ok)
				{
					if (numordersofscatter > 1 )
					{
						NXTRACE_ONCEONLY(firsttime,("SKTRANSO_Engine::CalculateRadiance, Got to check that we have the correct number of orders of scatter\n"));
						num1 = numordersofscatter - 2;						//  Changed from numordersofscatter-1 to numordersofscatter - 2 to get the correct number of orders of scatter
						for (idx = 0; idx < num1; idx++)
						{
							ignorehighalt =  (idx < (num1-1));
							ok = ok && CalculateIncomingRadianceAndScatter( ignorehighalt );								// and calculate the incoming rays
							ok = ok && ScatterIncomingRadianceAtGround ( true );											// Scatter the signal at ground, using the incoming rays from the last iteration
							orderofscatter++;
							if (diag != nullptr) diag->Diagnose(SKTRAN_DIAGNOSE_CALCULATEINCOMING, orderofscatter, ok );
						}
						if (diag != nullptr) diag->Diagnose(SKTRAN_DIAGNOSE_SCATTERINCOMING, orderofscatter, ok );
					}																			
				}
			}
			losradiance->resize( LinesOfSight()->NumRays() );
			buffer = &( (*losradiance)[0]);
			ok = ok && m_linesofsight.OpticalTableVar()->CalculateObserverIncomingRadiance( buffer, losradiance->size() );
			if (diag != nullptr) diag->Diagnose( SKTRAN_DIAGNOSE_FINISH, orderofscatter, ok );
		}

		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_EngineThreadController::ExecuteSasktranThread, Error performing SASKTRAN Calculation, results are untrustworthy");
		}
	}
	return ok;
}

