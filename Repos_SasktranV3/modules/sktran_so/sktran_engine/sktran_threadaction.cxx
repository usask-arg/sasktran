#include "../sasktranv21_internals.h"
#include <float.h>


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::DiffusePointsTable_ConfigureGeometryStage2		2010-3-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::DiffusePointsTable_ConfigureGeometryStage2( SKTRANSO_TableDiffusePoints* table )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_DiffusePointsTable_ConfigureGeometryStage2;
	m_diffusetable_geometry = table;
	numpoints               = table->NumDiffusePoints();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_diffusetable_geometry = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::TLS_DiffusePointsTable_ConfigureGeometryStage2		2010-3-23*/
/** The multi-threaded function configures points in the duffuse point table.
 *	Each thread will continue running until the diffuse table reports that
 *	there are no more points in the table to be processed.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_DiffusePointsTable_ConfigureGeometryStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool								ok;
	SKTRAN_DiffusePointGeometry_V21*	point;

	NXASSERT(( m_diffusetable_geometry != NULL ));
	ok    = m_diffusetable_geometry->LookupDiffusePoint( pointindex, &point );
	ok    = ok && point->ConfigureGeometry_Stage2MT( threadquadrature );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::DiffusePointsTable_LinkJIndices		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::DiffusePointsTable_CreateJIndices( SKTRANSO_TableDiffusePoints* table  )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_DiffusePointsTable_CreateJIndices;
	m_diffusetable_geometry = table;
	numpoints               = table->NumDiffusePoints();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_diffusetable_geometry = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_DiffusePointsTable_CreateJIndices		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_DiffusePointsTable_CreateJIndices( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool									ok;
	SKTRAN_DiffusePointGeometry_V21*		point;

	NXASSERT(( m_diffusetable_geometry != NULL ));
	ok  = m_diffusetable_geometry->LookupDiffusePoint( pointindex, &point );
	ok  = ok && point->CreateJIndexTables_MT( threadquadrature );
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::DiffusePointsTable_ConfigureOpticalStage2		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::DiffusePointsTable_ConfigureOpticalStage2( SKTRAN_TableDiffusePointsOptical_V21* table, const SKTRAN_TableOpticalProperties_V21*	optprop )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction       = &SKTRAN_ThreadManager::TLS_DiffusePointsTable_ConfigureOpticalStage2;
	m_diffusetable_optical = table;
	m_optprop              = optprop;
	numpoints              = table->NumDiffusePoints();
	ok                     = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_diffusetable_optical = NULL;
	m_optprop              = NULL;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_DiffusePointsTable_ConfigureOpticalStage2		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_DiffusePointsTable_ConfigureOpticalStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool									ok;
	SKTRAN_DiffusePointOptical_V21*			point;

	NXASSERT(( m_diffusetable_optical != NULL ));
	NXASSERT(( m_optprop              != NULL ));
	ok    =       m_diffusetable_optical->LookupDiffusePoint( pointindex, &point);
	ok    = ok && point->ConfigureOptical  ( m_optprop, threadquadrature );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::DiffusePointsTable_EvaluateIncomingRays		2010-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::DiffusePointsTable_EvaluateIncomingRays( SKTRAN_TableDiffusePointsOptical_V21* table , bool ignorehighalt)
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_DiffusePointsTable_EvaluateIncomingRays;
	m_diffusetable_optical  = table;
	m_ignorehighalt         = ignorehighalt;
	numpoints               = table->NumDiffusePoints();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_diffusetable_optical  = NULL;
	m_ignorehighalt         = false;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_DiffusePointsTable_ConfigureOpticalStage2		2010-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_DiffusePointsTable_EvaluateIncomingRays( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool									ok;
	SKTRAN_DiffusePointOptical_V21*			point;

	NXASSERT(( m_diffusetable_optical != NULL ));
	ok    =       m_diffusetable_optical->LookupDiffusePoint( pointindex, &point);
	ok    = ok && point->CalculateIncomingRadiances( m_ignorehighalt );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::DiffusePointsTable_EvaluateIncomingRays		2010-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::DiffusePointsTable_ScatterIncomingRadiance( SKTRAN_TableDiffusePointsOptical_V21* table , bool ignorehighalt)
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_DiffusePointsTable_ScatterIncomingRadiance;
	m_diffusetable_optical  = table;
	m_ignorehighalt         = ignorehighalt;
	numpoints               = table->NumDiffusePoints();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_diffusetable_optical  = NULL;
	m_ignorehighalt         = false;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_DiffusePointsTable_ScatterIncomingRadiance		2012-9-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_DiffusePointsTable_ScatterIncomingRadiance( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool									ok;
	SKTRAN_DiffusePointOptical_V21*			point;

	NXASSERT(( m_diffusetable_optical != NULL ));
	ok    =       m_diffusetable_optical->LookupDiffusePoint( pointindex, &point);
	ok    = ok && point->ScatterIncomingRadiance   ( m_ignorehighalt );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::DiffusePointsTable_FirstOrderInitialize		2010-6-8*/
/** Brings in the first order signal to initailize the incoming signal and
 *	scatters to the outbound directions
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::DiffusePointsTable_FirstOrderInitialize( SKTRAN_TableDiffusePointsOptical_V21* table)
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_DiffusePointsTable_FirstOrderInitialize;
	m_diffusetable_optical  = table;
	numpoints               = table->NumDiffusePoints();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_diffusetable_optical  = NULL;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_DiffusePointsTable_FirstOrderInitialize		2010-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_DiffusePointsTable_FirstOrderInitialize( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool									ok;
	SKTRAN_DiffusePointOptical_V21*			point;

	NXASSERT(( m_diffusetable_optical != NULL ));
	ok    =       m_diffusetable_optical->LookupDiffusePoint( pointindex, &point);
	ok    = ok && point->InitializeFirstOrderIncomingRadiances();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SolarTransmissionTable_ConfigureGeometryStage2		2010-3-25*/
/** Called by the main processing thread to launch multiple threads. Each thread
 *	will call TLS_SolarTransmissionTable_ConfigureGeometryStage2 and thus all
 *	the threads will process the geometry points in the solar transmission table
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::SolarTransmissionTable_ConfigureGeometryStage2( SKTRANSO_TableSolarTransmission* solartable, const SKTRAN_SpecsInternal_V21* specs  )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction      = &SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureGeometryStage2;
	m_solartable_geometry = solartable;
	m_specifications      = specs;
	numpoints             = solartable->NumPoints();
	ok                    = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_solartable_geometry = NULL;
	m_specifications      = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureGeometryStage2		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureGeometryStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool										ok;

	NXASSERT(( m_solartable_geometry != NULL ));
	NXASSERT(( m_specifications    != NULL ));

	ok = m_solartable_geometry->ConfigureGeometry_Stage2MT( pointindex, threadquadrature, m_specifications );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::SolarTransmissionTable_ConfigureGeometryStage2		2010-3-25*/
/** Called by the main processing thread to launch multiple threads. Each thread
 *	will call TLS_SolarTransmissionTable_ConfigureGeometryStage2 and thus all
 *	the threads will process the geometry points in the solar transmission table
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::SolarTransmissionTable_ConfigureOpticalStage2( SKTRANSO_TableSolarTransmission*  solartable, const SKTRAN_TableOpticalProperties_V21*	optprop )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction     = &SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureOpticalStage2;
	m_solartable_geometry = solartable;
	m_optprop            = optprop;
	numpoints            = solartable->NumPoints();
	ok                   = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_solartable_geometry = NULL;
	m_optprop            = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureGeometryStage2		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureOpticalStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool										ok;

	NXASSERT(( m_solartable_geometry != NULL ));
	NXASSERT(( m_optprop            != NULL ));

	ok    = m_solartable_geometry->ConfigureOptical_Stage2MT( pointindex, threadquadrature );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::GroundPointsTable_ConfigureGeometryStage2		2010-3-25*/
/** Called by the main processing thread to launch multiple threads. Each thread
 *	will call TLS_GroundPointsTable_ConfigureGeometryStage2 and thus all
 *	the threads will process the geometry points in the diffuse ground point table
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::GroundPointsTable_ConfigureGeometryStage2( SKTRANSO_TableGroundPointDiffuse* groundpointstable )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction             = &SKTRAN_ThreadManager::TLS_GroundPointsTable_ConfigureGeometryStage2;
	m_groundpointtable_geometry  = groundpointstable;
	numpoints                    = groundpointstable->NumPoints();
	ok                           = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_groundpointtable_geometry  = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_SolarTransmissionTable_ConfigureGeometryStage2		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_GroundPointsTable_ConfigureGeometryStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool										ok;

	NXASSERT(( m_groundpointtable_geometry != NULL ));
	ok = m_groundpointtable_geometry->ConfigureGeometry_Stage2MT( pointindex, threadquadrature );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::GroundPointsTable_ConfigureOpticalStage2		2010-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::GroundPointsTable_ConfigureOpticalStage2( SKTRANSO_TableGroundPointDiffuse*		    groundpointstable,
																	 const SKTRAN_TableOpticalProperties_V21*		optprop )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction             = &SKTRAN_ThreadManager::TLS_GroundPointsTable_ConfigureOpticalStage2;
	m_groundpointtable_geometry   = groundpointstable;
	m_optprop                    = optprop;

	numpoints                    = groundpointstable->NumPoints();
	ok                           = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_groundpointtable_geometry  = NULL;
	m_optprop                    = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_GroundPointsTable_ConfigureOpticalStage2		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_GroundPointsTable_ConfigureOpticalStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex  )
{
	bool										ok;

	NXASSERT(( m_groundpointtable_geometry != NULL ));
	ok    = m_groundpointtable_geometry->ConfigureOptical_Stage2MT( pointindex, m_optprop, threadquadrature );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::GroundPointsTable_ScatterIncomingRadiance		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::GroundPointsTable_ScatterIncomingRadiance( SKTRANSO_TableGroundPointDiffuse* groundpointstable )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction             = &SKTRAN_ThreadManager::TLS_GroundPointsTable_ScatterIncomingRadiance;
	m_groundpointtable_geometry   = groundpointstable;

	numpoints                    = groundpointstable->NumPoints();
	ok                           = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_groundpointtable_geometry  = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_GroundPointsTable_ScatterIncomingRadiance		2010-3-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_GroundPointsTable_ScatterIncomingRadiance( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool										ok;

	NXASSERT(( m_groundpointtable_geometry != NULL ));
	ok    = m_groundpointtable_geometry->ScatterIncomingRadiance_MT( pointindex );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::GroundPointsTable_CreateJIndices		2010-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::GroundPointsTable_CreateJIndexTables( SKTRANSO_TableGroundPointDiffuse* groundpointstable)
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction             = &SKTRAN_ThreadManager::TLS_GroundPointsTable_CreateJIndexTables;
	m_groundpointtable_geometry  = groundpointstable;
	numpoints                    = groundpointstable->NumPoints();
	ok                           = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_groundpointtable_geometry  = NULL;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_GroundPointsTable_CreateJIndices		2010-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_GroundPointsTable_CreateJIndexTables( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool										ok;
	SKTRAN_GroundPointDiffuseGeometry_V21*		point;

	NXASSERT(( m_groundpointtable_geometry != NULL ));
	point = m_groundpointtable_geometry->PointAtVar(pointindex);
	ok    = threadquadrature->CreateJIndexTable_GroundPointIncomingRaysUpwardFlux( point );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::DiffusePointsTable_ConfigureGeometryStage2		2010-3-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::LinesOfSightTable_ConfigureGeometryStage2( SKTRAN_TableLinesOfSight_V21*					 lostable,
																      const SKTRAN_LineOfSightArray_V21*			     observerspecs,
																      const SKTRAN_SpecsInternal_V21*				 specs,
																	  bool											 singlescatter
																    )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_LinesOfSightTable_ConfigureGeometryStage2;
	m_lostable_geometry     = lostable;
	m_observerspecs         = observerspecs;
	m_specifications        = specs;
	m_singlescatter         = singlescatter;
	numpoints               = lostable->NumRays();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_lostable_geometry     = NULL;
	m_observerspecs         = NULL;
	m_specifications        = NULL;
	m_singlescatter         = false;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadInstanceBoost::TLS_DiffusePointsTable_ConfigureGeometryStage2		2010-3-23*/
/** Multi-threaded code to create the geometry for the lines of sight. The Sasktran model will
 *	eventually calculate the radiance along these lines of sight. In this stage we simply
 *	create the rays and then create the Jindex tables for each of the contributing
 *	rays.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_LinesOfSightTable_ConfigureGeometryStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature, size_t pointindex )
{
	bool										ok;
	SKTRANSO_RayLOSGeometry_V21*				ray;
	HELIODETIC_VECTOR							observer;
	HELIODETIC_UNITVECTOR						look;
	nxVector									offsetobserver;
	const SKTRAN_LineOfSightEntry_V2*			entry;
	const SKTRAN_CoordinateTransform_V2*		coords; 
	const SKTRAN_SpecsInternal_RayTracing_V21*	raytracingspecs;



	NXASSERT(( m_lostable_geometry != NULL ));
	NXASSERT(( m_specifications    != NULL ));
	NXASSERT(( m_observerspecs     != NULL ));

	coords          = m_specifications->CoordinateSystemPtr();
	raytracingspecs = m_specifications->RayTracingSpecs();
	entry           = m_observerspecs->Entry( pointindex );
	ray             = m_lostable_geometry->RayAtVar( pointindex );
//	#if defined(NXDEBUG)
//		NXTRACE_ONCEONLY(firsttime, ("****** SKTRAN_ThreadManager::TLS_LinesOfSightTable_ConfigureGeometryStage2, REMOVE THE TEST CODE THAT CHEKS TANGENT POINTS **************\n"));
//		nxGeodetic	geoid = coords->OsculatingGeoid();
//		double		latitude;
//		double		longitude;
//		double		height;
//
//		geoid.FromTangentPointLocation( entry->m_observer, entry->m_look);
//		latitude = geoid.GeodeticLatitude();
//		longitude = geoid.GeodeticLongitude();
//		height    = geoid.Height();
//		NXTRACE(("Internal line of sight tangent point (lat = %12.7f, lng = %12.7f, height= %15.7f)\n", (double)latitude, (double)longitude, (double)height));
//#endif
	offsetobserver = coords->TranslateGeoidToOsculatingSphere(entry->Observer());
	//offsetobserver = entry->Observer();		// 2013-06-25, ndl303, The translation to osculating sphere is now done in initialization in SKTRAN_TableLinesOfSight_V21::SetLinesOfSight where look is also modified to give the same tangent altitude
	observer       = coords->GeographicToHelio( offsetobserver );
	look           = coords->GeographicToHelio( entry->Look()  ).UnitVector();
		
	ok =       ray->RayVar()->MoveObserver( observer, look );
	ok = ok && ray->RayVar()->TraceRay_NewMethod(  );
	ok = ok && ray->ConfigureInternalSolarTransmissionTableGeometry	( threadquadrature );				// Configure the internal solar scattering table if it is being used (normally for observer lines of sight).

	ok = ok && ray->IndicesVar()->CreateJIndexTables( threadquadrature, ray, m_singlescatter );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::LinesOfSightTable_ConfigureOpticalStage2		2010-4-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::LinesOfSightTable_ConfigureOpticalStage2( SKTRAN_TableLinesOfSightOptical_V21*      lostable,
																	 bool									   singlescatter
																	  )
{
	size_t	numpoints;
	bool	ok;

	m_threadfunction        = &SKTRAN_ThreadManager::TLS_LinesOfSightTable_ConfigureOpticalStage2;
	m_lostable_optical      = lostable;
	m_singlescatter		    = singlescatter,
	numpoints               = lostable->NumRays();
	ok                      = NotifyWorkerThreadsAndWaitForEnd( numpoints );
	m_lostable_optical     = NULL;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_ThreadManager::TLS_LinesOfSightTable_ConfigureOpticalStage2		2010-4-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_ThreadManager::TLS_LinesOfSightTable_ConfigureOpticalStage2( SKTRANSO_Quadrature_TLS_V21* threadquadrature,  size_t pointindex  )
{
	bool									ok;
	SKTRAN_RayInternalDiffuseOptical_V2*	ray;


	NXASSERT(( m_lostable_optical != NULL ));

	ray     = m_lostable_optical->RayAtVar( pointindex );
	ok      = ray->ConfigureOptical( threadquadrature, m_singlescatter, false,false);									// Configure the optical, diffuse points are never single scatter, reset the transmission and cell factor caches
	return ok;
}

