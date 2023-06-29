#include "../sasktranv21_internals.h"




/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::SKTRANSO_TableGroundPointDiffuse		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableGroundPointDiffuse::SKTRANSO_TableGroundPointDiffuse()
{
	m_groundpoints    = NULL;
	m_numgroundpoints = 0;
	m_opticaltable    = new SKTRAN_TableGroundPointDiffuseOptical_V21(this);
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::~SKTRANSO_TableGroundPointDiffuse		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_TableGroundPointDiffuse::~SKTRANSO_TableGroundPointDiffuse()
{
	ReleaseResources();
	delete m_opticaltable;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::ReleaseResources		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_TableGroundPointDiffuse::ReleaseResources()
{
	if (m_groundpoints  != NULL ) delete [] m_groundpoints;
	m_groundpoints       = NULL;
	m_numgroundpoints    = 0;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::ConfigureGeometry		2010-4-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::ConfigureGeometry( const SKTRAN_SpecsInternal_V21* specs, SKTRAN_ThreadManager* threadmanager )
{
	bool	ok;

	ok = ConfigureGeometry_Stage1( specs );
	ok = ok && threadmanager->GroundPointsTable_ConfigureGeometryStage2( this );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRANSO_TableGroundPointDiffuse::ConfigureGeometry, Error configuring geometry, thats a problem");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::AttachToGeometry		2010-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::AttachOpticalToGeometry( )
{
	return m_opticaltable->AttachToGeometry();
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::PointAtVar		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GroundPointDiffuseGeometry_V21* SKTRANSO_TableGroundPointDiffuse::PointAtVar(size_t idx)
{
	NXASSERT(( idx < m_numgroundpoints ));
	return m_groundpoints+idx;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::PointAt		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_GroundPointDiffuseGeometry_V21* SKTRANSO_TableGroundPointDiffuse::PointAt(size_t idx) const
{
	NXASSERT(( idx < m_numgroundpoints ));
	return m_groundpoints+idx;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::ConfigureGeometry_Stage1		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::ConfigureGeometry_Stage1( const SKTRAN_SpecsInternal_V21*	specs  )
{
	bool											ok;

	ReleaseResources();																		// Releease all of the existing sources
	m_numgroundpoints = NumberOfGroundPointsRequired(specs);										// Get the
	ok = (m_numgroundpoints ==0);
	if (!ok) 																				// If we have at least one ground point
	{																						// then 
		m_groundpoints    = new SKTRAN_GroundPointDiffuseGeometry_V21 [ m_numgroundpoints ];	// allocate the memory
		ok = (m_groundpoints != NULL );														// and make sure it worked ok.
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableGroundPointDiffuse_Colocated::ConfigureGeometry_Stage1, Error configuring ground point geometry");
		ReleaseResources();																	// Releease all of the existing sources
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::ConfigureOptical		2008-2-10*/
/** This code will configure the optical component o fthe diffuse ground
 *	points. The code uses the threadmanager to use the worker threads to process
 *	each ground point.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::ConfigureOptical( bool singlescatter, const SKTRAN_TableOpticalProperties_V21* opticalprops, SKTRAN_ThreadManager* threadmanager )
{
	bool ok;

	NXTRACE_ONCEONLY(firsttime,("SKTRANSO_TableGroundPointDiffuse::ConfigureOptical, Need to move call to AttachToGeometry into ConfigureGeometry call\n"));

//	ok  =       m_opticaltable->AttachToGeometry(singlescatter);
	ok  = threadmanager->GroundPointsTable_ConfigureOpticalStage2( this, opticalprops );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableGroundPointDiffuse::ConfigureOptical::ConfigureOptical, Error configuring ground point wavelength dependent terms");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointOptical_V2::ConfigureOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRANSO_TableGroundPointDiffuse::ConfigureOptical_Stage2MT( size_t pointindex, const SKTRAN_TableOpticalProperties_V21* optprop, SKTRANSO_Quadrature_TLS_V21* threadquadrature)
{
	bool	ok;

	ok = PointAtVar(pointindex)->OpticalPoint()->ConfigureOptical(optprop, threadquadrature);
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointOptical_V2::ConfigureOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::ScatterIncomingRadiance(SKTRAN_ThreadManager* threadmanager)
{
	bool			ok;

	ok  = threadmanager->GroundPointsTable_ScatterIncomingRadiance( this );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableGroundPointOptical_V2::ScatterIncomingRadiance, Error calculating upward flux from incoming radiances");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointOptical_V2::ScatterIncomingRadiance_MT		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::ScatterIncomingRadiance_MT(size_t pointindex )
{
	bool	ok;

	ok = PointAtVar(pointindex)->OpticalPoint()->UpdateDiffuseUpwardFlux();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::ConvertJIndexToRadiancePtr		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRANSO_TableGroundPointDiffuse::ConvertJIndexToRadiancePtr( const SKTRANSO_JIndex* entry, ENUM_SKTRAN_JSOURCE jsource ) const
{
	size_t						positionindex;
	const SKTRAN_StokesScalar*	Jptr;

	positionindex = entry->PositionIndex(); 
	Jptr = PointAt(positionindex)->OpticalPoint()->DiffuseUpwardFlux();
	return Jptr;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableGroundPointDiffuse::CreateJIndexTables_MT		2010-3-26*/
/** **/
/*---------------------------------------------------------------------------*/

/* bool SKTRANSO_TableGroundPointDiffuse::CreateJIndexTables_HemisphereIntegral_MT( size_t pointindex, SKTRANSO_Quadrature_TLS_V21* threadquadrature)
{
	bool	ok;

	ok = PointAtVar(pointindex)->CreateIncomingJIndexTableForDownwardFlux( );
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::CreateJIndexTables		2010-3-26*/
/**	Create the JIndexTables for the ground point table. This is multi-threaded
 *	although in reality it does not really need it until we get to lots of diffuse ground
 *	points
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_TableGroundPointDiffuse::CreateJIndexTables_HemisphereIntegral( SKTRAN_ThreadManager* threadmanager)
{
	bool								ok;

	ok = threadmanager->GroundPointsTable_CreateJIndexTables( this );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_TableGroundPointDiffuse::CreateJIndexTables_HemisphereIntegral, Error Linking downward flux to Incoming radiances");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::SKTRAN_TableGroundPointDiffuseOptical_V21		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableGroundPointDiffuseOptical_V21::SKTRAN_TableGroundPointDiffuseOptical_V21(SKTRANSO_TableGroundPointDiffuse* groundpointtable)
{
	m_geometrytable        = groundpointtable;
	m_groundpoints         = NULL;
	m_numgroundpointsresvd = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::~SKTRAN_TableGroundPointDiffuseOptical_V21		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableGroundPointDiffuseOptical_V21::~SKTRAN_TableGroundPointDiffuseOptical_V21(  )
{
	ReleaseResources();
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::ReleaseResources		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableGroundPointDiffuseOptical_V21::ReleaseResources()
{
	size_t idx;

	if (m_groundpoints != NULL )
	{
		for (idx = 0; idx < m_numgroundpointsresvd; idx++)
		{
			if  (m_groundpoints[idx] != NULL ) m_groundpoints[idx]->Release();
		}
		delete [] m_groundpoints;
	}
	m_groundpoints  = NULL;
	m_numgroundpointsresvd = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::AllocateStorage		2009-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableGroundPointDiffuseOptical_V21::AllocateStorage( size_t numpoints )
{
	bool									ok;
	size_t									idx;
	SKTRANSO_GroundPointDiffuseOptical*	opticalpoint;

	ok = (numpoints == m_numgroundpointsresvd);
	if (!ok)
	{
		ReleaseResources();
		m_groundpoints = new SKTRANSO_GroundPointDiffuseOptical* [numpoints];
		ok = (m_groundpoints != NULL );
		if (ok)
		{
			m_numgroundpointsresvd = numpoints;
			for (idx = 0; idx < numpoints; idx++)
			{
				opticalpoint = m_geometrytable->PointAtVar(idx)->OpticalPoint();
				opticalpoint->AddRef();
				m_groundpoints[idx] = opticalpoint;
			}
		}
		else
		{
			m_numgroundpointsresvd = 0;
			nxLog::Record(NXLOG_ERROR,"SKTRAN_TableGroundPointDiffuseOptical_V21::AllocateStorage, Error allocating memory for %u ground points", (size_t)numpoints);
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::AttachToGeometry		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableGroundPointDiffuseOptical_V21::AttachToGeometry()
{
	size_t									numpoints;
	bool									ok;
	bool									ok1;
	size_t									idx;


	numpoints = m_geometrytable->NumPoints();						
	ok = AllocateStorage( numpoints );
	if (ok)
	{
		for (idx = 0; idx < numpoints; idx++ )
		{
			ok1 = m_groundpoints[idx]->AttachToGeometry();
			ok = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableGroundPointDiffuseOptical_V21::AttachToGeometry, Error attaching to the geometry, This is a problem");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::ConfigureOptical_Stage2MT		2010-4-15*/
/** This function is called by the thread manager worker threads as part of the
 *	ConfigureOptical processing stage. It is intrinsically multi-threaded
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableGroundPointDiffuseOptical_V21::ConfigureOptical_Stage2MT( size_t									pointindex,
																		   const SKTRAN_TableOpticalProperties_V21*	optprop,
																		   SKTRANSO_Quadrature_TLS_V21*				 threadquadrature )
{
	bool	ok;

	ok  = m_groundpoints[pointindex]->ConfigureOptical( optprop, threadquadrature );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::ConfigureOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_TableGroundPointDiffuseOptical_V21::ScatterIncomingRadiance( )
{
	size_t			idx;
	bool			ok;

	ok           = true;
	for (idx = 0; idx < m_numgroundpointsresvd; idx++ )
	{
		ok  = ok && m_groundpoints[idx]->UpdateDiffuseUpwardFlux();
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableGroundPointDiffuseOptical_V21::UpdateOutboundJ, Error configuring the optical ground points");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::LookupJValue_DiffuseUpwardFlux		2010-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRAN_TableGroundPointDiffuseOptical_V21::LookupJValue_DiffuseUpwardFlux( const SKTRANSO_JIndex* index ) const
{
	size_t									pointindex;
	const SKTRAN_StokesScalar*				radianceptr = NULL;

	pointindex  = index->PositionIndex();
	radianceptr = m_groundpoints[pointindex]->DiffuseUpwardFlux(); 

	return radianceptr;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuseOptical_V21::LookupJValue_DiffuseIncomingRays		2010-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRAN_TableGroundPointDiffuseOptical_V21::LookupJValue_DiffuseIncomingRays( const SKTRANSO_JIndex* index ) const
{
	size_t									rayindex;
	size_t									pointindex;
	const SKTRAN_DiffusePointOptical_V21*	opticalpoint;
	const SKTRAN_StokesScalar*				radianceptr = NULL;

	pointindex   = index->PositionIndex();
	rayindex     = index->VertexIndex();
	opticalpoint = m_groundpoints[pointindex]->Geometry()->DiffuseGeometryPoint()->OpticalPoint();
	radianceptr  = opticalpoint->IncomingRadianceArray() + rayindex;
	return radianceptr;
}









