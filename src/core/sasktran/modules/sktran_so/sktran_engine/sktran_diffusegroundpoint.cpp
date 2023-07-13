#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::SKTRAN_GroundPointDiffuseGeometry_V21		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GroundPointDiffuseGeometry_V21::SKTRAN_GroundPointDiffuseGeometry_V21()
{
	m_opticalpoint   = NULL;
	m_diffusepoint   = NULL;
	m_pointindex     = (SKTRAN_GridIndex)(-1);
	m_opticalpoint   = new SKTRANSO_GroundPointDiffuseOptical(this);
	m_opticalpoint->AddRef();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::~SKTRAN_GroundPointDiffuseGeometry_V21		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GroundPointDiffuseGeometry_V21::~SKTRAN_GroundPointDiffuseGeometry_V21()
{
	ReleaseResources();
	m_opticalpoint->Release();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::ReleaseResources		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_GroundPointDiffuseGeometry_V21::ReleaseResources()
{
	m_diffusepoint   = NULL;
	m_pointindex     = (SKTRAN_GridIndex)(-1);
	m_diffusedownwardfluxtable.Clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::CreateIncomingJIndexTableForDownwardFlux		2009-1-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GroundPointDiffuseGeometry_V21::CreateIncomingJIndexTableForDownwardFlux()
{
	SKTRANSO_JIndex						jdescriptor;
	size_t									numrays;
	const SKTRANSO_RayInternalGeometry*		ray;
	HELIODETIC_UNITVECTOR					zenith;
	double									mu;
//	const double*							domega;
//	size_t									numdomega;
	double									omega;
	bool									ok;
	bool									ok1;
	size_t									idx;
	HELIODETIC_UNITVECTOR					look;

	numrays = m_diffusepoint->NumLOSIn();															// Get the number of rays coming into this ground point
	zenith  = m_diffusepoint->LocalZenith();														// Get the local zenith
//	ok      = m_diffusepoint->IncomingSolidAngleArray(&domega, &numdomega);							// Get the incoming solid angle array from the scattering matrix manager
//	ok      = ok && (numdomega == numrays );														// Make sure the array sizes still jive with each other
	ok      = m_diffusedownwardfluxtable.AllocateMaximumStorage(1,numrays,numrays);			// Allocate the storage for this.
	if (ok)
	{
		m_diffusedownwardfluxtable.ResetCounters();													// Reset the counters in this temporary storage area
		for (idx = 0; idx < numrays; idx++ )
		{
			ray   = &m_diffusepoint->LOSAt(idx);														// Get the next ray	
			look  = ray->Storage()->AverageLookVectorAwayFromObserver(0);										// Get the look vector (away from point)
			mu    = fabs(zenith & look);																// Get cos(angle) between zenith and ray. Note rays are shooting out of this point)
			NXASSERT((mu >= 0.0));
			omega = m_diffusepoint->IncomingUnitVectors()->CubatureWeightAt(idx);
			jdescriptor.ConfigureGroundPointTableIndex(m_pointindex, 1.0, idx, mu*omega );					// get the point indexing, weights and the vertex index points to the incoming ray
			ok1   = m_diffusedownwardfluxtable.InsertQuadraturePointEntries( &jdescriptor, 1);		// Each ground point only uses 
			ok = ok && ok1;																			// 
		}																							// and we are done
		m_diffusedownwardfluxtable.FinishCellEntries();
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GroundPointDiffuseGeometry_V21::CreateDiffuseIncomingJIndexTable, Error creating upward flux Jindex. this is a probblem that needs to be fixed!!");
		m_diffusedownwardfluxtable.Clear();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::CreateIncomingSolarJIndexTable		2009-1-28*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_GroundPointDiffuseGeometry_V21::CreateIncomingSolarJIndexTable( SKTRAN_JInterpolator_V2*	interpolator )
{
	HELIODETIC_UNITVECTOR		dummy;
	const SKTRANSO_JIndex*		index;
	size_t						numvertex;
	bool						ok;
	
	ok =       interpolator->InterpolateSourceFunctionTable( SKTRAN_JGROUNDSS, m_diffusepoint->Location(), dummy, &index, &numvertex, m_diffusepoint->Location().CosSZA() );
	ok = ok && m_incomingsolartable.AllocateMaximumStorage(1,1,numvertex);
	ok = ok && m_incomingsolartable.ResetCounters();													// Reset the counters in this temporary storage area
	ok = ok && m_incomingsolartable.InsertQuadraturePointEntries( index, numvertex );
	ok = ok && m_incomingsolartable.FinishCellEntries();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GroundPointDiffuseGeometry_V21::CreateIncomingSolarJIndexTable, Error creating table for single scatter upward flux. This is a probblem that needs to be fixed!!");
		m_incomingsolartable.Clear();
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::ConfigureGeometry		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GroundPointDiffuseGeometry_V21::ConfigureGeometry( SKTRAN_GridIndex pointindex, const SKTRAN_DiffusePointGeometry_V21* diffusepoint)
{
	bool	ok;

	ReleaseResources();
	m_diffusepoint = diffusepoint;
	m_pointindex   = pointindex;
	ok = (m_diffusepoint != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GroundPointDiffuseGeometry_V21::ConfigureGeometry, Error allocating AlbedoGeometry object");
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GroundPointDiffuseGeometry_V21::LocalZenith		2009-1-20*/
/** **/
/*---------------------------------------------------------------------------*/

HELIODETIC_UNITVECTOR SKTRAN_GroundPointDiffuseGeometry_V21::LocalZenith() const
{
	NXASSERT((m_diffusepoint != NULL));
	return m_diffusepoint->LocalZenith();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_GroundPointDiffuseOptical::SKTRANSO_GroundPointDiffuseOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_GroundPointDiffuseOptical::SKTRANSO_GroundPointDiffuseOptical( const SKTRAN_GroundPointDiffuseGeometry_V21* geometry  )
{
	m_FupMS         = -9999;
	m_geometrypoint = geometry;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_GroundPointDiffuseOptical::~SKTRANSO_GroundPointDiffuseOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_GroundPointDiffuseOptical::~SKTRANSO_GroundPointDiffuseOptical()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_GroundPointDiffuseOptical::ReleaseResources		2008-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_GroundPointDiffuseOptical::ReleaseResources()
{
	m_FupMS = -9999;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GrounPointOptical::ConfigureOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_GroundPointDiffuseOptical::AttachToGeometry()
{
	bool								ok;
	
	ReleaseResources();
	NXASSERT(( (m_geometrypoint != NULL) ));
	ok = m_diffuseupwardflux.AttachToGeometry( *(m_geometrypoint->DownwardFluxJIndexTable()) );
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"SKTRANSO_GroundPointDiffuseOptical::AttachToGeometry, Error attaching to ground point geometry, probably failed creating albedo object");
		ReleaseResources();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::ConvertJIndexToRadiancePtr		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar* SKTRANSO_GroundPointDiffuseOptical::ConvertJIndexToRadiancePtr( const SKTRANSO_JIndex* entry, ENUM_SKTRAN_JSOURCE jsource ) const
{
	size_t						incomingrayindex;
	const SKTRAN_StokesScalar*	JPtr;

	incomingrayindex  = entry->VertexIndex();
	JPtr              = m_geometrypoint->DiffuseGeometryPoint()->OpticalPoint()->IncomingRadianceArray() + incomingrayindex;
	return JPtr;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GrounPointOptical::ConfigureOptical		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_GroundPointDiffuseOptical::ConfigureOptical( const SKTRAN_TableOpticalProperties_V21* optProp, SKTRANSO_Quadrature_TLS_V21* quadrature)
{
	double										albedo;
	bool										ok;
	

	m_FupMS  = 0.0;																										// Initialize the diffuse upward flux to zero.
	ok       =       optProp->Get_AlbedoForDeprecatedLegacyCode( m_geometrypoint->DiffuseGeometryPoint()->Location(), &albedo );								// Get the albedo at this heliodetic location
	ok       = ok && m_diffuseupwardflux.SetWeightsAndRadiancePtrs( this, SKTRAN_JGROUND_UPFLUX ); 
	ok       = ok && m_diffuseupwardflux.AdjustWeightsByConstantFactorAndTrim( albedo );												// Convert downward flux to upward flux
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_GroundPointDiffuseOptical::ConfigureOptical, Error configuring the optical component of the ground point. Thats a problem. Setting initial ground signals to 0 but this may cause problems later on");
		m_FupMS = -999999.0;																					// and set up the outbound radiance pointer immediately
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_GroundPointDiffuseOptical::CalculateOutboundJ		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/


bool SKTRANSO_GroundPointDiffuseOptical::UpdateDiffuseUpwardFlux()
{
	m_FupMS = m_diffuseupwardflux.Evaluate();		// Total upward flux = solar single scattered plus diffuse
	NXASSERT((m_FupMS >= 0.0 ));
	return true;
}
