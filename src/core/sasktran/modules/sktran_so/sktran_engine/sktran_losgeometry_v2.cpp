#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorage_InternalJindex::SKTRANSO_RayStorage_InternalJindex		2007-11-9*/
/** Constructor for the line of sight geometry
**/
/*---------------------------------------------------------------------------*/


SKTRANSO_RayStorage_InternalJindex::SKTRANSO_RayStorage_InternalJindex()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorage_InternalJindex::~SKTRANSO_RayStorage_InternalJindex		2007-11-9*/
/** Destructor
 *	Release  the spatial grid and the array of line of sight elements.
**/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayStorage_InternalJindex::~SKTRANSO_RayStorage_InternalJindex()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorage_InternalJindex::CreateSingleScatterJIndexTables		2010-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayStorage_InternalJindex::CreateSingleScatterJIndexTables( SKTRANSO_Quadrature_TLS_V21*	quadrature, const SKTRANSO_RayInternalGeometry* ray)
{
	bool	ok1;
	bool	ok2;
	bool	ok3;
	bool	ok4;
	bool	groundishit;

	ok1 = !m_singlescatterJ.IsEmpty();
	if (!ok1) ok1 = quadrature->CreateJIndexTable_AtmosphericSingleScatter( ray, &m_singlescatterJ );
	

	groundishit = ray->Storage()->GroundIsHit();
	ok2 = !groundishit;
	if (!ok2)
	{
		ok2 = !m_groundpointSSJ.IsEmpty();
		if (!ok2) ok2 = quadrature->CreateJIndexTable_GroundSingleScatter( ray, &m_groundpointSSJ );
	}
	else
	{
		m_groundpointSSJ.Clear();
		ok2 = true;
	}

	ok3 = !m_emissionJ.IsEmpty();
	if (!ok3) ok3 = quadrature->CreateJIndexTable_AtmosphericEmissions( ray, &m_emissionJ);

	ok4 = !groundishit;
	if (!ok4)
	{
		ok4 = !m_groundppointemissionJ.IsEmpty();
		if (!ok4) ok4 = quadrature->CreateJIndexTable_GroundEmissions( ray, &m_groundppointemissionJ);
	}
	return ok1 && ok2 && ok3 && ok4;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorage_InternalJindex::CreateMultipleScatterJIndexTables		2010-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayStorage_InternalJindex::CreateMultipleScatterJIndexTables( SKTRANSO_Quadrature_TLS_V21*	quadrature,
																		    const SKTRANSO_RayInternalGeometry*   ray,
																			bool							singlescatter)
{
	bool	ok;
	bool	ok1;
	bool	ok2;
	
	NXTRACE_ONCEONLY(firsttime, ("SKTRANSO_RayStorage_InternalJindex::CreateMultipleScatterJIndexTables, need to implement high altitude truncation\n"));
	if (!singlescatter)
	{
		ok1 = !m_diffuseJ.IsEmpty();
		if (!ok1) ok1 = quadrature->CreateJIndexTable_AtmosphericDiffuseScatter( ray, &m_diffuseJ );

		if ( ray->Storage()->GroundIsHit() )
		{
			ok2 =  !m_groundpointMSJ.IsEmpty();
			if (!ok2) ok2 = quadrature->CreateJIndexTable_InterpolateGroundDiffuseScatter( ray, &m_groundpointMSJ);
		}
		else
		{
			m_groundpointMSJ.Clear();
			ok2 = true;
		}
		ok = ok1 && ok2;
	}
	else
	{
		ok = true;
//		m_diffuseJ.Clear();				// Only a call to configureModel should clear these guys once they are initializes
//		m_groundpointMSJ.Clear();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorage_InternalJindex::CreateJIndexTables		2008-2-8*/
/** Initialize the source function index tables for rays in the diffuse points
 *	table. This set of rays usually use the 
**/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_RayStorage_InternalJindex::CreateJIndexTables( SKTRANSO_Quadrature_TLS_V21*	quadrature, const SKTRANSO_RayInternalGeometry* ray, bool singlescatter)
{
	bool ok;

	ok  =       CreateSingleScatterJIndexTables  ( quadrature, ray );					// Configure the solar scattering table if required
	ok  = ok && CreateMultipleScatterJIndexTables( quadrature, ray, singlescatter );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_RayStorage_InternalJindex::LinkJIndexTables, Error creating the JIndex tablkes for the ray");
		m_singlescatterJ.Clear();
		m_groundpointSSJ.Clear();
		m_diffuseJ.Clear();
		m_groundpointMSJ.Clear();
	}
	return ok;
}




