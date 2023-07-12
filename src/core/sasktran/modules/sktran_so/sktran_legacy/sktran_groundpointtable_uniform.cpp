#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"
#include "sktran_legacy_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::SKTRAN_TableGroundPointDiffuse_Colocated		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableGroundPointDiffuse_Colocated::SKTRAN_TableGroundPointDiffuse_Colocated( const SKTRAN_TableDiffusePoints_2D_Height_SZA* diffusetable )

{
	m_diffusetable              = diffusetable;
	NXASSERT((m_diffusetable != NULL));
	m_diffusetable->AddRef();
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::~SKTRAN_TableGroundPointDiffuse_Colocated		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableGroundPointDiffuse_Colocated::~SKTRAN_TableGroundPointDiffuse_Colocated()
{
	m_diffusetable->Release();
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::NumberOfGroundPointsRequired		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_TableGroundPointDiffuse_Colocated::NumberOfGroundPointsRequired( const SKTRAN_SpecsInternal_V21*	specs  ) const
{
	return m_diffusetable->NumProfiles();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::ConfigureGeometry_Stage2MT		2008-2-8*/
/** This method is called by the multi-threaded worker threads from
 *	threadmanager->GroundPointsTable_ConfigureGeometryStage2. Please be aware that this
 *	is MULTI-THREADED and use caution when modifying, especially global or static
 *	variables.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableGroundPointDiffuse_Colocated::ConfigureGeometry_Stage2MT( size_t pointindex, SKTRANSO_Quadrature_TLS_V21* threadquadrature )
{
	bool											ok;
	SKTRAN_GridIndex								ptidx = (SKTRAN_GridIndex)pointindex;
	
	ok = PointAtVar(pointindex)->ConfigureGeometry( ptidx, m_diffusetable->PointAt( pointindex, 0 ));	// and configure this point
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableGroundPointDiffuse_Colocated::ConfigureGeometry_Stage2, Error configuring ground point geometry");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::FindingBoundingSZAIndices		2008-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

size_t  SKTRAN_TableGroundPointDiffuse_Colocated::FindingBoundingLocationIndices(  const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const
{
	return m_diffusetable->FindingBoundingLocationIndices( location, index, weight, maxvertices);
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableGroundPointDiffuse_Colocated::InterpolateTable		2010-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableGroundPointDiffuse_Colocated::InterpolateTable( const HELIODETIC_POINT&			location,
																   const HELIODETIC_UNITVECTOR&		look,
																   SKTRANSO_JIndex*				vertexdescriptortable, 
																   size_t							maxpts, 
																   size_t*							npts,
																   double							weight ) const
{
	SKTRAN_GridIndex								posindex	[12];	// The indices of the diffuse profiles that bound the solar zenith angle
	double											posweight   [12];	// The weights for each of the sza indices 
	size_t											numvertex;									
	bool											ok;
	size_t											posCtr;

	numvertex      = FindingBoundingLocationIndices( location, posindex,  posweight, N_ELEMENTS(posindex) );
	ok = (numvertex <= maxpts);																	// It is possible there are no bounding altitude vertices if the diffuse profile array
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableGroundPointDiffuse_Colocated::InterpolateTable, more points (%u) than passe din (%u) are required, truncating terms, thats a problem", (int)numvertex, (int)maxpts);
		numvertex = maxpts;
	}
	for (posCtr = 0; posCtr < numvertex; posCtr++ )												// For each of the potential solar zenith angles
	{																							// solar zenith angle
		vertexdescriptortable[posCtr].ConfigureGroundPointTableIndex( posindex[posCtr], posweight[posCtr], 0, weight );
	}
	*npts  = numvertex;
	return ok;
}


