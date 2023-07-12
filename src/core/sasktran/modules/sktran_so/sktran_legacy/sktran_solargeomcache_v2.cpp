#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"
#include "sktran_legacy_internals.h"
/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::SKTRAN_TableSolarTransmission_2D_Height_SZA		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmission_2D_Height_SZA::SKTRAN_TableSolarTransmission_2D_Height_SZA()
{
	m_altitudegrid   = NULL;
	m_szagrid        = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::~SKTRAN_TableSolarTransmission_2D_Height_SZA		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableSolarTransmission_2D_Height_SZA::~SKTRAN_TableSolarTransmission_2D_Height_SZA()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableSolarTransmission_2D_Height_SZA::ReleaseResources()
{
	m_altitudegrid   = NULL;
	m_szagrid        = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::FindingBoundingSZAIndices		2008-2-5*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_TableSolarTransmission_2D_Height_SZA::FindingBoundingLocationIndices	(  const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const
{
	NXASSERT((  (location.CosSZA() >= -1.0) && (location.CosSZA() <= 1.0) ));
	NXASSERT((  (m_szagrid != NULL)  ));
	NXASSERT(( maxvertices >= 2 ));
	return m_szagrid->FindingBoundingIndices( location.CosSZA(), SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, index, weight, maxvertices);
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::FindingBoundingAltitudeIndices		2008-1-15*/
/** Finds the indices and linear interpolating weights for the altitudes
 *	that bound the given radius.  If the radius is outside the altitude
 *	range of the grid it returns zero points contributing to the given radius.
 **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_TableSolarTransmission_2D_Height_SZA::FindingBoundingAltitudeIndices(  const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const
{
	NXASSERT((  (m_altitudegrid != NULL)  ));
	NXASSERT(( maxvertices >= 2 ));
	return m_altitudegrid->FindingBoundingIndices( location.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ZERO, index, weight, maxvertices);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::InterpolateSolarTransmissionTable		2010-2-19*/
/** Interpolates the solar transmission table and returns a set of indices that
 *	can be used to evaluate the solar transmission at arbitrary radius and
 *	angle.   The algorithm uses linear interpolation.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmission_2D_Height_SZA::InterpolateTable( const HELIODETIC_POINT& location, const HELIODETIC_UNITVECTOR& look, SKTRANSO_JIndex* vertexdescriptortable, size_t maxpts, size_t* npts, double weight ) const
{
	SKTRAN_GridIndex								szaindex	[12];	// The indices of the diffuse profiles that bound the solar zenith angle, 2 is sufficient for original Sasktran, 12 should be sufficient for 3D geometry
	SKTRAN_GridIndex								radiusindex	[12];	// The indices of the radii in each diffuse profile that bound the altitude
	double											szaweight   [12];	// The weights for each of the sza indices 
	double											radiusweight[12];	// The weights for each of the radii
	size_t											numszavertex;		// The number of solar zenith angle vertices in this interpolation
	size_t											numradvertex;		// The number of radius vertices in this interpolation
	size_t											numvertex;									
	bool											ok;
	size_t											szaCtr;
	size_t											radCtr;
	size_t											idx;

	NXASSERT(( vertexdescriptortable != NULL ));
	NXASSERT(( npts != NULL ));
	NXASSERT(( location.Altitude() < 250000.0 ));

	numszavertex  = FindingBoundingLocationIndices( location,  szaindex,    szaweight,    N_ELEMENTS(szaindex)    );
	numradvertex  = FindingBoundingAltitudeIndices( location,  radiusindex, radiusweight, N_ELEMENTS(radiusindex) );

	numvertex     = (numszavertex*numradvertex);										// Get the number of vertices requires
	ok            = ( numvertex <= maxpts );														// It is possible there are no bounding altitude vertices if the diffuse profile array
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableSolarTransmission_2D_Height_SZA::InterpolateSolarTransmissionTable, the user supplied buffer is too small (%d), Its needs to be at least %d vertices", (int)maxpts, (int)numvertex);
	}
	else
	{
		idx           = 0;
		for (szaCtr = 0; szaCtr < numszavertex; szaCtr++ )										// For each of the potential solar zenith angles
		{																						// solar zenith angle
			for (radCtr = 0; radCtr < numradvertex; ++radCtr)
			{
				vertexdescriptortable[idx++].ConfigureSolarTransmissionTableIndex( szaindex[szaCtr], radiusindex[radCtr], szaweight[szaCtr], radiusweight[radCtr], 0, weight );
			}
		}
		*npts = numvertex;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableSolarTransmission_2D_Height_SZA::InterpolateSolarTransmissionTable, error interpolating Solar Transmission Table");
		*npts     = 0;
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry		2010-4-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableSolarTransmission_2D_Height_SZA::ConfigureGeometry( const SKTRAN_SpecsInternal_V21* specs, SKTRAN_ThreadManager* threadmanager )
{
	bool									ok;

	ok = SKTRANSO_TableSolarTransmission::ConfigureGeometry( specs, threadmanager );
	m_altitudegrid = ProfileHeightsPtr(0);
	m_szagrid      = SZAGrid();
	return ok;
}


