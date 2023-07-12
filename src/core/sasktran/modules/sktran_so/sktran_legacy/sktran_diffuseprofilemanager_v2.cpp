#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"
#include "sktran_legacy_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::SKTRAN_TableDiffusePoints_2D_Height_SZA		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableDiffusePoints_2D_Height_SZA::SKTRAN_TableDiffusePoints_2D_Height_SZA()
{
	m_radii = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::~SKTRAN_TableDiffusePoints_2D_Height_SZA		2007-12-13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableDiffusePoints_2D_Height_SZA::~SKTRAN_TableDiffusePoints_2D_Height_SZA()
{
	if (m_radii != NULL) m_radii->Release();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::FindingBoundingSZAProfiles		2007-12-19*/
/** Retrieve the indexes and linear interpolation weights of the diffuse
 *	profiles that bound a given solar zenith angle. Solar zenith angles out
 *	of range of the table are truncated to the end points.  The code returns the
 *	number of diffuse profiles that contribute to the linear interpolation. This will be
 *	either 1 or 2.
 *
 **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_TableDiffusePoints_2D_Height_SZA::FindingBoundingLocationIndices	(  const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const
{
	size_t	numvertex;

	NXASSERT((  (location.CosSZA() >= -1.0) && (location.CosSZA() <= 1.0) ));
	NXASSERT(( maxvertices >= 2 ));
	numvertex = CosSZAGrid()->FindingBoundingIndices( location.CosSZA(), SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, index, weight, maxvertices);
	return numvertex;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::FindingBoundingAltitudeIndices		2008-1-15*/
/** Finds the indices and linear interpolating weights for the altitudes
 *	that bound the given radius.  If the radius is outside the altitude
 *	range of the grid it returns zero points contributing to the given radius.
 **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_TableDiffusePoints_2D_Height_SZA::FindingBoundingAltitudeIndices( const HELIODETIC_POINT& location, SKTRAN_GridIndex* index, double* weight, size_t maxvertices) const
{
	size_t	numvertex;
	NXASSERT((  (m_radii != NULL)  ));
	NXASSERT(( maxvertices >= 2 ));
	numvertex = m_radii->FindingBoundingIndices(  location.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ZERO, index, weight, maxvertices);
	return numvertex;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_TableDiffusePoints::InterpolateTable		2007-12-20*/
/** Finds the vertices in the diffuse table that can be used to evaluate the
 *	source function at a location, solar zenith angle, altitude and
 *	geographic direction look.
 *
 *	This object uses internal storage to store the vertices calculated from this
 *	subroutine and returns a pointer to the internal buffer. It is upto the user
 *	to copy them from the internal buffer into their own storage.  This code is
 *	intrinsically not multi-thread safe and must be used in a single, geometry
 *	configuration, thread of the SASKTRAN Model.
 *
 *	\param cossza
 *	The cosine of the solar zenith angle at the desired radiance location
 *
 *	\param radius
 *	The altitude of the desired location, expressed as a radius. If the radius is negative
 *	then this implies we are looking up a ground point and we use the first point in
 *	the 
 *
 *	\param look
 *	The unit vector direction of the desired radiance expressed in geographic
 *	coordinates (w.r.t. osculating sphere)
 *
 *	\param vertexdescriptortable
 *	The vertexdescriptortable.  If this is NULL no vertex descriptors are written
 *	otherwise descriptors are written.  The user must make sure there is sufficient space
 *	as indicated in parameter maxvertex.  It will throw an exception if there is not enough space.

 *	\returns
 *	The number of vertices that contribute to the interpolated radiance,
 *	this is typically 0,1,2 or 4.
 **/
/*---------------------------------------------------------------------------*/


bool SKTRAN_TableDiffusePoints_2D_Height_SZA::InterpolateTable( const HELIODETIC_POINT& location, const HELIODETIC_UNITVECTOR& look, SKTRANSO_JIndex* vertexdescriptortable, size_t maxpts, size_t* npts, double weight ) const
{
	SKTRAN_GridIndex							profileidx	[12];	// The indices of the diffuse profiles that bound the solar zenith angle, 2 is enough for normal Sasktran, 12 should be sufficient for more complex interpolation schemes.
	SKTRAN_GridIndex							radiusindex	[12];	// The indices of the radii in each diffuse profile that bound the altitude
	double										szaweight   [12];	// The weights for each of the sza indices 
	double										radiusweight[12];	// The weights for each of the radii
	SKTRAN_GridIndex							vertexindex	[3];	// The unit vector triangulation indices for the second solar zenith angle
	double										vertexweight[3];	// The weights for the triangulation indices of the first profile  (same used for upper and lower radii)
	size_t										numprofiles;		// The number of solar zenith angle vertices in this interpolation
	size_t										numradvertex;		// The number of radius vertices in this interpolation
	size_t										numvertex;									
	bool										ok;
	size_t										profileCtr;
	size_t										radCtr;
	HELIODETIC_UNITVECTOR						locallook;
	HELIODETIC_UNITVECTOR						localunitvectors[3];
	nxVector									locallookv;
	size_t										idx;
//	double										zen, azi;
	const SKTRAN_UnitSphere_V2*					unitsphere;

	NXASSERT(( weight == 1.0 ));
	NXASSERT(( vertexdescriptortable != NULL ));
	NXASSERT(( npts != NULL ));
//	NXASSERT(( location.Altitude() < 250000.0 ));

	numprofiles   = FindingBoundingLocationIndices( location,  profileidx,  szaweight,    N_ELEMENTS(profileidx)    );
	numradvertex  = FindingBoundingAltitudeIndices( location,  radiusindex, radiusweight, N_ELEMENTS(radiusindex) );
	numvertex     = (numprofiles*numradvertex);											// Calculate the number of vertices
	ok            = (maxpts >= numvertex*3);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableDiffusePoints_2D_Height_SZA::InterpolateTable, Insufficient storage space (%d) provided by user, %d is required. Thats a problem", (int)maxpts, (int)(numvertex*3));
	}
	else
	{

		ok        = location.LocalUnitVectors( localunitvectors, N_ELEMENTS(localunitvectors) );		// Look up the local unit vectors at this point
		locallook = location.TransformToLocalZenithCoords(look, localunitvectors );						// Transform look direction to local coordinates
		locallookv.SetCoords( locallook.X(), locallook.Y(), locallook.Z() );
		idx = 0;
		for (profileCtr = 0; profileCtr < numprofiles; profileCtr++ )																		// For each of the potential solar zenith angles
		{																																	// solar zenith angle
			for (radCtr = 0; radCtr < numradvertex; ++radCtr)
			{
				unitsphere = PointAt( profileidx[profileCtr], radiusindex[radCtr] )->ScatteringMatrix()->OutboundUnitSphere();
				ok = ok && unitsphere->Triangulate( locallookv, vertexindex, vertexweight, N_ELEMENTS(vertexindex) );
				if ( vertexweight[0] != 0.0 ) vertexdescriptortable[idx++].ConfigureDiffuseTableIndex( profileidx[profileCtr], radiusindex[radCtr], szaweight[profileCtr], radiusweight[radCtr], vertexindex[0], vertexweight[0] );//, weight);
				if ( vertexweight[1] != 0.0 ) vertexdescriptortable[idx++].ConfigureDiffuseTableIndex( profileidx[profileCtr], radiusindex[radCtr], szaweight[profileCtr], radiusweight[radCtr], vertexindex[1], vertexweight[1] );//, weight);
				if ( vertexweight[2] != 0.0 ) vertexdescriptortable[idx++].ConfigureDiffuseTableIndex( profileidx[profileCtr], radiusindex[radCtr], szaweight[profileCtr], radiusweight[radCtr], vertexindex[2], vertexweight[2] );//, weight);
			}
		}
		*npts  = idx;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRANSO_TableDiffusePoints::InterpolateDiffuseOutBoundSourceFunctions, error triangulating vertices");
		}
	}
	if (!ok)
	{
		*npts = 0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::ConfigureGeometry		2010-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableDiffusePoints_2D_Height_SZA::ConfigureGeometry( const SKTRAN_SpecsInternal_V21* specs, SKTRAN_ThreadManager* threadmanager )
{
	if (m_radii != NULL ) m_radii->Release();
	m_radii = specs->DiffuseSpecs()->DiffuseHeights(0);
	m_radii->AddRef();

	return SKTRANSO_TableDiffusePoints::ConfigureGeometry( specs, threadmanager );
}



