#include "../sktran_common.h"
#include <climits>


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::SKTRAN_UnitSphere_V2		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphere_V2::SKTRAN_UnitSphere_V2	()
{
	m_unitvectors    = NULL;
	m_numunitvectors = 0;

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::~SKTRAN_UnitSphere_V2		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphere_V2::~SKTRAN_UnitSphere_V2	()
{
	ReleaseVertices();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::ReleaseVertices		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_UnitSphere_V2::ReleaseVertices()
{
	if (m_unitvectors != NULL) delete [] m_unitvectors;
	m_weights.erase();
	m_unitvectors = NULL;
	m_numunitvectors = 0;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::AllocateVertices		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_V2::AllocateVertices( size_t numunitvectors)
{
	bool	ok;

	ok = (numunitvectors ==  m_numunitvectors);
	if (!ok)
	{
		ReleaseVertices();
		ok = (numunitvectors == 0);
		if (!ok)
		{
			m_unitvectors = new nxVector[numunitvectors];
			ok =        m_weights.SetSize(numunitvectors);
			ok = ok && (m_unitvectors != NULL);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphere_V2::AllocateVertices, Error allocating %Iu vertices", (size_t)numunitvectors);
			}
			else
			{
				m_numunitvectors = numunitvectors;
				m_weights.SetTo(0.0);
			}
		}
	}
	if (!ok) ReleaseVertices();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::LocalLookToAziZen		2011-1-5*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_V2::LocalLookToAziZen( const nxVector& unit, double* azi, double*zen ) const
{
	double	z;

	z     = unit.Z();
	if      ( z >  1.0) z = 1.0;
	else if ( z < -1.0) z = -1.0;
	*zen    = acos(z);																// Get the zenith angle
	*azi   = atan2(unit.Y(), unit.X());												// get the azimuth angle
	if (*azi < 0.0) *azi += nxmath::TWOPI;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::UnitVectorAt		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector& SKTRAN_UnitSphere_V2::UnitVectorAtVar( size_t idx )
{
	NXASSERT(( idx < m_numunitvectors ));
	return m_unitvectors[idx];
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::At		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

const nxVector& SKTRAN_UnitSphere_V2::UnitVectorAt( size_t idx ) const
{
	NXASSERT(( idx < m_numunitvectors ));
	return m_unitvectors[idx];
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::FindThreeClosestVertices		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_V2::FindThreeClosestVertices( const nxVector& unit, size_t* idxmax1, size_t* idxmax2, size_t* idxmax3, size_t* idxmax4 ) const 
{
	double			dot;
	double			max1 = -10;
	double			max2 = -10;
	double			max3 = -10;
	double			max4 = -10;
	size_t			idx;
	bool			ok;

	*idxmax1 = INT_MAX;
	*idxmax2 = INT_MAX;
	*idxmax3 = INT_MAX;
	*idxmax4 = INT_MAX;
    const size_t numUnitVectors =  NumUnitVectors(); // Loop invariant for MS auto vectorizer
	for( idx=0; idx < numUnitVectors; ++idx )
	{
		dot = UnitVectorAt(idx) & unit;
		if ( dot > max1 )
		{
			max4     = max3;
			max3     = max2;
			max2     = max1;
			max1     = dot;
			*idxmax4 = *idxmax3;
			*idxmax3 = *idxmax2;
			*idxmax2 = *idxmax1;
			*idxmax1 = idx;
		}
		else if ( dot > max2 )
		{
			max4     = max3;
			max3     = max2;
			max2     = dot;
			*idxmax4 = *idxmax3;
			*idxmax3 = *idxmax2;
			*idxmax2 = idx;
		}
		else if ( dot > max3 )
		{
			max4     = max3;
			max3     = dot;
			*idxmax4 = *idxmax3;
			*idxmax3 = idx;
		}
		else if ( dot > max4 )
		{
			max4     = dot;
			*idxmax4 = idx;
		}
	}
	ok = (*idxmax1 < NumUnitVectors()) && (*idxmax2 < NumUnitVectors()) && (*idxmax3 < NumUnitVectors()) && (*idxmax4 < NumUnitVectors());
	NXASSERT(( ok ));
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::InterpolateTrianlge		2008-1-20*/
/** 
 *	Let the trinagle have vertices A,B,C. Let P be the point inside the triangle.
 *	Let ab be the vector from A to B.
 *	let ac be the vector from A to C.
 *	Let ap be the vector from A to P.
 *	We can write the parameterized equation,
 *	ap = t *ab + u *ac
 *	where t and u are >= 0 and <= 1.
 *	If we break the vectors into in orthogonal coordinates x,y in the plane of the triangle
 *	then,
 *
 *	xp = xb*t + xc*u
 *	yp = yb*t + yc*u
 *
 *	We further define x parallel to ab so yb is always zero.
 *	Then the interpolated value v is given by
 *	vp = va + t *(vb - va) + u * (vc - va)
 *	or
 *	vp = (1-t-u)va + (t)vb + (u)vc
 *
 *	Note this analysis also works for all points above or below the plane of
 *	the triangle at point P.
 *
 *	\par Determining the cubature weights from triangular interpolation
 *	The addition of sharply peaked phase matrices has demanded that our
 *	code also supply cubature weights for each vertex to simplify
 *	integration over the unit sphere. I spent a couple of hours working out
 *	that the cubature weight of each vertex of the triangle is 1/3 the area
 *	of the triangle. I.E. we use the interpolation scheme vp = (1-t-u)va + (t)vb + (u)vc 
 *	to interpolate the function at any point inside the triangle and we integrate
 *	over the  area of the triangle. The total integral is then  1/3.A( v1 + v2 + v3) where
 *	the v's are the vertex values and A is the area of the triangle.
 *	
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_V2::InterpolateTriangle( const nxVector& P, size_t* unit_indexptr, double* unit_weightptr ) const
{
	nxVector	ab;
	nxVector	ac;
	nxVector	ap;
	nxVector	xunit;
    nxVector	yunit;
	nxVector	A;
	nxVector	B;
	nxVector	C;
	double		xb;
	double		xc;
	double		yc;
	double		xp;
	double		yp;
	double		u;
	double		t;

	A = m_unitvectors[unit_indexptr[0]];
	B = m_unitvectors[unit_indexptr[1]];
	C = m_unitvectors[unit_indexptr[2]];

	ab  = B-A;
	ac  = C-A;
	ap  = P-A;													// Vector from point A to point P. Note P is above plane of triangle but analysis is still good

	xunit = ab.UnitVector();									// x is parallel to AB
    yunit = ac.ComponentPerpendicularTo(xunit).UnitVector();	// y is perpendicular to AB

	xb = ab & xunit;											// X Coordinate of B (by definition B is along X axis
	xc = ac & xunit;											// X coordinate of C
//  yb = 0.0;													// Y coordinate of B is implicitly 0.0
	yc = ac & yunit;											// Y coordinate of C
	xp = ap & xunit;											// X coordinate of P
	yp = ap & yunit;											// Y coordinate of P
	u  = yp/yc;
	t  = (xp - u*xc)/xb;
	
//	NXTRACE_ONCEONLY(firsttimea,("SKTRAN_UnitSphere_V2::InterpolateTriangle, Might want to use Delaunay traingulation as current system does not get surrounding points\n"));
//	NXASSERT(( (u >= -0.01) && ( t >= -0.0000001) && (u+t < 1.000001) ));

	// Point is inside triangle if u and t are between 0 and 1

	unit_weightptr[0] = 1.0 - u - t;
	unit_weightptr[1] = t;
	unit_weightptr[2] = u;

#if defined(NXDEBUG)
	NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_UnitSphere_V2::InterpolateTriangle, Might want to include test to make sure point is inside a triangle\n"));
	nxVector	PTest;
	nxVector	diff;
	nxVector	zunit;
	nxVector	aps;

	zunit  = xunit ^ yunit;
	aps    = ap.ComponentPerpendicularTo(zunit);				// Get the point P in the plane of the triangle
	PTest  = t*ab + u*ac;
	diff   = (PTest-aps);
	NXASSERT(( diff.Magnitude() < 0.000001));
#endif

	return true;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::Triangulate		2007-12-17*/
/** Finds the three vertices that surrounds the unit sphere.
 *	This needs to be improved. There are lots of potential bugs at moment
 *	(eg we dont currently guarantee three points surrounding the point)
 *	We need a proper triangular interpolation scheme.
 *
 *	The current algorithm calculates the weights from the linear combination 
 *	of the 3 vertex unit vectors required to make the specified unit vector.
 *	Personally I dont see the rationale for this choice.
 *
 *	I would rather see a proper triangular interpolation
 *	eg,
 *
 * Ax0 + By0 + C = z0
 * Ax1 + By1 + C = z1
 * Ax2 + By2 + C = z2
 *
 *  (x0,y0) are the vertex coordinates in 2-D, z0 is the vertex value.
 *	A,BC are the linear interpolants.
 *	define x = 0  y = 0 to be the unit vector direction.
 *	then the interpolated value is  C, which is a linear combination of z0, z1 and z3
 *
 *	/param unit
 *	The outgoing unit vector in local coordinatess.
 *
 *	/param unit_indexptr
 *	An array of at least 3. Upon exit it will contain the indexes of the the 3
 *	rays that are the vertices of the triangle
 *
 *	/param unit_weightptr
 *	An array of at least 3 doubles. Upon exit will contain the weights to be applied
 *	to the radiances at the 3 vertices to estimate the radiance in direction given
 *	by parameter unit.
 *
 *	/param maxvertices
 *	The size of unit_indexptr and unit_weightptr.  Must be at least 3. The first 3 values
 *	are filled and the remainder are set to 0
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_V2::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const
{
	double		test;
	size_t		idxmax1;
	size_t		idxmax2;
	size_t		idxmax3;
	size_t		idxmax4;
	bool		ok;
	size_t		idx;

	for ( idx = 0; idx < maxvertices; idx++)									// First
	{																			// fill up the return
		unit_indexptr[idx]  = 0;												// arrays 
		unit_weightptr[idx] = 0.0;												// with zero as the caller will then ignore them
	}																			// and that is that
	ok = (maxvertices >= 3);													// make sure we have space for at least vertics
	if (!ok)																	// if we failed then
	{																			// log a message
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphere_V2::Triangulate, Insufficient vertices, you need to supply at least 3 vertices");
	}																					// and that is that
	else																				// otherwise
	{																					// we have space
		ok = FindThreeClosestVertices( unit, &idxmax1, &idxmax2, &idxmax3, &idxmax4 );	// find the 3 surrounding points
		if (ok)																			// if that worked as it should
		{																				// then
			test = (UnitVectorAt(idxmax1).Cross( UnitVectorAt(idxmax2) )).Dot( UnitVectorAt(idxmax3)) ;				// See if the first three points are coplanar
			ok   = ( fabs(test) >  1e-2 );												// if they are coplanar 
			if (!ok)																	// then
			{																			// replace the most distant of the first three
				NXTRACE_ONCEONLY(firsttime, ("**** 2008-12-23 ***** SKTRAN_UnitSphere_V2::Triangulate, There are some points on the unit sphere that are planar\n"));
				idxmax3 = idxmax4;														// with the fourth most distant
				test = (UnitVectorAt(idxmax1).Cross( UnitVectorAt(idxmax2) )).Dot( UnitVectorAt(idxmax3)) ;			// See if the the three points are coplanar
				ok   = ( fabs(test) >  1e-2 );											// if they are
				if (!ok)																// then we are up the creek.
				{
					nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphere_V2::Triangulate, four closest points are coplanar arghh...");
				}
			}
			if (ok)																		// we have three points
			{																			// now we need to do the interpolation
				unit_indexptr [0] = idxmax1;											// set up
				unit_indexptr [1] = idxmax2;											// the three indices on
				unit_indexptr [2] = idxmax3;											// our sphere
				ok = InterpolateTriangle( unit, unit_indexptr, unit_weightptr );		// and linearly interpolate
			}							
			if (!ok)																	// if that did not work
			{																			// then log a message
				nxLog::Record( NXLOG_WARNING, "SKTRAN_UnitSphere_V2::Triangulate, error triangulating sphere, setting interpolant to zero.");
				unit_indexptr [0] = 0;													// and return zeros.
				unit_indexptr [1] = 0;
				unit_indexptr [2] = 0;
				unit_weightptr[0] = 0.0;
				unit_weightptr[1] = 0.0;
				unit_weightptr[2] = 0.0;
			}
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::Triangulate		2013-09-18*/
/** Some lookups may be helped by returning an index to the user to be passed
 *  back to Triangulate on the next call, e.g. if two points are close together
 *  they might be bounded by the triangulation. Give the user the option to take
 *  advantage of this. 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_V2::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const
{
	return Triangulate(unit, unit_indexptr, unit_weightptr, maxvertices);
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::CubatureWeightAt		2010-12-14*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_UnitSphere_V2::CubatureWeightAt( size_t idx ) const
{
	return m_weights[idx];
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::CubatureWeightAtVar		2010-10-22*/
/** **/
/*---------------------------------------------------------------------------*/

double& SKTRAN_UnitSphere_V2::CubatureWeightAtVar( size_t idx )
{
	return m_weights[idx];
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphere_WithLookupTable_V2		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphere_WithLookupTable_V2	()
{
	m_numazimuth = 360;								// Num Azimuths resolution for lookup table;
	m_numzenith  = 180;								// Num Zeniths resolution for lookup table;
	m_deltazen   = nxmath::Pi/m_numzenith;			// Lookup table zenith resolution in radians
	m_deltaazi   = nxmath::TWOPI/m_numazimuth;		// Lookup table azimuth resolution in radians
	m_numlookup      = 0;								// The number of lookup entries.
	m_lookupentries  = NULL;							// The look up table array[numazimuth, numzen] = NULL;

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::~SKTRAN_UnitSphere_WithLookupTable_V2		2007-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphere_WithLookupTable_V2::~SKTRAN_UnitSphere_WithLookupTable_V2	()
{
	ReleaseLookupTable();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::AllocateAziZenTable		2008-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_WithLookupTable_V2::AllocateAziZenTable()
{
	bool		ok;
	bool		ok1;
	size_t		outidx;
	size_t		zenidx;
	size_t		aziidx;
	size_t		numlookup;
	SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphereInterpolationEntry*		entry;

	NXASSERT(( m_numazimuth > 12 ));				// Num Azimuths resolution for lookup table;
	NXASSERT(( m_numzenith  > 12 ));				// Num Zeniths resolution for lookup table;
	m_deltazen   = nxmath::Pi/m_numzenith;			// Lookup table zenith resolution in radians
	m_deltaazi   = nxmath::TWOPI/m_numazimuth;		// Lookup table azimuth resolution in radians
	numlookup    = (m_numzenith+1)*m_numazimuth;	// Add a little extra at the end to protect for found off
 
	ok = AllocateLookupTable( numlookup );
	if (ok)
	{
		for (outidx=0; outidx < numlookup; outidx++)
		{
			zenidx = outidx/m_numazimuth;
			aziidx = outidx%m_numazimuth;
			entry  = InterpolationEntryAtVar( outidx );
			ok1    = GenerateLookupEntry( entry,  zenidx, aziidx );
			ok     = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphere_WithLookupTable_V2::AllocateAziZenTable, Error allocating azimuth, zenith, lookup table");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::GenerateLookupEntry		2008-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_WithLookupTable_V2::GenerateLookupEntry( SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphereInterpolationEntry* entry, size_t zenidx, size_t aziidx )
{
	nxVector	unit; 
	double		zenith;
	double		coszen;
	double		sinzen;
	double		azimuth;
	bool		ok;

	NXTRACE_ONCEONLY( firsttime, ("**** Check if Multithreading protection required in SKTRAN_UnitSphere_WithLookupTable_V2::GenerateLookupEntry ***** \n "));
	zenith = (zenidx + 0.5)*m_deltazen;
	coszen = cos(zenith);
	sinzen = sin(zenith);
	azimuth = (aziidx + 0.5)*m_deltaazi;
	unit.SetCoords( cos(azimuth)*sinzen, sin(azimuth)*sinzen, coszen );
	ok      = SKTRAN_UnitSphere_V2::Triangulate( unit, entry->unit_indexptr, entry->unit_weightptr, 3 );
	entry->isvalid = ok;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::InitializeLookupTable		2008-4-28*/
/** Initializes the internal tables of the unit sphere. This class uses a
 *	fairly high resolution azimuth-zenith table to look up the surrounding points.
 *	This is fairly brute force and can probably be optimized by using cleverer
 *	search algorithms (eg corase grids down to fine grids.)
 *
 *	This code is not thread safe and should be invoked in a thread safe manner.
 *	We may have to change this in the future.
**/
/*---------------------------------------------------------------------------*/


bool SKTRAN_UnitSphere_WithLookupTable_V2::InitializeLookupTable()
{
	bool	ok;
	ok = (NumLookupEntries() > 0);
	if (!ok) ok = AllocateAziZenTable();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::Triangulate		2008-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_WithLookupTable_V2::Triangulate( const nxVector& locallook, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const
{
	const SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphereInterpolationEntry* entry;

	size_t	outidx;
	size_t	aziidx;
	size_t	zenidx;
	double	zenith;
	double	azimuth;
	bool	ok;

	LocalLookToAziZen( locallook, &azimuth, & zenith );						// Convert the local look vector into azimuth and zenith
	NXASSERT(( (m_numlookup > 0) && (m_lookupentries != NULL) ));			// Make sure we have a lookup table. otehrwise we need to call initialize (in a thread safe manner) 
	NXASSERT(( (azimuth >= 0.0) && (azimuth <= nxmath::TWOPI) ));
	NXASSERT(( (zenith >= 0.0) && (zenith <= nxmath::Pi) ));
	zenidx = (size_t)(zenith/m_deltazen);
	aziidx = (size_t)(azimuth/m_deltaazi);
	if( aziidx == m_numazimuth)										// djz828 added this code in Version 3 (2013-10-15). We need to wrap indexing around the world. There were pointers going out of bounds in the old code. 
	{																// This actually changes the testintallation answer by few thousandths of a percent
		aziidx = 0;
	}
	outidx = zenidx*m_numazimuth + aziidx;

	entry = InterpolationEntryAt(outidx);
	ok    = entry->isvalid;
	NXASSERT((maxvertices >= 3));
	if (ok)
	{
		unit_indexptr [0] = entry->unit_indexptr [0];
		unit_indexptr [1] = entry->unit_indexptr [1];
		unit_indexptr [2] = entry->unit_indexptr [2];
		unit_weightptr[0] = entry->unit_weightptr[0];
		unit_weightptr[1] = entry->unit_weightptr[1];
		unit_weightptr[2] = entry->unit_weightptr[2];
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::ReleaseaziZenTable		2008-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_UnitSphere_WithLookupTable_V2::ReleaseLookupTable()
{
	if ( m_lookupentries != NULL) delete [] m_lookupentries;
	m_lookupentries = NULL;
	m_numlookup    = 0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_WithLookupTable_V2::AllocateAziZenTable		2008-4-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphere_WithLookupTable_V2::AllocateLookupTable(size_t numlookup)
{
	bool		ok;
	size_t		outidx;

	ok = (numlookup == m_numlookup);
	if (!ok)
	{
		ReleaseLookupTable();
		if ( numlookup > 0 )
		{
			m_lookupentries = new SKTRAN_UnitSphereInterpolationEntry [numlookup];
			ok = (m_lookupentries != NULL);
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphere_V2::AllocateLookupTable, Error allocating lookup table memory");
			numlookup = 0;
		}
		m_numlookup = numlookup;
		
	}
	for (outidx=0; outidx < m_numlookup; outidx++)
	{
		m_lookupentries[outidx].isvalid = false;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::InterpolationEntryAt		2010-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphereInterpolationEntry*	SKTRAN_UnitSphere_WithLookupTable_V2::InterpolationEntryAt( size_t lookupindex ) const
{
	NXASSERT(( lookupindex < m_numlookup));
	return m_lookupentries+lookupindex;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2::InterpolationEntryAtVar		2010-12-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphereInterpolationEntry*	SKTRAN_UnitSphere_WithLookupTable_V2::InterpolationEntryAtVar( size_t lookupindex )
{
	NXASSERT(( lookupindex < m_numlookup));
	return m_lookupentries+lookupindex;
}
