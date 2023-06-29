#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefOutboundUnitSphereME100::SKTRAN_GridDefOutboundUnitSphereME100		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphereLatLonGrid::SKTRAN_UnitSphereLatLonGrid(  )
{
	m_zenith  = NULL;
	m_azieven = NULL;
	m_aziodd  = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::~SKTRAN_UnitSphereLatLonGrid		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphereLatLonGrid::~SKTRAN_UnitSphereLatLonGrid()
{
	ReleaseGrids();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::IsGroundPoint		2011-2-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::IsGroundPoint() const
{
	bool isground;

	isground = (m_zenith != NULL);
	isground = isground && m_zenith->IsGroundPoint();
	return isground;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::ReleaseGrids		2011-1-6*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_UnitSphereLatLonGrid::ReleaseGrids()
{
	if (m_zenith  != NULL) m_zenith ->Release();
	if (m_azieven != NULL) m_azieven->Release();
	if (m_aziodd  != NULL) m_aziodd ->Release();
	m_zenith  = NULL;
	m_azieven = NULL;
	m_aziodd  = NULL;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::IsPole		2011-1-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_UnitSphereLatLonGrid::IsPole( double angle) const
{
	double deps = 1.0E-10;
	return ( fabs(angle) < deps) || ( fabs(angle-nxmath::Pi) < deps);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::AllocateInternalZenithGrid		2011-1-6*/
/** Allocates the internal zenith angle grid for the unit sphere. The code
 *	handles ground point hemispheres as well as full sphere conditions.
 *	The biggest problem we need to deal with is that the user may send in
 *	points that exclude the pole points. This is okay for incoming
 *	rays but is problematic for interpolating the outgoing unit sphere.
 *	Thus I force the Lat/Lon Grid to place a point at both the upward and downward pole.
 *
 *	For ground points I force the code to include an upward pole but I also
 *	demand the disribution include a zenith angle exactly at the equator.
 *
 *	The inclusion of points at the poles (or equator for ground points) 
 *	helps guarantee that interpolations and integrals over the sphere are
 *	properly bounded and normalized.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::AllocateInternalZenithGrid(const SKTRAN_GridDefDiffuseIncomingZenith_V21* userzenith)
{
	double				minzen;
	double				maxzen;
	double				lastzenval;
	double				lastzenrad;
	size_t				numzen;
	bool				ok;
	bool				isgroundpoint;
	nx1dArray<double>	zenbuffer;
	size_t				idx;
	size_t				i;

	NXASSERT(( m_zenith == NULL));
	NXTRACE_ONCEONLY(firsttime, ("**** SKTRAN_UnitSphereLatLonGrid::AllocateInternalZenithGrid, check incoming angles are in degrees. Change accordingly\n"));

	isgroundpoint = userzenith->IsGroundPoint();
	minzen        = (userzenith->front());
	maxzen        = (userzenith->back());		
	lastzenval    = isgroundpoint ? 90 : 180;
	lastzenrad    = isgroundpoint ? (nxmath::Pi*0.5 -1.0E-07) : nxmath::Pi;	// FOr ground points dont place last zenith exactly at 90 degrees as it creates round off issues

	ok =    ( maxzen > 3.5 )										// check that we are not using radians (0 to pi)														
		 && ( minzen >= 0.0 )										// check that we have sensible range of values
		 && ( maxzen <= 180)										// for the min and max. Note that groundpoints will automatically filter out downward facing rays
		 && ( minzen < maxzen);										// check that the data are in ascending zenith angle. 
	ok = ok && zenbuffer.SetSize( userzenith->NumAngles() + 2);		// Create a buffer to temporarily hold  the zenith angles.
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AllocateInternalZenithGrid, The zenith angles are not defined proeprly. They must be in ascending order, in degrees 0-180");
	}	
	else
	{
		if (isgroundpoint)														// If we are dealing with a ground point upward facing hemisphere so filter out downward pointing rays.
		{																		// so
			idx = userzenith->NumAngles();										// Get the number of zenith angles sent in in by the user
			do																	// scan backwards
			{																	// until we have only the upward
				idx--;															// or horizontal points
			} while (userzenith->At(idx) > lastzenval);							// in the upper hemisphere
			numzen = idx + 1;													// adjust the number of useful zenith points.
		}																		// and that is that
		else																	// otherwise not a ground point so its a whole sphere
		{																		// simply
			numzen   = userzenith->NumAngles();									// Get the number of zenith angles sent in in by the user
		}																		// and that is that

		idx = 0;																// Copy the user zenith angles to our internal buffer
		if (userzenith->At(0) > 0.0) zenbuffer.At(idx++) = 0.0;					// Make sure we have the upward pole  always inserted
		for (i =0; i < numzen; i++)												// Now copy over the useful range of zenith angles
		{																		// that have been set by the user
			zenbuffer.At(idx++) = nxmath::DegreesToRadians( userzenith->At(i) );
		}
		if ( userzenith->back() < lastzenval) zenbuffer.At(idx++) = lastzenrad;	// make sure we have downward pole (or equator for ground points) the last bounding zenith point in place.
		zenbuffer.TrimSize(idx);
	}

	if (ok)																		// IF we are good so far
	{																			// then
		m_zenith = new SKTRAN_GridDefDiffuseIncomingZenith_V21;					// create the grid of zenith angles that we shal use internally
		ok = (m_zenith != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AllocateInternalZenithGrid, Error allocating grid object for internal zenith angles");
		}
		else
		{
			m_zenith->AddRef();
			m_zenith->SetIsGroundPoint(isgroundpoint);
			ok = m_zenith->CopyGridArray( zenbuffer.UnsafeArrayBasePtr(), zenbuffer.size() );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AllocateInternalZenithGrid, Error copying zeniths angles into zenith object");
				m_zenith->Release();
				m_zenith = NULL;
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::AssignCubatureWeights		2011-1-31*/
/** Assigns the cubature weights to each point. By this stage we have an internal
 *	array of zeniths that start at 0 radians and finishes exactly on pi or 2 pi
 *	depending upon whether the we have a ground point or an atmospheric point. The
 *	azimuth grid has been padded with an extra element at the beginning
 *	and end to ensure we completely cover the 2 pi radiancs of azimuth.
 *
 *	Most of the vertices are laid out as triangle strips although the points
 *	around the two pole positions are laid out as triangle fans. 
 *

 *       0----|-1----2-----3-----4-----5-----6----|-7--  EVEN ROW AZIMUTHS
 *        \####/\    /\    /\  A /\ C  /\    /\   |/\  
 *         \##/##\  /  \  /  \  /  \  /  \  /  \  /  \ 
 *          \/####\/    \/    \/  B \/	  \/    \/|   \
 *       ---0-|---1-----2-----3--*--4-----5-----6-|---7	 ODD ROW AZIMUTHS
 *
 *
 *
 *       ---0-|---1-----2-----3--*--4-----5-----6-|---7	 ODD ROW AZIMUTHS
 *         / \####/\    /\    /\  B /\    /\    /\|   /
 *        /   \##/##\  /  \  /  \  /  \  /  \  /  \  / 
 *       /    |\/####\/    \/ A  \/ C  \/    \/   |\/   
 *       0----|-1----2-----3-----4-----5-----6----|-7--	 EVEN ROW AZIMUTHS



 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::AssignCubatureWeights()
{
	size_t	numazi;
	size_t	vindextop;			// vertex index at the start of the top    row
	size_t	vindexbot;			// vertex index at the start of the bottom row 
	size_t	iazip1;
	size_t	izen;
	size_t	iazi;
	double	dazi;
	double	mutop;				// cos(zenith) on top row
	double	mubot;				// cos(zenith) on bottom row
	double	dmu;				// The delta mu from one azimuth row to the next
	double	h;					// Half height of triangle (used in area calculation).
	double  triarea;
	size_t	a1;
	size_t	a2;
	size_t	a3;
	size_t	b1;
	size_t	b2;
	size_t	b3;

	#if defined(NXDEBUG)
	double	totalarea = 0.0;
	#endif

	numazi = m_azieven->NumAngles() - 2;									// Get the number of actual vertices in each azimuth row, exclude the padded elements
	mubot  = cos(m_zenith->At(0));											// get the first
	NXASSERT((mubot > 0.9999999));											// mubot should be 1.0 exactly
	vindexbot = ZenithIndexToVertexIndex(0);
	for (izen = 0; izen < (m_zenith->NumAngles()-1); izen++)
	{
		mutop = mubot;
		mubot = cos( m_zenith->At(izen+1) );
		dmu   = (mutop - mubot);
		h     = 0.5*dmu;
		vindextop = izen;
		vindexbot = izen+1;

		for (iazi = 0; iazi < numazi; iazi++)
		{
			iazip1  = (iazi + 1);																// get the index of the next point
			dazi    = m_azieven->At(iazip1) - m_azieven->At(iazi);							// Get the change in azimuth, index over the extra entry at front of array
			triarea = (h*dazi)/3;																// Each triangle contributes exactly one third of its area to the weight of each vertex in our grid.
		 
			a1 = ZenAziIndexToVertexIndex( vindextop, iazi); 
			a2 = ZenAziIndexToVertexIndex( vindextop, iazip1);
			a3 = ZenAziIndexToVertexIndex( vindexbot, iazi) ;

			CubatureWeightAtVar( a1  ) += triarea;			// Set the weights for the top half triangular strip
			CubatureWeightAtVar( a2  ) += triarea;
			CubatureWeightAtVar( a3  ) += triarea;

			b1 = a3; //ZenAziIndexToVertexIndex( vindexbot, iazi);
			b2 = ZenAziIndexToVertexIndex( vindexbot, iazip1);
			b3 = a2; //ZenAziIndexToVertexIndex( vindextop, iazip1);

			CubatureWeightAtVar( b1 ) += triarea;			// Set the weights for the bottom half triangular strip
			CubatureWeightAtVar( b2 ) += triarea;
			CubatureWeightAtVar( b3 ) += triarea;

			#if defined(NXDEBUG)
			totalarea += triarea*6.0;
			#endif
		}
	}
	#if defined(NXDEBUG)
	NXASSERT(  ((mubot > -1.00001) && (mubot < -0.9999999) && (totalarea > 12.56637) && (totalarea < 12.56638) ) || ((mubot > -0.00001) && (mubot < 0.00001)&& (totalarea > 6.283) && (totalarea < 6.284)) );
	#endif
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::ZenithIndexToVertexIndex		2011-2-1*/
/** Returns the offset from the start of the vertex array to the start of a
 *	strip of vertices at the specifed zenithith index. This has the gotcha that
 *	The first zenith angle is forced to be the UP pole and there is only one
 *	vertex placed at the pole.
 **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_UnitSphereLatLonGrid::ZenithIndexToVertexIndex( size_t zenidx ) const
{
	size_t	index;

	if (zenidx > 0)
	{
		index = 1 + (zenidx - 1)*NumAzimuthVertex();
	}
	else
	{
		index = 0;
	}
	return index;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::NumZenith		2011-2-9*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_UnitSphereLatLonGrid::NumZenith() const
{
	return m_zenith->NumAngles();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::GetZenithVertexIndexAndNumAzimuth		2011-2-9*/
/** Get the starting vertex index for this zenith index as well as the number of
 *	azimuth points at this zenith angle. This is used by the diffuse optical point code
 *	to process all o fthe incoming rays at a given zenith angle, which can be optimized
 *	for homogenous atmospheres.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::GetZenithVertexIndexAndNumVertexAzimuth	( size_t zenidx, size_t* firstvertexidx, size_t* numazi ) const
{
	bool	ok;

	ok = (m_zenith != NULL) && ( zenidx < m_zenith->NumAngles() );
	if (ok)
	{

		*firstvertexidx = ZenithIndexToVertexIndex	( zenidx );
		if (IsPole(m_zenith->At(zenidx))) *numazi = 1;
		else                              *numazi = NumAzimuthVertex();
	}
	else
	{
		*firstvertexidx = 0;
		*numazi = 0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::AziIndexToVertexOffset		2011-2-1*/
/** Converts the index in the azimuth array to a vertex offset. The sneaky part
 *	is that the azimuth array has 2 extra points, one at either end, to help
 *	with wrap around issues.
 **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_UnitSphereLatLonGrid::AziIndexToVertexOffset( size_t aziidx ) const
{
	size_t	index;
	size_t	numaziv;

	numaziv =  NumAzimuthVertex();
	if ((aziidx > 0) && (aziidx <= numaziv))	// If we are not dealing with the two extra end points
	{											// then we can just do the conversion
		index = aziidx - 1;
	}
	else if (aziidx == 0)						// if we the first point in the azimuth angle array
	{											// then this maps
		index = numaziv - 1;					// to the last point in the vertex strip aroun dthe azimuth
	}											// and
	else if (aziidx == (numaziv+1))				// if we have the last point in the azimuth angle array
	{											// then
		index = 0;								// this maps to the first point in the vertex strip.
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AziIndexToVertexOffset, the incoming azimuth index (%u) is out of the range 0 to %u. Thats a problem!", (unsigned int)aziidx, (unsigned int)(numaziv+1) );
		index = 0;
	}
	return index;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::ZenAziIndexToVertexIndex		2011-2-1*/
/** **/
/*---------------------------------------------------------------------------*/

 size_t SKTRAN_UnitSphereLatLonGrid::ZenAziIndexToVertexIndex( size_t zenidx, size_t aziidx ) const
 {
	 size_t	vertex;
	 bool   ispole;

	 vertex  = ZenithIndexToVertexIndex( zenidx );
	 ispole  =   (zenidx == 0) || ( (zenidx == (m_zenith->NumAngles()-1) && !IsGroundPoint()));
	 if (!ispole)
	 {
		 vertex += AziIndexToVertexOffset( aziidx );
	 }
	 NXASSERT( vertex < NumUnitVectors() );
	 return vertex;
 }



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::InterpTriangleFan		2011-2-4*/
/**
 *	This is the code used to interpolate the uppe rand lower triangle fans that
 *	surround the top and bottom pole. The fans start out looking like a pizza pie
 *	which i sthen unrolled. The unrolling operation actually converts the triangle
 *	slices into rectangles with the pole vertex on one edge and the first row of
 *	azimuths on the other edge. We can then just write down the linear interpolation
 *	within the rectangle. The north pole and south pole use the same technique but
 *	have slightly different formulations that I have explicitly written out.
 *
 *	North Pole Fomulation
 *	----------------------
 *
 *                      VPOLE  VPOLE
 *	     POLE    0------1--------2-----3-----4-----5-----6------7--  NORTH POLE
 *		                |        |
 *	                    |        |
 *	                    |     V  |
 *	                    |        |
 *	    zenidx+1 0----- 1-----^--2-----3-----4-----5-----6-------7	 ODD ROW AZIMUTHS
 *                      V1    |  v2              
 *                            |             
 *                            V3
 *	V3 = (1-f)V1 + f.V2
 *  V  = (1-g)VPOLE + g(V3)
 *	V  = g(1-f)V1 + gf.V2 + (1-g)VPOLE
 * 
 *	South Pole formulation
 *  ----------------------
 *
 *
 *                     V1        V2
 *	     POLE    0------1----V3--2-----3-----4-----5-----6------7-- ODD/ OR EVEN ROW AZIMUTHS
 *		                |        |
 *	                    |        |
 *	                    |     V  |
 *	                    |        |
 *	    zenidx+1 0----- 1-----^--2-----3-----4-----5-----6-------7	SOUTH POLE
 *                      VPOLE    VPOLE              
 *                           
 *	V3 = (1-f)V1 + f.V2
 *  V  = (g)VPOLE + (1-g)(V3)
 *	V  = (1-g)(1-f)V1 + (1-g)f.V2 + (g)VPOLE

**/
/*---------------------------------------------------------------------------*/

 bool SKTRAN_UnitSphereLatLonGrid::InterpTriangleFan( double									zen,
													  double									azi,
													  size_t									zenidx1,
													  SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	azigrid,
													  trianglevertexstruct*						thegoodtriangle,
													  double*									weights) const
{
	size_t	aziidx1;
	bool	ok;
	size_t	numazi;
	double	x1;
	double	x2;
	double	y1;
	double	y2;
	double	f;
	double	g;

	numazi   = azigrid->NumAngles();
	ok       = azigrid->IndexOfPointBelowOrEqual( azi, &aziidx1 );							// interpolate the azimuth grid
	if (aziidx1 == (numazi-1)) aziidx1--;													// make sure we are in bounds
	if (ok)
	{
		x1 = azigrid->At(aziidx1);
		x2 = azigrid->At(aziidx1+1);
		y1 = m_zenith->At(zenidx1);
		y2 = m_zenith->At(zenidx1+1);

		f = (azi-x1)/(x2-x1);
		g = (zen-y1)/(y2-y1);


		if (zenidx1 == 0)										// If this is the "North pole" we are deal ing with 
		{
			weights[0] = (1-f)*g;
			weights[1] = f*g;
			weights[2] = 1.0-g;

			thegoodtriangle[0].zen     = zenidx1+1;
			thegoodtriangle[0].azi     = aziidx1;
			thegoodtriangle[0].azigrid = azigrid;

			thegoodtriangle[1].zen     = zenidx1+1;
			thegoodtriangle[1].azi     = aziidx1+1;
			thegoodtriangle[1].azigrid = azigrid;

			thegoodtriangle[2].zen     = zenidx1;
			thegoodtriangle[2].azi     = aziidx1;
			thegoodtriangle[2].azigrid = azigrid;
		}
		else													// Otherwise we are dealing with the South Pole
		{
			weights[0] = (1-f)*(1-g);
			weights[1] = f*(1-g);
			weights[2] = g;

			thegoodtriangle[0].zen     = zenidx1;
			thegoodtriangle[0].azi     = aziidx1;
			thegoodtriangle[0].azigrid = azigrid;

			thegoodtriangle[1].zen     = zenidx1;
			thegoodtriangle[1].azi     = aziidx1+1;
			thegoodtriangle[1].azigrid = azigrid;

			thegoodtriangle[2].zen     = zenidx1+1;
			thegoodtriangle[2].azi     = aziidx1;
			thegoodtriangle[2].azigrid = azigrid;
		}
	}
	else
	{
		weights[0] = 0;
		weights[1] = 0;
		weights[2] = 0;
	}
	return ok;
 }



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::IsInsideTriangle		2011-2-1*/
/** Returns true if the point (zen, azi) is inside the triangle defined by
 *	the 3 vertics in zenazivertices. Note that the code is optimized so that the first two
 *	points are assumed to be at the same Y coordinate.
 *      _     _     _
 *	Let r = a.v2 +b.v3			                                            ........eqn 1.
 *        _     
 *	where v2 = vector from point 1 to point 2   = ( x2-x1).i			    ........eqn 2. (Note y2 == y1)
 *	where v3 = vector from point 1 to point 3   = ( x3-x1).i + (y3-y1).j    ........eqn 3.
 *        r  = vector from point 1 to our point = (x-x1)i + (y-y1)j			........eqn 4.
 *
 *	Therefore it can be shown by equating eqn 1 to eqn 4 that
 *  b = (y-y1)/(y3-y1)
 *  a = ((x-x1) - b(x3-x1))/(x2-x1)
 *
 *	Let the three vertics have a value v1, v2 and v3 then we can say in a manner similar to eqn 1
 *
 *	v = v1 + a.(v2-v1) + b.(v3-v1)
 *
 *  or
 *
 *	v = (1-a-b).v1 + a.v2 + b.v3
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::IsInsideTriangle( double zen, double azi, const trianglevertexstruct* zenazivertices, trianglevertexstruct* thetriangle, double* weights) const
{
	double  x;
	double  y;
	double  x1;
	double  y1;
	double  x2;
	double  y2;
	double  x3;
	double  y3;
	double	a;
	double	b;
	double  ab;
	bool	isinside;

	x = azi;
	y = zen;

	x1 = zenazivertices[0].azigrid->At( zenazivertices[0].azi );				// Get the azimuth of vertex 1
	x2 = zenazivertices[1].azigrid->At( zenazivertices[1].azi );				// Get the azimuth of vertex 2
	x3 = zenazivertices[2].azigrid->At( zenazivertices[2].azi );				// Get the azimuth of vertex 3.

	y1 = m_zenith->At( zenazivertices[0].zen);									// Get the zenith of vertex 1
	y2 = m_zenith->At( zenazivertices[1].zen);									// Get the zenith of vertex 2
	y3 = m_zenith->At( zenazivertices[2].zen);									// Get the zenith of vertex 3

	NXASSERT(( y2 == y1 ));														// Make sure we have the first two point on the same zenith coordinate

	b        = (y-y1)/( y3-y1);													// Get the b and a coefficients
	a        = ((x-x1) - b*(x3-x1))/(x2-x1);									// for linear interpolation. Both should be between 0 and 1
	ab       = a + b;															// Get the sum it should also be in range 0 to 1
	isinside = (a >= -0.000001) && ( b >= -0.000001) && (ab <= 1.000001);		// See if we are inside or so close that its good enough
	if (isinside)																// if we are
	{																			// then get the weights to apply for linear interpolation
		weights[0] = 1.0 - ab;													// of values located at the vertices.
		weights[1] = a;
		weights[2] = b;

		thetriangle[0] = zenazivertices[0];
		thetriangle[1] = zenazivertices[1];
		thetriangle[2] = zenazivertices[2];
	}
	return isinside;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::FindInsideTriangle		2011-2-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::FindInsideTriangle( double zen, double azi, const trianglevertexstruct* triangleptr[3], trianglevertexstruct* thegoodtriangle, double* weights) const
{
	size_t	i = 0;
	bool	isinside = false;

	while ((i < 3) && (!isinside))
	{
		isinside = IsInsideTriangle( zen, azi, triangleptr[i], thegoodtriangle, weights );	// See if our point is inside this guy
		i++;
	}
	NXASSERT(isinside);
	if (!isinside)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::FindInsideTriangle, The point( zen = %e, azi = %e) was not in the 3 selected triangles, Thats not good", (double)zen, (double)azi);
	}
	return isinside;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::EvenUpperTriangleStrip		2011-2-2*/
/** This is the triangular interpolation used for triangular strips with the
 *	upper strip (at lower zenith value) uses the even azimuths and the
 *	lower strip (at higher zenith angles) uses the odd azimuths. zenidx and
 *	aziidx are the indexes in to the zenith and even azimuth arrays of the
 *	value les sthan equal to zen and azi.
 *
 *                       azimuth of point
 *                              |
 *                        aziidx|  aziidx+1
 *                           |  |  |
 *                           |VVVVV|
 *	     zenidx  0----|-1----2-----3-----4-----5-----6----|-7--  EVEN ROW AZIMUTHS
 *		        > \   |/\    /\#B##/\    /\    /\    /\   |/\  
 *	zenith of -->  \  /	 \  /A#\##/#C\  /  \  /  \  /  \  /  \ 
 *	point       >   \/|	  \/  ##\/##  \/    \/	  \/    \/|   \
 *	    zenidx+1 ---0-|---1-----2-----3-----4-----5-----6-|---7	 ODD ROW AZIMUTHS
 *                        |     |     |
 *                        |     |     |
 *                 aziidx-1   aziidx  aziidx+1                 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::InterpEvenUpperTriangleStrip( double zen, double azi, size_t zenidx, size_t aziidx, trianglevertexstruct* thegoodtriangle, double* weights) const
{
	bool ok;

	trianglevertexstruct triangleverticesB[3] = {  {zenidx,   aziidx,   m_azieven},			// There is a 50:50 chance the point is inside this triangle
												   {zenidx,   aziidx+1, m_azieven},
												   {zenidx+1, aziidx,   m_aziodd}
	                                            };
	
	trianglevertexstruct triangleverticesA[3] = {  {zenidx+1,   aziidx-1, m_aziodd},
												   {zenidx+1,   aziidx,   m_aziodd},
												   {zenidx,     aziidx,   m_azieven}
											    };

	trianglevertexstruct triangleverticesC[3] = { {zenidx+1,    aziidx,   m_aziodd},
												  {zenidx+1,    aziidx+1, m_aziodd},
												  {zenidx,      aziidx+1, m_azieven}
												};

	const trianglevertexstruct*	triangleptr[3]      =  { triangleverticesB, triangleverticesA, triangleverticesC};

	ok = FindInsideTriangle( zen, azi, triangleptr, thegoodtriangle, weights);
	return ok;
}

/*----------------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::InterpOddUpperTriangleStrip 2011-2-2*/
/** This is the triangular interpolation used for triangular strips with the
 *	upper strip (at lower zenith value) uses the odd azimuths and the
 *	lower strip (at higher zenith angles) uses the even azimuths. zenidx and
 *	aziidx are the indexes in the zenith and even azimuth arrays of the
 *	value less than equal to zen and azi.
 *
 *
 *                 aziidx-1   aziidx  aziidx+1                 
 *                        |     |     |
 *                        |     |     |
 *		 zenidx  ---0-|---1-----2-----3-----4-----5-----6-|---7	 ODD ROW AZIMUTHS
 *		        >  / \|   /\  ##/\##  /\    /\    /\    /\|   /
 *	zenith of --> /   \  /  \A#/##\#C/  \  /  \  /  \  /  \  / 
 *	point       >/    |\/    \/#B##\/    \/   \/    \/    |\/   
 *	   zenidx+1  0----|-1----2-----3-----4-----5-----6----|-7--  EVEN ROW AZIMUTHS
 *                           |^^|^^|
 *                           |  |  |
 *                        aziidx|  aziidx+1
 *                              |
 *                       azimuth of point

 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::InterpOddUpperTriangleStrip( double zen, double azi, size_t zenidx, size_t aziidx, trianglevertexstruct* thegoodtriangle, double* weights) const
{
	bool	ok;


	trianglevertexstruct triangleverticesB[3] = {  {zenidx+1,   aziidx,   m_azieven},			// There is a 50:50 chance the point is inside this triangle
												   {zenidx+1,   aziidx+1, m_azieven},
												   {zenidx,     aziidx,   m_aziodd}
	                                            };
	
	trianglevertexstruct triangleverticesA[3] = {  {zenidx,     aziidx-1, m_aziodd},
												   {zenidx,     aziidx,   m_aziodd},
												   {zenidx+1,   aziidx,   m_azieven}
											    };

	trianglevertexstruct triangleverticesC[3] = { {zenidx,    aziidx,     m_aziodd},
												  {zenidx,    aziidx+1,   m_aziodd},
												  {zenidx+1,  aziidx+1,   m_azieven}
												};

	const trianglevertexstruct*	triangleptr[3]    =  { triangleverticesB, triangleverticesA, triangleverticesC};


	ok = FindInsideTriangle( zen, azi, triangleptr, thegoodtriangle, weights);
	return ok;
}

 /*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::Triangulate		2011-1-31*/
/** Triangulates a unit vector onto our unit sphere. The entire unit sphere is
 *	completely covered in triangles so our point must be inside at one of them.
 *	The ordered structure of the triangles lets us quickly narrow the search
 *	down to 1 of 3 triangles.
 *
 *            0 radians                         2pi radians
 *            |					                  |
 *       0----|-1----2-----3-----4-----5-----6----|-7--  EVEN ROW AZIMUTHS
 *        \   |/\    /\    /\  A /\ C  /\    /\   |/\  
 *         \  /	 \  /  \  /  \  /  \  /  \  /  \  /  \ 
 *          \/|	  \/    \/    \/  B \/	  \/    \/|   \
 *       ---0-|---1-----2-----3--*--4-----5-----6-|---7	 ODD ROW AZIMUTHS
 *         / \|   /\    /\    /\  B /\    /\    /\|   /
 *        /   \  /  \  /  \  /  \  /  \  /  \  /  \  / 
 *       /    |\/    \/    \/ A  \/ C  \/    \/   |\/   
 *       0----|-1----2-----3-----4-----5-----6----|-7--	 EVEN ROW AZIMUTHS
 *
 *	The first stage decides whether the top row is an even azimuth or an odd azimuth,
 *	(we also check for the possibility of upper or lower poles). We then find out the
 *	lower bounding index in the even azimuth distribution, in our example the lower
 *	bounding azimuth is 3. 
 *
 *	There are three triangles, A, B and C which a point could lie in. The point could be
 *	in triangle B but also might be triangle A or C.  We check all three until we find
 *	one that works.  The triangle checking function also calculates the weights of
 *	the 3 vertices  if the point is inside.
 *	of
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const
{
	SKTRAN_GridIndex		zenidx1;
	SKTRAN_GridIndex		aziidx1;
	size_t					numzen;
	size_t					numazi;
	size_t					i;
	bool					isodd;
	double					zen;
	double					azi;
	bool					ok;
//	bool					isinside = false;
	trianglevertexstruct	thegoodtriangle[3];
	SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	azigrid;




	ok = ( (maxvertices >= 3) && (m_zenith != NULL));
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphereLatLonGrid::Triangulate, Undefined lat/lon grid or insufficient vertex buffer space passed in by caller. The incoming weights and index buffer must support at least 3 vertex interpolation points");
	}
	else
	{
		numzen = m_zenith ->NumAngles();
		ok =       LocalLookToAziZen( unit, &azi, &zen );
		ok = ok && m_zenith ->IndexOfPointBelowOrEqual( zen, &zenidx1 );
		if (zenidx1 == (numzen-1)) zenidx1--;												// Watch out for points user has right on the bottom zenith (e.g. at pole)

		if (ok)
		{
			if (zenidx1 == 0)																			// If the upper bound is the pole
			{																							// the we do a triangle fan
				ok = InterpTriangleFan( zen, azi, zenidx1, m_aziodd, thegoodtriangle, unit_weightptr);	// and interpolate the triangle fan
			}
			else if ( (zenidx1 == (numzen-2)) && !IsGroundPoint() )										// otherwise if this the "south Pole"
			{																							// then
				isodd    = (zenidx1%2 ) != 0;															// See if the last row befor ethe pole is odd or even
				azigrid  = isodd ? m_aziodd : m_azieven;												// and get the correct grid
				ok = ok && InterpTriangleFan( zen, azi, zenidx1, azigrid, thegoodtriangle, unit_weightptr);	// and interpolate the triangle fan
			}
			else																					// otherwise we are in the triangle strip region
			{																						// so
				numazi = m_azieven->NumAngles();
				isodd = (zenidx1%2 ) != 0;															// See if the first row is odd
				ok   = m_azieven->IndexOfPointBelowOrEqual( azi, &aziidx1 );						// and get the bound azimuths form the even coords
				if (aziidx1 == (numazi-1)) aziidx1--;
				if (ok)
				{
					if (isodd)																		// If the first row is odd thens econd row is even
					{																				// This point(zen, azi) can be in 1 of only 3 triangle
						ok = InterpOddUpperTriangleStrip ( zen, azi, zenidx1, aziidx1, thegoodtriangle, unit_weightptr);
					}
					else
					{
						ok = InterpEvenUpperTriangleStrip( zen, azi, zenidx1, aziidx1, thegoodtriangle, unit_weightptr);
					}
				}
			}
			if (ok)
			{
				for (i = 0; i < 3; i++)
				{
					unit_indexptr[i] = ZenAziIndexToVertexIndex( thegoodtriangle[i].zen, thegoodtriangle[i].azi );
				}
				for (i = 3; i < maxvertices ; i++)
				{
					unit_indexptr[i]  = 0;
					unit_weightptr[i] = 0.0;
				}
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::Triangulate, The LAT/Lon triangulation failed. That is not good as goofy things are likely to happen. Sasktran code needs debugging");
		for (i = 0; i < maxvertices ; i++)
		{
			unit_indexptr[i]  = 0;
			unit_weightptr[i] = 0.0;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::AziInRange0To360		2011-1-6*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_UnitSphereLatLonGrid::AziInRange0To360( double val )
{
	NXASSERT(( (val >= -180) && (val <= 360) ));

	if (val < 0) val += 360.0;
	return val;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::AllocateInternalAzimuthGrid		2011-1-6*/
/** Allocates the azimuth grid. This is composed of an odd and even grid so we
 *	interlace the points across the unit sphere. The last azimuth point is
 *	prepended to the list to make the table span the full two pi.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::AllocateInternalAzimuthGrid(const SKTRAN_GridDefDiffuseIncomingAzimuth_V21* userazi)
{
	double				minazi;
	double				maxazi;
	bool				ok;
	size_t				i;
	size_t				numazi;

	NXASSERT(( m_azieven == NULL));
	NXASSERT(( m_aziodd  == NULL));
	NXTRACE_ONCEONLY(firsttime, ("**** SKTRAN_UnitSphereLatLonGrid::AllocateInternalAzimuthGrid, check incoming angles are in degrees. Change accordingly\n"));

	minazi        = userazi->front();
	maxazi        = userazi->back();
	numazi        = userazi->NumAngles();

	ok =    ( maxazi > 6.3   )										// check that we are not using radians (0 to 2pi)														
		 && ( minazi >= -180 )										// check that we have sensible range of values
		 && ( maxazi <= 360  )										// for the min and max. Note that groundpoints will automatically filter out downward facing rays
		 && ( minazi < maxazi);										// check that the data are in ascending zenith angle. 
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AllocateInternalAzimuthGrid, The azimuth angles are not defined proeprly. They must be in ascending order, in degrees -180 to +360 ");
	}	
	else
	{
		m_azieven = new SKTRAN_GridDefDiffuseIncomingAzimuth_V21;					// create the grid of zenith angles that we shal use internally
		m_aziodd  = new SKTRAN_GridDefDiffuseIncomingAzimuth_V21;					// create the grid of zenith angles that we shal use internally
		ok = ( (m_azieven != NULL) && (m_aziodd != NULL));
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AllocateInternalAzimuthGrid, Error allocating grid object for internal azimuth angles");
		}
		else
		{
			m_azieven->AddRef();
			m_aziodd ->AddRef();
			ok =       m_azieven->AllocateGridArray( numazi + 2 );
			ok = ok && m_aziodd ->AllocateGridArray( numazi + 2 );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::AllocateInternalAzimuthGrid, Error allocating interal storage for azimuth angles");
			}
			else
			{
				for (i = 0; i <  numazi; i++ )																	// Copy over the  azimuth grid
				{																								// note that we  leave the first entry blank
					m_azieven->AtVar(i+1) = nxmath::DegreesToRadians( AziInRange0To360(userazi->At(i)));		// copy over the azimuths, convert to 0 ... 2 Pi
				}																								// do all of the measurements
				m_azieven->AtVar( 0 )        = m_azieven->At(numazi) - nxmath::TWOPI;									// prepend the front with the last entry to the first point
				m_azieven->AtVar( numazi+1 ) = m_azieven->At(1)      + nxmath::TWOPI;									// prepend the front with the last entry to the first point

				for (i = 0; i <  numazi+1; i++ )																	// now do the odd points
				{																								// by placing them halfway between
					m_aziodd->AtVar(i) = 0.5*( m_azieven->At(i) + m_azieven->At(i+1));							// the even points
				}																								// and
				m_aziodd->AtVar( numazi+1 ) = m_aziodd->At(1)      + nxmath::TWOPI;								// prepend the front with the last entry to the first point
			}
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::CreateVertexGrid		2011-1-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::CreateVertexGrid()
{
	size_t		numvertices;
	size_t		numazi;
	size_t		nazi;
	size_t		i;
	size_t		iazi;
	double		zen;
	double		azi;
	double		coszen;
	double		sinzen;
	size_t		idx;
	bool		ispole;
	bool		ok;
	bool		isodd;
	SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	azimuthgrid;

	numvertices = 0;
	numazi      = m_azieven->NumAngles() - 2;

	for (i = 0; i < m_zenith->NumAngles(); i++ )					// The first stage is simply to count the number of vertics
	{																// so loop over all of the zenith angles
		ispole       = IsPole( m_zenith->At(i) );					// If this is a pole
		nazi         = ispole ? 1 : numazi;							// then we only need 1 vertex otherwise we need all numazi vertices
		numvertices += nazi;										// sum up the total number of vertices
	}																// loop over all of the zeniths

	ok = AllocateVertices( numvertices );							// allocate the vertices
	if (ok)															// if that worked then
	{																// we are going to assign the directions  to eachvertex
		idx = 0;
		for (i = 0; i < m_zenith->NumAngles(); i++ )				// loop over all of the zenith angles
		{															// for each zenith angle
			zen          = m_zenith->At(i);							// Get the zenith angle in radians.
			coszen       = cos(zen);								// Get the cosine of the zenith angle
			sinzen       = sin(zen);								// and the size of the zenith angle.
			ispole       = IsPole( zen );							// find out if we are at the pole.
			isodd        = ((i%2) != 0);							// work out if this is an odd or even row
			azimuthgrid  = isodd  ? m_aziodd : m_azieven;			// get the odd or even azimuth angles appropriate for this row
			nazi         = ispole ? 1 : numazi;						// work out how many vertices are in this row.

			for (iazi = 0; iazi < nazi; iazi++)						// for all of the vertics
			{														// get the unit vector
				azi = azimuthgrid->At( iazi+1 );					// index the azimuth table, skip over the wrap around point at the front of the buffer
				UnitVectorAtVar(idx).SetCoords( cos(azi)*sinzen, sin(azi)*sinzen, coszen );
				idx++;
			}
		}
		NXASSERT(( idx == numvertices));
		ok = (idx == numvertices);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphereLatLonGrid::CreateVertexGrid, Error creating thevertex grid. Thats a problem");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid::DefineGrid		2008-4-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_UnitSphereLatLonGrid::DefineGrid( const SKTRAN_GridDefDiffuseIncomingZenith_V21* userzenith, const SKTRAN_GridDefDiffuseIncomingAzimuth_V21* userazimuth )
{
	bool	ok;

	ReleaseGrids();
	ok   =       AllocateInternalZenithGrid ( userzenith );
	ok   = ok && AllocateInternalAzimuthGrid( userazimuth);
	ok   = ok && CreateVertexGrid           ( );
	ok   = ok && AssignCubatureWeights		( );
	
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_UnitSphereLatLonGrid::Constructor, Error making latlon grids. Thats not good");
		ReleaseGrids();
	}
	return ok;
}
