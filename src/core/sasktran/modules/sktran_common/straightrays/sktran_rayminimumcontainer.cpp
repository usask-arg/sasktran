#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Base::SKTRAN_RayStorage_Base		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayStorage_Base::SKTRAN_RayStorage_Base( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )
{
	m_coords = coords;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Base::~SKTRAN_RayStorage_Base		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayStorage_Base::~SKTRAN_RayStorage_Base( )
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayStorage_Straight::SKTRAN_RayStorage_Straight	(  std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) 
	                       : SKTRAN_RayStorage_Base(coords) 
{
	m_Robs   = std::numeric_limits<double>::quiet_NaN();
	m_Tobs   = std::numeric_limits<double>::quiet_NaN();
	m_Rt	 = std::numeric_limits<double>::quiet_NaN();

	m_geoidradius = coords->AltitudeToRadius(0.0);
}


/*-----------------------------------------------------------------------------
 *					~SKTRAN_RayStorage_Straight		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayStorage_Straight:: ~SKTRAN_RayStorage_Straight()
{
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Base::CalculateStraightLineTangentPoint		2014-1-26*/
/** This is a helper function that calculates the radius of the observer,
 *	the straight line tangent point radius and straight line distance from
 *	the observer to the tangent point. The code is used by both the straight
 *	line ray classes and the piecewise linear class. The code may slightly
 *	adjust the observer position if he is close to the ground. 
 *
 *	\param groundaltitude
 *		The altitude of the ray tracing ground above true ground. The code 
 *		ensures the observer if above or at this altitude. Returns NaN if failure.
 *
 *	\param Robs
 *		returns the radius of the observer.
 *
 *	\param Tobs
 *		returns the straight line distance from the observer to the tangent point.
 *		Assumes that all rays are straight.  Returns NaN if failure.
 *
 *	\param Rt
 *		returns the radius of the tangent point. This can be below the ground. Returns NaN if failure.
 *
 *	\returns
 *		True if successful.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::CalculateTangentPointDetails()
{
	HELIODETIC_POINT	pt;	
	bool				ok;
	double				coszenith;
	double				minradius;
	double				diff;

	NXASSERT((GetCoordsPtr() != nullptr ));
	m_Robs = GetObserver().Magnitude();
	minradius    = GetCoordsPtr()->AltitudeToRadius( GetCoordsPtr()->GroundAltitude() );
	if (m_Robs <= minradius)
	{		
		HELIODETIC_VECTOR	observer = GetObserver();
		diff       = (minradius- (m_Robs));								// then look for minor floating point differences 
		if (diff > 0.002)												// if we are outside of a roundoff error domain
		{																// thenlog an error to the user
			nxLog::Record(NXLOG_WARNING,"SKTRAN_RayGeometry_Base::CalculateBaseLineTangentPointDetails, the observer at (%e, %e, %e) was substantially below the RT ground (%e meters). It should be with 0.0001 meters. This indicates a problem in configuration", GetObserver().X(), GetObserver().Y(), GetObserver().Z(), (double)diff);
		}
		minradius  += 0.001;											// add on a tenth of a millimeter.
		observer.SetCoords( observer.UnitVector(), minradius);			// Fix it so observer is just above the ground
		m_observer = observer;
		m_Robs = minradius;
		GetCoordsPtr()->HelioVectorToHelioPoint(observer, &pt);
	}

	ok    = GetCoordsPtr()->HelioVectorToHelioPoint( GetObserver(), &pt );								// get the observer heliodetic point
	coszenith = pt.CosZenithAngle( LookVector() );															// and its coszenith
	m_Rt   = (  coszenith*coszenith <= 0.99999999999999 ? (m_Robs)*sqrt( 1.0 - coszenith*coszenith ) : 0 );	// Get the radius of the tangent point
	m_Tobs =  -1.0*(m_Robs)*coszenith;																		// and the tangent distance to the observer
//	NXASSERT(( m_Tobs >= 0.0 ));
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::InitializeObserver		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::InitializeObserver( const HELIODETIC_VECTOR&   observer, const HELIODETIC_UNITVECTOR& look )
{
	m_observer = observer;
	m_lookaway = look;
	m_looktowards.SetCoords( -look.X(), -look.Y(), -look.Z() );
	CalculateTangentPointDetails();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::ReleaseResources		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_RayStorage_Straight::ReleaseResources()
{
	m_distancefromorigin.clear();
	m_distFromTan.clear();
	m_radii.clear();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::PushBack		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::PushBack( SKTRAN_Distance r, SKTRAN_Distance distFromTan, SKTRAN_Distance s)
{
	double	alt;

	//NXASSERT( (distFromTan >= 0.0) );							// distFromTan can be negative. Maybe need to check this out
	NXTRACE_ONCEONLY( firsta, ("SKTRAN_RayStorage_Straight::PushBack, need to check out that negative distance to tangent point is ok\n"));
	alt = GetCoordsPtr()->RadiusToAltitude( r );
	m_distancefromorigin.push_back(s);
	m_radii.           push_back( r );
	m_distFromTan.     push_back( distFromTan );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::Insert		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::Insert( SKTRAN_Distance r, SKTRAN_Distance distFromTan, SKTRAN_Distance s, size_t index)
{
	double	alt;

	NXASSERT( (distFromTan >= 0.0) );
	alt = GetCoordsPtr()->RadiusToAltitude( r );
	m_distancefromorigin.insert( m_distancefromorigin.begin()+index, s);
	m_radii.           insert( m_radii.begin() + index, r );
	m_distFromTan.     insert( m_distFromTan.begin() + index, distFromTan );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::Reserve		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::Reserve( size_t numquadraturepoints )
{

	m_distancefromorigin.reserve( numquadraturepoints);
	m_radii.           reserve( numquadraturepoints);
	m_distFromTan.     reserve( numquadraturepoints );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::Resize		 2014- 12- 2*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::Resize( size_t numquadraturepoints )
{

	m_distancefromorigin.resize( numquadraturepoints);
	m_radii.           resize  ( numquadraturepoints);
	m_distFromTan.     resize  ( numquadraturepoints);
	return true;
}

bool SKTRAN_RayStorage_Straight::SplitCell(size_t cellindex)
{
	double s0 = DistanceOfPointFromOrigin(cellindex);
	double s1 = DistanceOfPointFromOrigin(cellindex + 1);
	double s = (s0 + s1) / 2;

	double t0 = DistanceOfPointFromCellTangentPoint(cellindex, cellindex);
	double t1 = DistanceOfPointFromCellTangentPoint(cellindex + 1, cellindex);
	double t = (t0 + t1) / 2;

	double r1 = RadiusOfPoint(cellindex + 1);
	double r = std::sqrt(t*t - t1 * t1 + r1 * r1);

	return Insert(r, t, s, cellindex + 1);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::EndIntercept		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

//SKTRAN_Distance SKTRAN_RayStorage_Straight::EndIntercept() const
//{
//	SKTRAN_Distance	endpt;
//
//	if ( NumCells() > 0 ) endpt =  m_distancefromorigin[NumCells()];
//	else                  endpt = 0;
//	return endpt;
//}
/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::EndPoint		 2014- 11- 26*/
/** **/
/*---------------------------------------------------------------------------*/
//
//HELIODETIC_VECTOR SKTRAN_RayStorage_Straight::EndPoint() const
//{
//	HELIODETIC_VECTOR	endpt;
//
//	if ( NumCells() > 0 ) LocationAlongRayAsVector( m_distancefromorigin[NumCells()], &endpt );
//	else                  endpt = m_observer;
//	return endpt;
//}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Base::GetCellIntercepts		 2014- 11- 26*/
/** **/
/*---------------------------------------------------------------------------*/

//bool SKTRAN_RayStorage_Straight::GetCellIntercepts( size_t cellidx, SKTRAN_Distance* startintercept, SKTRAN_Distance* exitintercept ) const
//{
//	*startintercept = m_distancefromorigin[cellidx];
//	*exitintercept  = m_distancefromorigin[cellidx+1];
//	return true;
//}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::CellLength		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayStorage_Straight::CellLength( size_t cellindex ) const
{
	double	s0;
	double	s1;
	double	length;

	s0 = DistanceOfPointFromOrigin( cellindex );
	s1 = DistanceOfPointFromOrigin( cellindex + 1 );
	length = s1-s0;	
	return length;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::CellMidPoint		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::CellMidPoint( size_t cellindex, HELIODETIC_POINT* pt) const
{
	double	s0;
	double	s1;
	double	s;

	s0 = DistanceOfPointFromOrigin( cellindex );
	s1 = DistanceOfPointFromOrigin( cellindex + 1 );
	s = (s0 + s1)/2.0;
	
	HELIODETIC_VECTOR	v(m_lookaway, s );
	v += m_observer;
	pt->FromVector( v, GetCoordsPtr() );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayBaseGeometry_V21::LocationAlongRayAsVector		2008-1-29*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayStorage_Straight::LocationOfPoint( size_t quadraturepoint_index, HELIODETIC_POINT* pt ) const
{
	double					s;
	HELIODETIC_VECTOR		v;
	HELIODETIC_UNITVECTOR	unit;
	double					r;
	double					f;

	s = m_distancefromorigin.at( quadraturepoint_index);
	v.SetCoords( m_lookaway, s);
	v += m_observer;
	// r = v.Magnitude();
	r = RadiusOfPoint(quadraturepoint_index);
	f = (1.0/r);
	unit.SetCoords( v.X()*f, v.Y()*f, v.Z()*f );
	//pt->Initialize( unit, r, GetCoordsPtr() );
	pt->InitializeFromRaw(unit, r, GetCoordsPtr()->RadiusToAltitude(r));
//	NXASSERT( (pt->Altitude() < 150000.0) );
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::TruncateToNumElements		 2014- 11- 27*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_RayStorage_Straight::TruncateToNumElements  ( size_t numels )
{
	m_distancefromorigin.resize( std::min( numels, m_distancefromorigin.size() ) );
	m_radii.resize(           std::min( numels, m_radii.size() ) );
	m_distFromTan.resize(     std::min( numels, m_distFromTan.size() ) );
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::RadiusOfQuadraturePoint		 2014- 11- 28*/
/** **/
/*---------------------------------------------------------------------------*/

double	 SKTRAN_RayStorage_Straight::RadiusOfPoint( size_t quadraturepoint_index ) const
{
	return m_radii.at(quadraturepoint_index);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::ClearStorage		 2014- 11- 28*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_RayStorage_Straight::ClearStorage()
{
	m_distancefromorigin.clear();
	m_radii.clear();
	m_distFromTan.clear();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::DistanceOfPointFromOrigin		 2014- 12- 4*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayStorage_Straight::DistanceOfPointFromOrigin( size_t quadraturepoint_index ) const 
{
	return m_distancefromorigin[quadraturepoint_index];
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::DistanceFromCellTangentPoint		 2014- 11- 28*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayStorage_Straight::DistanceOfPointFromCellTangentPoint( size_t quadraturepoint_index, size_t raysegment_index ) const
{
	return m_distFromTan.at(quadraturepoint_index);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight::AltitudeOfPoint		 2014- 12- 2*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayStorage_Straight::AltitudeOfPoint( size_t quadraturepoint_index ) const
{
	return GetCoordsPtr()->RadiusToAltitude( m_radii.at(quadraturepoint_index) );
}


