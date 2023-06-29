#include "../sktran_common.h"


bool SKTRAN_GeometryObject_Sphere::EqualTo( const SKTRAN_GeometryObject& rhs ) const
{
	const SKTRAN_GeometryObject_Sphere* rhspointer = dynamic_cast<const SKTRAN_GeometryObject_Sphere*> ( &rhs );

	if( rhspointer )
	{
		// say that shells within 1mm are equal
		return ( abs( m_radius - rhspointer->m_radius ) < 1E-5 );
	}
	else
	{
		return false;
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Sphere::FindIntersections		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

std::array<double,2> SKTRAN_GeometryObject_Sphere::FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const 
{
	std::array<double,2> ret;
	
	NXASSERT( fabs( look.Magnitude() - 1.0 ) < 0.0001 );

	const double obsRadius = positionmagnitude;
	const double lookDotObs = look & position;
	const double discrim = lookDotObs*lookDotObs - obsRadius*obsRadius + m_radius*m_radius;

	if( discrim >= 0.0 )
	{
		ret[0] = -1.0*lookDotObs - sqrt( discrim );
		ret[1] = -1.0*lookDotObs + sqrt( discrim );
	}
	else
	{
		ret[0] = -1.0;
		ret[1] = -1.0;
	}
	
	return ret;
}

bool SKTRAN_GeometryObject_Plane::EqualTo( const SKTRAN_GeometryObject& rhs ) const
{
	const SKTRAN_GeometryObject_Plane* rhspointer = dynamic_cast<const SKTRAN_GeometryObject_Plane*> ( &rhs );

	if( rhspointer )
	{
		return abs( m_normal.ComponentPerpendicularTo( rhspointer->m_normal ).Magnitude() ) < 1E-8;
	}
	else
	{
		return false;
	}

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Plane::FindIntersections		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

std::array<double,2> SKTRAN_GeometryObject_Plane::FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const
{
	std::array<double,2> ret;

	const double normaldotlook = look & m_normal;

	if( abs(normaldotlook) < 1E-12 )
	{
		// no intersection
		ret[0] = -1.0;
		ret[1] = -1.0;
		return ret;
	}

	// there is an intersection
	const double d = -1.0 * ( position & m_normal );
	
	ret[0] = d / normaldotlook;
	ret[1] = -1.0;

	return ret;
}

bool SKTRAN_GeometryObject_Cylinder::EqualTo( const SKTRAN_GeometryObject& rhs ) const
{
	const SKTRAN_GeometryObject_Cylinder* rhspointer = dynamic_cast<const SKTRAN_GeometryObject_Cylinder*> ( &rhs );

	if( rhspointer )
	{
		return ( abs( m_radius - rhspointer->m_radius) < 1E-5 ) && ( m_z.ComponentPerpendicularTo( rhspointer->m_z ).Magnitude() < 1E-8 );
	}
	else
	{
		return false;
	}
}

std::array<double,2> SKTRAN_GeometryObject_Cylinder::FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const
{
	std::array<double,2> ret;

	const double obsrad = position.Magnitude();
	const double obsdotz = m_z & position;
	const double lookdotz = m_z & look;
	const double obsdotlook = position & look;

	const double c = obsrad*obsrad - obsdotz*obsdotz - m_radius*m_radius;
	const double b = -2*obsdotz*lookdotz + 2*obsdotlook;
	const double a = 1 - lookdotz*lookdotz;

	const double discrim = b*b - 4*a*c;

	if( discrim > 0 )
	{
		ret[0] = (-b + sqrt(discrim)) / (2*a);
		ret[1] = (-b - sqrt(discrim)) / (2*a);
	}
	else
	{
		ret[0] = -1;
		ret[1] = -1;
	}
	
	return ret;
}

bool SKTRAN_GeometryObject_Cone::EqualTo( const SKTRAN_GeometryObject& rhs ) const
{
	const SKTRAN_GeometryObject_Cone* rhspointer = dynamic_cast<const SKTRAN_GeometryObject_Cone*> ( &rhs );

	if( rhspointer )
	{
		return (abs( m_halfapexangle - rhspointer->m_halfapexangle) < 1E-8) && ( m_dir.ComponentPerpendicularTo( rhspointer->m_dir ).Magnitude() < 1E-8 );
	}
	else
	{
		return false;
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Cone::FindIntersections		2015-02-09*/
/** Cone is collection of points, c, where | c . u | = cos(theta) / |c|, where
 *  u is the unit vector of the cone and theta is the half apex angle
 **/
/*---------------------------------------------------------------------------*/

std::array<double,2> SKTRAN_GeometryObject_Cone::FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const
{
	std::array<double,2> ret;

	const double obsdotunit  = position.Dot( m_dir );
	const double lookdotunit = look.Dot( m_dir );
	const double obsmag      = positionmagnitude;
	const double costheta    = nxmath::cosd( m_halfapexangle );
	const double obsdotlook  = position.Dot( look );

	const double a = lookdotunit*lookdotunit - costheta*costheta;
	const double b = 2.0*lookdotunit*obsdotunit - 2.0*costheta*costheta*obsdotlook;
	const double c = obsdotunit*obsdotunit - costheta*costheta*obsmag*obsmag;

	const double discrim = b*b - 4*a*c;

	if( discrim < 0 )
	{
		// no intersections
		ret[0] = -1.0;
		ret[1] = -1.0;
		return ret;
	}
	/* else */
	ret[0] = (-1.0*b - sqrt(discrim))/(2.0*a);
	ret[1] = (-1.0*b + sqrt(discrim))/(2.0*a);
	return ret;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Straight_Generic::SKTRAN_RayTracer_Straight_Generic		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracer_Straight_Generic::SKTRAN_RayTracer_Straight_Generic(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords)
	                              :SKTRAN_RayTracer_Shells(coords)
{
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Straight_Generic::Push_BackDistance		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Straight_Generic::PushBack_Distance( double s, SKTRAN_RayStorage_Straight* storage) const
{
	double							tobs;
	double							rt;
	double							r;
	double							t;

	rt   = storage->Rt();
	tobs = storage->Tobs();
	t    = fabs(tobs - s);
	r    = sqrt( rt*rt + t*t);
	storage->PushBack( r,t,s);
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Straight_Generic::TraceRay		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Straight_Generic::TraceStraightRay( SKTRAN_RayOptical_Straight*	straightray ) const
{
	//TODO: use storage from the ray or thread manager to avoid reallocs

	bool							ok;
	SKTRAN_RayStorage_Straight*		storage        = straightray->StraightStorageVar();
	HELIODETIC_UNITVECTOR			lookhelio;
	HELIODETIC_VECTOR				obshelio;
	std::vector<double>*			trajdistances;
	std::array<double,2>			groundintersect; 
	std::array<double,2>			atmointersect;
	double							atmodistancetoend;
	double							atmodistancetostart;
	bool							groundhit;
	double							grounddistance;
	double							obsradius;

//	straightray = dynamic_cast<SKTRAN_RayOptical_Straight*>( aray );					// Make sure the ray coming in looks like a structure we can handle
	ok   =     (straightray != nullptr);
	ok = ok && m_threadtrajectory.LookupUpThreadData( &trajdistances );
	if (ok)
	{
		trajdistances->resize(0);
		trajdistances->reserve( m_container.size() * 2 );
		lookhelio = storage->AverageLookVectorAwayFromObserver( 0);
		obshelio  = storage->GetObserver();

		nxVector look( lookhelio.X(), lookhelio.Y(), lookhelio.Z() );
		nxVector obs( obshelio.X(), obshelio.Y(), obshelio.Z() );
		obsradius = obs.Magnitude();

		groundintersect = m_earth.FindIntersections( look, obs, obsradius );
		atmointersect   = m_upperatmo.FindIntersections( look, obs, obsradius );
		atmodistancetoend    = std::max( atmointersect[0], atmointersect[1] );
		atmodistancetostart  = std::min( atmointersect[0], atmointersect[1] ) - 1E-6;
//		if( atmodistancetoend < 0 ) atmodistancetoend = 1E20;
		
		groundhit      = false;
		grounddistance = 9999999999;

		if( groundintersect[0] > 0.0 || groundintersect[1] > 0.0 )
		{
			groundhit = true;
			std::pair<double,double> t = std::minmax( groundintersect[0], groundintersect[1] );
			if( t.first < 0 )
				grounddistance = t.second;
			else
				grounddistance = t.first;
			grounddistance = grounddistance - 1E-6;	 // subtract off a small number to help with rounding
			trajdistances->push_back( grounddistance );
		}
	
		if( obshelio.Magnitude() < m_upperatmo.Radius() + 1.0)
		{
			// in atmosphere
			trajdistances->push_back( 0.0 );
		}
		for( const auto& ele : m_container )
		{
			std::array<double,2> intersect = ele->FindIntersections( look, obs, obsradius );
			for( const auto& s : intersect )
			{
				if( (s > 0.0) && (!groundhit || s < grounddistance) && s <= atmodistancetoend && s >= atmodistancetostart )
				{
					trajdistances->push_back( s );
				}
			}
		}
		if( storage->Tobs() > 0.0 && storage->Tobs() < grounddistance )
		{
			trajdistances->push_back( storage->Tobs() );
		}

		std::sort( trajdistances->begin(), trajdistances->end() );
		
		storage->Resize(0);
		storage->Reserve( trajdistances->size() );
		storage->SetGroundIsHit( groundhit );
		for( const auto& s : *trajdistances )
		{
			PushBack_Distance(s, storage);
		}

		// set up the geometry endpoints
//		if( storage->NumCells() > 0 )
//		{
//			straightray->GRAY()->SetDistanceToStart( storage->DistanceOfPointFromOrigin( 0 ) );
//			straightray->GRAY()->SetDistanceToEnd  ( storage->DistanceOfPointFromOrigin( storage->NumCells() ) );
//		}
//		else
//		{
//			straightray->GRAY()->SetDistanceToStart( 0 );
//			straightray->GRAY()->SetDistanceToEnd  ( 0 );
//		}	

	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Straight_Generic::CreateRay		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/
//bool SKTRAN_RayTracer_Straight_Generic::TraceRay( const HELIODETIC_VECTOR&	        origin,
//												  const HELIODETIC_UNITVECTOR&		look,
//												  SKTRAN_RayOptical_Base*			ray) const
//
//// bool SKTRAN_RayTracer_Straight_Generic::CreateRay( const SKTRAN_CoordinateTransform_V2*			coords,
////													const HELIODETIC_VECTOR&						pt, 
////													const HELIODETIC_UNITVECTOR&					look, 
////													std::unique_ptr<SKTRAN_RayStorage_Straight>		storage, 
////													std::unique_ptr<SKTRAN_RayOptical_Base>*		rayptr  ) const
//{
//	bool ok;
//	double robs, tobs, rt;
//	SKTRAN_RayGeometry_Straight*	straightgeoray;
//	SKTRAN_RayOptical_Straight
//	std::unique_ptr<SKTRAN_RayOptical_Straight> straightray( new SKTRAN_RayOptical_Straight(coords, std::move(storage), this) );
//
//	ok =      (straightray.get() != nullptr);
//	ok = ok && straightray->TraceRay( pt, look);
//
//	if (ok)
//	{
//		straightgeoray = straightray->GRAY();
//
//		ok = straightgeoray->CalculateBaseLineTangentPointDetails( 0.0, &robs, &tobs, &rt );
//		straightgeoray->SetRobs( robs );
//		straightgeoray->SetTobs( tobs );
//		straightgeoray->SetRt( rt );
//
//		ok = ok && TraceRay( straightgeoray );
//		if( ok )
//		{
//		// set up the geometry endpoints
//			if( straightray->GetNumQuadraturePoints() >= 2 )
//			{
//				straightgeoray->SetDistanceToStart( straightray->GRAY()->TrajectoryDistance().front());
//				straightgeoray->SetDistanceToEnd  ( straightray->GRAY()->TrajectoryDistance().back());
//			}
//			else
//			{
//				straightgeoray->SetDistanceToStart( 0 );
//				straightgeoray->SetDistanceToEnd  ( 0 );
//			}	
//		}
//	}
//	*rayptr = std::move(straightray);
//	return ok;
//}
//

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Straight_Generic::AddGeometryObject		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Straight_Generic::AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject> obj )
{
	bool ok = (obj != nullptr);
	if( ok )
	{
		// we dont add duplicate elements to the container, this can happen when the wf grid and
		// optical properties grids overlap for example
		if( std::find_if( std::begin( m_container ), std::end( m_container ), [&obj] ( std::unique_ptr<SKTRAN_GeometryObject>& s ) { return *s == *obj; } ) == std::end( m_container ) )
		{
			// object is not already not in the container
			m_container.emplace_back( std::move(obj) );
		}
	}
	return ok;
}

bool SKTRAN_RayTracer_Straight_Generic::SetEarthRadius( double radius )
{
	m_earth.RadiusVar() = radius;
	return radius > 0.0;
}