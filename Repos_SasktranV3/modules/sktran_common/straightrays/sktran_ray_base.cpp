#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Base::SKTRAN_RayGeometry_Base		2014-1-24*/
/** Constructs the ray geometry and provides it with its own instance of appropriate storage. The storage
 *	will last for the duration of this ray. 
 * 
 *	@param trajectorystorage
 *		The storage object to be used for this ray trajectory.										
**/
/*---------------------------------------------------------------------------*/

//SKTRAN_RayGeometry_Base::SKTRAN_RayGeometry_Base( SKTRAN_RayStorage_Base* trajectorystorage)
//
//
//{
//	double	nan = std::numeric_limits<double>::quiet_NaN();
//
//	m_observer.SetCoords( nan, nan, nan); ;					//!< The location of the obserevr
//	m_look.SetCoords(nan, nan, nan);						//!< The observers look direction
//	m_trajectorystorage = trajectorystorage;
//}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Base::~SKTRAN_RayGeometry_Base		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

//SKTRAN_RayGeometry_Base::~SKTRAN_RayGeometry_Base()
//{
//}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Base::SetTrajectoryStorage		 2014- 11- 28*/
/** **/
/*---------------------------------------------------------------------------*/

//void SKTRAN_RayGeometry_Base::SetTrajectoryStorage( sSKTRAN_RayStorage_Base>	trajectorystorage)
//{
//	m_trajectorystorage = std::move(trajectorystorage);
//}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Base::SetCoordinates		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

//bool SKTRAN_RayGeometry_Base::SetCoordinates(const SKTRAN_CoordinateTransform_V2* coords)
//{
//	if (coords   != nullptr)   coords->AddRef();
//	if (m_coords != nullptr) m_coords->Release();
//	m_coords = coords;
//	return true;
//}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base::MoveObserver		 2016- 7- 15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayOptical_Base::MoveObserver( const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look )
{
	m_observer = observer;
	m_look     = look;
	m_trajectorystorage->InitializeObserver( observer, look );
	NotifyDerived_RayInvalid();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base::CalculateStraightLineTangentPoint		2014-1-26*/
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

bool SKTRAN_RayOptical_Base::CalculateBaseLineTangentPointDetails( double groundaltitude, double* Robs, double* Tobs, double* Rt)
{
	double				minRadius;
	HELIODETIC_POINT	pt;	
	bool				ok;
	double				coszenith;
	double				diff;

	NXASSERT((Coordinates() != nullptr ));
	minRadius    = Coordinates()->AltitudeToRadius( groundaltitude  );
	*Robs = GetObserver().Magnitude();
	if (*Robs <= minRadius)
	{		
		HELIODETIC_VECTOR	observer = GetObserver();
		diff       = (minRadius- (*Robs));								// then look for minor floating point differences 
		if (diff > 0.001)												// if we are outside of a roundoff error domain
		{																// thenlog an error to the user
			nxLog::Record(NXLOG_WARNING,"SKTRAN_RayGeometry_Base::CalculateBaseLineTangentPointDetails, the observer at (%e, %e, %e) was substantially below the RT ground (%e meters). It should be with 0.0001 meters. This indicates a problem in configuration", GetObserver().X(), GetObserver().Y(), GetObserver().Z(), (double)diff);
		}
		minRadius  += 0.001;											// add on a tenth of a millimeter.
		observer.SetCoords( observer.UnitVector(), minRadius);			// Fix it so observer is just above the ground
		m_observer = observer;
		*Robs = minRadius;
		Coordinates()->HelioVectorToHelioPoint(observer, &pt);

	}
	ok    = Coordinates()->HelioVectorToHelioPoint( GetObserver(), &pt );								// get the observer heliodetic point
	coszenith = pt.CosZenithAngle( LookVector() );															// and its coszenith
	*Rt   = (  coszenith*coszenith <= 0.99999999999999 ? (*Robs)*sqrt( 1.0 - coszenith*coszenith ) : 0 );	// Get the radius of the tangent point
	*Tobs =  -1.0*(*Robs)*coszenith;																		// and the tangent distance to the observer
	
	return ok;
}
	 


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base::SKTRAN_RayOptical_Base		 2014- 12- 2*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayOptical_Base::SKTRAN_RayOptical_Base()	
{
	double	nan = std::numeric_limits<double>::quiet_NaN();

	m_observer.SetCoords( nan, nan, nan); ;					//!< The location of the obserevr
	m_look.SetCoords(nan, nan, nan);						//!< The observers look direction
	m_trajectorystorage = nullptr;					
	m_coords = nullptr;

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base::InitializeStorage		 2016- 7- 15*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_RayOptical_Base::InitializeStorage(SKTRAN_RayStorage_Base* trajectorystorage )
{
	m_trajectorystorage = trajectorystorage; 
	m_coords = m_trajectorystorage->GetCoordsPtr(); 
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base::~SKTRAN_RayOptical_Base		 2014- 12- 2*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayOptical_Base::~SKTRAN_RayOptical_Base()
{
};




/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Base::GetQuadratureInterpParamsFromPointsAndLook		2014-1-24*/
/** Static helper function that calculates the geocentric radius and distance to the
 *	straight line tangent point of both the start and end points. This is used by the
 *	optical depth calculation as it linearly interpolates through a shell.
 *
 *	@param r0
 *		returns the radius of the start point in meters from the center of the Earth
 *	@param r1
 *		returns the radius of the end point in meters from the center of the Earth
 *	@param t0
 *		returns the distance of the startpoint from the rays tangent point
 *	@param t1
 *		returns the distance of the endpoint from the rays tangent point
 *	@param rt
 *		returns the radius of the ray's tangentpoint in meters from the center of the Earth.
 *	@param startpoint
 *		The location at which we are goint to start the quadrature
 *	@param endpoint
 *		The location at which we are goint to end quadrature
 *	@param look
 *		The look unit vector (away from observer) of the ray.
 *	@return 
 *	true if the trajectoryu is defined.
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayOptical_Base::GetQuadratureInterpParamsFromPointsAndLook( double* r0, double* r1, double* t0, double* t1, double* rt, 
																		  const HELIODETIC_POINT& startpoint, 
																		  const HELIODETIC_POINT& endpoint, 
																		  const HELIODETIC_UNITVECTOR& look )
{
	bool ok = true;

	double					rtsq;
	HELIODETIC_VECTOR		vec;

	// do the end point first
	vec = endpoint.Vector();
	*r1 = endpoint.Radius();
	*t1 = fabs(look.X()*vec.X() + look.Y()*vec.Y() + look.Z()*vec.Z());

	// now the start point
	vec = startpoint.Vector();
	*r0 = startpoint.Radius();
	*t0 = fabs(look.X()*vec.X() + look.Y()*vec.Y() + look.Z()*vec.Z());

	rtsq = (*r0)*(*r0) - (*t0)*(*t0);
	*rt = rtsq < 0 ? 0 : sqrt(rtsq);

	return ok;
}

