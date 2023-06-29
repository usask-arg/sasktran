#include "../sktran_common.h"


/**-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Straight::SKTRAN_RayOptical_Straight		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayOptical_Straight::SKTRAN_RayOptical_Straight( std::unique_ptr<SKTRAN_RayStorage_Straight> trajectorystorage, std::shared_ptr< SKTRAN_RayTracer_Shells> raytracer)
	                       : m_trajectorystorageobject( std::move(trajectorystorage) ),			// Note that the order here is important
	                         m_raytracer(raytracer)
{
	m_trajectorystorage = m_trajectorystorageobject.get();
	InitializeStorage	(m_trajectorystorage);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayOptical_Straight::~SKTRAN_RayOptical_Straight		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayOptical_Straight::~SKTRAN_RayOptical_Straight()
{
}

/*-----------------------------------------------------------------------------
	*					TraceRay		 2014- 11- 20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayOptical_Straight::TraceRay_NewMethod() 
{
	bool	ok;

	ok = m_raytracer->TraceStraightRay( this );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Straight::LocationAlongRayAsVector		2013-06-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayOptical_Straight::LocationAlongRayAsVector( const double& distancealongray, HELIODETIC_VECTOR* pt ) const
{
	pt->SetCoords( LookVector(), distancealongray );
	(*pt) += GetObserver();

	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Straight::GetQuadratureInterpParams_startPoint		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayOptical_Straight::GetQuadratureInterpParams_startPoint( size_t quadraturepoint, size_t cellidx, double* r0, double* t0, double* rt, HELIODETIC_POINT* startpoint ) const 
{
	bool ok = true;

	*r0 = m_trajectorystorage->RadiusOfPoint( quadraturepoint );
	*t0 = m_trajectorystorage->DistanceOfPointFromCellTangentPoint( quadraturepoint, cellidx );
	*rt = m_trajectorystorage->RadiusOfCellTangentPoint( cellidx);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Straight::GetQuadratureInterpParams		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayOptical_Straight::GetQuadratureInterpParams( size_t cellidx, double* r0, double* r1, double* t0, double* t1, double* rt, HELIODETIC_POINT* startpoint, HELIODETIC_POINT* endpoint ) const
{
	bool ok = true;

	// do the end point first (equal to start point of next cell)
	ok = ok && GetQuadratureInterpParams_startPoint( cellidx+1, cellidx, r1, t1, rt, endpoint);

	// now the start point
	ok = ok && GetQuadratureInterpParams_startPoint( cellidx, cellidx,   r0, t0, rt, startpoint);


	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Straight::NotifyDerived_RayInvalid		2014-1-27*/
/** The user has changed a parameter, (usually MoveObserver) that invalidates
 *	the ray. This code enforces that by destroying all of the relevant ray
 *	information. Note that we use resize for the array operation as this does not
 *	re-allocate memory.
 **/
/*---------------------------------------------------------------------------*/

void SKTRAN_RayOptical_Straight::NotifyDerived_RayInvalid()
{
}


