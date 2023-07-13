#include "skoccultation.h"


/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableLinesOfSight::SKOCCULT_TableLinesOfSight		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_TableLinesOfSight::SKOCCULT_TableLinesOfSight()
{
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableLinesOfSight::~SKOCCULT_TableLinesOfSight		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_TableLinesOfSight::~SKOCCULT_TableLinesOfSight()
{
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableLinesOfSight::ClearRayTrajectories		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

void SKOCCULT_TableLinesOfSight::ClearRayTrajectories()
{
	std::for_each( m_raytrajectory.begin(), m_raytrajectory.end(), []( SKOCCULT_RayGeometry_Curved_Piecewise& entry){ entry.ClearRay(); } );
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableLinesOfSight::SetLinesOfSight		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_TableLinesOfSight::SetLinesOfSight(	const SKTRAN_LineOfSightArray_V21&						observerlinesofsight,
													std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	coordinates,
													SKTRAN_RayTracingRegionManager*							raytracingregionmanager)
{
	HELIODETIC_VECTOR				observer;
	HELIODETIC_UNITVECTOR			look;
	nxVector						obsoff;
	bool							ok1;
	bool							ok = true;
#if defined(NXDEBUG)
	nxGeodetic						geoid( coordinates->OsculatingGeoid() );
#endif

	m_raytrajectory.resize( observerlinesofsight.NumRays() );
	for ( size_t idx = 0; idx < NumRays(); idx++)
	{
		obsoff   = coordinates->TranslateGeoidToOsculatingSphere( observerlinesofsight.Entry(idx)->Observer() );
		observer = coordinates->GeographicToHelio( obsoff );
		look     = coordinates->GeographicToHelioUnitVector( observerlinesofsight.Entry(idx)->Look() );

#if defined(NXDEBUG)
		nxVector	obs( observer.X(), observer.Y(), observer.Z() );
		nxVector	l(   look.X(), look.Y(), look.Z()) ;
		geoid.FromTangentPointLocation( obs, l);
		double h = geoid.Height();
		double lat = geoid.GeodeticLatitude();
		double lng = geoid.GeodeticLongitude();
		NXTRACE(( "Line of sight [%3d] in osculating sphere system tangent point at (lat = %8.3f, lng = %8.3f, h = %8.3f)\n", (int)idx, (double) lat, (double)lng, (double)h));
//		geoid.FromGeocentric( obs);
//		h = geoid.Height();
//		lat = geoid.GeodeticLatitude();
//		lng = geoid.GeodeticLongitude();
//		NXTRACE(( "Line of sight [%3d] in osculating sphere system observer      at (lat = %8.3f, lng = %8.3f, h = %8.3f)\n", (int)idx, (double) lat, (double)lng, (double)h));
#endif
		ok1 = m_raytrajectory.at(idx).Initialize( coordinates, observer, look);
		ok = ok && ok1;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_TableLinesOfSight::SetLinesOfSight, there were errors configuring the lines of sight for the occultation engine. Results may be untrustworthy.");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_TableLinesOfSight::AddLinesOfSightFromTangentAltitude		 2014- 5- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_TableLinesOfSight::AddLinesOfSightFromTangentAltitude(	const std::list<SKOCCULT_LOSFromTangentPoint>&	tanpoint_los,
																		SKTRAN_RayRefracted_ThomPepSim&					raytracer,
																		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>			coordinates)
{
	nxVector						obsgeo;
	nxVector						lookgeo;
	HELIODETIC_VECTOR				observer;
	HELIODETIC_UNITVECTOR			look;
	nxVector						obsoff;
	bool							ok1;
	bool							ok = true;
	auto							iter = tanpoint_los.begin();
	size_t							nrays;
	size_t							idx;

#if defined(NXDEBUG)
	nxGeodetic						geoid( coordinates->OsculatingGeoid() );
#endif

	if ( tanpoint_los.size() > 0)
	{
		nrays = m_raytrajectory.size();
		m_raytrajectory.resize( nrays + tanpoint_los.size() );
		idx = 0;
		for ( ; !(iter == tanpoint_los.end()); ++iter )
		{
			ok1 = raytracer.FindObserverAndLookFromTangentAltitude( (*iter).TangentAltitude, (*iter).ObserverAltitude, coordinates->GetSunUnitVector(), &obsgeo, &lookgeo );
			obsoff   = coordinates->TranslateGeoidToOsculatingSphere( obsgeo );
			observer = coordinates->GeographicToHelio( obsoff );
			look     = coordinates->GeographicToHelioUnitVector( lookgeo );

	//#if defined(NXDEBUG)
	//		nxVector	obs( obsgeo.X(), obsgeo.Y(), obsgeo.Z() );
	//		nxVector	l(   lookgeo.X(), lookgeo.Y(), lookgeo.Z()) ;
	//		geoid.FromTangentPointLocation( obs, l);
	//		double h = geoid.Height();
	//		double lat = geoid.GeodeticLatitude();
	//		double lng = geoid.GeodeticLongitude();
	//		NXTRACE(( "Line of sight [%3d] from tangent point alt in osculating sphere system, straight line tangent point at (lat = %8.3f, lng = %8.3f, h = %8.3f)\n", (int)idx, (double) lat, (double)lng, (double)h));
	//#endif
			ok1 = m_raytrajectory.at(nrays + idx).Initialize( coordinates, observer, look);
			ok = ok && ok1;
			++idx;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_TableLinesOfSight::SetLinesOfSight, there were errors configuring the lines of sight for the occultation engine. Results may be untrustworthy.");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::SKOCCULT_RayGeometry_Curved_Piecewise		2014-1-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_RayGeometry_Curved_Piecewise::SKOCCULT_RayGeometry_Curved_Piecewise()
{
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::~SKOCCULT_RayGeometry_Curved_Piecewise		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_RayGeometry_Curved_Piecewise::~SKOCCULT_RayGeometry_Curved_Piecewise()
{
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::Initialize		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_RayGeometry_Curved_Piecewise::Initialize( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords, const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look)
{
	ClearRay();
	m_coords = coords;
	m_observer = observer;
	m_look     = look;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::ReserveTrajectoryStorage		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_RayGeometry_Curved_Piecewise::ReserveTrajectoryStorage( size_t n )
{
	ClearRay();
	m_curvedpathdistance.reserve(n);
	m_path.reserve(n);
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::GetCellLength		2014-1-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_RayGeometry_Curved_Piecewise::GetCellLength( size_t idx, double* length ) const
{
	*length = m_curvedpathdistance[idx+1];   //m_curvedpathdistance[idx+1] - 	
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::GetCellMidPoint		2014-2-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_RayGeometry_Curved_Piecewise::GetCellMidPoint( size_t cellidx, HELIODETIC_POINT* pt        ) const
{
	HELIODETIC_VECTOR v = m_path.at(cellidx).Vector();			// Get the start point

	v += m_path.at(cellidx+1).Vector();							// add on the end point
	v *= 0.5;													// multiply by 0.5
	m_coords->HelioVectorToHelioPoint( v, pt );
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::GetCellStartLocation		2014-2-1*/
/** **/
/*---------------------------------------------------------------------------*/


bool SKOCCULT_RayGeometry_Curved_Piecewise::GetCellStartLocation( size_t cellidx, HELIODETIC_POINT* pt        ) const
{
	m_coords->HelioVectorToHelioPoint( m_path.at(cellidx).Vector(), pt );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::GetCellRayUnitVector		2014-2-1*/
/** Get the **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_RayGeometry_Curved_Piecewise::GetCellRayUnitVector( size_t cellidx, HELIODETIC_UNITVECTOR* look ) const
{
	HELIODETIC_VECTOR	v1 = m_path.at(cellidx  ).Vector();
	HELIODETIC_VECTOR	v2 = m_path.at(cellidx+1).Vector();

	v2 -= v1;							// Get v2 - v1, ie vector away from observer, its the straight line from point to point
	v2.UnitVector();					// Convert it to a unit vector			
	*look = v2.UnitVector();			// and assign
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::Push_Point		2014-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_RayGeometry_Curved_Piecewise::Push_Point( const HELIODETIC_POINT& pt, double s)
{
	m_path.push_back(pt);
	m_curvedpathdistance.push_back(s);
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKOCCULT_RayGeometry_Curved_Piecewise::ClearRay		2014-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

void SKOCCULT_RayGeometry_Curved_Piecewise::ClearRay()
{
	m_path.resize(0);
	m_curvedpathdistance.resize(0);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Curved_ThomPepSim::SKTRAN_RayTracer_Curved_ThomPepSim		2014-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracer_Curved_ThomPepSim::SKTRAN_RayTracer_Curved_ThomPepSim( )
{
	NXTRACE(("SKTRAN_RayTracer_Curved_ThomPepSim, This class is still in beta. It does not handle observer in atmosphere\n"));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Curved_ThomPepSim::~SKTRAN_RayTracer_Curved_ThomPepSim		2014-2-1*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracer_Curved_ThomPepSim::~SKTRAN_RayTracer_Curved_ThomPepSim( )
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Curved_ThomPepSim::Initialize		2014-2-1*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Curved_ThomPepSim::Initialize ( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>			coords,
													  std::shared_ptr<const SKTRAN_GridDefRayTracingShells_V21>		raytracinggrid,
													  skClimatology*								raytracing_atmosphericstate,
													  const GEODETIC_INSTANT&						riprofile_location,		
													  double										raytracing_wavenumber)
{
	bool	ok;

	ok =       m_curvedraytracer.Configure( raytracinggrid, coords);
	ok = ok && m_curvedraytracer.SetRIAtmosphericState( raytracing_atmosphericstate);
	ok = ok && m_curvedraytracer.SetRILocationAndWavenumber( riprofile_location, raytracing_wavenumber);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Curved_ThomPepSim::CreateRay		2014-2-4*/
/** A curved ray tracer that traces ray through a spherical shell atmosphere.
 *	The method draws the ray through a spherical  atmosphere as a series
 *	radii, angles and curved lengths. It then converts the radii and angles into
 *	actual points in the atmosphere.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Curved_ThomPepSim::TraceRay( SKOCCULT_RayGeometry_Curved_Piecewise* georaycurved ) const
{
	double												zenithangle;
	double												observerheight;
	HELIODETIC_POINT									pt;
	HELIODETIC_POINT									trajectorypt;
	bool												ok;
	SKTRAN_RayRefracted_TrajectoryData::const_iterator	lastiter;
	SKTRAN_RayRefracted_TrajectoryData::const_iterator	iter;
	HELIODETIC_VECTOR									v;
	double												s;
	const HELIODETIC_VECTOR&							observer     = georaycurved->Observer();
	const HELIODETIC_UNITVECTOR&						look         = georaycurved->LookVectorAtObserver();
	SKTRAN_RayRefracted_TrajectoryData					trajectory;

	NXTRACE_ONCEONLY(firstime,("2014-02-04 ndl303, SKTRAN_RayTracer_Curved_ThomPepSim::CreateRay, It would make sense to have variable 'trajectory' as a member variable but this does not work if we have to have const attributes for multi-thread safety\n"));
	georaycurved->ClearRay();
	georaycurved->CoordinatesPtr()->HelioVectorToHelioPoint( observer, &pt );
	observerheight = pt.Altitude();
	zenithangle    = nxmath::acosd( look & pt.LocalZenith() );
	ok             = m_curvedraytracer.REFRAC( observerheight, zenithangle, &trajectory );		// Get the ray trajectory as a series of radii, angles and curved lengths
	if (ok && (trajectory.size() > 0))																						// if that worked
	{																							// then 
		iter     = trajectory.begin();															// convert the trajectory points defined in that ray tracing
		lastiter = iter;																		// into actual 3d locations in HELIODETIC coordinates
		++iter;
		while ( !(iter == trajectory.end()))
		{
			s = (*iter).length - (*lastiter).length;
			v = trajectory.ConvertTo3DLocation( observer, pt, look, iter );
			georaycurved->CoordinatesPtr()->HelioVectorToHelioPoint( v, &trajectorypt);
			georaycurved->Push_Point( trajectorypt,s);
			lastiter = iter;
			++iter;
		}
	}
	return ok;
}

