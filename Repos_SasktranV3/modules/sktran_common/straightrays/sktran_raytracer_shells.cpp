#include "../sktran_common.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::SKTRAN_RayTracer_Shells		 2014- 11- 6*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracer_Shells::SKTRAN_RayTracer_Shells(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords)
	: SKTRAN_RayTracer_Base(coords)
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::~SKTRAN_RayTracer_Shells		 2014- 11- 6*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracer_Shells::~SKTRAN_RayTracer_Shells()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::Initialize		2013-06-07*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::Initialize( std::shared_ptr< const SKTRAN_GridDefRayTracingShells_V21> raytracinggrid )
{
	m_parentgrid = raytracinggrid;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::CreateRay		2013-05-28*/
/** **/
/*---------------------------------------------------------------------------*/

//bool SKTRAN_RayTracer_Shells::CreateRayOLD( std::unique_ptr<SKTRAN_RayStorage_Straight>	storage,
//										 std::unique_ptr<SKTRAN_RayOptical_Base>*	    rayptr  ) const
//{
//	bool							ok;
//	SKTRAN_RayOptical_Straight*		straightray;
//
//	straightray = new SKTRAN_RayOptical_Straight( std::move(storage), RayFactoryObject() );
//	ok   =  (straightray != nullptr);
//	if (!ok)
//	{
//		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracer_Shells::CreateRay, Error creating ray; should allow user to handle error.");
//	}
//
//	*rayptr = std::move(std::unique_ptr<SKTRAN_RayOptical_Straight>(straightray));
//	return ok;
//}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::TraceRay		 2014- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceStraightRay( SKTRAN_RayOptical_Straight*			straightray) const
{
	bool							ok;


	ok   =  (straightray != nullptr && m_parentgrid != nullptr);
	if (ok)
	{
		ok = ok && TraceRayInternal( straightray );
//		if(ok)
//		{
//			if( straightray->GetNumQuadraturePoints() > 1 )
//			{
//				straightgeoray->SetDistanceToStart( straightray->GRAY()->StraightStorageVar()->DistanceOfPointFromOrigin(0) );
//				straightgeoray->SetDistanceToEnd  ( straightray->GRAY()->StraightStorageVar()->DistanceOfPointFromOrigin(straightray->GetNumQuadraturePoints()-1) );
//			}
//			else
//			{
//				straightgeoray->SetDistanceToStart( 0 );
//				straightgeoray->SetDistanceToEnd  ( 0 );
//			}
//		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracer_Shells::CreateRay, Error tracing ray; Hmmm. should allow user to handle error.");
	}
	return ok;
}



//bool SKTRAN_RayTracer_Shells::JoinTwoVectors ( const HELIODETIC_VECTOR& start, const HELIODETIC_VECTOR& end, SKTRAN_RayOptical_Base& ray, bool truncateToEnd ) const 
//{
//	bool ok = true;
//
//	HELIODETIC_VECTOR lookdir;
//
//	lookdir = end-start;
//	if( 1e-6 < lookdir.Magnitude() ){
//		// Trace ray and truncate to #end location
//		ok = ok && ray.Initialize( *ray.GeometryRay()->Coordinates(), start, lookdir.UnitVector() );
//		ok = ok && CreateRay( ray );
//		SKTRAN_RayOptical_Straight*    straightray		= reinterpret_cast<SKTRAN_RayOptical_Straight*> (&ray);
//		SKTRAN_RayGeometry_Straight*   straightgeoray	= straightray->GRAY();
//		if(truncateToEnd) straightgeoray->TruncateToDistance( lookdir.Magnitude() );
//	} else{
//		// Points are too close together -- make ray containing only the observer
//		SKTRAN_RayOptical_Straight*		straightray		= reinterpret_cast<SKTRAN_RayOptical_Straight*> (&ray);
//		SKTRAN_RayGeometry_Straight*	straightgeoray	= straightray->GRAY(); 
//		ok = ok && ray.Initialize( *ray.GeometryRay()->Coordinates(), start, ray.GeometryRay()->LookVector() );
//		AllocatePathElements( 1, straightray );
//		straightgeoray->Push_BackPoint( 0.0, straightgeoray->Robs(), fabs(straightgeoray->Tobs()) );
//	}
//
//	return ok;
//}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::ALlocatePathElements		2013-05-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::AllocatePathElements( size_t numelements, SKTRAN_RayOptical_Straight* straightray ) const
{
	straightray->StorageVar()->Resize(0);
	if( 0 == numelements )
	{
		// allocate at least one element to be safe
		//straightray->m_quaddistances.resize( 1 );
		straightray->StorageVar()->Reserve(1);
	}
	else
	{
		//straightray->m_quaddistances.resize( numelements );
		straightray->StorageVar()->Reserve(numelements);
	}
	//straightray->m_numquadraturepoints = straightray->m_quaddistances.size();
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells::TraceRayInternal		2013-05-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceRayInternal( SKTRAN_RayOptical_Straight* ray ) const
{
	double		highestshell;							// The radius of the highest shell from the spatial grid, Upper bound for atmosphere
	double		lowestshell;							// The radius of the lowest shell from the spatial grid, Lower bound for atmosphere
	bool		isUp;
	bool		ok;
	bool		above;
	double		Rt;
	double		Robs;
	double		Tobs;

	SKTRAN_RayStorage_Straight*	storage = ray->StraightStorageVar();
	
	storage->Resize(0);
	Rt   = storage->Rt();
	Tobs = storage->Tobs();
	Robs = storage->Robs();

#if defined(NXDEBUG)
	double tangentalt = ray->Coordinates()->RadiusToAltitude(Rt);
#endif

	isUp				   = (ray->GetObserver().UnitVector() & ray->LookVector()) >= 0.0;			// True if the observer is looking upwards (at local zenith of observer)
	highestshell           =  ray->Coordinates()->AltitudeToRadius( m_parentgrid->HighestShell() );
	lowestshell            =  ray->Coordinates()->AltitudeToRadius( m_parentgrid->LowestShell()  );

	above = (Robs > highestshell );						// See if the observer is outside the atmosphere (Cases 1 and 2). 
	if (above)												// if the user is ABOVE the atmosphe
	{														// then
		ok = (!isUp && ( Rt < highestshell) );			// See if the LOS goes through the top shell (if not, the geometry is invalid).
		if (ok)												// as long as rays goe through the atmosphere
		{													// then
			if ( Rt <= lowestshell ) ok = TraceObserverABOVE_LOSHitsGround ( ray );
			else                     ok = TraceObserverABOVE_LOSPassesThrough( ray );
		}
		else
		{
			storage->SetGroundIsHit( false );
			ok  =  AllocatePathElements(1, ray);				// Observer looking up, ray never hits atmosphere., Allocating at least 1 element avoids a lot of "out-of-range" exceptions from std::vector
			storage->PushBack(  Robs, Tobs, 0.0 );				// Get the distance of the observer from the shell;
			NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_RayGeometry_Shells::TraceRaysThroughShells, got to check that 0 shell rays work ok\n"));
//			nxLog::Record(NXLOG_WARNING,"SKTRAN_RayGeometry_Shells::TraceRaysThroughShells, observer is above atmosphere looking upwards, Path is set to NULL");
		}
	}
	else													// The user is not above the atmosphere
	{														// so make sure he is inside
		ok = (Robs >= lowestshell) ;						// See if the observer is inside the atmosphere
		if (!ok)
		{
			if ((lowestshell-Robs) < 0.000001)				// Watch out for numerical roundoff
			{
				Robs = lowestshell;
				ok = true;
			}
		}
		if (ok)												// then handle
		{													// the three cases
			if (isUp)						        ok = TraceObserverINSIDE_LookingUp            ( ray );
			else if ( Rt <= lowestshell )			ok = TraceObserverINSIDE_LookingDownHitsGround( ray );
			else									ok = TraceObserverINSIDE_LookingDownPassesThru( ray );
		}													// otherwise
		else												// observer is below, 
		{													// which is valid but log a verbose message for debugging purposes as it often indicates a bug
			nxLog::Record(NXLOG_WARNING,"SKTRAN_RayGeometry_Shells::TraceRaysThroughShells, observer is below lowest shell by (%e) meters, Path is set to NULL", (double)(lowestshell-Robs) );
		}
	}

	//if(ok) 
	//{
	//	std::stringstream	stream;
	//	std::ofstream		outfile;
	//	stream << "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2013/curvedRays/" << "shellpath_" << ray->GeometryRay()->GetTangentRadius() << ".txt";
	//	outfile.open(stream.str());
	//	HELIODETIC_VECTOR vec;
	//	HELIODETIC_POINT pt;
	//	for( size_t idx = 0; idx < ray->GetNumQuadraturePoints(); idx++ )
	//	{
	//		ray->GetCellQuadraturePoint(idx,pt);
	//		vec = pt.Vector();
	//		outfile << vec.X() << " " << vec.Y() << " " << vec.Z() << std::endl;
	//	}
	//	outfile.close();
	//}


	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Shells::DistanceToTangentPoint_fromTrig		2007-11-10*/
/** Get the distance of a shell, pointed to by shelliter, from the tangent point as defined by
 *	tangentradiussquared.
 **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayTracer_Shells::DistanceToTangentPoint_fromTrig( double rshell, double tangentradiussquared) const
{
	double sqrdist;

	NXASSERT(( rshell > 6300000.0 ));								// Make sure we have a sensible Earth like radius (And not an altitude)'
	sqrdist = rshell*rshell - tangentradiussquared;					// we occassionally get round off issues
	return (sqrdist > 1.0E-12) ? sqrt(sqrdist) : 0.0;				// any distance less than 1 micron squared is really 0.0
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Shells::ObserverABOVELOShitsGround		2007-11-10*/
/** Setup the pathlength and cell indices for an observer outside of the
 *	atmosphere with a ray that hits the ground. This geometry guarantees
 *	all grid cells are used in the pathlength elements.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceObserverABOVE_LOSHitsGround( SKTRAN_RayOptical_Straight* ray ) const
{
	size_t										numshells;
	SKTRAN_GridDefBase_V2::const_iterator		shelliter;
	double										tangentradiussquared;
	double										t;
	size_t										i;
	bool										ok;
	double										r;
	double										s;
	SKTRAN_RayStorage_Straight*					storage;

	storage               = ray->StraightStorageVar();
	NXASSERT( (storage->Tobs() >= 0.0) );
	tangentradiussquared  = nxmath::sqr(storage->Rt());
	numshells             = m_parentgrid->NumShells();
	ok                    = AllocatePathElements( numshells, ray );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayGeometry_Shells::TraceObserverABOVE_LOSHitsGround, error allocating memory for data arrays");
	}
	else
	{
		shelliter  = m_parentgrid->end();														// get the end of the grid shells array
		for (i =0; i < numshells; i++ )															// For all of the shellintercepts
		{																						// calculate distance of observer from shell intercepts
			--shelliter;																		// point to the highest (last) altitude on the grid
			r = storage->GetCoordsPtr()->AltitudeToRadius(*shelliter);
			t = DistanceToTangentPoint_fromTrig( r,tangentradiussquared);						// Get the distance of shell intercept from the tangent point
			s = storage->Tobs() - t;
			storage->PushBack(  r, t, s );							// Get the distance of the observer from the shell;
		}
	}
	storage->SetGroundIsHit( true );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Shells::TraceObserverINSIDE_LookingUp		2007-11-10*/
/** Setup the pathlength and cell indices for an observer inside the
 *	atmosphere looking up
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceObserverINSIDE_LookingUp( SKTRAN_RayOptical_Straight* ray ) const
{
	SKTRAN_GridDefBase_V2::const_iterator			shelliter;
	SKTRAN_GridDefBase_V2::const_iterator			shellfinish;
	SKTRAN_GridDefBase_V2::const_iterator			shellstart;
	size_t											shellMin;
	size_t											numshells;
	size_t											numelements;
	size_t											i;
	bool											ok;
	double											t;
	double											tangentradiussquared;
	double											r;
	double											obsalt;
	SKTRAN_RayStorage_Straight*						storage;

	storage                 = ray->StraightStorageVar();
	storage->SetGroundIsHit( false );
	tangentradiussquared	= nxmath::sqr( storage->Rt() );
	numshells				= m_parentgrid->NumShells();
	shellstart				= m_parentgrid->begin();
	shellfinish				= m_parentgrid->end();
	obsalt					= ray->Coordinates()->RadiusToAltitude( storage->Robs());
	shelliter				= std::upper_bound( shellstart, shellfinish, obsalt);	// Find the lowest shell in the ray tracing grid that is greater than the observer's  radius

	shellMin      = (shelliter - shellstart);										// Convert the pointer to integer index
	numelements   = numshells - shellMin + 1;									// Get the number of shell intercepts, include 1 point for the observers position.
	ok            = AllocatePathElements( numelements, ray );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayGeometry_Shells::TraceObserverINSIDE_LookingUp, error allocating arrays for line of sight");
	}
	else
	{
		storage->PushBack( storage->Robs(), fabs(storage->Tobs()), 0.0  );			// Get the distance of the observer from the shell;
		NXASSERT( (storage->Tobs() <= 1.0E-05) );									// Make sure Tobs is less than 0 with slop for numerical roundoff
		shelliter  = m_parentgrid->begin() + shellMin;								// get the first ray tracing end of the grid shells array
		for (i =1; i < numelements; i++ )											// For all of the shellintercepts
		{																			// calculate distance of observer from shell intercepts
			r = storage->GetCoordsPtr()->AltitudeToRadius(*shelliter);
			t = DistanceToTangentPoint_fromTrig( r,tangentradiussquared);			// Get the distance of shell intercept from the tangent point
			storage->PushBack( r, t, t + storage->Tobs());							// Get the distance of the observer from the shell, Note: m_Tobs is intrinsically negative so we are really subtracting
			++shelliter;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Shells::TraceObserverINSIDE_LookingDownHitsGround		2007-11-14*/
/** Setup the pathlength and cell indices for an observer inside the
 *	atmosphere looking down where the ray hit the gorund
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceObserverINSIDE_LookingDownHitsGround( SKTRAN_RayOptical_Straight* ray ) const
{
	SKTRAN_GridDefBase_V2::const_iterator			shelliter;
	SKTRAN_GridDefBase_V2::const_iterator			shellfinish;
	SKTRAN_GridDefBase_V2::const_iterator			shellstart;
	size_t											numelements;
	size_t											i;
	bool											ok;
	double											t;
	double											tangentradiussquared;
	double											obsalt;
	double											r;
	SKTRAN_RayStorage_Straight*						storage;

	storage                 = ray->StraightStorageVar();
	storage->SetGroundIsHit( true );
	tangentradiussquared = nxmath::sqr( storage->Rt());
	obsalt        = storage->GetCoordsPtr()->RadiusToAltitude( storage->Robs() );				// Get the (rounded) altitude of the observer					
	shellstart    = m_parentgrid->begin();
	shellfinish   = m_parentgrid->end();
	shelliter     = std::lower_bound( shellstart, shellfinish, obsalt )	;	// Find the lowest shell in the grid that is greater than or equal to observer
	numelements   = (shelliter - shellstart) + 1;										// Convert the pointer to integer index
	ok            = AllocatePathElements( numelements, ray );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayGeometry_Shells::TraceObserverINSIDE_LookingDownHitsGround, error allocating arrays for line of sight");
	}
	else
	{
		storage->PushBack( storage->Robs(), fabs(storage->Tobs()), 0.0 );				// The first ray intercept is at the observer
		for (i =1; i < numelements; i++ )												// For all of the shellintercepts
		{																				// calculate distance of observer from shell intercepts
			--shelliter;																// point to the first shell below the observer
			r = storage->GetCoordsPtr()->AltitudeToRadius( *shelliter );
			t = DistanceToTangentPoint_fromTrig( r, tangentradiussquared);				// Get the distance of shell intercept from the tangent point
			storage->PushBack(  r, t, storage->Tobs() - t  );							// Get the distance of the observer from the shell, Note: m_Tobs is intrinsically positive
		}
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Shells::TraceObserverABOVE_LOSPassesThrough		2007-11-10*/
/** Setup the pathlength and cell indices for an observer outside of the
 *	atmosphere with a ray that passes through.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceObserverABOVE_LOSPassesThrough( SKTRAN_RayOptical_Straight* ray ) const
{
	SKTRAN_GridDefBase_V2::const_iterator			shelliter;
	SKTRAN_GridDefBase_V2::const_iterator			shellfinish;
	SKTRAN_GridDefBase_V2::const_iterator			shellstart;
	size_t											shellMin;
	size_t											shellInts;
	size_t											numshells;
	size_t											numelements;
	size_t											i;
	double											t;
	bool											ok;
	double											rtalt;
	double											r;
	double											tangentradiussquared;
	SKTRAN_RayStorage_Straight*						storage;

	storage                 = ray->StraightStorageVar();
	tangentradiussquared    = nxmath::sqr( storage->Rt() );
	storage->SetGroundIsHit( false );
	numshells     = m_parentgrid->NumShells();
	shellstart    = m_parentgrid->begin();
	shellfinish   = m_parentgrid->end();
	rtalt         = storage->GetCoordsPtr()->RadiusToAltitude(storage->Rt());
	shelliter     = std::upper_bound( shellstart, shellfinish, rtalt );			// Find the lowest shell in the grid that is greater than tangent radius of this ray
	shellMin      = (shelliter - shellstart);									// Convert the pointer to integer index
	shellInts     = numshells - shellMin;										// Get the number of cells this ray passes through
	numelements   = 2*shellInts + 1;											// Include 1 point for the tangent point itself
	
	ok =  AllocatePathElements(  numelements, ray );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayGeometry_Shells::TraceObserverABOVE_LOSPassesThrough, error allocating arrays for line of sight");
	}
	else
	{
		shelliter  = m_parentgrid->end();
		for (i =0; i < shellInts; i++ )											// For all of the shellintercepts
		{																		// calculate distance of observer from shell intercepts
			--shelliter;														// point to next shell going downwards
			r = storage->GetCoordsPtr()->AltitudeToRadius(*shelliter);
			t = DistanceToTangentPoint_fromTrig( r,tangentradiussquared);		// Get the distance of shell intercept from the tangent point
			storage->PushBack(  r, t, (storage->Tobs() - t) );									// Get the distance of the observer from the shell, Note: m_Tobs is intrinsically positive
		}
		storage->PushBack( storage->RadiusOfCellTangentPoint(0), 0.0, storage->Tobs()  );		// Insert the tangent point as one of the shells

		for (i =0; i < shellInts; i++ )											// Now step back out of the atmosphere
		{																		// from the far side of the tangent point
			r = storage->GetCoordsPtr()->AltitudeToRadius(*shelliter);
			t = DistanceToTangentPoint_fromTrig( r, tangentradiussquared);		// calculate the distance of shell from tangent point
			storage->PushBack( r,t, ( storage->Tobs() + t ) );							// insert its distance from the observer
			++shelliter;														// point to next shell going upwards
		}
		NXASSERT(( storage->NumQuadraturePoints() == numelements ));
	}
	return ok;
}	

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Shells::TraceObserverINSIDE_LookingDownPassesThru		2007-11-14*/
/** Setup the pathlength and cell indices for an observer inside the
 *	atmosphere looking down where the ray passes through the atmospherehit the gorund
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracer_Shells::TraceObserverINSIDE_LookingDownPassesThru( SKTRAN_RayOptical_Straight* ray ) const
{
	SKTRAN_GridDefBase_V2::const_iterator			titer;			// Pointer to the shell just above the tangent point
	SKTRAN_GridDefBase_V2::const_iterator			obsiter;			// Pointer to the shell just above the observer
	SKTRAN_GridDefBase_V2::const_iterator			shellfinish;
	SKTRAN_GridDefBase_V2::const_iterator			shellstart;
	SKTRAN_GridDefBase_V2::const_iterator			shelliter;
	size_t					numshells;
	size_t					numelements;
	size_t					numnearelements;
	size_t					numfarelements;
	bool					ok;
	double					t;
	double					r;
	double					rtalt;
	double					obsalt;
	double					tangentradiussquared;
	SKTRAN_RayStorage_Straight*						storage;

	storage              = ray->StraightStorageVar();
	tangentradiussquared = nxmath::sqr( storage->Rt() );
	storage->SetGroundIsHit( false );
	numshells      = m_parentgrid->NumShells();
	shellstart     = m_parentgrid->begin();
	shellfinish    = m_parentgrid->end();
	rtalt          = storage->GetCoordsPtr()->RadiusToAltitude( storage->Rt() );
	obsalt         = storage->GetCoordsPtr()->RadiusToAltitude( storage->Robs() );
	titer          = std::upper_bound( shellstart, shellfinish, rtalt  );			// Find the lowest shell in the grid greater than the tangent point
	obsiter        = std::lower_bound( shellstart, shellfinish, obsalt );			// Find the lowest shell in the grid greater than or equal to the observer
	if ( titer > obsiter )															// THis happens if observer and tangent point are at same altitude , and its exactly on a shell boundary
	{																				// if this has happened
		NXASSERT(( rtalt == obsalt ));												// then, make sure that is that case
		titer = obsiter;															// and make them equal
	}

	numfarelements   = (shellfinish - titer);										// Number of path elements on far side of tangent point
	numnearelements  = (obsiter     - titer)+1;										// Number of path elements on near side of tangent point (Note that lowest shell is common to near and far side)
	numelements      = numfarelements + numnearelements + 1;						// Get the number of cells this ray passes through
	ok               = AllocatePathElements( numelements, ray );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayGeometry_Shells::TraceObserverINSIDE_LookingDownPassesThru, error allocating arrays for line of sight");
	}
	else
	{
		shelliter  = obsiter;
		storage->PushBack( storage->Robs(), fabs(storage->Tobs()), 0.0);
		for (size_t i =1; i < numnearelements; i++ )									// For all of the shellintercepts
		{																		// calculate distance of observer from shell intercepts
			--shelliter;														// point to next shell going downwards
			r = storage->GetCoordsPtr()->AltitudeToRadius( *shelliter );
			t = DistanceToTangentPoint_fromTrig( r ,tangentradiussquared);				// Get the distance of shell intercept from the tangent point
			storage->PushBack( r, t, ( storage->Tobs() - t ) );							// Get the distance of the observer from the shell, Note: m_Tobs is intrinsically positive
		}
		storage->PushBack( storage->RadiusOfCellTangentPoint(0), 0.0, storage->Tobs()  );									// Insert the tangent point as one of the shells

		NXASSERT(( shelliter == titer ));
		for (size_t i =0; i < numfarelements; i++ )									// Now step back out of the atmosphere
		{																		// from the far side of the tangent point
			r = storage->GetCoordsPtr()->AltitudeToRadius( *shelliter );
			t = DistanceToTangentPoint_fromTrig( r, tangentradiussquared);					// calculate the distance of shell from tangent point
			storage->PushBack( r, t, ( storage->Tobs() + t ));	// insert its distance from the observer
			++shelliter;														// point to next shell going upwards
		}
		NXASSERT(( storage->NumQuadraturePoints() == numelements ));
		NXASSERT(( shelliter == shellfinish ));
	}
	return ok;
}

