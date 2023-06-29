#include "skoccultation.h"
#include "omp.h"

/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::SKOCCULT_OCC_Engine		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_OCC_Engine::SKOCCULT_OCC_Engine()
{
	m_raytracingwavenumber   = 2400.0;
	m_maxpathlength          = 10000.0;
	m_geometrychanged		 = true;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::~SKOCCULT_OCC_Engine		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_OCC_Engine::~SKOCCULT_OCC_Engine()
{
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::SetRayTracingWavenumber		 2014- 4- 30*/
/** Sets the wavenumber used for the refractive index of the atmosphere which are
 *	are used to trace curved rays through the atmosphere. This value will take
 *	effect in subsequent calls to CalculateRadiance tracing **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OCC_Engine::SetRayTracingWavenumber( double wavenumber )
{
	m_raytracingwavenumber = wavenumber;
	return true;
}

 /*-----------------------------------------------------------------------------
 *					SKTRANSO_Engine::ConfigureModel		2007-12-12*/
/** Configures the radiative transfer model. Each call to this function
 *	resets all of the internal structures and tables.  This will force the
 *	all internal tables to recreate themselves upon the next call to
 *	CalculateRadiance.
 *
 *	The code also sets a request for the climatologies
 **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OCC_Engine::ConfigureModel( SKTRAN_SpecsUser_Base& specifications, const SKTRAN_LineOfSightArray_V21& linesofsight, size_t numthreads )
{
	bool	ok = true;
	
    SKOCCULT_Specs_User* userspecifications = dynamic_cast< SKOCCULT_Specs_User* >( &specifications);

    ok = ok && nullptr!=userspecifications;

	if ( linesofsight.NumRays()                               > 0) ok = ok && (userspecifications->TangentPointLinesOfSight().size() == 0);
	if ( userspecifications->TangentPointLinesOfSight().size() > 0) ok = ok && (linesofsight.NumRays() ==0);
	if (!ok) nxLog::Record(NXLOG_WARNING,"SKOCCULT_OCC_Engine::ConfigureModel, You cannot define lines of sight using both the tangent altitude technique and the observer/look technique. Use one or the other");
	ok = ok && userspecifications->UpdateUndefinedParametersFromLinesOfSight ( linesofsight );			// Set the unit vector from Earth to Sun (Note: only one sun is used for all lines of sight)
	ok = ok && m_modelspecifications.Initialize ( userspecifications );								// setup the new model specifications
	ok = ok && m_linesofsight.SetLinesOfSight   ( linesofsight, m_modelspecifications.CoordinateSystemObject(), userspecifications->RayTracingRegionManagerVar() );							// setup the new lines of sight
	ok = ok && m_opticalpropertiestable.ConfigureGeometry( &m_modelspecifications );
	if (!ok)																			// see if it worked
	{																					// if it did not then log an error
		nxLog::Record(NXLOG_WARNING," SKTRANSO_Engine::ConfigureModel, there was an error configuring the lines of sight, thats a problem");
	}
	if (numthreads > 0) omp_set_num_threads( (int)numthreads );

	m_geometrychanged = true;
	return ok;																			// and return the status;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::ReferencePoint		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

GEODETIC_INSTANT SKOCCULT_OCC_Engine::ReferencePoint() const
{
	GEODETIC_INSTANT point;

	point.latitude  = m_modelspecifications.CoordinateSystemPtr()->ReferencePtLatitude();
	point.longitude = m_modelspecifications.CoordinateSystemPtr()->ReferencePtLongitude();
	point.heightm   = 0.0;
	point.mjd       = m_modelspecifications.CoordinateSystemPtr()->ReferencePointMJD();
	return point;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::TraceLineOfSightRays		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OCC_Engine::TraceLineOfSightRays( SKTRAN_AtmosphericOpticalState_V21*	opticalstate)
{
	skClimatology*			atmosphericstate;
	bool					ok;
	bool					ok1;

	ok =       opticalstate->SetTimeAndLocation( ReferencePoint(), true);									// Update the internal cache of the climatologies.
	ok = ok && opticalstate->GetAtmosphericStateModel(&atmosphericstate);									// Fetch the the atmospheric state
	ok = ok && m_raytracer.Initialize( m_modelspecifications.CoordinateSystemObject(),							// Initialize the ray tracing software
									   m_modelspecifications.RayTracingSpecs()->RayTracingShells(),
									   atmosphericstate,
									   ReferencePoint(), 
									   m_raytracingwavenumber);

	ok = ok && m_linesofsight.AddLinesOfSightFromTangentAltitude(	m_modelspecifications.TangentPointLinesOfSight(), m_raytracer.CurvedRayTracer() ,m_modelspecifications.CoordinateSystemObject());

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_OCC_Engine::TraceLineOfSightRays, There were errors initializing the ray tracing software. Thats not good");
	}
	else																									// Otherwise we are good to go
	{																										// so
		for (size_t i = 0; i < m_linesofsight.NumRays(); i++)												// for each ray
		{																									// lets trace
			ok1 = m_raytracer.TraceRay( &m_linesofsight.RayAtVar(i) );										// the ray through the atmosphere
			ok = ok && ok1 ;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKOCCULT_OCC_Engine::TraceLineOfSightRays, Error tracing rays through atmosphere. Thats not good");
		}
	}

	if (!ok)
	{
		m_linesofsight.ClearRayTrajectories();
	}
	return ok;
};

static double StraightLength( const HELIODETIC_POINT& startpt, const HELIODETIC_POINT& endpt)												// And get the straight length between the end points
{
	HELIODETIC_VECTOR	v0(startpt.UnitVector(), startpt.Radius() );
	HELIODETIC_VECTOR	v1(endpt.UnitVector(),   -endpt.Radius() );

	v0 += v1;
	return v0.Magnitude();
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::IntegrateExtinctionAlongRay		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OCC_Engine::IntegrateExtinctionAlongRay( const SKOCCULT_RayGeometry_Curved_Piecewise& ray, const std::vector<double>& wavenumber, std::vector<double>* extinction)
{
	SKTRAN_OpticalDepthCalculator_LinearWithHeight	od;
	size_t											numwave;
	size_t											numcells;
	HELIODETIC_POINT								startpt;
	HELIODETIC_POINT								endpt;
	double											curvedlen;
	double											straightlen;
	bool											ok1;
	double											sigma0;
	double											sigma1;
	double											adjustment;
	std::vector<double>								lastsigma;
	size_t											iw;
	bool											ok = true;

	numwave  = wavenumber.size();																	// Get the number of wavenumbers to process
	numcells = ray.NumCells();																		// Get the number of ray cells to process
	extinction->assign( numwave, 0.0);																// reset the size to store extinction for each wavenumber


	if ( (numwave > 0) && (numcells > 0))
	{
		ray.GetCellStartLocation( 0, &endpt);															// Get the start location of the ray

		lastsigma.resize(numwave);																		// ow create an intermediate storage buffer
		for (iw = 0; iw < numwave; iw++)																// and
		{																								// initialize it to the
			lastsigma.at(iw) = m_opticalpropertiestable.TotalExtinctionPerCM( endpt, iw );				// extinction of the first point along the ray
		}																								// do it for wall wavenumbers

		for (size_t idx = 0; idx < numcells; idx++)														// for each cell along the ray
		{																								// calculate total extinction
			startpt = endpt;																			// so the start point of this cell is end point of last cell
			ok1 =        ray.GetCellLength( idx, &curvedlen );											// Get the curved length of this ray
			ok1 = ok1 && ray.GetCellStartLocation( idx+1, &endpt );										// Get the start location of the next cell. This is the end of thsi cell
			ok1 = ok1 && od.ConfigureQuadratureCoefficients( startpt, endpt );							// Initialize the optical depth calculator with the end points
			if (ok1)																					// we have configured the geometry part of the optical depth calculation
			{																							// so
				straightlen = StraightLength(startpt, endpt);												// And get the straight length between the end points
				adjustment = (straightlen > 0.001) ? curvedlen/straightlen : 1.0;						// now get the curved ray adjustment factor. Wath out for divide by zero
				for (iw = 0; iw < numwave; iw++)														// and iterate over wavenumber 
				{																						// for each wavenumber
					sigma0 = lastsigma.at(iw);															// get the extinction per cm at the start point, which is the last end point
					sigma1 = m_opticalpropertiestable.TotalExtinctionPerCM( endpt,   iw );				// Get the extinction per cm at the end point
					extinction->at(iw) += od.OpticalDepthFromStartToEnd( sigma0, sigma1)*adjustment;	// Get optical depth along segment (sigma  varies linearly with radius). Adjust for curvature of ray, compared to straight
					lastsigma.at(iw) = sigma1;															// Save the extinction at the end point for this wavenumber. We'll use on the next cell.
				}																						// do all the wavenumebrs
			}																							// do all the cells
			ok = ok && ok1;																				// and check to make sure it works
		}																								// and that is that
	}
	if (!ok)																							// Check for errors and pump out message if problem.
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_OCC_Engine::IntegrateExtinctionAlongRay, error calculating extinction along the ray. Thats not good. Setting values to NaN");
		extinction->assign( numwave, 0 );																// reset the size to store extinction for each wavenumber

	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_OCC_Engine::CalculateRadiance		2010-4-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_OCC_Engine::CalculateMultiWavelengthRadiance( std::vector< std::vector<SKTRAN_StokesScalar> >*	losradiance,
															const std::vector<double>&							wavenumber,
															size_t												unused_numordersofscatter,
															SKTRAN_AtmosphericOpticalState_V21*					opticalstate,
															std::vector< std::vector<skRTStokesVector> > *      unused_losvector, 
															bool												userupdateclimatology,
															SKTRAN_DiagnosticInterface*							diag)
{
	bool					ok;
//	SKTRAN_CodeTimer		s1;

	ok = ContainerIsAscendingOrder( wavenumber.begin(), wavenumber.end() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_OCC_Engine::CalculateRadiance, The wavelength array is not in descending order. It should be in descending order to make the code work efficiently");
	}
	

	ok = !m_geometrychanged;
	if (!ok)
	{
		ok = TraceLineOfSightRays( opticalstate);
		m_geometrychanged = !ok;
	}
	ok = ok && m_opticalpropertiestable.ConfigureOptical ( wavenumber, *opticalstate );						// Configure the optical properties table. This can be numerically intensive if it includes HITRAN in the IR
	if (ok)
	{

		losradiance->resize( wavenumber.size() );
		for (size_t iw = 0; iw < wavenumber.size(); iw++) losradiance->at(iw).resize(m_linesofsight.NumRays() );

		std::vector<double>	ext;
		#pragma omp parallel for schedule(dynamic) shared(ok) private( ext )
		for ( int ih = 0; ih < (int)m_linesofsight.NumRays(); ih++)
		{
			bool ok1;

			ok1 = IntegrateExtinctionAlongRay( m_linesofsight.RayAt(ih), wavenumber, &ext );
			if (ok1) for (size_t iw = 0; iw < wavenumber.size(); iw++) losradiance->at(iw).at(ih) = ext.at(iw);
			#pragma omp critical
			{
				ok = ok && ok1;
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_EngineThreadController::ExecuteSasktranThread, Error performing SASKTRAN Calculation, results are untrustworthy");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_Internal::SKOCCULT_Specs_Internal		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_Specs_Internal::SKOCCULT_Specs_Internal	()
{
	m_raytracingspecs   = nullptr;
	m_opticaltablespecs = nullptr;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_Internal::~SKOCCULT_Specs_Internal		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_Specs_Internal::~SKOCCULT_Specs_Internal()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_Internal::ReleaseResources		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

void SKOCCULT_Specs_Internal::ReleaseResources()
{
	if ( m_raytracingspecs   != nullptr ) m_raytracingspecs  ->Release();
	if ( m_opticaltablespecs != nullptr ) m_opticaltablespecs->Release();

	m_raytracingspecs   = nullptr;
	m_opticaltablespecs = nullptr;
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_Internal::Initialize		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_Internal::Initialize( const SKOCCULT_Specs_User* userspecifications )
{
	bool	ok;

	ReleaseResources();
	m_tanpoint_los = userspecifications->TangentPointLinesOfSight();
	ok =       userspecifications->CreateCoordinateSystem( &m_coordinatesystem );
	ok = ok && userspecifications->RayTracingSpecs ()->CreateInternalSpecs( &m_raytracingspecs   );
	ok = ok && userspecifications->OpticalGridSpecs()->CreateInternalSpecs( &m_opticaltablespecs );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_Specs_Internal::Initialize, Error fetching internal specifications from user specifications. Thats a problem");
		ReleaseResources();
	}
	return ok;

}



/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::SKOCCULT_Specs_User		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_Specs_User::SKOCCULT_Specs_User()
{
	m_manualrefpt.latitude  = std::numeric_limits<double>::quiet_NaN();
	m_manualrefpt.longitude = std::numeric_limits<double>::quiet_NaN();
	m_manualrefpt.heightm   = std::numeric_limits<double>::quiet_NaN();
	m_manualrefpt.mjd       = std::numeric_limits<double>::quiet_NaN();
	m_manualsun.SetInvalid();
	ConfigureEvenSpacedShells( 0.0, 1000.0, 100000.0 );
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::~SKOCCULT_Specs_User		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

SKOCCULT_Specs_User::~SKOCCULT_Specs_User()
{
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::ConfigureEvenSpacedShells		 2014- 4- 28*/
/**  Configures even spaced spherical shells for the radiative transfer calculations.
 *	This actually configures the ray tracing shell boundaries, the diffuse point
 *	locations and the location of the optical properties
 *
 *	The diffuse points are located in the mid-point of the ray tracing shells
  *	and the optical properties are located at both the mid point and shell boundaries.
**/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::ConfigureEvenSpacedShells( double minshellheight_meters, double shellwidth, double maxshell)
{
	bool						ok;
	nx1dArray<double>			shellalts;
	nx1dArray<double>			diffusealts;
	nx1dArray<double>			optpropheights;
	size_t						npts;

	npts                           =  (size_t) ((maxshell - minshellheight_meters)/shellwidth + 1);
	shellalts.Indgen(npts)        *= shellwidth;							// Ray Shell Altitudes at 1 km boundaries to 100 km, 0.0, 1.0, 2.0, 3.0 ... 100.0
	optpropheights.Indgen(npts+1) *= shellwidth;							// Optical properties Shell altitudes on same boundary as shells 

	shellalts		+= minshellheight_meters;
	optpropheights  += minshellheight_meters;

	ok    =       m_rayregionmanager.SetGroundAltitude(minshellheight_meters);
	ok    = ok && m_raytracingspecs.ConfigureRayTracingShellAlts     ( shellalts.UnsafeArrayBasePtr(),        shellalts.size()   );
	ok    = ok && m_opticalpropspecs.ConfigureOpticalPropertyShells	 ( optpropheights.UnsafeArrayBasePtr(),   optpropheights.size() );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_Specs_User::ConfigureEvenSpacedShells, Error configuring the default specifications. This is a problem");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_SpecificationsUser_Legacy::ConfigureUserDefinedShells		2013-6-20*/
/** Configures the shell specifications using the shells defined by the user. The
 *	shellalts parameter must pass in the altitudes of the shells in ascending order
 *	in meters.
 *	The diffuse points are placed in the middle of the shells and the optical
 *	properties are placed in the middle and on the shell
 */
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::ConfigureUserDefinedShells( const std::vector<double>& shellalts )
{
	size_t					numshells;
	double					minshellheight;
	std::vector<double>		diffusealts;
	std::vector<double>		optpropheights;		
	size_t					idx;
	bool					ok;

	NXTRACE_ONCEONLY(firsttime,("****** NEW CODE: NEEDS TO BE CHECKED OUT  ***** SKOCCULT_Specs_User::ConfigureUserDefinedShells\n"));
	numshells      = shellalts.size();
	minshellheight = shellalts.front();

	optpropheights.assign ( numshells+1, 0.0 );

	for (idx = 0; idx < numshells; idx++)
	{
		optpropheights.at(idx)   = shellalts.at(idx);
	}
	optpropheights.at(numshells) = shellalts.back() + 1.0;


	ok    =       m_rayregionmanager.SetGroundAltitude(minshellheight);
	ok    = ok && m_raytracingspecs.ConfigureRayTracingShellAlts     ( &shellalts.front(),        shellalts.size()   );
	ok    = ok && m_opticalpropspecs.ConfigureOpticalPropertyShells	 ( &optpropheights.front(),    optpropheights.size() );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKOCCULT_Specs_User::ConfigureUserDefinedShells, Error configuring the user defined shell specifications. This is a problem");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::CreateCoordinateSystem		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::CreateCoordinateSystem( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* usercoords ) const
{
	bool	ok;

	ok = m_rayregionmanager.MakeCoordinateSystem( usercoords, m_raytracingspecs.GroundAltitude(), m_raytracingspecs.TOAAltitude() );
	return ok;;

}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::MnualRefPtIsDefined		 2014- 5- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::ManualRefPtIsDefined() const
{
	return NXFINITE( m_manualrefpt.latitude) && NXFINITE( m_manualrefpt.longitude) && NXFINITE( m_manualrefpt.mjd);
}

/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::UpdateUndefinedParametersFromLinesOfSight		 2014- 4- 28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::UpdateUndefinedParametersFromLinesOfSight	( const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool ok = true;

	if (m_tanpoint_los.size() > 0)							// If the user has defined some tangent altitude lines of sight
	{														// then 
		ok = ManualRefPtIsDefined();						// we must make sure the reference point has been manually defined
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKOCCULT_Specs_User::UpdateUndefinedParametersFromLinesOfSight, You must define the latitude, longitude and mjd of the Reference Point if you use lines of sight defined from tangent altitude");
		}
	}

	if (ManualRefPtIsDefined()) ok = ok && m_rayregionmanager.SetReferencePoint( m_manualrefpt.latitude, m_manualrefpt.longitude, m_manualrefpt.heightm, m_manualrefpt.mjd );
	if (m_manualsun.IsValid() ) ok = ok && m_rayregionmanager.SetSun           ( m_manualsun );
	ok = ok && m_rayregionmanager.UpdateUndefinedParametersFromLinesOfSight(linesofsight);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::SetReferencePointManually		 2014- 5- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::SetReferencePointManually( double latitude, double longitude, double mjd )
{
	m_manualrefpt.latitude  = latitude;
	m_manualrefpt.longitude = longitude;
	m_manualrefpt.heightm   = 0.0;
	m_manualrefpt.mjd       = mjd;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::SetSunManually		 2014- 5- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::SetSunManually( const nxVector& sun )
{
	m_manualsun = sun;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKOCCULT_Specs_User::AddLineOfSightFromTangentAlt		 2014- 5- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKOCCULT_Specs_User::AddLineOfSightFromTangentAlt( double tangentaltitude, double observeraltitude)
{
	SKOCCULT_LOSFromTangentPoint entry;

	entry.TangentAltitude = tangentaltitude;
	entry.ObserverAltitude = observeraltitude;
	m_tanpoint_los.push_back( entry);
	return true;
}
