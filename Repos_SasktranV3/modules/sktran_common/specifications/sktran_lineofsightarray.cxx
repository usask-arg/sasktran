#include "../sktran_common.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightEntry_V2::Configure		2014-4-3*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_LineOfSightEntry_V2::Configure( const nxVector& observer, const nxVector& lookunit, double mjd )
{
	m_observer = observer;
	m_look     = lookunit.UnitVector();
	m_mjd      = mjd;
	m_viewingtype = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_UNDEFINED;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightEntry_V2::IntegerToViewType		2014-4-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE SKTRAN_LineOfSightEntry_V2::IntegerToViewType( size_t idx)
{
	VIEWING_TYPE	viewtype;

	switch (idx)
	{
	case 1:  viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE; break;
	case 2:  viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATLIMB; break;
	case 3:  viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR; break;
	case 4:  viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE; break;
	case 5:  viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB; break;
	case 6:  viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR; break;
	case 7:	 viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_NEARGROUND; break;
	default: viewtype  = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_UNDEFINED; break;
	};
	return viewtype;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightEntry_V2::ViewTypeToInteger		2014-4-4*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_LineOfSightEntry_V2::ViewTypeToInteger( VIEWING_TYPE viewtype)
{
	return (size_t)( viewtype) -  (size_t)SKTRAN_VIEWING_TYPE_UNDEFINED;;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightEntry_V2::DefaultViewingType		2014-4-3*/
/** Returns the type of viewing geometry. This is used to select how the 
 *	models reference point is chosen for the osculating sphere and 
 *	diffuse profiles.
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE SKTRAN_LineOfSightEntry_V2::DefaultViewingType( nxGeodetic& geoid, double toa_altitude ) const
{
	double								obslat;
	double								obslng;
	double								obsh;
	double								Re;
	bool								inspace;							// True if inspace
	bool								inatmos;							// True if floating in the atmosphere (above 5 km)
	bool								onground;							// True if observer is close to the ground.
	double								maxzenang;	
	double								minzenang;	
	double								zenang;
	bool								lookingatspace;
	bool								lookingatlimb;
	bool								lookingatnadir;																	
	VIEWING_TYPE						index;

	geoid.FromGeocentricVector(m_observer);														// Get the location of the observer
	obslat = geoid.GeodeticLatitude();
	obslng = geoid.GeodeticLongitude();
	obsh   = geoid.Height();
	if (obsh < 0.0)  obsh = 0.0;
	Re     = geoid.Location().Magnitude() - obsh; 

	inspace   = (obsh > toa_altitude);																// Is the observer in space
	inatmos   = !inspace && (obsh > 5000.0);														// In the atmosphere floating on a high altitude balloon) 
	onground  = obsh <= 5000.0;																		// or on or near the ground

	if ( inspace)																				// if observer is in space then
	{																							// Then calculate maximum zenith angles for limb viewing geometery
		minzenang  = 180.0 - nxmath::asind( (Re+toa_altitude)/(Re+obsh) )-0.5;				// Min limb zenith angle is when we look at the top of the atmosphere
		maxzenang  = 180.0 - nxmath::asind( Re/(Re+obsh) ) + 0.5;								// Max limb zenith is when we look at the surface of the Earth
	}
	if ( inatmos)																				// if observer is on a high altitude balloon 
	{																							// Then calculate maximum zenith angles for limb viewing geometery
		minzenang  = 100.0;																		// Min zenith is when we looking just above horizontal
		maxzenang  = 180.0 - nxmath::asind( Re/(Re+obsh) ) + 5.0;								// Max limb zenith is when we looking just below the surface of the Earth (about 60 kms below for a 40 km balloon) 
	}
	if (onground)
	{
		minzenang = 180.0;
		maxzenang = 180.0;
	}

	zenang         = ZenithAngleOfLook();
	lookingatspace = (zenang <= minzenang);
	lookingatlimb  = (zenang >  minzenang) && ( zenang < maxzenang);
	lookingatnadir = (zenang >= maxzenang );

	index = SKTRAN_VIEWING_TYPE_UNDEFINED;
	if (inspace)
	{
		if (lookingatspace) index = SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE;
		if (lookingatlimb ) index = SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATLIMB;
		if (lookingatnadir) index = SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR;
	}
	else if (inatmos)
	{
		if (lookingatspace) index = SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE;
		if (lookingatlimb ) index = SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB;
		if (lookingatnadir) index = SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR;
	}
	else if (onground)
	{
		index = SKTRAN_VIEWING_TYPE_NEARGROUND;
	}
	NXASSERT(( index != SKTRAN_VIEWING_TYPE_UNDEFINED));			// This should never happen if we  have done it correctly
	return index;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::SKTRAN_LineOfSightArray_V21		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_LineOfSightArray_V21::SKTRAN_LineOfSightArray_V21()
{
	m_linesofsight.reserve(200);
	Clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::~SKTRAN_GridSpecsLineOfSightEntry_V2		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_LineOfSightArray_V21::~SKTRAN_LineOfSightArray_V21()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::Clear		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_LineOfSightArray_V21::Clear()
{

	m_linesofsight.clear();
//	SetMjdDirty();

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::AddLineOfSightFromTangentPoint		2008-7-9*/
/** Adds a line of sight using a specified tangent point location and altitude as well as the
 *	geographic bearing of the observer from the tangent point,   The height of the observer can be included. Typically an altitude of 600,000 can be used
 *	for satellite work.
 *
 *	\param tangentpoint
 *	A geodetic instant defining the location of the tangent point.
 *
 *	\param geographicbearingofobserver_degrees
 *	The geographic bearing of the (satellite) observer from the tangent point, due North = 0, East = 90
 *	etc. 
 *
 *	\param heightofobserver_meters
 *	The height of the observer in meters. This must be greater than the height of the tangent pint.
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::AddLineOfSightFromTangentPoint( const GEODETIC_INSTANT& tangentpoint, double geographicbearingofobserver_degrees, double heightofobserver_meters, nxGeodetic* geoid)
{
	nxVector	north;
	nxVector	east;
	nxVector	up;
	nxVector	look;
	nxVector	tp;
	nxVector	entrypoint;
	nxVector	exitpoint;
	nxVector	observer;
	nxVector	dv;
	bool		ok;
	double		theta;
	

	ok = (heightofobserver_meters > tangentpoint.heightm);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::AddLineOfSightFromTangentPoint, The observer must be above the tangent point");
	}
	else
	{
		geoid->FromGeodetic( tangentpoint.latitude, tangentpoint.longitude, tangentpoint.heightm );			// Set the location o fthe tangent point
		tp = geoid->Location();																				// Get the geocentric location
		geoid->GetGeodeticWestSouthUp( &east, &north, &up );													// Get west, south and up
		east  = -1.0*east;																						// Convert west to east
		north = -1.0*north;																						// Convert south to North
		look  = nxmath::cosd(geographicbearingofobserver_degrees)*north + nxmath::sind(geographicbearingofobserver_degrees)*east;	// Get look vector from tangent point towards satellite
		ok    = geoid->GetShellHeightLocation( heightofobserver_meters, tp, look, &entrypoint, &exitpoint );						// Find out where this intercepts the observers height shell
		if (ok)																									
		{
			dv = (exitpoint-tp);
			if (( dv & look) > 0) observer = exitpoint;
			else                  observer = entrypoint;
			dv    = (observer - tp );
			theta = dv.AngleTo( look) ;
			ok = (theta > -0.01) && (theta < 0.01);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "The calculated observer does not lie along the line of sight. Thats not good, I'm going to ignore this ray definition!");
			}
			else
			{
				look = -1.0*look;
				ok = AddLineOfSight( observer, look, tangentpoint.mjd );
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_LineOfSightArray_V21::AddLineOfSightFromTangentPoint,There was an error defining the requested line of sight");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::SetRaysFromTangentHeightArray		2009-4-24*/
/** Sets the rays for this calculation from an array of tangent heights. Returns true if successful
 *	otherwise returns false. This method is designed for observers who are located outside the atmosphere.
 *	It may not work properly for observers located inside the atmopshere.
 *
 *	\param mjd
 *	The Universal Time of these rays expressed as a Modified Julian Date. 2009-01-01 0:0:0 is represented as 54832.0
 *	This value is only used by climatologies to lookup up particle number densities.
 *
 *	\param lat
 *	The latitude of the tangent point of the set of rays. Expressed in degrees (-90 to +90). All rays are assumed to have a tangent point at
 *	the same latitude and longitude
 *
 *	\param lng
 *	The longitude of the tangent point of the rays. Expressed in degrees, positive east. (-180 to +360.0)
 *
 *	\param sza
 *	The solar zenith angle at the tangent point. Expressed in degrees (0 to 180.0)
 *
 *	\param saa
 *	The geographic azimuth of the sun at the tangent point in local compass coordinates. Expressed in degrees, 0 is North, 90 is East, 180 is South, 270 is West.
 *	
 *	\param double rayazi
 *	The geographic azimuth of all of the rays at their tangent points in local compass coordinates. Expressed in degrees, 0 is North, 90 is East, 180 is South, 270 is West.
 *	
 *	\param tangentheight_meters
 *	An array of #numheights elements that specifies the tangent height of the rays on this one profile. All heights are specified in meters above the surface of the Earth.
 *
 *	\param numheights
 *	The number of elements in array #tangentheight_meters. 	
 *
 *	\param nominal_observerheight_meters
 *	The nominal height of the observer. This height does not need to be particularly accurate. This specific method is intended for observers
 *	located outside/above the atmosphere and this variable should be set big enough to ensure that (eg a value of 600000 will guarantee that the
 *	observer is always located above an atmosphere that stops at 100 km).
 *
 *	\param sunvector
 *	returns the unit vector toward the sun. This should be passed to the user configuration to manually define the sun position
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::SetRaysFromTangentHeightArray( double  mjd, double lat, double lng,
																	    double  sza, double saa,
																		double  rayazi,
																		double* tangentheight_meters, int numheights,
																		double  nominal_observerheight_meters,
																		nxVector*	sunvector)
{
	nxGeodetic	geoid;
	nxVector	west;
	nxVector	south;
	nxVector	up;
	nxVector	sunh;
	nxVector	sun;
	nxVector	rayh;
	nxVector	point;
	nxVector	base;
	nxVector	observer;
	double		R0;
	double		RL;
	bool		ok = true;
	int			i;
	nxVector	looktoobserver;
	nxVector	entrypoint;
	nxVector	exitpoint;
	nxVector	lostangentpoint;
	bool		useexitpoint;


	Clear();
	geoid.FromGeodetic( lat,lng,0.0);
	geoid.GetGeodeticWestSouthUp ( &west, &south, &up );
	sunh  = -nxmath::cosd(saa)*south - nxmath::sind(saa)*west;			// Horizontal component towards the sun
	sun   = nxmath::cosd(sza)*up + nxmath::sind(sza)*sunh;				// Unit vector towards the sun

	rayh = -nxmath::cosd(rayazi)*south - nxmath::sind(rayazi)*west;			// Unit vector parallel to ray 
	looktoobserver = -rayh;
	base  = geoid.Location();												// vector to surface of oblate spheroid
	R0    = base.Magnitude();												// Distance to oblate spheroid
	RL    = (R0 + nominal_observerheight_meters);							// Distance of observer from center of spherical Earth (no correction for oblateness)

	/* Need to get nominal observer location. Dont Assume
	 * a spherical Earth and use spherical trig */

	for (i=0; i < numheights; i++ )											// for each tangent point
	{
		point = base + tangentheight_meters[i]*up;							// Get the exact location of the tangent point
		geoid.GetShellHeightLocation(nominal_observerheight_meters, point, looktoobserver, &entrypoint, &exitpoint, &point, tangentheight_meters[i]);	//!< Calculate where a look direction enters and exits  a specific shell height.

		useexitpoint = (exitpoint - point).Dot(looktoobserver) > 0.0;
		if (!useexitpoint)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::SetRaysFromTangentHeightArray, Hmmm we should be using the exit point, lets try the entrance point");
			exitpoint = entrypoint;
			useexitpoint = (exitpoint - point).Dot(looktoobserver) > 0.0;
			if (!useexitpoint) 
			{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::SetRaysFromTangentHeightArray, Entry was no good as well. Thats not good");
				ok = false;
			}
		}
#if defined(NXDEBUG)
		double		obsh;
		double		th;
		geoid.FromGeocentricVector(exitpoint);
		obsh = geoid.Height();
		geoid.FromTangentPointLocation( exitpoint, rayh );
		th   = geoid.Height();
#endif
		ok = ok && AddLineOfSight(  exitpoint, rayh, mjd );
	}
	if (sunvector != NULL) *sunvector = sun;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::SetRaysFromTangentHeightArray, Error configuring specs, results are undefined");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::AddEquatorialLineOfSight		 2014- 4- 25*/
/**	Assumes the sun is in the direction (1,0,0) and places tangent point on the equator. 
 *  Should only be used if it is the only function used to add LOS's to the LOSArray; 
 *	other methods to add lines of sight do not assume the same sun position. 
 *	Written for quick testing of the monte carlo engine. Seth Dueck, 2012/06/14 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::AddEquatorialLineOfSight(double sza, double saa, double tanheight, double satAltitude, double mjd){

	bool		ok;
	nxGeodetic	geoid;
	nxVector	west;
	nxVector	south;
	nxVector	up;
	nxVector	tanpoint;
	nxVector	look;
	nxVector	observer;
	double		ell;
	double		r_earth;

	geoid.FromGeodetic(0.0, sza, 0.0);
	geoid.GetGeodeticWestSouthUp ( &west, &south, &up );

	tanpoint = geoid.Location() + tanheight*up;

	look = -south*nxmath::sind(saa) + west*nxmath::cosd(saa);
	if(nxmath::sind(sza)< 0.0) look *= -1.0;

	r_earth = geoid.Location().Magnitude();
	ell = (r_earth+satAltitude)*(r_earth+satAltitude) - tanpoint.Magnitude()*tanpoint.Magnitude();
	NXASSERT( (ell > 0) );			// make sure the satelite altitude is high enough
	ell = sqrt(ell);
	
	observer = tanpoint - ell*look;

	ok = AddLineOfSight( observer, look, mjd );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::AddEquatorialLineOfSight, Error configuring specs, results are undefined");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::AddLineOfSight		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::AddLineOfSight( const nxVector& observer, const nxVector& lookunit, double mjd)
{
	SKTRAN_LineOfSightEntry_V2	entry;

	entry.Configure( observer, lookunit, mjd );
	m_linesofsight.push_back(entry);
//	SetMjdDirty();
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::DeepCopy		2010-5-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::DeepCopy( const SKTRAN_LineOfSightArray_V21& other )
{
//	m_mjd = other.m_mjd;
	m_linesofsight.resize(other.m_linesofsight.size());
	std::copy( other.m_linesofsight.begin(), other.m_linesofsight.end(), m_linesofsight.begin() );
	return (m_linesofsight.size() == other.m_linesofsight.size());
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::GetRay		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::GetRay( size_t idx, const SKTRAN_LineOfSightEntry_V2** entry ) const
{
	bool	ok;
	iterator	iter;

	ok = (idx < NumRays());
	if (ok)
	{
		*entry = &m_linesofsight.at(idx);
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_LineOfSightArray_V21::GetRay, Requested ray is out of bounds or the configuraion is dirty, call update or adjust index bounds");
		*entry = NULL;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::GetRayVar		2013-6-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_LineOfSightArray_V21::GetRayVar( size_t idx, SKTRAN_LineOfSightEntry_V2** entry )
{
	bool	ok;
	iterator	iter;

	ok = (idx < NumRays());
	if (ok)
	{
		*entry = &m_linesofsight.at(idx);
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_LineOfSightArray_V21::GetRayVar, Requested ray is out of bounds or the configutaion is dirty, call update or adjust index bounds");
		*entry = NULL;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::MeanMJDNoCache		2011-5-25*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_LineOfSightArray_V21::MeanMJD() const
{
	size_t	npts;
	size_t	idx;
	double mjd;

	mjd = 0.0;												// from the supplied lines of sight.
	npts  = m_linesofsight.size();
	if (npts > 0)
	{
		for (idx = 0; idx < m_linesofsight.size(); idx++)
		{
			mjd += m_linesofsight[idx].Mjd();
		}
		mjd /= npts;
	}
	else
	{
		mjd = std::numeric_limits<double>::quiet_NaN();
	}
	return mjd;
}

