#include "../sktran_common.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::SKTRAN_RayTracingReferencePointEstimator		2014-4-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracingReferencePointEstimator::SKTRAN_RayTracingReferencePointEstimator( SKTRAN_RayTracingRegionManager* manager, double targetaltitude, double targetrange)
{
	m_regionmanager   = manager;
	m_primaryscantype = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_UNDEFINED;
	ConfigureTargetAltitude( targetaltitude, targetrange, 5.0);

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::~SKTRAN_RayTracingReferencePointEstimator		2014-4-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracingReferencePointEstimator::~SKTRAN_RayTracingReferencePointEstimator()
{
	m_primaryscantype = SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_UNDEFINED;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::ConfigureTargetAltitude		2014-4-4*/
/** Configures the target altitude and range for limb measurements. A parabola centred on the target altitude is created
 *	whose value at the target altitude is "maxweightattargetaltitude", typically 3 or 5. The
 *	width of the parabol is adjusted so its weight at the end of the range is 1.0. 
 *	A max weight of  weight of 1.0 
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingReferencePointEstimator::ConfigureTargetAltitude( double targetaltitude_meters, double targetrange_meters, double maxweightattargetaltitude)
{
	bool ok;

	ok = (targetrange_meters > 0.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::ConfigureTargetAltitude, The target range (%f) must be greater than zero. setting it to 10,000 (ie 10 km)", (double)targetrange_meters );
		targetrange_meters = 10000;
	}
	else
	{
		ok = (maxweightattargetaltitude >= 1.0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::ConfigureTargetAltitude, The maxweightattargetaltiude (%f) should be unity or greater. setting it to 1.0", (double)targetrange_meters );
			maxweightattargetaltitude = 1.0;
		}
	}
	m_targetaltitude    = targetaltitude_meters;
	m_targetrange       = targetrange_meters;								//!< A rough range indicating the range of interesting height eg 15000, ie. desired range is from 10000m to 40000m (25000-15000 to 25000+15000)
	m_targetquadraticweights[0] = -(maxweightattargetaltitude-1.0)/(m_targetrange*m_targetrange);
	m_targetquadraticweights[1] = maxweightattargetaltitude;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::LimbAltitudeWeight		2014-4-4*/
/** Returns the weight that will be applied to a limb location when trying to find
 *	the "mean" reference point. Altitudes in the target altitude range are
 *	weighted more heavily than lines of sight outside of the region. Lines are 
 *	never weighted less than 1 and lines 
 **/
/*---------------------------------------------------------------------------*/

double SKTRAN_RayTracingReferencePointEstimator::LimbAltitudeWeight( double altitude_meters )
{
	double	x;
	double	w;

	x = (m_targetaltitude - altitude_meters);
	w = m_targetquadraticweights[0]*x*x + m_targetquadraticweights[1];
	if ( w < 1.0) w = 1.0;
	return w;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GuessViewingGeometryMode		2014-4-3*/
/** Guesses the best "line of sight" configuration from the lines of sight.
 *	The configuration follows the VIEWING_TYPES defined in class SKTRAN_LineOfSightEntry_V2
 *	Each line of sight is categorized into one of the 7 viewing types. We then calculate
 *	which viewing type is most popular, ignoring observers in space looking at space, and
 *	use this as the VIEWING TYPE of the set of lines of sight.
 *
 *	\par View type details
 *	the view types stats are stored in local variable geometryviewingstats with the following
 *	indexing
 *		-Entry 0 = Num of lines of sight, Undefined look direction. "Should never happen"
 *		-Entry 1 = Num of lines of sight, "In Space looking at space"
 *		-Entry 2 = Num of lines of sight, "In Space looking at limb"
 *		-Entry 3 = Num of lines of sight, "In Space looking at nadir"
 *		-Entry 4 = Num of lines of sight, "On balloon looking at space"
 *		-Entry 5 = Num of lines of sight, "On balloon looking at limb"
 *		-Entry 6 = Num of lines of sight, "On Balloon looking at nadir"
 *		-Entry 7 = Num of lines of sight, "On or near Ground"
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingReferencePointEstimator::GuessViewingGeometryMode(const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE& primaryscantype )
{
	size_t										idx;
	const SKTRAN_LineOfSightEntry_V2*			entry;	
	size_t										index;
	bool										ok;
	SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE	viewtype;
	std::vector<size_t>							geometryviewingstats;
	std::vector<size_t>::iterator				iter;

	geometryviewingstats.assign( SKTRAN_LineOfSightEntry_V2::NumViewingTypes(), 0);
	for (idx = 0; idx < linesofsight.NumRays(); idx++ )											// For each ray specification
	{
		entry    = linesofsight.Entry(idx);
		viewtype = entry->DefaultViewingType( m_regionmanager->Geoid(),  m_regionmanager->MaxAltitude() );
		index    = SKTRAN_LineOfSightEntry_V2::ViewTypeToInteger( viewtype);
		geometryviewingstats[index]++;
	}
	geometryviewingstats[ (size_t)SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE] = 0;			// Ignore lines of sight looking at space
	iter               = std::max_element( geometryviewingstats.begin(), geometryviewingstats.end() );					// Find the most common line of sight
	index              = iter - geometryviewingstats.begin();															// Get the index of the mostt common view type
	primaryscantype  = SKTRAN_LineOfSightEntry_V2::IntegerToViewType( index );										// and Convert this index to a view type
	ok                 =		(*iter > 0)																				// Make sure our most common has at least one member
							&&  (primaryscantype != SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_UNDEFINED);		// and is not undefined
	return ok;

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::EstimateReferencePoint		2014-4-4*/
/** Uses the lines of sight  to define the location of the reference point. The
 *	reference point is used to fit the osculating sphere and is often the location
 *	of the  primary diffuse profile. The osculating sphere simplifies much of the
 *	radiative transfer calculation but the reference point is one place where it complicates
 *	matters. The issue is that any tangent point or observer not located at the reference
 *	point latitude,longitude will be slightly shifted in altitude when translated to the osculating sphere.
 *	It is possible to change rays in the osculating sphere system so this inconsistency
 *	is removed but it is more difficult to move the observer, keep the tangent point the same 
 *	and not change the solar scattering angle.
  **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingReferencePointEstimator::EstimateReferencePoint(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd )
{
	bool ok = true;

	bool onground;
	ok = ok && m_regionmanager->GetNadirReferencePointOnGround(&onground);

	if (m_primaryscantype == SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_UNDEFINED)
		ok = ok && GuessViewingGeometryMode( linesofsight, m_primaryscantype );												// Guess the most appropriate viewing mode 

	if (ok)
	{
		switch ( m_primaryscantype )
		{
		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE :		// Obververs in space looking up (should never happen) but if it does
				ok = AverageObserverLocation( linesofsight, referencepoint, mjd);				// use the average observers position
				break;
			
		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATLIMB :			// Observer in space looking at limb
				ok = AverageWeightedLimbLocation( linesofsight, referencepoint, mjd);			// Use the average limb location of the limb looking rays
				break;

		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR:			// Observer in space looking down
				if (onground)																		// if it is requested by the user
				{																					// then	
					m_targetaltitude = m_regionmanager->GroundAlt();								// set the target altitude at the ground
					ok = AverageTargetHeightLocation( linesofsight, referencepoint, mjd);			// and make the reference point the average intersection with the ground
				}
				else																				// if it is not requested by the user
				{																					// we explicitly decided not to use average location of intersection with target altitude
					ok = AverageObserverLocation(linesofsight, referencepoint, mjd);				// but use use the average observer location.
				}
				break;
		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE:			// Observer in atmosphere (eg balloon) looking up
//				ok = AverageTargetHeightLocation( linesofsight, referencepoint, mjd);			// We explicitly decided not to use average location of intersection with target altitude
				ok = AverageObserverLocation( linesofsight, referencepoint, mjd);				// but use use the average observer location.
				break;

		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB:			// Observer in atmosphere (eg balloon) looking at limb
				ok = AverageWeightedLimbLocation( linesofsight, referencepoint,mjd);			// use the average limb location of the limb looking rays 
				break;

		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR:			// Observer in atmosphere looking down at nadir
				if (onground)																		// if it is requested by the user
				{																					// then	
					m_targetaltitude = m_regionmanager->GroundAlt();								// set the target altitude at the ground
					ok = AverageTargetHeightLocation( linesofsight, referencepoint, mjd);			// and make the reference point the average intersection with the ground
				}
				else																				// if it is not requested by the user
				{																					// we explicitly decided not to use average location of intersection with target altitude
					ok = AverageObserverLocation(linesofsight, referencepoint, mjd);				// but use use the average observer location.
				}
				break;
		case	SKTRAN_LineOfSightEntry_V2::SKTRAN_VIEWING_TYPE_NEARGROUND:						// Observer is with 5 km of gound 
//				ok = AverageTargetHeightLocation( linesofsight, referencepoint,mjd);			// We explicitly decided not to use average location of intersection with target altitude
				ok = AverageObserverLocation( linesofsight, referencepoint,mjd);				// but use use the average observer location.
				break;

		default:nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::EstimateReferencePoint, The viewing geometry of this set of lines of sight is un-determined. That should not happen!");
				ok = AverageObserverLocation( linesofsight, referencepoint,mjd);
				break;
		};
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::EstimateReferencePoint, There were errors using the specialty algorithms for determining the reference point. You might have to re-adjust the default target height properties. I am just using the average observer's location to get the reference point.");
		AverageObserverLocation( linesofsight, referencepoint, mjd);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::AverageWeightedLimbLocation		2014-4-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_RayTracingReferencePointEstimator::AverageWeightedLimbLocation(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd)
{
	size_t										idx;
	const SKTRAN_LineOfSightEntry_V2*			entry;	
	SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE	viewtype;
	bool										ok;
	double										h;
	double										w;
	double										sum;
	double										sumh;
	nxVector									l;

	referencepoint->SetCoords( 0,0,0);							
	sum    = 0.0;
	sumh   = 0.0;
	(*mjd) = 0.0;
	for (idx = 0; idx < linesofsight.NumRays(); idx++ )											// For each ray specification
	{
		entry    = linesofsight.Entry(idx);																	// Get the line of sight
		viewtype = entry->DefaultViewingType( m_regionmanager->Geoid(),  m_regionmanager->MaxAltitude() );	// get the type of line of sight
		if (viewtype == m_primaryscantype)																	// if this ray is the same type as the primary type
		{																									// then
			m_regionmanager->Geoid().FromTangentPointLocation( entry->Observer(), entry->Look() );			// Get the tangent pont locations
			h = m_regionmanager->Geoid().Height();															// get the tangent height in meters
			w = LimbAltitudeWeight(h);																		// get the reference point "altitude" weight for this height
			l = m_regionmanager->Geoid().Location()*w;														// get the weighted tangent point location
			sum    += w;																						// get the weighted sum for normalization
			sumh   += h*w;
			(*mjd) += entry->Mjd()*w;
			(*referencepoint) += l;																			// get the weighted sum of tangent point locations
		}
	}
	ok = (sum > 0);
	if (ok)
	{
		(*referencepoint) /= sum;		// normalize the reference point calculation
		(*mjd)            /= sum;
		sumh              /= sum;
		NXTRACE(("SKTRAN_RayTracingReferencePointEstimator::AverageWeightedLimbLocation, The weighted reference point height is %10.3f kilometers \n", (double)sumh/1000.0));
	}
	else nxLog::Record( NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::AverageWeightedLimbLocation, no lines of sight met the limb viewing criterion,. Thats odd.");	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::AverageObserverLocation		2014-4-4*/
/** Calculate the reference point from the average of the users location. This
 *	is usually a fallback, guaranteed to work, solution. It just averages all of
 *	the observer positions.
 **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_RayTracingReferencePointEstimator::AverageObserverLocation(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd)
{
	size_t										idx;
	const SKTRAN_LineOfSightEntry_V2*			entry;	
	bool										ok;

	referencepoint->SetCoords( 0,0,0);
	*mjd = 0.0;
	ok = linesofsight.NumRays() > 0;
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::AverageObserverLocation, No lines of sight are defined. We cannto defeine a refercen point without lines of sight");
	}
	else
	{
		for (idx = 0; idx < linesofsight.NumRays(); idx++ )											// For each ray specification
		{
			entry    = linesofsight.Entry(idx);
			(*referencepoint) += entry->Observer();
			(*mjd)            += entry->Mjd();
		}
		(*referencepoint) /= (double)linesofsight.NumRays();
		(*mjd)            /= (double)linesofsight.NumRays();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator::AverageTargetHeightLocation		2014-4-4*/
/** If the user is looking up or looking down in the nadir (not at the limb) then
 *	we should find the average location of where the rays intersect the target altitude.
 **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_RayTracingReferencePointEstimator::AverageTargetHeightLocation(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd)
{
	size_t										idx;
	const SKTRAN_LineOfSightEntry_V2*			entry;	
	SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE	viewtype;
	bool										ok;
	double										w;
	double										sum;
	nxVector									entrypoint;
	nxVector									exitpoint;
	nxVector									l;
	nxVector									obs;
	double										zenang;
	double										hobs;
	bool										good;
	double										lat,lng;
	double										Re;
	double										maxdist;

	referencepoint->SetCoords( 0,0,0);							
	*mjd = 0.0;
	sum = 0.0;

	for (idx = 0; idx < linesofsight.NumRays(); idx++ )														// For each ray specification
	{
		entry    = linesofsight.Entry(idx);																	// Get the line of sight
		viewtype = entry->DefaultViewingType( m_regionmanager->Geoid(),  m_regionmanager->MaxAltitude() );	// get the type of line of sight
		if (viewtype == m_primaryscantype)																	// if this ray is the same type as the primary type
		{																									// then
			zenang = entry->ZenithAngleOfLook();															// Get the nominal zenith angle (actual would use south west and up unit vectors)
			m_regionmanager->Geoid().FromGeocentricVector( entry->Observer());									// Get the location of the observer
			hobs = m_regionmanager->Geoid().Height();														// and get his height
			lat  = m_regionmanager->Geoid().GeodeticLatitude();
			lng  = m_regionmanager->Geoid().GeodeticLongitude();
			m_regionmanager->Geoid().FromGeodetic( lat, lng, 0.0);
			Re   = m_regionmanager->Geoid().Location().Magnitude();
			if (hobs <= m_targetaltitude)
			{
				maxdist = sqrt( nxmath::sqr(Re+m_targetaltitude) - nxmath::sqr(Re+hobs) );
				good    =  (zenang <= 90.0); 																	// If the observer is below the target altitude then he must be looking up.
				if (good)
				{
					good    =  good && m_regionmanager->Geoid().GetShellHeightLocation( m_targetaltitude, entry->Observer(), entry->Look(), &entrypoint, &l);			// Get the tangent pont locations
					good    =  good && (1.1*maxdist > (l - entry->Observer()).Magnitude() );
					if (!good)
					{
						nxLog::Record( NXLOG_WARNING, "SKTRAN_RayTracingReferencePointEstimator::AverageTargetHeightLocation, There was some inconsistency looking up as the calculated reference point is too far away. I'm going to ignore it");
					}
				}
			}
			else
			{																								// else the observer is above the target height
				maxdist = sqrt( nxmath::sqr(Re+hobs) - nxmath::sqr(Re) );
				good =  (zenang > 90.0); 																	// then the observer must be looking down.
				if (good)
				{
					good =  good && m_regionmanager->Geoid().GetShellHeightLocation( m_targetaltitude, entry->Observer(), entry->Look(), &l, &exitpoint);			// Get the tangent pont locations
					good =  good && (1.1*maxdist > (l - entry->Observer()).Magnitude() );
					if (!good)
					{
						nxLog::Record( NXLOG_WARNING, "SKTRAN_RayTracingReferencePointEstimator::AverageTargetHeightLocation, There was some inconsistency looking down as the calculated reference point is too far away. I'm going to ignore it");
					}
				}
			}
			if (good)
			{
				w  = LimbAltitudeWeight(hobs);																	// get the reference point "altitude" weight for this height
				l *= w;																							// apply the weight to the intercept location
				sum += w;																						// get the weighted sum for normalization
				(*mjd)            += entry->Mjd()*w;
				(*referencepoint) += l;																			// get the weighted sum of intercept location 
			}
		}
	}
	ok = (sum > 0);
	if (ok)
	{
		(*referencepoint) /= sum;		// normalize the reference point calculation
		(*mjd)            /= sum;
	}
	else  nxLog::Record( NXLOG_WARNING,"SKTRAN_RayTracingReferencePointEstimator::AverageTargetHeightLocation, no lines of sight intyersected the target height. You should change the ReferencePoint Target Height property.");	
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SKTRAN_RayTracingRegionManager		2010-4-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracingRegionManager::SKTRAN_RayTracingRegionManager()
{
	m_referencepoint.SetInvalid();
	m_referencepoint_targetaltitude = 25000.0;			// The target altitude when determining the reference altitude ie either limb, nadir or zenith measurements, default 25000 meters
	m_referencepoint_targetrange    = 15000.0;			// The altitude range above and below the traget altitude for enhanced weighting of limb viewing lines of sight, default 15000 meters.
	m_nadir_referencepoint_onground = false;

	m_sza    = -9999;
	m_minsza = -9999;
	m_maxsza = - 9999;
	m_mjd    = -1.0;
	m_sun.SetInvalid();
	m_maxalt    = 100000.0;													//
	m_minalt    = 0.0;														//!< the minium altitude in meters used when consoidering reference point
	m_groundalt = 0.0;														//!< The ground altitude in meters. Used to determine rays that hit ground
	m_isdirty      = true;
	m_fixedearthradius = std::numeric_limits<double>::quiet_NaN();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::~SKTRAN_RayTracingRegionManager		2010-4-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_RayTracingRegionManager::~SKTRAN_RayTracingRegionManager()
{
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetSun		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetSun( const nxVector& sun )
{
	SetDirty();
	m_sun = sun.UnitVector();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::Clear		2010-4-27*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_RayTracingRegionManager::Clear()
{
	m_sza    = -9999;
	m_minsza = -9999;
	m_maxsza = - 9999;
	m_mjd    = -1.0;
	m_sun.SetInvalid();
	m_maxalt    = 100000.0;													//
	m_minalt    = 0.0;														//!< the minium altitude in meters used when consoidering reference point
	m_groundalt = 0.0;														//!< The ground altitude in meters. Used to determine rays that hit ground
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::CheckParameters		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::CheckParameters() const
{
	bool	ok;
	bool	ok1;
	bool	ok2;

	ok = (m_maxalt > 1000) && ( m_maxalt > m_minalt ) && (m_maxalt > m_groundalt);
	if (!ok)                                      nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::CheckParameters, The user specified altitude ranges are not acceptable");
	ok2 = (m_mjd >= 1000);							if (!ok2) nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::CheckParameters, The user has not defined the mjd");
	ok1 = (m_sun.IsValid());	if (!ok1) nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::CheckParameters, The sun has not been defined");
	ok = ok && ok1 && ok2;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetReferencePoint		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetReferencePoint( double latitude, double longitude, double height_meters, double mjd )
{
	SetDirty();
	m_mjd  = mjd;
	m_geoid.FromGeodetic( latitude, longitude, height_meters);
	m_referencepoint = m_geoid.Location();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetEarthRadius		 2017- 2- 2*/
/** Set the Earth radikus to a manual value. Handy for comparing with codes
 *	that explicitly set the Earth radius.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetEarthRadius( double earthradius_meters)
{
	bool ok;

	ok = (earthradius_meters > 100000.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::SetEarthRadius, it looks like your Earth radius (%f) is not in meters, you should enter a value similar to 6371000.0 for Earth. Nothing has been changed.", (double)earthradius_meters);
	}
	else
	{
		m_fixedearthradius = earthradius_meters;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetReferencePointTargetAltitude		2014-4-8*/
/** Sets the target altitude to be used when determining the reference point. The
 *	target altitude is used in zenith and nadir observations to find the location 
 *	where a ray intersects this altitude. This altitude is then used to find an
 *	average reference point location. The target altitude is used with the
 *	target range variable in limb viewing geometries to apply an increased weighting to 
 *	lines of sight tangential in the vicinity of the target altitude. This encourages
 *	the reference point to be closer to the lines of sight which are tangential in the
 *	region of the reference point.
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetReferencePointTargetAltitude( double targetaltitude_meters)
{
	m_referencepoint_targetaltitude = targetaltitude_meters;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetReferencePointTargetRange		2014-4-8*/
/** Sets the reference point target range parameter. This variable is only used
 *	when calculating the reference point from limb viewing geometries. The parameter
 *	is sued with the 
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetReferencePointTargetRange( double range_meters)
{
	m_referencepoint_targetrange    = range_meters;				//!< The altitude range above and below the traget altitude for enhanced weighting of limb viewing lines of sight, default 15000 meters.
	return true;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetSZARange		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetSZARange( double sza, double minsza, double maxsza )
{
	SetDirty();
	m_sza = sza;
	m_minsza = minsza;
	m_maxsza = maxsza;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetGroundAltitude		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetGroundAltitude( double groundheight_meters )
{
	m_groundalt = groundheight_meters;
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetUpperBoundAltitude		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetUpperBoundAltitude( double upperheight_meters )
{
	m_maxalt = upperheight_meters;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetLowerBoundAltitude		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetLowerBoundAltitude( double lowerheight_meters )
{
	m_minalt = lowerheight_meters;
	return true;
}

bool SKTRAN_RayTracingRegionManager::SetNadirReferencePointOnGround(bool onground)
{ 
	m_nadir_referencepoint_onground = onground;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::SetSunFromTangentPointZenAzi		2009-6-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::SetSunFromTangentPointZenAzi( GEODETIC_INSTANT pt, double solar_zenith, double solar_azimuth )
{
	nxVector	north;
	nxVector	east;
	nxVector	up;
	nxVector	horiz;
	nxVector	tp;
	nxVector	sun;
	

	m_geoid.FromGeodetic( pt.latitude, pt.longitude, pt.heightm);											// Set the location o fthe tangent point
	m_geoid.GetGeodeticWestSouthUp( &east, &north, &up );													// Get west, south and up
	east  = -1.0*east;																						// Convert west to east
	north = -1.0*north;																						// Convert south to North
	horiz = nxmath::cosd(solar_azimuth)*north + nxmath::sind(solar_azimuth)*east;							// Get look vector from tangent point towards satellite
	sun   = nxmath::sind(solar_zenith)*horiz  + nxmath::cosd(solar_zenith)*up;
	return SetSun( sun );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::ConfigureModelGeographicRegion		2008-1-14*/
/**	Configures the coordinates. This sets the following variables,
 *	#- m_sun;	
 *	#- m_latitude;
 *	#- m_longitude;
 *	#- m_mincossza;
 *	#- m_maxcossza;
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::ConfigureGeographicRegionManually( double latitude, double longitude, double mjd, nxVector sun, double deltaszaDegrees)
{
	nxVector	west;
	nxVector	south;
	nxVector	up;
	double		sza;
	double		minsza;
	double		maxsza;
	bool		ok;

	SetReferencePoint( latitude, longitude, 0.0, mjd );
	SetSun( sun );

	m_geoid.FromGeodetic( latitude,longitude,0.0);
	m_geoid.GetGeodeticWestSouthUp ( &west, &south, &up );
	sza    = up.AngleTo( sun );
	minsza = sza - deltaszaDegrees;
	maxsza = sza + deltaszaDegrees;
	if (minsza < 0    ) minsza = 0.0;
	if (maxsza > 180.0) maxsza = 180.0;
	ok =       SetSZARange( sza, minsza, maxsza );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::ConfigureGeographicRegionManually, There was an error with the manual configuration. Thats not good");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21::MoveInsideObserverOutsidee		2009-2-5*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector SKTRAN_RayTracingRegionManager::MoveInsideObserverOutside(const nxVector& observer, const nxVector& look)
{
	nxVector	ofsobs;
	double		dist;

	dist         = 2.1*(m_geoid.SemiMajorAxis()+m_maxalt);											// get diameter of planet plus upper alt
	ofsobs       = observer - dist*look;																// Step back that distance, takes the observer outside the atmosphere
	return ofsobs;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetRayEndpointsObserverInside		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetRayEndpointsObserverInside(const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt)
{
	nxVector	entrance;
	nxVector	exit;
	nxVector	dv;
	nxVector	ofsobs;
	bool		ok;
	bool		hitsground;
	bool		lookingattp;

	ofsobs       = MoveInsideObserverOutside( observer, look);											// Step back that distance, takes the observer outside the atmosphere
	m_geoid.FromTangentPointLocation( ofsobs, look);													// Calculate the tangent point
	dv           = observer - m_geoid.Location();														// Get vector from tangent point to observer
	lookingattp  = (dv & look) < 0;																		// is the observer looking at the tangent point
	hitsground   = (m_geoid.Height() <= m_groundalt);													// Is the tangent ray below ground.

	if (!lookingattp || !hitsground)																	// If we are looking away from tangent point, or the tangent ray does not hit the ground
	{																									// then we use the observer and the exit point
		ok       = m_geoid.GetShellHeightLocation( m_maxalt, ofsobs, look, &entrance, &exit );			// find the entrance and exit points
		*startpt = observer;																			// Ray starts at observer
		*endpt   = exit;																				// and finishes at exit point
	}																									// and that is that
	else																								// if we are looking at tangent point and rays hit ground
	{																									// then
		ok       = m_geoid.GetShellHeightLocation( m_groundalt, ofsobs, look, &entrance, &exit );	// find the entrance and exit points of ground intersection
		*startpt = observer;																			// Ray starts at observer
		*endpt   = entrance;																			// and finishes where tangent ray first intersects ground
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetRayEndpointsObserverOutside		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetRayEndpointsObserverOutside( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt)
{
	double		tpheight;
	nxVector	tplocation;
	nxVector	dummy;
	bool		ok;

	m_geoid.FromTangentPointLocation( observer, look );													// Get the tangent point	
	tplocation = m_geoid.Location();																	// Get the location
	tpheight = m_geoid.Height();																		// Get its altitude
	if (tpheight < m_groundalt)																			// If the observer is looking at the ground
	{																									// Then get the two intercepts, one with upper shell one with ground
		ok =       m_geoid.GetShellHeightLocation( m_maxalt,       observer, look, startpt, &dummy, &tplocation, tpheight );	// start of ray is at upper shell
		ok = m_geoid.GetShellHeightLocation( m_groundalt,    observer, look, endpt,   &dummy, &tplocation, tpheight );
	}																									// and that is that
	else																								// the observer is looking above ground
	{																									// so get the entry and exit points.
		ok =  m_geoid.GetShellHeightLocation( m_maxalt, observer, look, startpt, endpt, &tplocation, tpheight );
	}
//#if defined(NXDEBUG)
//		m_geoid.FromGeocentric(*startpt );
//		NXTRACE(("Start point of ray at (%f, %f, %f)\n", (double)m_geoid.GeodeticLatitude(), (double) m_geoid.GeodeticLongitude(), m_geoid.Height() ));
//		m_geoid.FromGeocentric(*endpt);
//		NXTRACE(("End   point of ray at (%f, %f, %f)\n", (double)m_geoid.GeodeticLatitude(), (double) m_geoid.GeodeticLongitude(), m_geoid.Height() ));
//#endif
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::UpdateSun		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::UpdateSun( )
{
	bool	ok;

	ok = (m_sun.IsValid());
	if (!ok )																		// If the sun is not defined
	{																				// then
		PlanetSun	sun;														// we can calculate the position of the sun
		nxTimeStamp	mjd(m_mjd);													

		ok = (m_mjd > 1000);
		if (ok)
		{
			sun.UpdateECIPosition( mjd );											// Get the ECI position of the sun, accurate to about 1 arc minute
			m_sun = sun.Location().EquatorialToGeographic( mjd ).UnitVector();		// Get the unit vector from Earth to SUN
		}
		else
		{
			nxLog::Record(NXLOG_WARNING," SKTRAN_LineOfSightArray_V21::UpdateSun, Cannot update sun as the user has not defined an mjd");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::UpdateReferencePointSZA		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::UpdateReferencePointSZA( )
{
	nxVector	west;
	nxVector	south;
	nxVector	up;
	m_geoid.FromGeocentricVector( m_referencepoint );
	if (m_sza < -1)																	// if the sza is not yet defined												
	{																				// then 
		m_geoid.GetGeodeticWestSouthUp ( &west, &south, &up );						// Get local up
		m_sza = up.AngleTo(m_sun);													// and set the solar zenith angle
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::UpdateReferencePointAndMJD		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::UpdateReferencePointAndMJD(const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool		ok;
	SKTRAN_RayTracingReferencePointEstimator	refpoint_estimator(this, 25000, 15000);

	ok = (linesofsight.NumRays() > 0);
	if (!ok)
	{
		ok =		( m_mjd > 1000)
				&&  ( m_referencepoint.IsValid());
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_RayTracingRegionManager::FindReferencePoint, No lines of sight are defined and geographic region is not manually defined.  We cannot find a sensible reference point");
		}
	}
	else
	{
		ok = ok && refpoint_estimator.GuessViewingGeometryMode(linesofsight, m_primaryscantype); // still useful to classify the lines of sight, even if a manual reference point is set
		if (!m_referencepoint.IsValid())
		{
			ok = ok && refpoint_estimator.EstimateReferencePoint(linesofsight, &m_referencepoint, &m_mjd);
		}
	}
#if defined(NXDEBUG)
		m_geoid.FromGeocentricVector(m_referencepoint );
		NXTRACE(("Final, average reference point of ray at (%f, %f, %f)\n", (double)m_geoid.GeodeticLatitude(), (double) m_geoid.GeodeticLongitude(), m_geoid.Height() ));
#endif

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::UpdateMinMaxSZA		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::UpdateMinMaxSZA( const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool		ok;
	bool		ok1;
	nxVector	point[3];			//  0 = entry point, 1 = exit point, 2 = tangent point
	double		sza;
	nxVector	west;
	nxVector	south;
	nxVector	up;
	size_t		i;
	bool		insideatmos;
	nxVector	observer;
	nxVector	look;
	size_t		idx;
	const SKTRAN_LineOfSightEntry_V2* entry;

	NXTRACE_ONCEONLY(firsttimne, ("SKTRAN_RayTracingRegionManager::UpdateMinMaxSZA, We need to account for rays that cross the solar pole as they will under-estimate the zenith angle\n"));
//	nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::UpdateMinMaxSZA, Please be careful with this function. The code is old and needs upgrading");
	ok = (m_minsza >= 0) && (m_maxsza <=180) && (m_minsza <= m_maxsza);
	if (!ok)
	{
		m_minsza = m_sza;
		m_maxsza = m_sza;
		ok = (linesofsight.NumRays() == 0);
		for (idx = 0; idx < linesofsight.NumRays(); idx++)						// For each ray specification
		{																									// in our list of rays
			entry = linesofsight.Entry(idx);
			observer = entry->Observer();
			look     = entry->Look();

			m_geoid.FromGeocentricVector( observer);																// Get the location of the observer
			insideatmos = (m_geoid.Height() <= m_maxalt);												// make sure he is above the atmosphere
			if (insideatmos)																				// if he is inside the atmosphere
			{																								// then step the observer back sufficient distance so hei outside
				ok1 = GetRayEndpointsObserverInside(  observer, look, &point[0], &point[1]);
				point[2] = 0.5*( point[0] + point[1] );
			}																								// and that is that
			else																							// otherwise observer is above
			{
				ok1 = GetRayEndpointsObserverOutside(  observer, look, &point[0], &point[1]);
				point[2] = 0.5*( point[0] + point[1] );
			}
			if (ok1)
			{
				ok = true;
				
				for (i = 0; i < 3; i++)
				{
					m_geoid.FromGeocentricVector( point[i] );
					m_geoid.GetGeodeticWestSouthUp ( &west, &south, &up );
					sza = up.AngleTo( m_sun );
					if (sza < m_minsza)
					{
						m_minsza = sza;
					}
					if (sza > m_maxsza)
					{
						m_maxsza = sza;
					}
				}
				
			}
		}
	}
	return ok;
}

bool SKTRAN_RayTracingRegionManager::GetRayEndpoints( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt )
{
	bool			ok = true;
	bool			insideatmos;

	m_geoid.FromGeocentricVector( observer);																// Get the location of the observer
	insideatmos = (m_geoid.Height() <= m_maxalt);												// make sure he is above the atmosphere
	if (insideatmos)																				// if he is inside the atmosphere
	{																								// then step the observer back sufficient distance so hei outside
		ok = ok && GetRayEndpointsObserverInside(  observer, look, startpt, endpt);
	}																								// and that is that
	else																							// otherwise observer is above
	{
		ok = ok && GetRayEndpointsObserverOutside(  observer, look, startpt, endpt);
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::UpdateUsingRayRegion		2008-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight(  const SKTRAN_LineOfSightArray_V21& linesofsight)
{
	bool	ok;

	ok =	   UpdateReferencePointAndMJD( linesofsight );
	ok = ok && UpdateSun();
	ok = ok && UpdateReferencePointSZA	( );
	ok = ok && UpdateMinMaxSZA			( linesofsight );
	ok = ok && CheckParameters();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight, Error updating the bounds defined by the rays. Thats a problem");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::ReferencePoint		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetReferencePoint(double* latitude, double* longitude) const
{
	bool	ok;
	nxGeodetic	geoid( m_geoid );

	ok = m_referencepoint.IsValid();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::GetReferencePoint, The reference point is not defined. Thats a problem");
		*latitude  = std::numeric_limits<double>::quiet_NaN();
		*longitude = *latitude;
	}
	else
	{
		geoid.FromGeocentricVector( m_referencepoint );
		*latitude = geoid.GeodeticLatitude();			// Changed from m_geoid 2012-10-31, Found by djz828
		*longitude = geoid.GeodeticLongitude();
	}
	return ok;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetSZA		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetSZA( double* refptsza, double* minsza, double* maxsza) const
{
	bool	ok;

	*refptsza = m_sza;
	*minsza   = m_minsza;
	*maxsza   = m_maxsza;

	ok =		(m_minsza >= 0.0)
			&&	(m_maxsza >= 0.0)
			&&  (m_minsza <= 180.0)
			&&  (m_maxsza <= 180.0)
			&&  (m_minsza <= m_maxsza)
			&&  (m_sza    >= 0.0)
			&&  (m_sza    <= 180.0);

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::GetSZA, the solar zenith angles are not in the correct range. Thats a problem");
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetSun		2008-4-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetSun( nxVector* sun) const
{
	bool	ok;

	ok   = m_sun.IsValid();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_RayTracingRegionManager::GetSun, The solar vector is not defined, thats a problem");
	}
	*sun      = m_sun;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetMJD		2008-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetMJD( double* mjd ) const
{
	bool ok;

	NXASSERT(( m_mjd > 10000 ));
	ok = (m_mjd > 10000);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_LineOfSightArray_V21::GetMJD, The object is dirty, please call Update before fetching the coordinate object");
	}
	*mjd = m_mjd;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetNadirReferencePointOnGround		2021-12-09*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetNadirReferencePointOnGround(bool* onground) const
{
	*onground = m_nadir_referencepoint_onground;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::GetPrimaryScanType		2021-12-10*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::GetPrimaryScanType(SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE* primaryscantype) const
{
	*primaryscantype = m_primaryscantype;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::MakeCoordinateSystem		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::MakeCoordinateSystem( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>* usercoordinates, double groundaltitude, double toa_altitude ) const
{
	bool							ok;
	double							latitude;
	double							longitude;
	nxVector						sun;
	std::unique_ptr<SKTRAN_CoordinateTransform_V2>	coordinates(new SKTRAN_CoordinateTransform_V2); 

	ok =       IsProperlyDefined();
	ok = ok && CheckParameters  ();

	if (!ok)
	{
		coordinates.release();
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::MakeCoordinateSystem, Cannot make coordinates as the ray tracing region is nor properly defined. Thats a problem");
	}
	else
	{
		ok =       GetReferencePoint( &latitude, &longitude );
		ok = ok && GetSun( &sun );
		ok = ok && coordinates->ConfigureCoordinates( latitude, longitude, m_mjd, sun );
		ok = ok && coordinates->SetAtmosphereAltitudeBounds( groundaltitude, toa_altitude);
		if (NXFINITE(m_fixedearthradius)) ok = ok && coordinates->ManuallySetOsculatingSphereRadius( m_fixedearthradius );
		coordinates->SetStatic();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::MakeCoordinateSystem, Error making coordinate system");
		}
	}
	*usercoordinates = std::move(coordinates);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager::IsProperlyDefined		2010-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_RayTracingRegionManager::IsProperlyDefined() const
{
	bool	ok;

	ok =    (m_referencepoint.IsValid())
		 && (m_sun.IsValid())
		 && ( m_sza     >= 0 ) && (m_sza    <= 180.0 )
		 && ( m_minsza  >= 0 ) && (m_minsza <= 180.0 )
		 && ( m_maxsza  >= 0 ) && (m_maxsza <= 180.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_RayTracingRegionManager::IsProperlyDefined, the ray tracing region is not properly defined");
	}
	return ok;
}
