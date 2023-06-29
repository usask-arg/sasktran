#include <skclimatology21.h>
#include <algorithm>

/*
template< class ITERTYPE>
bool ContainerIsAscendingOrder( ITERTYPE first, ITERTYPE last)
{
	bool		ok = true;
	ITERTYPE	next(first);

	++next;
	while (ok && !(next == last))
	{
		ok = ok && (*next >= *first );
		++first;
		++next;
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table::UserDefined_LatLon_Table		2013-09-11*/
/** **/
/*---------------------------------------------------------------------------*/

UserDefined_LatLon_Table::UserDefined_LatLon_Table()
{

}

/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table::~UserDefined_LatLon_Table		2013-09-11*/
/** **/
/*---------------------------------------------------------------------------*/

UserDefined_LatLon_Table::~UserDefined_LatLon_Table()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table::LinearInterpWeights		2013-09-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool UserDefined_LatLon_Table::LinearInterpWeights( double val, 
																  const std::vector<double>& profile, 
																  double* weights, 
																  size_t* indices, 
																  size_t* numindex ) const
{
	bool ok = true;
	
	std::vector<double>::const_iterator	low,up;
	up = std::upper_bound( profile.begin(), profile.end(), val );
	if( up == profile.end() )
	{
		up = profile.end() - 1;
	}
	if( up == profile.begin() )
	{
		low = profile.begin();
	}
	else
	{
		low = up - 1;
	}
	if( fabs(*up - *low) < 1E-8 )
	{
		weights[0] = 1;
		indices[0] = low - profile.begin();
		*numindex = 1;
	}
	else if ( val > *up )
	{
		weights[0] = 1;
		indices[0] = up - profile.begin();
		*numindex = 1;
	}
	else
	{
		weights[0] = ((val - *low) / (*up - *low));
		weights[1] = ((*up - val) / ( *up - *low));
		indices[0] = up - profile.begin();
		indices[1] = low - profile.begin();
		*numindex = 2;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table::ReleaseResources		2013-09-11*/
/** **/
/*---------------------------------------------------------------------------*/

void UserDefined_LatLon_Table::ReleaseResources()
{
	m_alts.clear();
	m_lats.clear();
	m_lons.clear();
	m_profile.erase();
}


/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table::InterpTable		2014-4-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool UserDefined_LatLon_Table::InterpTable(double alt, double lat, double lon, double* value, double badvalue ) const
{
	bool ok = true;

	double altweights[2];
	size_t altindices[2];
	double latweights[2];
	size_t latindices[2];
	double lonweights[2];
	size_t lonindices[2];

	size_t numalt, numlat, numlon;

	double inrangelon = lon;
	if( inrangelon < 0 ) inrangelon += 360;

	ok = ok && LinearInterpWeights( alt, m_alts, altweights, altindices, &numalt );
	ok = ok && LinearInterpWeights( lat, m_lats, latweights, latindices, &numlat );
	ok = ok && LinearInterpWeights( inrangelon, m_lons, lonweights, lonindices, &numlon );

	NXASSERT(( numlon <= 2));
	NXASSERT(( numlat <= 2));
	NXASSERT(( numalt <= 2));
	if (ok)
	{
		*value = 0;
		for( size_t latidx = 0; latidx < numlat; latidx++ )
		{
			for( size_t lonidx = 0; lonidx < numlon; lonidx++)
			{
				for( size_t altidx = 0; altidx < numalt; altidx++ )
				{
					(*value) += m_profile.At( altindices[altidx] , lonindices[lonidx], latindices[latidx] ) * altweights[altidx] * latweights[latidx] * lonweights[lonidx];
				}
			}
		}
	}
	else
	{
		*value  = badvalue;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					UserDefined_LatLon_Table::LoadProfileFromData		2013-09-11*/
/** Adds a profile to the latitude/longitude mesh.  All species must be on the same
 *  grid, each call to this function will overwrite the alts/lats/lons values.
 *
 *  \param profile
 *		3d array (alts, lons, lats) of profile values.  Total number of values is
 *		alts.size()*lons.size()*lats.size()
 *  \param alts
 *		vector of altitudes in m, must be monotomically increasing
 *  \param lons
 *		vector of longitudes in degrees, must be monotomically increasing in the range 0-360
 *  \param lats
 *		vector of latitudes in degrees, must be monotomically increasing in the range -90 - 90
 **/
/*---------------------------------------------------------------------------*/

bool UserDefined_LatLon_Table::LoadProfileFromData( const std::vector<double>& alts,
																  const std::vector<double>& lons,
																  const std::vector<double>& lats,
																  const nx3dArray<double>& profile
																 )
{
	bool ok = true;

	m_profile        = profile;
	m_alts			 = alts;
	m_lats			 = lats;
	m_lons			 = lons;
	ok =     (m_alts.size() == m_profile.XSize())
		  && (m_lats.size() == m_profile.ZSize())
		  && (m_lons.size() == m_profile.YSize());

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"UserDefined_LatLon_Table::LoadProfileFromData, Data and dimension sizesare inconsistent, Data = (%z,%z,%z), alts = %z, lons=%z, last = %z",
			          (size_t) m_profile.XSize(),(size_t) m_profile.YSize(),(size_t) m_profile.ZSize(),
					  (size_t)m_alts.size(), (size_t)m_lons.size(), (size_t)m_lats.size());
	}
	else
	{
		ok =  (m_lons.size() == 0) || ((m_lons.front() == 0.0) && (m_lons.back() == 360.0));
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"UserDefined_LatLon_Table::LoadProfileFromData, The first longitude (%g) must be 0.0 and the last longitude (%g) must be 360.0", (double)m_lons.front(), (double)m_lons.back());
		}
		if (ok)
		{
			ok = ok && ContainerIsAscendingOrder(m_lons.begin(), m_lons.end());
			ok = ok && ContainerIsAscendingOrder(m_lats.begin(), m_lats.end());
			ok = ok && ContainerIsAscendingOrder(m_alts.begin(), m_alts.end());
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"UserDefined_LatLon_Table::LoadProfileFromData, Alts, lons and lats must be specified in ascending order");
			}
		}
	}
	if (!ok ) ReleaseResources();
	return ok;
}




/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight::skClimatology_UserDefined3D_LatLonHeight		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserDefined3D_LatLonHeight::skClimatology_UserDefined3D_LatLonHeight( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight::~skClimatology_UserDefined3D_LatLonHeight		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserDefined3D_LatLonHeight::~skClimatology_UserDefined3D_LatLonHeight( )
{
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight::LoadProfile		2014-1-30*/
/** Add or replace a species profile in this climatology. Sets up the spline 
*	interpolation of this species. Note that different species in the same
*	climatology do not share the same grid and can be expressed on different
*	grids. Note that the current version demands that longitudes span 0 to 360 and that
*	the first element be 0.0 and the last element be 360.0
*
*	\param species
*		The id code of the species. If this id code already exists in the climatology 
*		it will be replaced otherwise a new profile is created.
*
*	\param alts
*		The heights in meters above sea level at which the variable "profile" is specified.
*		Teh array must be in ascending order (this is checked for by the code).
*
*	\param lons
*		The geodetic longitude of the grid points. Note the current implementation demands that the
*		user provide values for all longitudes in the range (0 to 360.0). The first element must be 0.0
*		and the last element must be 360.0.
*
*	\param lats
*		The geodetic latitude of the grid points. The array must be in ascending order (this is checked for by the code).

*	\param profile
*		The profile of the desired species atthe grid points. It is a 3-D array profile( als, lons, lats).
*
**/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefined3D_LatLonHeight::LoadProfile( const CLIMATOLOGY_HANDLE&		species,
														    const std::vector<double>&		alts,
															const std::vector<double>&		lons,
															const std::vector<double>&		lats,
															const nx3dArray<double>&		profile,
															double							badval)
{

	value_type							v(species, UserDefined3D_LatLonHeightEntry() );			// get a blank entry
	UserDefined3D_LatLonHeightEntry*	entryptr;
	std::pair<iterator, bool>			status;
	iterator							iter;
	bool								ok;
	

	status = m_species.insert(v);				// Try to do the insertion
	ok = status.second;							// See if it worked
	if (!ok)									// If it failed then the component may already exist
	{											// So delete
		m_species.erase(species);				// any existing entry
		status = m_species.insert(v);			// and try again
		ok = status.second;
	}
	if (ok)
	{

		iter = status.first;
		entryptr = &((*iter).second);
		ok = entryptr->SetBadValue(badval);
		ok = entryptr->TableVar()->LoadProfileFromData( alts, lons,lats,profile);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight::UpdateCache		2014-1-30*/
/** Nothing to do as our cache is not place dependent**/
/*---------------------------------------------------------------------------*/


bool skClimatology_UserDefined3D_LatLonHeight::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	bool	ok = true;

	std::for_each( m_species.begin(), m_species.end(), [&ok]( value_type& entry)			// iterate through all entries using a lambda function to check that each member is defined
	{
		ok = ok && entry.second.Table().IsDefined();
	}
	);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skClimatology_UserDefined3D_LatLonHeight::UpdateCache, Some of the user defined entries are blank. Thats not good"); 
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight::GetParameter		2014-1-30*/
/** Fetch the value of the desired species at the specified location.
 *
 *	\param species
 *		ID of the species of interest
 *
 *	\param placeandtime
 *		The location at whichthe value is required.  The heightm, longitude and
 *		latitude fields of placeandtime are used.
 *
 *	\param value
 *		Returns the value of the species at the designated height, latitude and longitude. 
 *		The value will return NaN if the requested species does not exist in this climatology.
 *
 *	\param updatecache
 *		This parameter is ignored as we do not have any time and (ground) location dependency.
 **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefined3D_LatLonHeight::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	iterator	iter;
	bool		ok;		
	double		badval;
	
	iter =  m_species.find( species );
	ok   = !(iter == m_species.end());
	if (ok)
	{
		badval = (*iter).second.BadValue();
		ok = (*iter).second.Table().InterpTable( placeandtime.heightm, placeandtime.latitude, placeandtime.longitude, value, badval);
	}
	else
	{
		*value = std::numeric_limits<double>::quiet_NaN();
		nxLog::Record(NXLOG_WARNING,"skClimatology_UserDefined3D_LatLonHeight::GetParameter, requested species is not supported by this object");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefined3D_LatLonHeight::IsSupportedSpecies		2014-1-30*/
/** Returns true if the requested species is supported by this climatology.
**/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefined3D_LatLonHeight::IsSupportedSpecies( const CLIMATOLOGY_HANDLE& species )
{
	iterator					iter;
	bool						ok;

	iter =  m_species.find( species );
	ok   = !(iter == m_species.end());
	return ok;
}


