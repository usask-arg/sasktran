#include <skclimatology21.h>



/*-----------------------------------------------------------------------------
 *					UserTableSplineEntry::UserTableSplineEntry		2014-1-30*/
/** **/
/*---------------------------------------------------------------------------*/

UserTableSplineEntry::UserTableSplineEntry()
{
	m_badvalue    = 0.0;
	m_dolog       = false;
	m_dopiecewiselinear = false;
}


/*-----------------------------------------------------------------------------
 *					UserTableSplineEntry::~UserTableSplineEntry		2014-1-30*/
/** **/
/*---------------------------------------------------------------------------*/

UserTableSplineEntry::~UserTableSplineEntry()
{
}


/*-----------------------------------------------------------------------------
 *					UserTableSplineEntry::CheckHeightsAreAscending		 2014- 10- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool UserTableSplineEntry::CheckHeightsAreAscending(  const std::vector<double>& h_meters ) const
{
	double lasth = -9.0E30;
	bool	ok = true;

	for (auto iter = h_meters.begin(); !(iter == h_meters.end()); ++iter )
	{
		ok = ok && (*iter >= lasth );
		lasth = *iter;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"bool UserTableSplineEntry::CheckHeightsAreAscending, The input heights are not in ascending order. Please fix as that will create problems");
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					UserTableSplineEntry::CreateProfile		2014-1-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool UserTableSplineEntry::CreateProfile( const std::vector<double>& h_meters, const std::vector<double>& profile, bool dologinterpolation, bool dopiecewiselinear, double badvalue)
{
	bool	ok;
	bool	ok1 = true;

	m_badvalue          = badvalue;
	m_dolog             = dologinterpolation;
	m_dopiecewiselinear = dopiecewiselinear;
	ok                  = CheckHeightsAreAscending(  h_meters );
	if (ok)
	{
		m_heights           = h_meters;
		m_profile           = profile;

		if (dologinterpolation)
		{
			size_t				numbad = 0;
			double				v;

			for (size_t i = 0; i< profile.size(); i++)
			{
				v = profile.at(i);
				m_profile.at(i) = (v > 0)? log(v) : std::numeric_limits<double>::quiet_NaN();
				ok1 = ok1 && ( v > 0);
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"UserTableSplineEntry::CreateProfile, Error creating log interpolation tables, Check for zero or negative values in your profile.");
			}
		}
		if (!dopiecewiselinear) ok = ok && m_spline.Configure( h_meters, m_profile, m_badvalue);
	}
	if (!ok)
	{
		m_spline.Clear();
		m_profile.clear();
		m_heights.clear();
		nxLog::Record(NXLOG_WARNING,"bool UserTableSplineEntry::CreateProfile, There were errors and the profile has been set to a blank entry.");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					UserTableSplineEntry::Interpolate		2014-1-30*/
/** **/
/*---------------------------------------------------------------------------*/

double UserTableSplineEntry::Interpolate( double h_meters )
{
	double v;
	static size_t numWarnings = 0;

	if ((h_meters > m_heights.back() || h_meters < m_heights.front()) && numWarnings < 5 && (m_badvalue != m_badvalue))
	{
		if (h_meters > m_heights.back()) {
			nxLog::Record(NXLOG_WARNING, "double UserTableSplineEntry::Interpolate, The requested altitude (%.4e m) is above the highest altitude (%.4e m) in a user defined climatology. Make sure all user defined climatologies reach the top of atmosphere and remove any manualopticalheights that are higher than the top of atmosphere.", h_meters, m_heights.back());
		}
		if (h_meters < m_heights.front()) {
			nxLog::Record(NXLOG_WARNING, "double UserTableSplineEntry::Interpolate, The requested altitude (%.4e m) is below the lowest altitude (%.4e m) in a user defined climatology. Make sure all user defined climatologies reach the top of atmosphere and remove any manualopticalheights that are higher than the top of atmosphere.", h_meters, m_heights.front());
		}
		numWarnings++;
	}

	v = m_dopiecewiselinear ?  nxLinearInterpolate::EvaluateYatX( h_meters, m_heights, m_profile, nxLinearInterpolate::ENUM_MISSINGVALUE, m_badvalue)
		                    : m_spline.Interpolate( h_meters);

	if ((m_dolog) && ( v != m_badvalue) && NXFINITE(v) )
	{
		v = exp(v);
	}
	return v;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::skClimatology_UserTableSpline		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserTableSpline::skClimatology_UserTableSpline( )
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::~skClimatology_UserTableSpline		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserTableSpline::~skClimatology_UserTableSpline( )
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::ReleaseResources		2008-2-29*/
/** **/
/*---------------------------------------------------------------------------*/

void skClimatology_UserTableSpline::ReleaseResources()
{
	m_species.clear();
}

/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::LoadProfile		2014-1-30*/
/** Add or replace a species profile in this climatology. Sets up the spline 
*	interpolation of this species. Note that different species in the same
*	climatology do not share the same height grid and can be expressed on different
*	grids.
*
*	\param species
*		The id code of the species. If this id code already exists in the climatology 
*		it will be replaced otherwise a new profile is created.
*
*	\param h_meters
*		The heights in meters above sea level at which the variable "profile" is specified.
*		Teh array must be in ascending order (this is assumed but not checked for by the code).
*		The spline interpolation is only valid between the first and last value of h_meters.
*		The "badvalue" is returned for interpolation requests outside this range.
*
*	\param profile
*		The profile of the desired species. This array must be the same same as h_meters and should not
*		include NaN values as this will break the spline interpolation code. If log interpolation is
*		requested all values must be greater than 0.
*
*	\param dologinterpolation
*		True if spline interpolation of the log of the profile is required. False if spline interpolation of the profile
*		is required. The code will automatically disable log interpolation if any element of the profile is not greater than 0.0;
*
*	\param badvalue
*		The value returned by GetParameter if the requested height is outside the min/max values of h_meters.
*

**/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserTableSpline::AddProfile( const CLIMATOLOGY_HANDLE& species, const std::vector<double>& h_meters, const std::vector<double>& profile, bool dologinterpolation, bool dopiecewiselinear, double badvalue )
{
	value_type					v(species, UserTableSplineEntry() );
	std::pair<iterator, bool>	status;
	UserTableSplineEntry*		entryptr;
	iterator					iter;
	bool						ok;
	bool						ispositive;

	if (dologinterpolation)
	{
		ispositive = true;
		for (auto& v: profile)
		{
			ispositive = ispositive && ( v > 0.0 );
		}
		dologinterpolation = dologinterpolation && ispositive;
	}
	
	status = m_species.insert(v);				// Try to do the insertion
	ok = status.second;							// See if it worked
	if (!ok)									// If it failed then the component may already exist
	{											// So delete
		m_species.erase(species);				// any existing ientry
		status = m_species.insert(v);			// and try again
		ok = status.second;
	}
	if (ok)
	{

		iter = status.first;
		entryptr = &((*iter).second);
		ok = entryptr->CreateProfile( h_meters, profile, dologinterpolation, dopiecewiselinear, badvalue );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::UpdateCache		2014-1-30*/
/** Nothing to do as our cache is not place dependent**/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserTableSpline::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::GetParameter		2014-1-30*/
/** Fetch the value of the desired species at the specified location.
 *
 *	\param species
 *		ID of the species of interest
 *
 *	\param placeandtime
 *		The location at whichthe value is required.  Only the "heightm" field of
 *		the placeandtime is used.
 *
 *	\param value
 *		Returns the value of the species at the desited height. returns the result 
 *		of the spline interpolation. The spline interpolation will return the badvalue
 *		set on the previous call to AddProfile if the height is outside the spline range. The
 *		value will return NaN if the requested species does not exist in this climatology.
 *
 *	\param updatecache
 *		This parameter is ignored as we do not have any time and (ground) location dependency.
 **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserTableSpline::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	iterator	iter;
	bool		ok;		
	
	iter =  m_species.find( species );
	ok   = !(iter == m_species.end());
	if (ok)
	{
		*value = (*iter).second.Interpolate( placeandtime.heightm );
	}
	else
	{
		*value = std::numeric_limits<double>::quiet_NaN();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserTableSpline::IsSupportedSpecies		2014-1-30*/
/** Returns true if the requested species is supported by this climatology.
**/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserTableSpline::IsSupportedSpecies( const CLIMATOLOGY_HANDLE& species )
{
	iterator					iter;
	bool						ok;

	iter =  m_species.find( species );
	ok   = !(iter == m_species.end());
	return ok;
}

