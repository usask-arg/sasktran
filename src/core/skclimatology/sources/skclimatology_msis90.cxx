#include <skclimatology21.h>

/*-----------------------------------------------------------------------------
 *					skClimatology_MSIS90::skClimatology_MSIS90		2005-7-27*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_MSIS90::skClimatology_MSIS90()
{
}


/*-----------------------------------------------------------------------------
 *					skClimatology_MSIS90::~skClimatology_MSIS90		2005-7-27*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_MSIS90::~skClimatology_MSIS90()
{
}

/*-----------------------------------------------------------------------------
 *					skClimatology_MSIS90::IsSupportedSpecies		2005-6-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_MSIS90::IsSupportedSpecies ( const CLIMATOLOGY_HANDLE& species )
{
	bool	ok;

	ok  =   m_msis.IsSupportedSpecies(species);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_MSIS90::UpdateCache		2005-7-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_MSIS90::UpdateCache( const GEODETIC_INSTANT& placeandtime)
{
	m_msis.UpdateCachedProfiles	( placeandtime.mjd, placeandtime.latitude, placeandtime.longitude );
	SetCacheIsLoaded( true );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_MSIS90::GetParameter		2005-6-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_MSIS90::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	bool	ok;
	double	z;

	if (updatecache || CachedMSIS().CacheIsDirty() || !CacheIsLoaded() ) ok = UpdateCache( placeandtime );
	else             ok = CheckCache();

	if (!ok)
	{
		*value = MissingValue();
	}
	else
	{
		z  = placeandtime.heightm/1000.0;							// Get the height in kilometers
		ok = m_msis.InterpolateToHeight( species, z, value );
		if (ok && species == SKCLIMATOLOGY_PRESSURE_PA)
		{
			*value /= 10.0;		// Convert Dynes to Pascals
		}
	}
	return ok;
}
