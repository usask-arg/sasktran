#include <skclimatology21.h>


/*-----------------------------------------------------------------------------
 *					skClimatology::skClimatology		2005-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology::skClimatology()
{
	m_cacheisloaded = nxFALSE;
	m_missingvalue = -9999999.0;
}


/*-----------------------------------------------------------------------------
 *					skClimatology::DeepCopy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool skClimatology::DeepCopy( const skClimatology& other )
//{
//	m_cacheisloaded = other.m_cacheisloaded;
//	m_missingvalue  = other.m_missingvalue;
//	return true;
//}

/*-----------------------------------------------------------------------------
 *					skClimatology::CheckCache		2005-7-29*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology::CheckCache( /*const GEODETIC_INSTANT& placeandtime*/ )
{
	if (!m_cacheisloaded)
	{
		nxLog::Record(NXLOG_WARNING, "skClimatology::CheckCache, There is no cache loaded into this climatology. Use a succesful callto UpdateCache()");
	}
	return m_cacheisloaded;
}

/*-----------------------------------------------------------------------------
 *					skClimatology::GetHeightProfile		2005-7-22*/
/** Implements a default implementation that repeatedly gets the value at a
 *	specific altitude
 **/
/*---------------------------------------------------------------------------*/

bool skClimatology::GetHeightProfile( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, const double* geometricheightm, int numheights, double* value, bool updatecache, size_t* usernumbad)
{
	bool	ok;
	bool	ok1;
	double	v;
	GEODETIC_INSTANT	pt;
	size_t				numbad;

	ok = (numheights > 0);
	if (ok)
	{
			if (updatecache) ok = UpdateCache( placeandtime);
			else             ok = CheckCache();

			if (ok)
			{
				pt = placeandtime;
				numbad = 0;
				for (size_t i = 0; i < numheights; i++)
				{
					pt.heightm = geometricheightm[i];
					ok1 = GetParameter( species,  pt, &v, nxFALSE );
					if (!ok1)
					{
						value[i] = MissingValue();
						numbad++;
					}
					else value[i] = v;
				}
			}
	}
	if (!ok)
	{
		numbad = numheights;
		for ( size_t i = 0; i < numheights; i++) value[i] = MissingValue();
	}
	*usernumbad = numbad;
	return ok;
}

