#include <skopticalproperties21.h>
#include "hitran_xs_cache.h"

/*---------------------------------------------------------------------------
 *          hitran_geodetic_point::hitran_geodetic_point          2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/
hitran_geodetic_point::hitran_geodetic_point()
{
	m_height   = -99999.0;
	m_latitude = -99999.0;
	m_longitude= -99999.0;
}


/*---------------------------------------------------------------------------
 *          hitran_geodetic_point::hitran_geodetic_point          2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

hitran_geodetic_point::hitran_geodetic_point(const GEODETIC_INSTANT& point)
{
	FromGeodeticInstant(point);
}

/*---------------------------------------------------------------------------
 *           hitran_geodetic_point::FromGeodeticInstant           2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

void hitran_geodetic_point::FromGeodeticInstant( const GEODETIC_INSTANT& point)
{
	m_height = point.heightm;
	m_latitude = point.latitude;
	m_longitude = point.longitude;
	if (m_longitude < 0.0) m_longitude += 360.0;
}

/*---------------------------------------------------------------------------
 *				hitran_geodetic_point::operator<				  2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool hitran_geodetic_point::operator<(const hitran_geodetic_point& other ) const
{
	bool less;

	less =  m_height < other.m_height;
	if (!less)
	{
		if (m_height == other.m_height)
		{
			less =  m_latitude < other.m_latitude;
			if (!less)
			{
				if (m_latitude == other.m_latitude)
				{
					less = m_longitude < other.m_longitude;
				}
			}
		}
	}
	return less;
}

/*---------------------------------------------------------------------------
 *                          hitran_geodetic_point::operator==     2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool hitran_geodetic_point::operator==(const hitran_geodetic_point& other ) const
{
	bool equal;

	equal =    ( m_height    == other.m_height) 
			&& ( m_latitude  == other.m_latitude)
		    && ( m_longitude == other.m_longitude);
	return equal;
}
/*---------------------------------------------------------------------------
 *      Hitran_CrossSection_Cache::Hitran_CrossSection_Cache      2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/
Hitran_CrossSection_Cache::Hitran_CrossSection_Cache(skOpticalProperties_HitranChemical* parent)
{
	m_parent = parent;
	m_wbegin = nullptr;
	m_wend   = nullptr;
}

/*---------------------------------------------------------------------------
 *             Hitran_CrossSection_Cache::SetLocation             2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_CrossSection_Cache::SetLocation( const GEODETIC_INSTANT& geo_pt )
{
	bool					ok;
	hitran_geodetic_point	pt(geo_pt);
	iterator				iter = m_cached_entries.find( pt);

	ok = !(iter == m_cached_entries.end() );
	if (!ok)
	{
		m_parent->HitranChemical()->UpdateLocation(geo_pt, m_parent->AtmosphericStateClimatology());
		ok  = CreateNewEntry( geo_pt, &iter );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"Hitran_CrossSection_Cache::SetLocation. There were errors creating a cached entry for the requested location");
		}
	}
	m_current_absxs = ok ? &iter->second : &m_blank_entry;
	return ok;
}

/*---------------------------------------------------------------------------
 *           Hitran_CrossSection_Cache::CreateNewEntry            2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_CrossSection_Cache::CreateNewEntry( const GEODETIC_INSTANT& geo_pt, iterator* iter )
{
	bool					ok;
	hitran_geodetic_point	pt(geo_pt);
	nx1dArray<double>		dummy;
	nx1dArray<double>*		absxs;

	auto status = m_cached_entries.insert( value_type(pt,dummy) );
	ok  = status.second;
	if (ok)
	{
		*iter = status.first;
		absxs = &(*iter)->second;
		ok = absxs->SetSize( m_wavenum.size() );
		ok = ok && m_parent->CalculateCrossSectionsArray( m_wavenum.ArrayBasePtr(), (int)m_wavenum.size(), absxs->UnsafeArrayBasePtr(), m_workerextxs.UnsafeArrayBasePtr(), m_workerscatxs.UnsafeArrayBasePtr() );
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *     Hitran_CrossSection_Cache::HasWavenumberAlreadyInCache     2020-09-16 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_CrossSection_Cache::HasWavenumberAlreadyInCache( double wavenum) const
{
	bool ok;

	ok = std::find( m_wbegin, m_wend, wavenum) != m_wend;
	return ok;
}


/*---------------------------------------------------------------------------
 *        Hitran_CrossSection_Cache::SetCachedWavenumbers         2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_CrossSection_Cache::SetCachedWavenumbers( const nx1dArray<double>& wavenumbers )
{
	bool ok;
	double*	start;
	double* end;

	ok = m_wavenum.DeepCopy(wavenumbers);					// Copy over the wavenumbers from the user
	start    = m_wavenum.UnsafeArrayBasePtr();				// Reset the pointer to the start
	end      = start + m_wavenum.size();					// and end of the array
	m_wbegin = start;										// Reset the pointer to the start
	m_wend   = end;											// and end of the array
	std::sort( start, end );								// and sort the wavenumbers into ascending order
	m_workerextxs.SetSize( m_wavenum.size() );
	m_workerscatxs.SetSize( m_wavenum.size() );
	m_cached_entries.clear();								// clear all the location cached entries
	m_current_absxs = &m_blank_entry;						// and reset our current point to be a blacnk entry.
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"Hitran_CrossSection_Cache::SetCachedWavenumbers, There were errors caching the wavenumbers.");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *       Hitran_CrossSection_Cache::CalculateCrossSections        2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_CrossSection_Cache::CalculateCrossSections( double wavenumber, double *absxs, double* extxs, double* scattxs )
{
	const double*		iter;
	bool				ok;
	size_t				idx;

	iter = std::lower_bound( m_wbegin, m_wend, wavenumber);
	ok = !(iter == m_wend) && (*iter == wavenumber);
	if (ok)
	{
		idx = iter - m_wbegin;
		*absxs   = m_current_absxs->at(idx);
		*extxs   = *absxs;
		*scattxs = 0.0;
	}
	else
	{
		*absxs   = std::numeric_limits<double>::quiet_NaN();
		*extxs   = std::numeric_limits<double>::quiet_NaN();
		*scattxs = std::numeric_limits<double>::quiet_NaN();
		nxLog::Record(NXLOG_WARNING,"Hitran_CrossSection_Cache::CalculateCrossSections. Error looking up the requested wavenumber from the cache. Thats indicates a serious problem.");
	}
	return ok;
}

