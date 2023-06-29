#include <skopticalproperties21.h>






/*-----------------------------------------------------------------------------
 *					class Albedo_Entry		2013-6-13*/
/** \internal **/
/*---------------------------------------------------------------------------*/

class Albedo_Entry
{
	public:
		double	albedo;
		double	wavelennm;

	public:
		bool		operator < ( const Albedo_Entry& other ) const { return wavelennm < other.wavelennm;}
};

/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable_DEPRECATE::SetAlbedo		2010-3-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBRDF_VariableAlbedo::SetAlbedo( const double* albedo, const double* wavelennm, size_t npts )
{


	std::vector<Albedo_Entry>	entries;
	Albedo_Entry				entry;
	bool						ok;
	size_t						idx;

	
	entries.reserve( npts );							// We first
	for (idx = 0; idx < npts; idx++)					// Copy the user data
	{													// into
		entry.albedo   = albedo[idx];					// a local storage
		entry.wavelennm = wavelennm[idx];				// array
		entries.push_back(entry);						// and then we 
	}													// sort them
	std::sort( entries.begin(), entries.end() ) ;		// into ascending order

	ok =       m_wavelengths.SetSize( npts );			// and then we copy into the actual class 
	ok = ok && m_albedo.SetSize( npts );				// members

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skRTAlbedo_Variable_DEPRECATE::SetAlbedo, Error Setting the albedo. Thats a problem");
		m_wavelengths.erase();
		m_albedo.erase();
	}
	else
	{
		for (idx = 0; idx < npts; idx++)					// in sorted order.
		{
			m_wavelengths.At(idx) = entries[idx].wavelennm;
			m_albedo.At(idx)      = entries[idx].albedo;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable_DEPRECATE::GetAlbedo		2010-3-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skBRDF_VariableAlbedo::BRDF( double wavelennm, const GEODETIC_INSTANT& /*pt*/, double /*mu_in*/, double /*mu_out*/, double /*dcosphi*/, double* brdf) const
{
	double albedo;

	albedo = nxLinearInterpolate::EvaluateYatX ( wavelennm , m_wavelengths.ArrayBasePtr(), m_albedo.ArrayBasePtr(), m_wavelengths.size(), nxLinearInterpolate::ENUM_TRUNCATE , 0.0 );
	*brdf = albedo/nxmath::Pi;
	return true;
}



/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable_DEPRECATE::SetAlbedo		2010-3-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTAlbedo_Variable_DEPRECATE::SetAlbedo( const double* albedo, const double* wavelennm, size_t npts )
{


	std::vector<Albedo_Entry>	entries;
	Albedo_Entry				entry;
	bool						ok;
	size_t						idx;

	
	entries.reserve( npts );							// We first
	for (idx = 0; idx < npts; idx++)					// Copy the user data
	{													// into
		entry.albedo   = albedo[idx];					// a local storage
		entry.wavelennm = wavelennm[idx];				// array
		entries.push_back(entry);						// and then we 
	}													// sort them
	std::sort( entries.begin(), entries.end() ) ;		// into ascending order

	ok =       m_wavelengths.SetSize( npts );			// and then we copy into the actual class 
	ok = ok && m_albedo.SetSize( npts );				// members

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skRTAlbedo_Variable_DEPRECATE::SetAlbedo, Error Setting the albedo. Thats a problem");
		m_wavelengths.erase();
		m_albedo.erase();
	}
	else
	{
		for (idx = 0; idx < npts; idx++)					// in sorted order.
		{
			m_wavelengths.At(idx) = entries[idx].wavelennm;
			m_albedo.At(idx)      = entries[idx].albedo;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable_DEPRECATE::GetAlbedo		2010-3-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTAlbedo_Variable_DEPRECATE::GetAlbedo( double wavelennm, const GEODETIC_INSTANT& /*pt*/, double* albedo) const
{
	*albedo = nxLinearInterpolate::EvaluateYatX ( wavelennm , m_wavelengths.ArrayBasePtr(), m_albedo.ArrayBasePtr(), m_wavelengths.size(), nxLinearInterpolate::ENUM_TRUNCATE , 0.0 );
	return true;
}



/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Variable_DEPRECATE::CreateClone		2010-3-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTAlbedo_Variable_DEPRECATE::CreateClone( skRTAlbedo_DEPRECATED** clone ) const
{
	skRTAlbedo_Variable_DEPRECATE* copy;
	bool				 ok;


	copy = new skRTAlbedo_Variable_DEPRECATE;

	ok = (copy != NULL);;
	if (ok)
	{
		copy->AddRef();
		copy->m_albedo      = m_albedo;
		copy->m_wavelengths = m_wavelengths;
	}
	*clone = copy;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTAlbedo_Constant_DEPRECATE::CreateClone		2008-11-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTAlbedo_Constant_DEPRECATE::CreateClone( skRTAlbedo_DEPRECATED** clone ) const
{
	skRTAlbedo_Constant_DEPRECATE* copy;
	bool				 ok;


	copy = new skRTAlbedo_Constant_DEPRECATE;

	ok = (copy != NULL);;
	if (ok)
	{
		copy->AddRef();
		copy->m_albedo = m_albedo;
	}
	*clone = copy;
	return ok;
}



