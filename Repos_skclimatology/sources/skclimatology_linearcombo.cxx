#include <skclimatology21.h>


/*---------------------------------------------------------------------------
 *  skClimatologyLinearCombination::skOpticalPropertiesLinearCombination2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

 skClimatologyLinearCombination:: skClimatologyLinearCombination()
{
	m_climatologies.resize(2, nullptr);
	m_coefficients.resize(2, 0.0);
}


/*---------------------------------------------------------------------------
 *  skClimatologyLinearCombination::~skOpticalPropertiesLinearCombination 2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

 skClimatologyLinearCombination::~ skClimatologyLinearCombination()
{
	ReleaseResources();
}

/*---------------------------------------------------------------------------
 *    skClimatologyLinearCombination::ReleaseResources     2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

void  skClimatologyLinearCombination::ReleaseResources()
{
	for ( auto iter = m_climatologies.begin(); iter != m_climatologies.end(); ++iter)
	{
		if ((*iter ) != nullptr) (*iter)->Release();
	}
	m_climatologies.clear();
	m_coefficients.clear();
}

/*---------------------------------------------------------------------------
 *      skClimatologyLinearCombination::SetFirstClimatology       2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/
bool skClimatologyLinearCombination::SetFirstClimatology( skClimatology* climate)
{
	if (climate != nullptr) climate->AddRef();
	if (m_climatologies[0] != nullptr) m_climatologies[0]->Release();
	m_climatologies[0] = climate;
	return true;
}

/*---------------------------------------------------------------------------
 *      skClimatologyLinearCombination::SetSecondClimatology      2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/
bool skClimatologyLinearCombination::SetSecondClimatology( skClimatology* climate)
{
	if (climate != nullptr) climate->AddRef();
	if (m_climatologies[1] != nullptr) m_climatologies[1]->Release();
	m_climatologies[1] = climate;
	return true;
}


/*---------------------------------------------------------------------------
 *       skClimatologyLinearCombination::CheckHeightProfile       2020-01-09 */
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatologyLinearCombination::CheckHeightProfile(const nx1dArray<double>& heightm, const nx1dArray<double>& f  )
{
	size_t	npts = heightm.size();
	bool	ok   = npts > 0;
	double	maxh  = -9999.0;
	double	maxf  =  -9999.0;
	double	minf  = 1.0E10;
	double  minh  = 1.0E10;

	for (size_t i = 0; i < npts; i++)
	{
		maxh = std::max( maxh, heightm.at(i) );
		minh = std::min( minh, heightm.at(i) );
		maxf = std::max( maxf, f.at(i) );
		minf = std::min( minf, f.at(i) );
	}
	ok = ok && (maxh > 10.0);
	ok = ok && (minf > -0.1);
	ok = ok && (maxf < 1.1);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skClimatologyLinearCombination::CheckHeightProfile, the height profile coeffiecients dont look correct are you sure you have not mixed up the coeffs for h and f. minh = %f, maxh = %f, minf = %f, maxf = %f ", (double)minh, (double)maxh, (double)minf, (double)maxf);
	}
	return ok;

}
/*---------------------------------------------------------------------------
 *  skClimatologyLinearCombination::SetHeightProfileCoeffsOfFirstProperty 2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/
bool  skClimatologyLinearCombination::SetHeightProfileCoeffsOfFirstProperty( const nx1dArray<double>& heightm, const nx1dArray<double>& f )
{
	bool	ok;
	
	CheckHeightProfile(heightm, f  );

	m_f = f.STLVector();
	m_h = heightm.STLVector();
	ok = m_altitudecoeffs.Configure( m_h, m_f);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skClimatologyLinearCombination::SetHeightProfileCoeffsOfFirstProperty, there were errors setting the altitude coefficients. Thats a problem");
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *  skClimatologyLinearCombination::SetComboCoefficients   2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool  skClimatologyLinearCombination::SetComboCoefficients( const GEODETIC_INSTANT& pt)
{
	double h = pt.heightm;
	double f;
	bool	ok;

	f = m_altitudecoeffs.Interpolate(h, 0.0);
	ok = (f >= -0.0000001) && (f < 1.0000001);
	m_coefficients[0] = (m_climatologies[0] != nullptr) ? f   : 0.0;
	m_coefficients[1] = (m_climatologies[1] != nullptr) ? 1-f : 0.0;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," skClimatologyLinearCombination::SetComboCoefficients, The interpolated coefficient at height %e should be between 0 and 1 (we got %e)",(double)h, (double)f);
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *  skClimatologyLinearCombination::SetAtmosphericState    2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatologyLinearCombination::UpdateCache( const GEODETIC_INSTANT& placeandtime)
{
	bool	ok = true;
	bool	ok1;

	for ( auto iter = m_climatologies.begin(); iter != m_climatologies.end(); ++iter)
	{
		if (*iter != nullptr)
		{
			ok1 = (*iter)->UpdateCache(placeandtime);
			ok = ok && ok1;
		}
		else
		{
			if (ok) nxLog::Record(NXLOG_WARNING," skClimatologyLinearCombination::UpdateCache. At least one of the climatologies in the LINEAR_COMBO climatology is not set. Thats a problem and may give subtle bad results.");
			ok = false;
		}

	}
	return ok;
}

/*---------------------------------------------------------------------------
 *          skClimatologyLinearCombination::GetParameter          2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatologyLinearCombination::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	bool	ok = true;
	bool	ok1;
	double	f;
	double	v;

	*value = 0.0;
	ok = (m_climatologies.size() > 0) && SetComboCoefficients(placeandtime);
	if (ok)
	{
		for ( size_t i = 0; i < m_climatologies.size(); ++i)
		{
			f = m_coefficients[i];
			if (f != 0.0)
			{
				ok1 = m_climatologies[i]->GetParameter(species, placeandtime, &v, updatecache);
				if (ok1)
				{
					*value += f*v;
				}
				ok = ok && ok1;
			}
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *       skClimatologyLinearCombination::IsSupportedSpecies       2020-01-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool  skClimatologyLinearCombination::IsSupportedSpecies( const CLIMATOLOGY_HANDLE& species )
{
	bool issupported = (m_climatologies.size() > 0);

	for ( auto iter = m_climatologies.begin(); iter != m_climatologies.end(); ++iter)
	{
		if (*iter != nullptr)
		{
			issupported = issupported && (*iter)->IsSupportedSpecies( species);
		}
	}
	return issupported;
}
