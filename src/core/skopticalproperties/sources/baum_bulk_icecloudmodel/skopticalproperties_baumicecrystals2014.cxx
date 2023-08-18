#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::skOpticalProperties_BaumIceCrystals2014		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_BaumIceCrystals2014::skOpticalProperties_BaumIceCrystals2014()
{
	m_effectivesize = nullptr;
	ResetCurrentValues(0.0);
	m_useEddingtonApproximation = true;
	m_backgroundatmosphere = nullptr;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::~skOpticalProperties_BaumIceCrystals2014		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_BaumIceCrystals2014::~skOpticalProperties_BaumIceCrystals2014()
{
	if (m_effectivesize != nullptr ) m_effectivesize->Release();
	if (m_backgroundatmosphere != nullptr) m_backgroundatmosphere->Release();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::SetEffectiveSizeClimatology		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::SetEffectiveSizeClimatology( skClimatology* effectivesize)
{
	bool	ok;

	if (  effectivesize != nullptr) effectivesize->AddRef();
	if (m_effectivesize != nullptr) m_effectivesize->Release();
	m_effectivesize = effectivesize;
	ResetCurrentValues(0.0);
	
	ok = (m_effectivesize == nullptr);																// IF we have NULL passed in then we are good to go
	if (!ok)																						// otherwise
	{																								// check to make  sure the
		ok = m_effectivesize->IsSupportedSpecies( SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS );			// object supprts the EFFECTIVE SIZE parameter.
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::SetEffectiveSizeClimatology, The climatology passed in does not support SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS.");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::ResetCurrentValues		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_BaumIceCrystals2014::ResetCurrentValues( double de)
{
	double	nan = std::numeric_limits<double>::quiet_NaN();
	m_current_De             = de;
	m_current_wavenumber     = nan;
	m_current_absxs          = nan;
	m_current_extxs          = nan;
	m_current_scattxs        = nan;
	m_current_forwardscatter = nan;
	m_p11cachevalid = false;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::SetAtmosphericState		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::SetAtmosphericState( skClimatology* neutralatmosphere)
{
	if (neutralatmosphere      != NULL) neutralatmosphere->AddRef();
	if (m_backgroundatmosphere != NULL) m_backgroundatmosphere->Release();
	m_backgroundatmosphere = neutralatmosphere;
	return true;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::SetAtmosphericState		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )
{
	bool	ok;
	double	De = -9999.0;
	
	ok =       (m_effectivesize != nullptr );
	ok = ok &&  m_effectivesize->GetParameter( SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS, pt, &De, false);
	if (!ok)
	{
		if (m_effectivesize == nullptr)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::SetAtmosphericState, Cannot get effective size, De, as the climatology is not set, see SetEffectiveSizeClimatology");
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::SetAtmosphericState, Error fetching effective size for (lat, lng, height) = %g, %g, %g", (double)pt.latitude,(double)pt.longitude, (double)pt.heightm);
		}
	}
	*crosssectionschanged = !ok || (m_current_De != De);
	if (*crosssectionschanged) ResetCurrentValues(De);
	ok = ok && InternalClimatology_UpdateCache(pt);
	return  ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::InternalClimatology_UpdateCache		2014-4-15*/
/** Updates the internal effective size climatology. It will also load the
 *	Baum Ice crystal database at this time (this can take a few seconds).
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool	ok;
	double	De = 0.0;

	ok =       (m_effectivesize != nullptr );
	ok = ok &&  m_effectivesize->UpdateCache( pt );
	if (!ok)
	{
		if (m_effectivesize == nullptr)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::InternalClimatology_UpdateCache, Cannot update effective size climatology, De, as the climatology is not set, call SetEffectiveSizeClimatology as part of the intialization");
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::InternalClimatology_UpdateCache, Error updating effective size climatology for (lat, lng, height) = %g, %g, %g", (double)pt.latitude,(double)pt.longitude, (double)pt.heightm);
		}
	}
	else
	{
		ok = m_icecrystalsdb.IsLoaded();
		if (!ok) ok = m_icecrystalsdb.LoadDatabase( m_useEddingtonApproximation );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::InternalClimatology_UpdateCache, Error loading the Baum Ice CVrystal Databaseg");
		}
	}
	ok = ok &&  m_effectivesize->GetParameter( SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS, pt, &De, false);
	if( ok )
	{
		if (De != m_current_De)
		{
			ResetCurrentValues(De);
		}
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_BaumIceCrystals2014::InternalClimatology_UpdateCache, Error calculating De" );
		ResetCurrentValues(0.0);
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::CalculateCrossSections		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool	ok;

	static std::mutex	lock;

	lock.lock();
	ok = (wavenumber == m_current_wavenumber);
	if (!ok)
	{
		double	w = 1.0E7/wavenumber;
		ok =       m_icecrystalsdb.InterpolateCrossSections  ( w, m_current_De, &m_current_absxs, &m_current_extxs, &m_current_scattxs);
		ok = ok && m_icecrystalsdb.InterpolateForwardScatter ( w, m_current_De, &m_current_forwardscatter );
		if (ok) m_current_wavenumber = wavenumber;
		else
		{
			ResetCurrentValues(m_current_De);
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::CalculateCrossSections, error fetching crosssection data for wavelenghth = %g, De = %g", (double)w, (double)m_current_De);
		}
	}
	*absxs   = m_current_absxs;
	*extxs   = m_current_extxs;
	*scattxs = m_current_scattxs;
	lock.unlock();
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::CalculatePhaseMatrix		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::CalculatePhaseMatrix( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	double	w = 1.0E7/wavenumber;
	double	angle;

	bool	ok;

	angle = nxmath::acosd( cosscatterangle );
	ok    = m_icecrystalsdb.InterpolatePhaseMatrix( w, m_current_De, angle, phasematrix);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::CalculatePhaseMatrix, error fetching phase matrix data for wavelenghth = %g, De = %g", (double)w, (double)m_current_De);
	}
	return ok;
}

bool skOpticalProperties_BaumIceCrystals2014::CalculateP11(double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11)
{
	double	w = 1.0E7 / wavenumber;
	double	angle;

	bool	ok;

	if (wavenumber != m_current_wavenumber)
	{
		m_p11cachevalid = false;
	}

	if (!m_p11cachevalid && m_phasegridhint->size() > 0)
	{
		m_p11cache.resize(m_phasegridhint->size());
		for (int idx = 0; idx < m_phasegridhint->size(); idx++)
		{
			m_icecrystalsdb.InterpolatePhaseScalar(w, m_current_De, nxmath::acosd(m_phasegridhint->at(idx)), m_p11cache[idx]);
		}
		m_p11cachevalid = true;
	}

	if (m_p11cachevalid)
	{
		p11 = m_p11cache[cosscatterandindex.second];
		return true;
	}

	angle = nxmath::acosd(cosscatterandindex.first);
	ok = m_icecrystalsdb.InterpolatePhaseScalar(w, m_current_De, angle, p11);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_BaumIceCrystals2014::CalculateP1, error fetching phase matrix data for wavelenghth = %g, De = %g", (double)w, (double)m_current_De);
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::DeltaFunctionForwardScatterFraction		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

double skOpticalProperties_BaumIceCrystals2014::DeltaFunctionForwardScatterFraction() const
{
	NXASSERT(( NXFINITE(m_current_forwardscatter) ));
	return m_current_forwardscatter;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::IsScatterer		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::IsScatterer() const
{
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::IsAbsorber		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::IsAbsorber() const
{
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014::SetForwardScatterCutoffAngle		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_BaumIceCrystals2014::SetForwardScatterCutoffAngle( double cutoff_degrees ) 
{
	bool	ok;

	ok = !m_icecrystalsdb.IsLoaded();
	ok = ok && m_icecrystalsdb.SetForwardScatterCutoffAngle( m_useEddingtonApproximation ? cutoff_degrees : 0.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_BaumIceCrystals2014::SetForwardScatterCutoffAngle, You cannot set the forward scatter cutoff angle once the Ice database is loaded. Call this eralier in the initialization, (before  InternalClimatology_UpdateCache)");
	}
	return ok;
}


bool skOpticalProperties_BaumIceCrystals2014::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff)
{
	if (m_icecrystalsdb.IsLegendreCached())
	{
		std::vector<double> internalmoments;

		m_icecrystalsdb.InterpolateLegendre(wavenumber, m_current_De, internalmoments);

		opticalmaxcoeff = std::min(internalmoments.size(), (size_t)usermaxcoeff);
		for (size_t idx = 0; idx < opticalmaxcoeff; idx++)
		{
			coeff[idx] = internalmoments[idx] / (2*idx + 1);
		}
		return true;
	}
	else
	{
		return skOpticalProperties::LegendreCoefficientsP11(wavenumber, coeff, usermaxcoeff, opticalmaxcoeff);
	}


}
