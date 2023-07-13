#include <skopticalproperties21.h>

//#include "pch.h"
//#include <iostream>
//#include "hitranline_upperstates.h"


/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::skEmission_HitranChemical		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

skEmission_HitranChemical::skEmission_HitranChemical()
{
	init();
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::skEmission_HitranChemical		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

skEmission_HitranChemical::skEmission_HitranChemical	(const char* chemicalname, double lowerwavenumber, double upperwavenumber)
{

	init();
	SetChemicalName( chemicalname);
	SetWavenumberRange( lowerwavenumber, upperwavenumber);
}


/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::init		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

void skEmission_HitranChemical::init()
{
	m_xs_optimizer                    = nullptr;
	m_hapicompliant                   = true;
	m_isdirty                         = true;
	m_lowwavenum                      = 0.0;
	m_hihwavenum                      = 9999999.0;
	m_hitranchemical                  = nullptr;
	m_selfbroadeningclimatology       = nullptr;
	m_selfbroadeningclimatologyhandle = SKCLIMATOLOGY_UNDEFINED;
	m_upperstate_numberdensity        = nullptr;
	m_upperstate_numberdensity_handle = SKCLIMATOLOGY_UNDEFINED;	// The GUID of the species for the number density of the excited upper state molecules.
	m_upperstate_numberdensity_value  = 0.0;
	m_lineshapeobject                 = nullptr;
	m_atmospheric_state               = nullptr;
	m_isotopefilterid                 = 0;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::~skEmission_HitranChemical		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

skEmission_HitranChemical::~skEmission_HitranChemical()
{
	if ( m_xs_optimizer              != nullptr) delete m_xs_optimizer;
	if ( m_selfbroadeningclimatology != nullptr) m_selfbroadeningclimatology->Release();
	if ( m_lineshapeobject           != nullptr) m_lineshapeobject->Release();
	if (m_hitranchemical             != nullptr) delete m_hitranchemical;
	if ( m_atmospheric_state         != nullptr) m_atmospheric_state->Release();
	if (m_upperstate_numberdensity   != nullptr) m_upperstate_numberdensity->Release();
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetDirty		2013-3-20*/
/** Flag the obejct as dirty, ie the user has changed some configuration info.
 *	Configuration is normally done during the early stages of the object life. If
 *	for some reasona the user trie sto change the configuration after all of the 
 *	internal objects are created we issue a warning saying that that is not the
 *	right way to do things as we get large computational hits. 
**/
/*---------------------------------------------------------------------------*/

void skEmission_HitranChemical::SetDirty()
{
	m_isdirty = true;
	if ( m_hitranchemical != NULL)
	{
		delete m_hitranchemical;
		m_hitranchemical = NULL;
		nxLog::Record(NXLOG_INFO,"skEmission_HitranChemical::SetDirty, It is very inefficient to reset parameters of an skEmission_HitranChemical object after the HITRAN object is created. It makes more sense to create a brand new instance");
	}
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::CheckDirtyAndUpdate		2013-3-20*/
/**	We use CheckDirtyandUpdate to delay creation of the Hitran object until
 *	the user actually needs to do something with the obejct. This way we avoid
 *	quite large coputational overheads if the object is created but never
 *	used.
 *
 *	Typically the underlying Hitran objects are not created until the user
 *	makes a call to SetAtmopshericState or InternalClimatology_UpdateCache, both of which
 *	must be called before a call to CalculateCrossSections.
**/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::CheckDirtyAndUpdate( const GEODETIC_INSTANT& pt )
{
	bool	ok;
	bool	ok1, ok2;

	ok = !m_isdirty;
	if (!ok)
	{
		NXASSERT( m_hitranchemical == NULL );
		ok = !m_chemicalname.IsEmpty();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::CheckDirtyAndUpdate, No hitran chemical species is defined. Try calling SetChemicalName()");
		}
		else
		{
			if (m_lineshapeobject == NULL)											// IF we have no line shape object passed in by the suer
			{																		// then
				m_lineshapeobject = new skSpectralLineShape_VoigtKuntz;			// Use the tabulated voigt as the default
				m_lineshapeobject->AddRef();
			}
			m_hitranchemical = new skSpectralEmissionCollection_HitranChemical( m_chemicalname, m_lowwavenum, m_hihwavenum, m_hapicompliant, m_isotopefilterid, m_lowerstateglobalquantafilter.c_str(), m_upperstateglobalquantafilter.c_str() );
			
			ok1 = (m_hitranchemical != NULL);
			ok1 = ok1 && m_hitranchemical->SetLineShapeObject ( m_lineshapeobject );
			ok1 = ok1 && m_hitranchemical->SetSelfBroadeningClimatology( m_selfbroadeningclimatologyhandle, m_selfbroadeningclimatology);
			ok1 = ok1 && m_hitranchemical->UpdateCache( pt );
			if (!ok1)
			{
				nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::CheckDirtyAndUpdate, There were errors creating and updating the Hitran Chemical Instance");
			}

			ok2 = ( m_upperstate_numberdensity != nullptr);
			ok2 = ok2 && m_upperstate_numberdensity->UpdateCache(pt);
			if (!ok2)
			{
				if ( m_upperstate_numberdensity == nullptr) nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::CheckDirtyAndUpdate, No climatology is defined for the number of excited upper state molecules. It will default to 0.0");
				else                                        nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::CheckDirtyAndUpdate, error updating cache for the climatology of excited upper state molecules.");
				m_upperstate_numberdensity_value = 0.0;
			}
			ok = ok1 && ok2;
		}
		m_isdirty = !ok;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetAtmosphericState		2013-3-20*/
/** The skOpticalProperties method to SetAtmosphericState. Our spectral line
 *	object will typically need pressure (Pascals) and temperature (K) if they
 *	are based on regular, Doppler, Lorenz or Voigt profiles.
 *
 *	A separate climatology interface (which may or may not use the same actual
 *	skClimatology object) is used to describe chemical number desnity for partial
 *	pressure calculations
 *
 *	The user will typically use a standard "atmospheric" skClimatology class for calls
 *	to SetAtmosphericState, eg, ECMWF, NCEP, MSIS or user defined tables. If you choose to use
 *	a LineShapeObject that needs more info that pressure and temperature then
 *	the atmosphericstate climatology must support that information and the LineShapeObject
 *	must extract that information during this call (and the asscoiated function calls) 
 **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetAtmosphericState( skClimatology* atmosphericstate)
{
	if ( atmosphericstate != m_atmospheric_state)
	{
		if ( atmosphericstate != nullptr) atmosphericstate->AddRef();
		if ( m_atmospheric_state != nullptr) m_atmospheric_state->Release();
		m_atmospheric_state = atmosphericstate;
		SetDirty();
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::UpdateCache( const GEODETIC_INSTANT& pt )
{
	bool	ok;

	SetDirty();
	ok = UpdateLocation( pt, false );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::UpdateLocation( const GEODETIC_INSTANT& pt, bool isground )
{
	bool	ok;
	bool	ok1,ok2;

	ok = CheckDirtyAndUpdate(pt);
	if (m_xs_optimizer != nullptr)
	{
		ok = ok && m_xs_optimizer->SetLocation(pt);
	}
	else
	{
		ok1 = (m_atmospheric_state != nullptr);
		ok1 = ok1 && m_hitranchemical->UpdateLocation(pt, m_atmospheric_state);
		if (!ok1)
		{
			if (m_atmospheric_state == nullptr) nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::UpdateLocation, cannot update location as the atmospheric state object is not defined. Call SetAtmopshericState to fix this");
			else                                nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::UpdateLocation, error calling UpdateLocation of the Hitran Chemical object");
		}

		ok2 = (m_upperstate_numberdensity != nullptr);
		ok2 = ok2 && m_upperstate_numberdensity->GetParameter( m_upperstate_numberdensity_handle, pt, &m_upperstate_numberdensity_value, false);
		if (!ok2) m_upperstate_numberdensity_value = 0.0;

	}
	return ok;
}


/*---------------------------------------------------------------------------
 * skEmission_HitranChemical::FactorToConvertTo_PhotonsPerCM2PerSecPerSteradianPerNM 2020-08-21 */
/**	The HITRAN calculation calculates the fractional number of photons emitted over a 4 pi sphere  
 *	at a given wavenumber for one excited state in photons/sec/cms/wavenumber. We need to
 *	convert that to photons/cm2/sec/steradian/nm for all of the excited upper state
 *	molecules.
 **/
/*---------------------------------------------------------------------------*/

double skEmission_HitranChemical::FactorToConvertTo_PhotonsPerCM2PerSecPerSteradianPerNM( double wavenumber ) const
{
	double factor;
	
	factor = m_upperstate_numberdensity_value*wavenumber*wavenumber*1.0E-07/( 4.0*nxmath::Pi);
	return factor;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::CalculateCrossSectionsInternal		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::CalculateEmissionInternal( double wavenumber, double* signal ) const
{
	double		emission = 0.0;
	bool		ok;

	ok =       (m_hitranchemical != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skEmission_HitranChemical::CalculateEmissionInternal, The internal hitran object is not yet loaded. Try calling SetAtmopshericState first, otherwise get the debugger going");
	}
	ok = ok && m_hitranchemical->Emission(wavenumber, &emission);
	if (!ok)
	{
		emission = std::numeric_limits<double>::quiet_NaN();
	}
	*signal  = emission*FactorToConvertTo_PhotonsPerCM2PerSecPerSteradianPerNM( wavenumber );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::IsotropicEmission( double wavenumber, double* isotropicradiance)
{
	bool ok = false;

	if (m_xs_optimizer != nullptr) ok = m_xs_optimizer->CalculateEmissions( wavenumber, isotropicradiance);		// Are we using the optimized pre-cached cross-section configuration
	else                           ok = CalculateEmissionInternal         ( wavenumber, isotropicradiance);		// Guarantees thread safety as we call the const implementation
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::CheckWavenumberIsAscending		 2014- 10- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::CheckWavenumberIsAscending( const std::vector<double>&	wavenumber) const
{
	double lastval = wavenumber.front() - 1.0;
	bool	ok = true;

	for ( auto& entry: wavenumber)
	{
		ok = ok && entry >= lastval;
		lastval = entry;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::CalculateCrossSectionsArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::IsotropicEmissionArray(const std::vector<double>& wavenumber, std::vector<double>* isotropicradiance )
{
	bool	ok;

	ok = (m_hitranchemical != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skEmission_HitranChemical::IsotropicEmissionArray, The internal hitran object is not yet loaded. Try calling SetAtmopshericState first, otherwise get the debugger going");
	}
	else
	{
		ok = (wavenumber.size() == 0) || ( CheckWavenumberIsAscending(wavenumber) );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::IsotropicEmissionArray, The incoming wavenumber array must be in ascending order for CalculateCrossSectionsArray. Yours is not. Please correct and try again");
		}
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skEmission_HitranChemical::IsotropicEmissionArray, Error calculating cross-sections. Returning a zeroed array");
		isotropicradiance->assign( wavenumber.size(), 0.0);
	}
	else
	{

		isotropicradiance->assign( wavenumber.size(), 0.0);
		ok = m_hitranchemical->EmissionArray( wavenumber, isotropicradiance);
		if (ok)
		{
			for ( size_t i = 0; i < wavenumber.size(); i++)
			{
				isotropicradiance->at(i) *= FactorToConvertTo_PhotonsPerCM2PerSecPerSteradianPerNM( wavenumber[i] );
			}
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *     skEmission_HitranChemical::SetUpperStateNumberDensity      2020-08-21 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetUpperStateNumberDensity( skClimatology* upperstate_numberdensity)
{

	if ( upperstate_numberdensity != m_upperstate_numberdensity)
	{
		if ( upperstate_numberdensity   != nullptr) upperstate_numberdensity->AddRef();
		if ( m_upperstate_numberdensity != nullptr) m_upperstate_numberdensity->Release();
		m_upperstate_numberdensity = upperstate_numberdensity;
		SetDirty();
	}
	return true;
}


/*---------------------------------------------------------------------------
 *  skEmission_HitranChemical::SetUpperStateNumberDensityHandle   2020-08-21 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetUpperStateNumberDensityHandle	( const CLIMATOLOGY_HANDLE& parameterguid )
{
	m_upperstate_numberdensity_handle =  parameterguid;
	SetDirty();
	return true;

}

/*---------------------------------------------------------------------------
 *    skEmission_HitranChemical::SetSelfBroadeningClimatology     2020-08-21 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetSelfBroadeningClimatology( skClimatology* numberdensityclimatology  )
{
	if ( numberdensityclimatology != NULL) numberdensityclimatology->AddRef();
	if ( m_selfbroadeningclimatology    != NULL) m_selfbroadeningclimatology->Release();
	m_selfbroadeningclimatology = numberdensityclimatology;
	SetDirty();
	return true;
}


/*---------------------------------------------------------------------------
 * skEmission_HitranChemical::SetSelfBroadeningClimatologyHandle  2020-08-19 */
/**	Sets the GUID used by the m_selfbroadeningclimatology to calculate the number density  
 *	of this molecule in molecules/cm3. It is iused to calculate partial pressure as part of the
 *	self broadening calculation.
 **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetSelfBroadeningClimatologyHandle( const CLIMATOLOGY_HANDLE& parameterguid  )
{
	m_selfbroadeningclimatologyhandle =  parameterguid;
	SetDirty();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetLineShapeObject		2013-3-21*/
/** Sets the line shape object to be used in calculation of optical properties.
 *	By deafult, if this function is not called (or if NULL is passed in), then
 *	the optical properties will be calculated with skSpectralLineShape_VoigtTabulated
 *	which is pretty good
 *
 *	/param lineshapeobject
 *		The line shape object to be used in spectral line calculations. Several references are
 *		placed on the object but are released when the object is no longer required.
 **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetLineShapeObject( skSpectralLineShape* lineshapeobject )
{
	if (lineshapeobject != NULL) lineshapeobject->AddRef();
	if (m_lineshapeobject != NULL) m_lineshapeobject->Release();
	m_lineshapeobject = lineshapeobject;
	SetDirty();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetChemicalName		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetChemicalName( const char* chemicalname )
{
	m_chemicalname = chemicalname;
	SetDirty();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skEmission_HitranChemical::SetWavenumberRange		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetWavenumberRange( double lowwavenum, double highwavenum )
{
	bool	ok;

	SetDirty();
	if (lowwavenum < highwavenum)
	{
		m_lowwavenum = lowwavenum;
		m_hihwavenum = highwavenum;
	}
	else
	{
		m_lowwavenum = highwavenum;
		m_hihwavenum = lowwavenum;
	}
	ok = (lowwavenum != highwavenum);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::SetWavenumberRange, The low wavenumber (%g) is bigger than the high wavenumber (%g). Thats not good.", (double)m_lowwavenum, (double)m_hihwavenum);
	}
	SetDirty();
	return ok;
}


/*---------------------------------------------------------------------------
 *         skEmission_HitranChemical::SetIsotopeIdFilter          2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetIsotopeIdFilter					( int isotopefilterid )
{
	m_isotopefilterid = isotopefilterid;
	SetDirty();
	return true;
}

/*---------------------------------------------------------------------------
 *   skEmission_HitranChemical::SetLowerStateGlobalQuantaFilter   2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetLowerStateGlobalQuantaFilter		( std::string lowerstateglobalquantfilter)
{
	m_lowerstateglobalquantafilter = lowerstateglobalquantfilter;
	SetDirty();
	return true;
}

/*---------------------------------------------------------------------------
 *   skEmission_HitranChemical::SetUpperStateGlobalQuantaFilter   2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::SetUpperStateGlobalQuantaFilter		( std::string upperstateglobalquantfilter)
{
	m_upperstateglobalquantafilter = upperstateglobalquantfilter;
	SetDirty();
	return true;
}


/*---------------------------------------------------------------------------
 * skEmission_HitranChemical::EnableCachedCrossSections  2019-11-08 */
/** Enables the Hitran_Emissions_Cache class which calculates and caches
 *	cross-sections for the given wavenumbers at every call to SetLocation
 *	after this point.
 **/
/*---------------------------------------------------------------------------*/

bool skEmission_HitranChemical::EnableCachedEmissions( double* wavenumbers, size_t numwave )
{
	bool ok = true;

	SetDirty();
	if ( m_xs_optimizer != nullptr) 
	{ 
		delete m_xs_optimizer; 
		m_xs_optimizer = nullptr;
	}
	if (numwave > 0)
	{
		std::vector<double>		wavenum; //( numwave, wavenumbers );
		double minval = 1.0E20;
		double maxval = -99999.0;
		
		wavenum.assign( wavenumbers, wavenumbers+numwave);

		for (size_t i = 0; i < numwave; i++)
		{
			minval = nxmin( minval, wavenumbers[i]);
			maxval = nxmax( maxval, wavenumbers[i]);
		}
		ok = (maxval > minval) && (minval > 0) && (maxval < 100000.0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::EnableCachedEmissions, The range of wavelengths from %e to %e is invalid.", (double)minval, (double)maxval);
		}
		else
		{
			ok = SetWavenumberRange( minval, maxval);		
			if (ok)
			{
				m_xs_optimizer = new Hitran_Emission_Cache(this);
				ok =       m_xs_optimizer != nullptr;
				ok = ok && m_xs_optimizer->SetCachedWavenumbers( wavenum );
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"skEmission_HitranChemical::EnableCachedEmissions, There were errors enabling cross-section caching. This will need debugging.");
			}
		}
	}
	return ok;

}


/*---------------------------------------------------------------------------
 *      Hitran_Emission_Cache::Hitran_Emission_Cache      2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

Hitran_Emission_Cache::Hitran_Emission_Cache(skEmission_HitranChemical* parent)
{
	m_parent = parent;
}

/*---------------------------------------------------------------------------
 *             Hitran_Emission_Cache::SetLocation             2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_Emission_Cache::SetLocation( const GEODETIC_INSTANT& geo_pt )
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
			nxLog::Record(NXLOG_WARNING,"Hitran_Emission_Cache::SetLocation. There were errors creating a cached entry for the requested location");
		}
	}
	m_current_emission = ok ? &iter->second : &m_blank_entry;
	return ok;
}

/*---------------------------------------------------------------------------
 *           Hitran_Emission_Cache::CreateNewEntry            2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_Emission_Cache::CreateNewEntry( const GEODETIC_INSTANT& geo_pt, iterator* iter )
{
	bool					ok;
	hitran_geodetic_point	pt(geo_pt);
	std::vector<double>		dummy;
	std::vector<double>*	emission;

	auto status = m_cached_entries.insert( value_type(pt,dummy) );
	ok  = status.second;
	if (ok)
	{
		*iter = status.first;
		emission = &(*iter)->second;
		emission->resize( m_wavenum.size(), 0.0 );
		ok = m_parent->IsotropicEmissionArray( m_wavenum, emission );
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *        Hitran_Emission_Cache::SetCachedWavenumbers         2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_Emission_Cache::SetCachedWavenumbers( const std::vector<double>& wavenumbers )
{
	bool ok = true;

	m_wavenum = wavenumbers;								// Copy over the wavenumbers from the user
	std::sort( m_wavenum.begin(), m_wavenum.end());			// and sort the wavenumbers into ascending order
	m_cached_entries.clear();								// clear all the location cached entries
	m_current_emission = &m_blank_entry;					// and reset our current point to be a blacnk entry.
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"Hitran_Emission_Cache::SetCachedWavenumbers, There were errors caching the wavenumbers.");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *       Hitran_Emission_Cache::CalculateCrossSections        2019-11-06 */
/** **/
/*---------------------------------------------------------------------------*/

bool Hitran_Emission_Cache::CalculateEmissions( double wavenumber, double *signal )
{
	std::vector<double>::iterator		iter;
	bool								ok;
	size_t								idx;

	iter = std::lower_bound( m_wavenum.begin(), m_wavenum.end(), wavenumber);
	ok = !(iter == m_wavenum.end()) && (*iter == wavenumber);
	if (ok)
	{
		idx     = iter - m_wavenum.begin();
		*signal = m_current_emission->at(idx);
	}
	else
	{
		*signal   = std::numeric_limits<double>::quiet_NaN();
		nxLog::Record(NXLOG_WARNING,"Hitran_Emission_Cache::CalculateEmissions. Error looking up the requested wavenumber from the cache. Thats indicates a serious problem.");
	}
	return ok;
}






