#include <skopticalproperties21.h>

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::SKTRAN_AtmosphericEmissionEntry		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericEmissionEntry::SKTRAN_AtmosphericEmissionEntry()
{
	m_emission      = NULL;
	m_species       = SKCLIMATOLOGY_UNDEFINED;
	m_radiance		= 0.0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::SKTRAN_AtmosphericEmissionEntry		2008-11-18*/
/** Special constructor reserved for quickly making blank entries that we can use to
 *	search lists
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericEmissionEntry::SKTRAN_AtmosphericEmissionEntry( const CLIMATOLOGY_HANDLE& species)
{
	m_emission      = NULL;
	m_species       = species;
	m_radiance      = 0.0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::SKTRAN_AtmosphericEmissionEntry		2014-3-5*/
/** The copy constructor is only meant to be used for copying blank objects although
 *	it should work with regular objects. Its main original purpose was to ensure that entries 
 *	copied onto the "list of entries" in the atmospheric optical state are properly copied.
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericEmissionEntry::SKTRAN_AtmosphericEmissionEntry ( const SKTRAN_AtmosphericEmissionEntry& other )
{
	m_emission = NULL;
	Configure( other.m_species, other.m_emission);
	m_radiance = other.m_radiance;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::Configure		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmissionEntry::Configure( CLIMATOLOGY_HANDLE species, skEmission* emission)
{
	if (  emission  != NULL) emission->AddRef();
	if ( m_emission != NULL) m_emission->Release();;
	m_radiance = 0;
	m_species  = species;
	m_emission = emission;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::~SKTRAN_AtmosphericEmissionEntry		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericEmissionEntry::~SKTRAN_AtmosphericEmissionEntry()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::ReleaseResources		2008-3-5*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_AtmosphericEmissionEntry::ReleaseResources()
{
	if (m_emission != NULL ) m_emission->Release();
	m_emission = NULL;
	m_species       = SKCLIMATOLOGY_UNDEFINED;
	m_radiance = 0.0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::CalculateCrossSections		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmissionEntry::CalculateEmission( double wavenumber, /*skClimatology* neutralatmosphere,*/ const GEODETIC_INSTANT& placeandtime, bool isground)
{
	bool	ok;

//	ok =       m_emission->SetAtmosphericState( neutralatmosphere);
	
	ok =       m_emission->UpdateLocation     ( placeandtime, isground );
	ok = ok && m_emission->IsotropicEmission  ( wavenumber, &m_radiance);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmissionEntry::CalculateEmission, Error updating emission for given place and time");
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::CalculateMultiWaveExtinctionsPerCM		 2014- 4- 24*/
/** Calculates the extinction, absorption and scattering per cm for an
 *	array of wavelengths for one species. The optical properties classes, especially HITRAN,
 *	permit considerable speed enhancement if we calculate the cross-sections in ascending
 *	wavenumber order.
 **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_AtmosphericEmissionEntry::CalculateMultiWaveEmission (	const std::vector<double>&	wavenumber, 
																	/*skClimatology*				neutralatmosphere,*/
																	const GEODETIC_INSTANT&		placeandtime,
																	bool						isground,
																	std::vector<double>*		emission)
{
	bool	ok;

	//ok =       m_emission->SetAtmosphericState      ( neutralatmosphere);
	ok =       m_emission->UpdateLocation			 ( placeandtime, isground);
	ok = ok && m_emission->IsotropicEmissionArray    ( wavenumber, emission);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmissionEntry::CalculateMultiWaveEmission, Error updating emission radiance for given place and time");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmissionEntry::UpdateNumberDensityPerCM3		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmissionEntry::UpdateInternalClimatologies( const GEODETIC_INSTANT& placeandtime, bool /*updatecache*/ )
{
	bool ok;

	ok = (m_emission == NULL) || (m_emission->UpdateCache(placeandtime));										// update the internal cache of the optical propes (used for climatologies of mode radius and mode width for example)
	return ok;
}


/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericEmission::SKTRAN_AtmosphericEmission		2003-12-5
 *-------------------------------------------------------------------------*/

SKTRAN_AtmosphericEmission::SKTRAN_AtmosphericEmission()
{
	m_wavenumber             = 0.0;
	m_placeandtime.heightm   = 0.0;
	m_placeandtime.latitude  = 0.0;
	m_placeandtime.longitude = 0.0;
	m_placeandtime.mjd       = -99999.0;
	m_isground               = false;

	m_isdirty                = true;
	m_updateclimatologycache = true;
	m_radiance               = 0.0;
//	m_atmosphericstate       = nullptr;
	m_solarspectrum          = new skSolarSpectrum_SAO2010;
	m_solarspectrum->AddRef();
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericEmission::ReleaseResources		2003-12-5
 *-------------------------------------------------------------------------*/

void SKTRAN_AtmosphericEmission::ReleaseResources()
{
	m_species.erase( m_species.begin(), m_species.end() );
	m_wavenumber             = 0.0;
	m_placeandtime.heightm   = -9999999.0;
	m_placeandtime.latitude  = -9999999.0;
	m_placeandtime.longitude = -9999999.0;
	m_placeandtime.mjd       = -9999999.0;
	m_isground               = false;

	m_isdirty                = true;					// Flags that the cross-sections are out of date
	m_updateclimatologycache = true;					// Flags that a climatology cache update is pending (time and location has changed)
	m_radiance               = 0.0;
	SetDirty();
}


/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericEmission::SKTRAN_AtmosphericEmission		2003-12-5
 *-------------------------------------------------------------------------*/

SKTRAN_AtmosphericEmission::~SKTRAN_AtmosphericEmission()
{
	ReleaseResources();
//	if (m_atmosphericstate != nullptr) m_atmosphericstate->Release();
	if (m_solarspectrum    != nullptr) m_solarspectrum->Release();

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::SetTimeAndLocation		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::SetTimeAndLocation(const GEODETIC_INSTANT& point, bool isgroundpoint, bool updateclimatologycache)
{
	m_placeandtime = point;
	m_isground     = isgroundpoint;
	SetDirty();
	if (m_solarspectrum != nullptr) m_solarspectrum->SetSolarDistanceFromMjd( point.mjd );
	if (updateclimatologycache) SetPendingCacheUpdate ();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::SetWavelength		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::SetWavelength( double wavelen_nm )
{
	m_wavenumber = 1.0E7/wavelen_nm;
	SetDirty();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::GetSpeciesClimatology		2008-8-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::GetSpeciesEmissionObject( const CLIMATOLOGY_HANDLE& speciesinlist, skEmission** emission)
{
	iterator	iter;
	bool		ok = false;
	SKTRAN_AtmosphericEmissionEntry*	entry;

	*emission = NULL;
	CheckClimatologyCacheIsValid(false);										// Update the climatology (and internal optical property) caches, but dont warn of invalid dates and times as it may not have been set yet
	for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										
	{																										
		entry = &(*iter);
		if ( entry->GetSpecies( ) == speciesinlist )
		{
			ok = true;
			*emission = entry->EmissionObject();
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::CheckDirtyAndUpdate		2008-3-3*/
/** Check the dirty flag to see if the cross-section data is out of date
 *	and needs to be updated.
 **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_AtmosphericEmission::CheckDirtyAndUpdate()
{
	iterator						iter;
	bool							ok;

	ok = (!m_isdirty );
	if (!ok)
	{
		ok = CalculateEmissions();
		m_isdirty = !ok;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::CheckClimatologyCacheIsValid		2013-6-11*/
/** Checks to see if there is a pending request to update the climatology
 *	caches. We need to update not only the explicit particle climatologies but also
 *	the internal climatologies inside the optical properties.
 *	The cache update can only happen if the user has previously set the Time and Location
 *	and we provide a switch to warn about this issue depending upon the callers context.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::CheckClimatologyCacheIsValid( bool warnaboutbadtime)
{	
	iterator										iter;
	SKTRAN_AtmosphericEmissionEntry*		entry;
	bool											ok;
	bool											ok1;
	bool											ok2;
	
	ok = !m_updateclimatologycache;																	// OK if there no climatology cache updates pending
	if (!ok)																						// if not ok then cache updates are pending
	{																								// so 
		ok = TimeAndPlaceIsValid();																	// make sure the time and place has been set to a valid value
		if (!ok)																					// If the time and place is bad
		{																							// then 
			ok = !(warnaboutbadtime);																// ok is good if user does not want warning
			if (!ok)																				// if he does want warning then ok is bad and we shall send the warning
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::CheckClimatologyCacheIsValid, Cannot update the caches as the internal Time and Location is invalid. call SetTimeAndLocation");
			}
		}
		else
		{
//			ok = m_atmosphericstate->UpdateCache(m_placeandtime);													// If a cache update is requested then update the atmospheric state
//			if (!ok)
//			{
//				nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::CalculateCrossSections, Error updating atmospheric state");
//			}
			ok1 = true;
			for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
			{																										// iterate over all of the species
				entry  = &(*iter);																					// get the entry
				ok2    = entry->EmissionObject()->UpdateCache( m_placeandtime );			// Update the internal climatology of any optical properties
				ok1    = ok1 && ok2;
			}
			if (!ok1)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::UpdateClimatologyCache, Error updating atmospheric species caches. Thats not good");
			}
			ok = ok && ok1;
			m_updateclimatologycache  = !ok;
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericEmission::CalculateCrossSections		2003-12-5
 *-------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::CalculateEmissions()
{
	iterator							iter;
	SKTRAN_AtmosphericEmissionEntry*	entry;
	double								rad;
	bool								ok = true;
	bool								ok1;
	double								toa_solarirradiance;


	m_radiance	= 0;
	ok = (m_species.size() == 0);
	if (!ok)
	{
		ok = CheckClimatologyCacheIsValid(true);															// Check and execute any pending climatology cache updates, warn about bad time settings
		toa_solarirradiance = m_solarspectrum->Irradiance(1.0E7/m_wavenumber );
		ok = ok && (toa_solarirradiance > 0.0) && NXFINITE(toa_solarirradiance);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::CalculateEmissions, There were errors either checking the climatology cache or the toa solar irradiance for wavenumber %15.6f cm-1 (%15.6f nm)", (double)m_wavenumber, (double)(1.0E7/m_wavenumber));
		}
		else
		{
			for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
			{																										// iterate over all of the species
				entry = &(*iter);																					// get the entry
				ok1 = entry->CalculateEmission(  m_wavenumber, /*m_atmosphericstate,*/ m_placeandtime, m_isground );		// so calculate cross-sections
				if (ok1)
				{
					rad		    = entry->IsotropicRadiance();
					m_radiance  += rad;
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::CalculateEmissions, Error calculating emission for individual entry");
				}
				ok = ok && ok1;
			}
			m_radiance /= toa_solarirradiance;
		}
		m_updateclimatologycache = m_updateclimatologycache && !ok;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::CalculateMultiWaveCrossSections		 2014- 4- 24*/
/** Calculates the extinction, absoprtion and scattering per cm at the
 *	current location
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::CalculateMultiWaveEmissions(const std::vector<double>& wavenumber, std::vector<double>* radiance)
{
	iterator							iter;
	SKTRAN_AtmosphericEmissionEntry*	entry;
	bool								ok = true;
	bool								ok1;
	std::vector<double>					rad;
	double								toa_solarirradiance;
	
	radiance->assign ( wavenumber.size(), 0.0);	
	ok = (m_species.size() == 0);
	if (!ok)
	{
		ok = CheckClimatologyCacheIsValid(true);															// Check and execute any pending climatology cache updates, warn about bad time settings
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::CalculateMultiWaveEmissions, There were errors checking the climatology cache");
		}
		else
		{
			for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
			{																										// iterate over all of the species
				entry = &(*iter);																					// get the entry
				ok1 = entry->CalculateMultiWaveEmission( wavenumber, /*m_atmosphericstate,*/ m_placeandtime, m_isground, &rad);		// so calculate cross-sections
				if (ok1)
				{
					for (size_t iw = 0; iw < wavenumber.size(); iw++)
					{
						radiance->at(iw)  += rad.at(iw);
					}
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericEmission::CalculateMultiWaveEmissions, Error calculating emissions");
				}
				ok = ok && ok1;
			}

			for (size_t iw = 0; iw < wavenumber.size(); iw++)
			{
				toa_solarirradiance = m_solarspectrum->Irradiance(1.0E7/wavenumber.at(iw) );
				radiance->at(iw)  /= toa_solarirradiance;
			}

		}
		m_updateclimatologycache = m_updateclimatologycache && !ok;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::AddSpecies		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::AddEmission( const CLIMATOLOGY_HANDLE& species, skEmission* emission)
{
	iterator							iter;
	SKTRAN_AtmosphericEmissionEntry		dummy(species);
	SKTRAN_AtmosphericEmissionEntry*	entry;
	bool								ok;

	iter = std::find( m_species.begin(), m_species.end(), dummy );						// See if this entry already exists
	if (iter == m_species.end())														// If it does not
	{																					// then
		m_species.push_back( dummy );													// create a new entry on the back
		entry = &(m_species.back());													// and get the pointer to this entry
	}																					// otherwise
	else																				// this entry does exist
	{																					// so
		entry = &(*iter);																// get its address
	}																					// done that
	ok = entry->Configure( species, emission );											// now configure the entry
	SetDirty();																			// set the list as dirty
	SetPendingCacheUpdate();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_AtmosphericEmission::AddSpecies, Error adding the requested species  to the list");
	}
	return ok;																			// and we are done
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::RemoveSpecies		2013-6-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::RemoveEmission( const CLIMATOLOGY_HANDLE& species)
{
	iterator							iter;
	SKTRAN_AtmosphericEmissionEntry	dummy(species);

	iter = std::find( m_species.begin(), m_species.end(), dummy );						// See if this entry already exists
	if (!(iter == m_species.end()))														// If it does not
	{																					// then
		m_species.erase( iter );														// erase this entry 
		SetDirty();																		// and flag that the cross-sections need updating
	}
	return true;																			// and we are done
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::SetAtmosphericStateModel		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_AtmosphericEmission::SetAtmosphericStateModel( skClimatology* atmosphericstate )
{
	if (atmosphericstate   != NULL) atmosphericstate->AddRef();
	if (m_atmosphericstate != NULL) m_atmosphericstate->Release();
	m_atmosphericstate = atmosphericstate;
	SetDirty();
	SetPendingCacheUpdate();

	return (atmosphericstate != NULL);
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::GetAtmosphericStateModel		2009-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_AtmosphericEmission::GetAtmosphericStateModel( skClimatology** statemodel )
{
	bool ok = false;
	if ( m_atmosphericstate != NULL ) 
	{
		CheckClimatologyCacheIsValid(false);										// Update cache if one is pending
		*statemodel = m_atmosphericstate;
		ok = true;
	}  
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericEmission::SetSolarSpectrum		 2015- 3- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericEmission::SetSolarSpectrum( skSolarSpectrum* solarspectrum)
{
	if (solarspectrum   != nullptr) solarspectrum->AddRef();
	if (m_solarspectrum != nullptr) m_solarspectrum->Release();

	m_solarspectrum = solarspectrum;
	return true;
}

