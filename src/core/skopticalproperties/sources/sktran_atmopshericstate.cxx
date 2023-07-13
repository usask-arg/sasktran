#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::SKTRAN_AtmosphericOpticalStateEntry_V21		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericOpticalStateEntry_V21::SKTRAN_AtmosphericOpticalStateEntry_V21()
{
	m_climatology   = NULL;
	m_particleprops = NULL;
	m_species       = SKCLIMATOLOGY_UNDEFINED;
	m_numberdensity = 0.0;
	m_absxs         = 0.0;
	m_extxs         = 0.0;
	m_scattxs       = 0.0;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::SKTRAN_AtmosphericOpticalStateEntry_V21		2008-11-18*/
/** Special constructor reserved for quickly making blank entries that we can use to
 *	search lists
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericOpticalStateEntry_V21::SKTRAN_AtmosphericOpticalStateEntry_V21( const CLIMATOLOGY_HANDLE& species)
{
	m_climatology   = NULL;
	m_particleprops = NULL;
	m_species       = species;
	m_numberdensity = 0.0;
	m_absxs         = 0.0;
	m_extxs         = 0.0;
	m_scattxs       = 0.0;

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::SKTRAN_AtmosphericOpticalStateEntry_V21		2014-3-5*/
/** The copy constructor is only meant to be used for copying blank objects although
 *	it should work with regular objects. Its main original purpose was to ensure that entries
 *	copied onto the "list of entries" in the atmospheric optical state are properly copied.
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericOpticalStateEntry_V21::SKTRAN_AtmosphericOpticalStateEntry_V21 ( const SKTRAN_AtmosphericOpticalStateEntry_V21& other )
{
	m_climatology   = NULL;
	m_particleprops = NULL;
	Configure( other.m_species, other.m_climatology, other.m_particleprops);
	m_numberdensity = other.m_numberdensity;
	m_absxs         = other.m_absxs;
	m_extxs         = other.m_extxs;
	m_scattxs       = other.m_scattxs;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::Configure		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalStateEntry_V21::Configure( CLIMATOLOGY_HANDLE species, skClimatology* numberdensityclimatology, skOpticalProperties* particleopticalprops)
{
	if ( numberdensityclimatology != NULL) numberdensityclimatology->AddRef();
	if ( particleopticalprops     != NULL) particleopticalprops->AddRef();
    if ( m_climatology            != NULL) m_climatology->Release();
	if ( m_particleprops          != NULL) m_particleprops->Release();;
	m_numberdensity = 0;
	m_species       = species;
	m_climatology   = numberdensityclimatology;
	m_particleprops = particleopticalprops;
	m_absxs         = 0.0;
	m_extxs         = 0.0;
	m_scattxs       = 0.0;
	return  !((m_climatology != NULL) ^ (m_particleprops!= NULL));	// return true if (both are NULL) or (both are Not NULL), thats an XOR
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::~SKTRAN_AtmosphericOpticalStateEntry_V21		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_AtmosphericOpticalStateEntry_V21::~SKTRAN_AtmosphericOpticalStateEntry_V21()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::ReleaseResources		2008-3-5*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_AtmosphericOpticalStateEntry_V21::ReleaseResources()
{
	if (m_climatology   != NULL ) m_climatology->Release();
	if (m_particleprops != NULL ) m_particleprops->Release();
	m_climatology   = NULL;
	m_particleprops = NULL;
	m_species       = SKCLIMATOLOGY_UNDEFINED;
	m_numberdensity = 0.0;
	m_absxs         = 0.0;
	m_extxs         = 0.0;
	m_scattxs       = 0.0;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateCrossSections		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateCrossSections( double wavenumber, skClimatology* neutralatmosphere, const GEODETIC_INSTANT& placeandtime)
{
	bool	ok;
	bool	haschanged;

	NXTRACE_ONCEONLY( firsttime, ("SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateCrossSections, needs re-working so the climatology update is in one loop and the cross-sections is in another\n"));

	ok =       m_particleprops->SetAtmosphericState( neutralatmosphere);
	ok =       m_particleprops->SetLocation       ( placeandtime, &haschanged );
	ok = ok && m_particleprops->CalculateCrossSections( wavenumber, &m_absxs, &m_extxs, &m_scattxs);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateCrossSections, Error updating cross-section for given place and time");
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateMultiWaveExtinctionsPerCM		 2014- 4- 24*/
/** Calculates the extinction, absorption and scattering per cm for an
 *	array of wavelengths for one species. The optical properties classes, especially HITRAN,
 *	permit considerable speed enhancement if we calculate the cross-sections in ascending
 *	wavenumber order.
 **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateMultiWaveExtinctionsPerCM (	const std::vector<double>&	wavenumber,
																					skClimatology*				neutralatmosphere,
																					const GEODETIC_INSTANT&		placeandtime,
																					std::vector<double>*		absxs,
																					std::vector<double>*		extxs,
																					std::vector<double>*		scattxs )
{
	bool	ok;
	bool	haschanged;
	double	numberdensity;

	NXTRACE_ONCEONLY( firsttime, ("SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateMultiWaveCrossSections, needs re-working so the climatology update is in one loop and the cross-sections is in another\n"));

	absxs->resize( wavenumber.size() );
	extxs->resize( wavenumber.size() );
	scattxs->resize( wavenumber.size() );

	ok =       m_particleprops->SetAtmosphericState( neutralatmosphere);
	ok =       m_particleprops->SetLocation        ( placeandtime, &haschanged );
	ok = ok && UpdateNumberDensityPerCM3( placeandtime, false);
	ok = ok && m_particleprops->CalculateCrossSectionsArray( &wavenumber.front(), (int)wavenumber.size(), &absxs->front(), &extxs->front(), &scattxs->front() );
	numberdensity = ok ? m_numberdensity : 0.0;
	if (ok)
	{
		for (size_t i = 0; i < absxs->size(); i++)
		{
			absxs->at(i)   *= numberdensity;
			extxs->at(i)   *= numberdensity;
			scattxs->at(i) *= numberdensity;
			NXASSERT(( extxs->at(i) < 1000));
		}
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalStateEntry_V21::CalculateMultiWaveCrossSections, Error updating cross-section for given place and time");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::UpdateNumberDensityPerCM3		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalStateEntry_V21::UpdateNumberDensityPerCM3( const GEODETIC_INSTANT& placeandtime, bool updatecache )
{
	bool ok;

	ok = m_climatology->GetParameter( m_species, placeandtime, &m_numberdensity, updatecache );			// Get the effective number density for optical props
	if (updatecache)																					// and if we nede to update the cache
	{																									// then
		ok = ok && m_particleprops->InternalClimatology_UpdateCache(placeandtime );										// update the internal cache of the optical propes (used for climatologies of mode radius and mode width for example)
	}
	ok = ok && NXFINITE(m_numberdensity);
																									// and that is that
	if (!ok) m_numberdensity =0.0;																		// if there was an error then clear the numebr density.
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::UpdateClimatology           		2008-8-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalStateEntry_V21::UpdateClimatology( skClimatology* numberdensityclimatology )
{

	if (numberdensityclimatology != NULL ) numberdensityclimatology->AddRef();
	if ( m_climatology           != NULL ) m_climatology->Release();
	m_numberdensity = 0.0;
	m_climatology   = numberdensityclimatology;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalStateEntry_V21::DeepCopy		2008-3-5*/
/** Creates a complet, new copy of the other entry.  This is used to make
 *	separate copies of the species list for multi-threaded execution.
 **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_AtmosphericOpticalStateEntry_V21::DeepCopy( const SKTRAN_AtmosphericOpticalStateEntry_V21& other )
{
	bool	ok1;
	bool	ok2;
	bool	ok;

	ReleaseResources();
	ok1 = (other.m_climatology == NULL );
	if (!ok1 ) ok1  = other.m_climatology->CreateClone(&m_climatology);

	ok2 = (other.m_particleprops == NULL);
	if (!ok2 ) ok2  = other.m_particleprops->CreateClone(&m_particleprops);

	m_species       = other.m_species;
	m_numberdensity = other.m_numberdensity;
	ok = m_climatology->IsSupportedSpecies(m_species);
	NXASSERT(( ok ));
	ok = ok && (ok1 && ok2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalStateEntry_V21::DeepCopy, Error copying object  into this instance");
		ReleaseResources();
	}
	return ok;
}
*/


/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericOpticalState_V21::SKTRAN_AtmosphericOpticalState_V21		2003-12-5
 *-------------------------------------------------------------------------*/

SKTRAN_AtmosphericOpticalState_V21::SKTRAN_AtmosphericOpticalState_V21()
{
	m_wavenumber             = 0.0;
	m_placeandtime.heightm   = 0.0;
	m_placeandtime.latitude  = 0.0;
	m_placeandtime.longitude = 0.0;
	m_placeandtime.mjd       = -99999.0;

	m_isdirty                = true;
	m_updateclimatologycache = true;
	m_kabs                   = 0.0;
	m_kext                   = 0.0;
	m_kscat                  = 0.0;
	m_kdelta                 = 0.0;

	m_constantalbedo = new SKTRAN_BRDF_Lambertian(0.0);
	if (m_constantalbedo == NULL) nxLog::Record(NXLOG_WARNING, "SKTRAN_AtmosphericOpticalState_V21::Constructor, error allocaing memory for albedo constant object, thats not good");
	m_constantalbedo->AddRef();
	m_albedo = NULL;
	SetAlbedo( 0.0 );

	m_defaultatmosphericstate= new skClimatology_MSIS90;
	if (m_defaultatmosphericstate == NULL) nxLog::Record(NXLOG_WARNING, "SKTRAN_AtmosphericOpticalState_V21::Constructor, error allocaing memory for default atmosphere  object, thats not good");
	m_defaultatmosphericstate->AddRef();

	m_atmosphericstate = m_defaultatmosphericstate;
	m_atmosphericstate->AddRef();
	//m_atmosphericemission.SetAtmosphericStateModel( m_atmosphericstate );

}

/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericOpticalState_V21::ReleaseResources		2003-12-5
 *-------------------------------------------------------------------------*/

void SKTRAN_AtmosphericOpticalState_V21::ReleaseResources()
{
	m_species.erase( m_species.begin(), m_species.end() );
	m_wavenumber             = 0.0;
	m_placeandtime.heightm   = -9999999.0;
	m_placeandtime.latitude  = -9999999.0;
	m_placeandtime.longitude = -9999999.0;
	m_placeandtime.mjd       = -9999999.0;

	m_isdirty                = true;					// Flags that the cross-sections are out of date
	m_updateclimatologycache = true;					// Flags that a climatology cache update is pending (time and location has changed)
	m_kabs                   = 0.0;
	m_kext                   = 0.0;
	m_kscat                  = 0.0;
	m_kdelta                 = 0.0;

	SetDirty();
	m_atmosphericemission.ReleaseResources();
}


/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericOpticalState_V21::SKTRAN_AtmosphericOpticalState_V21		2003-12-5
 *-------------------------------------------------------------------------*/

SKTRAN_AtmosphericOpticalState_V21::~SKTRAN_AtmosphericOpticalState_V21()
{
	ReleaseResources();
	if (m_atmosphericstate != NULL ) m_atmosphericstate->Release();
	if (m_albedo           != NULL ) m_albedo->Release();
	if (m_constantalbedo   != NULL ) m_constantalbedo->Release();
	if (m_defaultatmosphericstate != NULL) m_defaultatmosphericstate->Release();
//	--m_numinstances;

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::SetAlbedoObject		2010-6-28*/
/** Sets the albedo object that will be used by the Sasktran engine. The
 *	default is to use the internal constant albedo obect. This method allows
 *	the end user supply an albedo that varies with geographical position
 *	and wavelength.
 *
 *	The object passed in will normally be allocated on the heap as this class
 *	will place a reference count on the object.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::SetAlbedoObject(skBRDF* albedo)
{
	if (albedo == NULL) albedo = m_constantalbedo;
	albedo->AddRef();
	if (m_albedo != NULL) m_albedo->Release();
	m_albedo = albedo;
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::GetAlbedoObject		2010-6-28*/
/** Returnss the current albedo object being used. No reference count is added **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::GetAlbedoObject(skBRDF** albedo)
{
	NXASSERT(( m_albedo != NULL ));
	*albedo = m_albedo;
	//m_albedo->AddRef();
	return (*albedo != NULL);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::SetAlbedo		2010-6-28*/
/** A legacy method that allow the end user to set the albedo as a constant
 *	value. This is probably the most common way of setting the albedo in the
 *	radiative transfer engine.  More sophisticated albedo models that vary with
 *	position, time and wavelength can be implemented via method SetAlbedoObject.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::SetAlbedo(double albedo)
{
	bool	ok;

	ok =       m_constantalbedo->SetAlbedo ( albedo );
	ok = ok && SetAlbedoObject( m_constantalbedo );
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_AtmosphericOpticalState_V21::SetAlbedo, There was an error setting the albedo within the atmospheric optical state");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::SetTimeAndLocation		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::SetTimeAndLocation(const GEODETIC_INSTANT& point, bool updateclimatologycache)
{
	m_placeandtime = point;
	SetDirty();
	if (updateclimatologycache) SetPendingCacheUpdate ();
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::SetWavelength		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::SetWavelength( double wavelen_nm )
{
	m_wavenumber = 1.0E7/wavelen_nm;
	SetDirty();
	m_atmosphericemission.SetWavelength(wavelen_nm);
	return true;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::UpdateSpeciesClimatology		2008-3-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::UpdateSpeciesClimatology( const CLIMATOLOGY_HANDLE& speciesinlist, skClimatology* numberdensityclimatology )
{
	iterator	iter;
	bool		ok = false;
	SKTRAN_AtmosphericOpticalStateEntry_V21		dummy( speciesinlist );
	SKTRAN_AtmosphericOpticalStateEntry_V21*	entry;

	iter = std::find( m_species.begin(), m_species.end(), dummy );
	ok = !(iter == m_species.end());
	if (ok)
	{
		entry = &(*iter);
		ok = numberdensityclimatology->IsSupportedSpecies(speciesinlist);
		NXASSERT((ok));
		ok = ok & entry->UpdateClimatology( numberdensityclimatology );
		SetDirty			  ();												// Flag that cross-sections are dirty
		SetPendingCacheUpdate ();												// Flag that climatologies will need re-caching
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::UpdateSpeciesClimatology, The requested species does not exist in the list or the new species is not supported by the new climatology. Cannot change its climatology");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::GetSpeciesClimatology		2008-8-20*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_AtmosphericOpticalState_V21::GetSpeciesClimatology( const CLIMATOLOGY_HANDLE& speciesinlist, skClimatology** numberdensityclimatology )
{
	iterator	iter;
	bool		ok = false;
	SKTRAN_AtmosphericOpticalStateEntry_V21*	entry;

	CheckClimatologyCacheIsValid(false);										// Update the climatology cache, but downt warn of invalid dates and times as it may not have been set yet
	for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)
	{
		entry = &(*iter);
		if ( entry->GetSpecies( ) == speciesinlist )
		{
			ok = true;
			*numberdensityclimatology = entry->GetClimatology( );
		}
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::GetSpeciesClimatology		2008-8-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::GetSpeciesOpticalProperties( const CLIMATOLOGY_HANDLE& speciesinlist, skOpticalProperties** opticalprops)
{
	iterator	iter;
	bool		ok = false;
	SKTRAN_AtmosphericOpticalStateEntry_V21*	entry;

	*opticalprops = NULL;
	CheckClimatologyCacheIsValid(false);										// Update the climatology (and internal optical property) caches, but dont warn of invalid dates and times as it may not have been set yet
	for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)
	{
		entry = &(*iter);
		if ( entry->GetSpecies( ) == speciesinlist )
		{
			ok = true;
			*opticalprops = entry->ParticleOpticalProps();
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::CheckDirtyAndUpdate		2008-3-3*/
/** Check the dirty flag to see if the cross-section data is out of date
 *	and needs to be updated.
 **/
/*---------------------------------------------------------------------------*/

bool  SKTRAN_AtmosphericOpticalState_V21::CheckDirtyAndUpdate()
{
	iterator						iter;
	bool							ok;

	ok = (!m_isdirty );
	if (!ok)
	{
		ok = CalculateCrossSections();
		m_isdirty = !ok;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::CheckClimatologyCacheIsValid		2013-6-11*/
/** Checks to see if there is a pending request to update the climatology
 *	caches. We need to update not only the explicit particle climatologies but also
 *	the internal climatologies inside the optical properties.
 *	The cache update can only happen if the user has previously set the Time and Location
 *	and we provide a switch to warn about this issue depending upon the callers context.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::CheckClimatologyCacheIsValid( bool warnaboutbadtime)
{
	iterator										iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
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
				nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::CheckClimatologyCacheIsValid, Cannot update the caches as the internal Time and Location is invalid. call SetTimeAndLocation");
			}
		}
		else
		{
			ok = m_atmosphericstate->UpdateCache(m_placeandtime);													// If a cache update is requested then update the atmospheric state
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::CalculateCrossSections, Error updating atmospheric state");
			}
			ok1 = true;
			for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
			{																										// iterate over all of the species
				entry  = &(*iter);																					// get the entry
				ok2    =        entry->GetClimatology()->UpdateCache(m_placeandtime);								// Update the species "number density" climatology
				ok2    = ok2 && entry->ParticleOpticalProps()->InternalClimatology_UpdateCache( m_placeandtime );			// Update the internal climatology of any optical properties
				ok1    = ok1 && ok2;
			}
			if (!ok1)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::UpdateClimatologyCache, Error updating atmospheric species caches. Thats not good");
			}
			ok = ok && ok1;
			m_updateclimatologycache  = !ok;
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericOpticalState_V21::CalculateCrossSections		2003-12-5
 *-------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::CalculateCrossSections()
{
	iterator						iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
	double							n;
	double							kabs;
	double							kext;
	double							ksca;
	bool							ok = true;
	bool							ok1;
	double							f;
	double							omega;

	m_kabs   = 0.0;
	m_kext   = 0.0;
	m_kscat  = 0.0;
	m_kdelta = 0.0;

	ok = CheckClimatologyCacheIsValid(true);															// Check and execute any pending climatology cache updates, warn about bad time settings
	if (ok)
	{
		for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
		{																										// iterate over all of the species
			entry = &(*iter);																					// get the entry
			ok1    = entry->UpdateNumberDensityPerCM3( m_placeandtime, m_updateclimatologycache );				// Update the number desnity for the climatology
			n      = entry->CurrentNumberDensityPerCM3();														// only update climatologies
			if (ok1 && (n > 0) && (m_wavenumber != 0.0) )														// if the wavenumber is zero
			{																									// otherwise wavenumber not zero
				ok1 = entry->CalculateCrossSections (  m_wavenumber, m_atmosphericstate, m_placeandtime );		// so calculate cross-sections
				if (ok1)
				{
					kabs		= entry->AbsorptionCrossSection();
					kext		= entry->ExtinctionCrossSection();
					ksca		= entry->ScatteringCrossSection();
					f			= entry->ParticleOpticalProps()->DeltaFunctionForwardScatterFraction();
					NXASSERT(( (f >= 0) && (f <= 1.0) ));
					omega		= kext > 0 ? ksca/kext : 0.0;
					m_kabs		+= kabs*n;
					m_kext		+= (1-omega*f)*kext*n;
					m_kscat		+= (1-f)*ksca*n;
					m_kdelta    += f*ksca*n;
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::CalculateCrossSections, Error calculating cross-sections");
				}
			}
			ok = ok && ok1;
		}
	}
	m_updateclimatologycache = m_updateclimatologycache && !ok;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::CalculateMultiWaveCrossSections		 2014- 4- 24*/
/** Calculates the extinction, absoprtion and scattering per cm at the
 *	current location
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::CalculateMultiWaveCrossSections(const std::vector<double>& wavenumber, std::vector<double>* kabs, std::vector<double>* kext, std::vector<double>* kscat )
{
	iterator						iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
	double							n;
	bool							ok = true;
	bool							ok1;
//	double							f;
//	double							omega;
	std::vector<double>				absxs;
	std::vector<double>				extxs;
	std::vector<double>				scattxs;


	kabs->assign ( wavenumber.size(), 0.0);
	kext->assign ( wavenumber.size(), 0.0);
	kscat->assign( wavenumber.size(), 0.0);

	ok = CheckClimatologyCacheIsValid(true);															// Check and execute any pending climatology cache updates, warn about bad time settings
	if (ok)
	{
		for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
		{																										// iterate over all of the species

			entry = &(*iter);																					// get the entry
			ok1    = entry->UpdateNumberDensityPerCM3( m_placeandtime, m_updateclimatologycache );				// Update the number desnity for the climatology
			n      = entry->CurrentNumberDensityPerCM3();														// only update climatologies
			if (ok1 && (n > 0) )																				// if twe have a nobn zero number density
			{																									// otherwise wavenumber not zero
				if (entry->ParticleOpticalProps()->IsDeltaFunctionForwardScatter())
				{
					// Loop over the wavelengths classically since delta function forward scatter is not supported by the array methods
					for (int iw = 0; iw < wavenumber.size(); iw++)
					{
						entry->CalculateCrossSections(wavenumber[iw], m_atmosphericstate, m_placeandtime);

						double f = entry->ParticleOpticalProps()->DeltaFunctionForwardScatterFraction();
						double kabstemp;
						double kexttemp;
						double kscattemp;

						kabstemp = entry->AbsorptionCrossSection();
						kexttemp = entry->ExtinctionCrossSection();
						kscattemp = entry->ScatteringCrossSection();

						double omega = kexttemp > 0 ? kscattemp / kexttemp : 0.0;

						kabs->at(iw) += kabstemp;
						kext->at(iw) += (1 - omega * f)*kexttemp;
						kscat->at(iw) += (1 - f)*kscattemp;
					}
				}
				else
				{
					ok1 = entry->CalculateMultiWaveExtinctionsPerCM( wavenumber, m_atmosphericstate, m_placeandtime, &absxs, &extxs, &scattxs );		// so calculate cross-sections
					if (ok1)
					{
						for (size_t iw = 0; iw < wavenumber.size(); iw++)
						{
							kabs->at(iw) += absxs.at(iw);
							kext->at(iw) += extxs.at(iw);
							kscat->at(iw) += scattxs.at(iw);
						}
					}
				}
			}
			else
			{
				if (!ok1)
				{
					nxLog::Record(NXLOG_WARNING, "SKTRAN_AtmosphericOpticalState_V21::CalculateMultiWaveCrossSections, Error fetching number density at location (lat =%8.4f, lng =%9.4f, height=%12.3f)", (double)m_placeandtime.latitude, (double)m_placeandtime.longitude, (double)m_placeandtime.heightm);
				}
			}
			ok = ok && ok1;
		}
	}
	m_updateclimatologycache = m_updateclimatologycache && !ok;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::CheckCosineRange		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_AtmosphericOpticalState_V21::CheckCosineRange(double * cosscatteringangle)
{
	if (*cosscatteringangle < -1.0)
	{
		NXASSERT( *cosscatteringangle > -1.0001 );
		*cosscatteringangle = -1.0;
	}
	if (*cosscatteringangle < -1.0)
	{
		NXASSERT( *cosscatteringangle < 1.0001 );
		*cosscatteringangle = 1.0;
	}
}

bool SKTRAN_AtmosphericOpticalState_V21::PhaseGridHint(const std::vector<double>& cosscatterangles)
{
	iterator						iter;
	for (iter = m_species.begin(); !(iter == m_species.end()); ++iter)
	{
		iter->ParticleOpticalProps()->PhaseGridHint(cosscatterangles);
	}
	return true;
}


bool SKTRAN_AtmosphericOpticalState_V21::CalculateMultiWaveCrossSectionsAndPhaseMatrix(const std::vector<double>& wavenumber,
																					   std::vector<double>* kabs,
																					   std::vector<double>* kext,
																					   std::vector<double>* kscat,
																					   const std::vector<double>& cosangles,
																					   nx2dArray<skRTPhaseMatrix>* P )
{
	iterator						iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
	double							n;
	bool							ok = true;
	bool							ok1;
	double							f;
	double							omega;
	std::vector<double>				absxs;
	std::vector<double>				extxs;
	std::vector<double>				scattxs ;

	skRTPhaseMatrix					dummy;

	double norm;

	kabs->assign ( wavenumber.size(), 0.0);
	kext->assign ( wavenumber.size(), 0.0);
	kscat->assign( wavenumber.size(), 0.0);

	ok = CheckClimatologyCacheIsValid(true);															// Check and execute any pending climatology cache updates, warn about bad time settings
	ok = CheckDirtyAndUpdate();
	for (size_t iw = 0; iw < wavenumber.size(); iw++ )
	{
		for( size_t angleidx = 0; angleidx < cosangles.size(); angleidx++ )
		{
			P->At( iw, angleidx ).SetTo(0.0);
		}
	}

	if (ok)
	{
		for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)										// now update all o fthe species
		{																										// iterate over all of the species

			entry = &(*iter);																					// get the entry
			ok1    = entry->UpdateNumberDensityPerCM3( m_placeandtime, m_updateclimatologycache );				// Update the number desnity for the climatology
			n      = entry->CurrentNumberDensityPerCM3();														// only update climatologies
			if (ok1 && (n > 0) )																				// if twe have a nobn zero number density
			{																									// otherwise wavenumber not zero
				ok1 = entry->CalculateMultiWaveExtinctionsPerCM( wavenumber, m_atmosphericstate, m_placeandtime, &absxs, &extxs, &scattxs );		// so calculate cross-sections
				if (ok1)
				{
					for (size_t iw = 0; iw < wavenumber.size(); iw++)
					{

						double ksca;
						if( entry->ParticleOpticalProps()->IsDeltaFunctionForwardScatter() )
						{
							entry->CalculateCrossSections(wavenumber[iw], m_atmosphericstate, m_placeandtime);
							f			= entry->ParticleOpticalProps()->DeltaFunctionForwardScatterFraction();
							NXASSERT(( (f >= 0) && (f <= 1.0) ));
							omega		= kext->at(iw) > 0 ? kscat->at(iw)/kext->at(iw) : 0.0;
							kabs->at(iw)		+= kabs->at(iw)*n;
							kext->at(iw)		+= (1-omega*f)*kext->at(iw)*n;
							kscat->at(iw)		+= (1-f)*kscat->at(iw)*n;

							norm = nxmath::Pi * 4.0 * (1-f);

							ksca = scattxs[iw] * (1-omega*f);
						}
						else
						{
							NXASSERT(( extxs.at(iw) < 1000));
							kabs->at(iw)  += absxs.at(iw);
							kext->at(iw)  += extxs.at(iw);
							kscat->at(iw) += scattxs.at(iw);
							norm = nxmath::Pi * 4.0;
							ksca = scattxs[iw];
						}

						for( size_t angleidx = 0; angleidx < cosangles.size(); angleidx++ )
						{
							if( scattxs.at(iw) > 0 )
							{
								entry->ParticleOpticalProps()->CalculatePhaseMatrix( wavenumber[iw], cosangles[angleidx], &dummy );
								double factor = ksca / norm;
								dummy *= factor;
								P->At( iw, angleidx ) += dummy;
							}
						}
					}




				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::CalculateMultiWaveCrossSections, Error calculating cross-sections");
				}
			}
			ok = ok && ok1;
		}
	}
	m_updateclimatologycache = m_updateclimatologycache && !ok;
	return ok;;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::PhaseMatrix		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::VectorPhaseMatrix( double cosscatteringangle, skRTPhaseMatrix* P)
{
	iterator						iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
	double							n,norm;
	double							ksca,kext,f,omega,ksca0;
	double							ntotal = 0;
	SKRTFLOAT						factor;
	skRTPhaseMatrix					dummy;
	bool							ok;


	ok = CheckDirtyAndUpdate();
	CheckCosineRange(&cosscatteringangle);
	P->SetTo(0.0);

	if (ok)
	{
		for ( iter = m_species.begin(); !(iter == m_species.end()); ++iter)
		{
			entry = &(*iter);
			n = entry->CurrentNumberDensityPerCM3();
			if (n > 0)
			{
				if (entry->ParticleOpticalProps()->IsDeltaFunctionForwardScatter())
				{
					ksca0 = entry->ScatteringCrossSection();
					kext = entry->ExtinctionCrossSection();
					f = entry->ParticleOpticalProps()->DeltaFunctionForwardScatterFraction();
					omega = (kext > 0) ? ksca0 / kext : 1.0;
					ksca = (1.0 - omega * f)*ksca0;					// 'effective' scat cross-section
					norm = 4.0*nxmath::Pi*(1.0 - f);				// normalization for 'smoothed' phase funct
				}
				else
				{
					ksca = entry->ScatteringCrossSection();
					norm = 4.0*nxmath::Pi;
				}
				if (ksca > 0)
				{
					factor = (SKRTFLOAT)(n*ksca / norm);
					if (factor > 0.0)
					{
						ntotal += factor;
						entry->ParticleOpticalProps()->CalculatePhaseMatrix(m_wavenumber, cosscatteringangle, &dummy);
						dummy *= factor;
						*P += dummy;
					}
				}
			}
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::ScalarPhaseMatrix		2021-08-25*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::ScalarPhaseMatrix(std::pair<double, size_t> cosscatterandindex, double& p11)
{
	iterator						iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
	double							n, norm;
	double							ksca, kext, f, omega, ksca0;
	double							ntotal = 0;
	SKRTFLOAT						factor;
	double							dummy;
	bool							ok;

	double cosscatteringangle = cosscatterandindex.first;

	ok = CheckDirtyAndUpdate();
	CheckCosineRange(&cosscatteringangle);
	p11 = 0.0;

	if (ok)
	{
		for (iter = m_species.begin(); !(iter == m_species.end()); ++iter)
		{
			entry = &(*iter);
			n = entry->CurrentNumberDensityPerCM3();
			if (n > 0)
			{
				if (entry->ParticleOpticalProps()->IsDeltaFunctionForwardScatter())
				{
					ksca0 = entry->ScatteringCrossSection();
					kext = entry->ExtinctionCrossSection();
					f = entry->ParticleOpticalProps()->DeltaFunctionForwardScatterFraction();
					omega = (kext > 0) ? ksca0 / kext : 1.0;
					ksca = (1.0 - omega * f)*ksca0;					// 'effective' scat cross-section
					norm = 4.0*nxmath::Pi*(1.0 - f);				// normalization for 'smoothed' phase funct
				}
				else
				{
					ksca = entry->ScatteringCrossSection();
					norm = 4.0*nxmath::Pi;
				}
				if (ksca > 0)
				{
					factor = (SKRTFLOAT)(n*ksca / norm);
					if (factor > 0.0)
					{
						ntotal += factor;
						entry->ParticleOpticalProps()->CalculateP11(m_wavenumber, cosscatterandindex, dummy);
						dummy *= factor;
						p11 += dummy;
					}
				}
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::ScalarScatteringCoefficient		2008-3-3*/
/** Calculates the scattering coefficient from one angle to another.
 *	This is given by $f$coeff = \frac{P(1,1)\sigma_{scat}}{4\pi}$f$
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::ScalarScatteringCoefficient( double cosscatteringangle, double* scatteringvalue)
{
	skRTPhaseMatrix		P;
	bool				ok;

	ok               = VectorPhaseMatrix( cosscatteringangle, &P);
	*scatteringvalue = P.At(1,1);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::AddEmission		 2015- 3- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::AddEmission( const CLIMATOLOGY_HANDLE&  species, skEmission* emissionobject)
{
	return m_atmosphericemission.AddEmission( species, emissionobject);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::AddSpecies		2008-3-3*/
/** Adds a species to the list of species used to describe the optical
 *	properties of the atmosphere. This method can be used to either add a new species,
 *	replace an existing species or update just the climatology  of an existing species.
 *
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::AddSpecies( const CLIMATOLOGY_HANDLE& species, skClimatology* numberdensityclimatology, skOpticalProperties* particleopticalprops)
{
	iterator							iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21		dummy(species);
	SKTRAN_AtmosphericOpticalStateEntry_V21*	entry = nullptr;
	bool								ok;

	ok =  (!(species == SKCLIMATOLOGY_PRESSURE_PA) && !(species == SKCLIMATOLOGY_TEMPERATURE_K));
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::AddSpecies, do not use SKCLIMATOLOGY_PRESSURE_PA or SKCLIMATOLOGY_TEMPERATURE_K within AddSpecies. Temperature and Pressure do not have absorption and scattering cross-sections Only use species that have cross-sections");
	}
	else
	{
		ok = numberdensityclimatology->IsSupportedSpecies(species);
		NXASSERT((ok));
		if (ok)
		{
			iter = std::find( m_species.begin(), m_species.end(), dummy );						// See if this entry already exists
			if (iter == m_species.end())														// If it does not
			{																					// then
				ok = (particleopticalprops != nullptr);											// The user must define optical properties for new entries.
				if (!ok)
				{
					nxLog::Record( NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::AddSpecies, You cannot add a new entry if optical property is empty");
				}
				else
				{
					m_species.push_back( dummy );													// create a new entry on the back
					entry = &(m_species.back());													// and get the pointer to this entry
				}
			}																					// otherwise
			else																				// this entry does exist
			{																					// so
	//			nxLog::Record(NXLOG_INFO, "SKTRAN_AtmosphericOpticalState_V21::AddSpecies, You are replacing an existing species in the list with the instance just passed in");
				entry = &(*iter);																// get its address
				if ( particleopticalprops     == nullptr) particleopticalprops  = entry->ParticleOpticalProps();	// Use existing properties if the user does not give us one.
			}																					// done that
			ok = ok && entry->Configure( species, numberdensityclimatology, particleopticalprops );	// now configure the entry
		}
		SetDirty();																			// set the list as dirty
		SetPendingCacheUpdate();
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_AtmosphericOpticalState_V21::AddSpecies, Error adding the requested species  to the list");
	}
	return ok;																			// and we are done
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::RemoveSpecies		2013-6-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::RemoveSpecies( const CLIMATOLOGY_HANDLE& species)
{
    bool ok(true);

	iterator							iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21	dummy(species);

	iter = std::find( m_species.begin(), m_species.end(), dummy );						// See if this entry already exists
	if (!(iter == m_species.end()))														// If it does not
	{																					// then
		m_species.erase( iter );														// erase this entry
		SetDirty();																		// and flag that the cross-sections need updating
//		SetPendingCacheUpdate();														// but the remaining climatologies caches are unaffectedso dont update them
	}

    ok = ok && m_atmosphericemission.RemoveEmission( species );                         // Also remove emission species

	return ok;																			// and we are done
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::SetAtmosphericStateModel		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::SetAtmosphericStateModel( skClimatology* atmosphericstate )
{
	if (atmosphericstate   != NULL) atmosphericstate->AddRef();
	if (m_atmosphericstate != NULL) m_atmosphericstate->Release();
	m_atmosphericstate = atmosphericstate;
	SetDirty();
	SetPendingCacheUpdate();
	//m_atmosphericemission.SetAtmosphericStateModel(atmosphericstate);
	return (atmosphericstate != NULL);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::GetAtmosphericStateModel		2009-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_AtmosphericOpticalState_V21::GetAtmosphericStateModel( skClimatology** statemodel )
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

/*-----------------------------------------------------------------------------
 *					SKTRAN_AtmosphericOpticalState_V21::DeepCopy		2008-3-5*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_AtmosphericOpticalState_V21::DeepCopy( const SKTRAN_AtmosphericOpticalState_V21& other )
{
	bool						ok;
	bool						ok1;
	const_iterator				iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21	blank;


	NXTRACE_ONCEONLY(firsttime,("SKTRAN_AtmosphericOpticalState_V21::DeepCopy. DeepCopy should look for entries with identical climatologies and properties and only clone them one. At the moment common climatologies may get cloned a few times\n"));

	ReleaseResources();
	if (m_atmosphericstate != NULL ) m_atmosphericstate->Release();
    m_atmosphericstate = NULL;

	m_defaultatmosphericstate->DeepCopy( *(other.m_defaultatmosphericstate) );

	ok = (other.m_atmosphericstate == other.m_defaultatmosphericstate) ;
	if (ok)
	{
		m_atmosphericstate = m_defaultatmosphericstate;
		m_atmosphericstate->AddRef();
	}
	else
	{
		ok = other.m_atmosphericstate->CreateClone( &m_atmosphericstate );
	}
	m_wavenumber	         = other.m_wavenumber;
	m_placeandtime		     = other.m_placeandtime;
	m_updateclimatologycache = other.m_updateclimatologycache;
	m_kabs                   = other.m_kabs;
	m_kext                   = other.m_kext;
	m_kscat                  = other.m_kscat;
	m_isdirty                = true;

	for (iter = other.m_species.begin(); !(iter == other.m_species.end()); ++iter)
	{
		m_species.push_back(blank);
		ok1 = m_species.back().DeepCopy( *iter );
		ok = ok && ok1;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_AtmosphericOpticalState_V21::DeepCopy, Error copying species list to this instance");
		ReleaseResources();
	}
	return ok;
}
*/

/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericOpticalState_V21::Get_RotatedPhaseMatrix		2003-12-9
 *-------------------------------------------------------------------------*/

/*
nxBOOL SKTRAN_AtmosphericOpticalState_V21::Get_RotatedPhaseMatrix( double mu, double muprime, double dphi, nxVRTEPhaseMatrix* rotatedmatrix )
{
	double				costheta;
	nxVRTEPhaseMatrix	phasematrix;
	nxBOOL				ok;

	CheckCosineRange(&mu);
	CheckCosineRange(&muprime);
	costheta =  phasematrix.GetScatteringAngle( mu, muprime, dphi );
	ok = Get_PhaseMatrix( costheta, &phasematrix );
	if (ok) phasematrix.ApplyStokesRotation( mu, muprime, dphi, rotatedmatrix );
	return ok;
}
*/

/*---------------------------------------------------------------------------
 *'					SKTRAN_AtmosphericOpticalState_V21::Get_RotatedPhaseMatrix		2004-4-27
 *-------------------------------------------------------------------------*/

/*nxBOOL SKTRAN_AtmosphericOpticalState_V21::Get_RotatedPhaseMatrix( nxVector& observer, nxGeodetic& geoid, nxVector& incoming, nxVector& outgoing, nxVRTEPhaseMatrix* phase)
{
	double				mu;
	double				muprime;
	double				phi;
	nxVector			inh;
	nxVector			outh;
	nxVRTEPhaseMatrix*	rotatedmatrix;
	nxBOOL				ok;

	geoid.FromGeocentric ( observer );
	geoid.GetGeodeticWestSouthUp( &west, &south, &up );
	muprime = up & lookdir;				// Get cosine of incoming zenith angle
	mu      = up & outdir;				// Get cosine of outgoing zenith angle
	inh     = lookdir.ComponentPerpendicularTo(up).UnitVector();		// Horizontal component of sun unit vector
	inh     = outdir.ComponentPerpendicularTo(up).UnitVector();		//Horizonatl component of outgoing direction
	phi     = acos( outh &  inh );
	ok		= species.Get_RotatedPhaseMatrix( muprime, mu, phi, phase );
	return ok;
}
*/


