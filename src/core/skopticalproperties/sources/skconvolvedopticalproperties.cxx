#include <skopticalproperties21.h>
#include <nxbase_threads.h>


/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfo_TemperatureDependent::KeyedIndexFromAtmosphericState		2013-6-24*/
/** Returns a unique index for this atmospheric state. This is used for 
 *	cross-sections that depend only upon atmospheric temperature (and wavelength)
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperty_AdditionalStateInfo_TemperatureDependent::KeyedIndexFromAtmosphericState( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, skOpticalProperty_AdditionalStateInfoKey* index ) 
{
	bool	ok;
	double	T[1];

	ok     =       neutralatmosphere->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, pt, &T[0], false );
	ok     = ok && index->SetKeyStateParameters( T, 1);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperty_AdditionalStateInfo_TemperatureDependent::KeyedIndexFromAtmosphericState, There were errors retrieving temperature for given atmospheric state");
		index->Clear();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfo_PressTemperatureDependent::KeyedIndexFromAtmosphericState		2013-6-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperty_AdditionalStateInfo_PressTemperatureDependent::KeyedIndexFromAtmosphericState( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, skOpticalProperty_AdditionalStateInfoKey* index ) 
{
	bool	ok;
	double  param[2];

	ok  =       neutralatmosphere->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, pt, &param[0], false );
	ok  = ok && neutralatmosphere->GetParameter( SKCLIMATOLOGY_PRESSURE_PA,   pt, &param[1], false );
	ok =  ok && index->SetKeyStateParameters( param, 2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperty_AdditionalStateInfo_PressTemperatureDependent::KeyedIndexFromAtmosphericState, There were errors retrieving pressure and temperature for given atmospheric state");
		index->Clear();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *				skOpticalProperties_ConvolvedDiscreteWavelenCachedState::skOpticalProperties_ConvolvedDiscreteWavelenCachedState		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenCachedState::skOpticalProperties_ConvolvedDiscreteWavelenCachedState()
{
	m_highresopticalproperties = NULL;
	m_highresdetails           = NULL;
	m_atmosphericstateinfo     = NULL;
	m_entriestable             = NULL;
	m_backgroundatmosphere     = NULL;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::~skOpticalProperties_ConvolvedDiscreteWavelenCachedState		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenCachedState::~skOpticalProperties_ConvolvedDiscreteWavelenCachedState()
{
	if (m_highresopticalproperties != NULL) m_highresopticalproperties->Release();
	if (m_backgroundatmosphere     != NULL) m_backgroundatmosphere->Release();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CreateClone		2012-7-31*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CreateClone( skOpticalProperties** clone) const
{
	nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState DOES NOT support CreateClone. This should not be needed if you are using sasktran V2.1 or higher");
	*clone = NULL;
	return false;
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetHighResOpticalProperties		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetHighResOpticalProperties( skOpticalProperties* highresopticalproperties,
																			   skWavelengthToPSF_Table* measurementdetails,
																			   skOpticalProperty_AdditionalStateInfo* atmosphericstateinfo)
{
	if (  highresopticalproperties != NULL)   highresopticalproperties->AddRef();
	if (m_highresopticalproperties != NULL) m_highresopticalproperties->Release();
	if (m_highresopticalproperties != NULL && (m_backgroundatmosphere != NULL)) m_highresopticalproperties->SetAtmosphericState(m_backgroundatmosphere);
	m_highresopticalproperties = highresopticalproperties;
    m_highresdetails           = measurementdetails;
	m_atmosphericstateinfo     = atmosphericstateinfo;
	m_cachedentriestable.clear();											// Clear all of the 
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetAtmosphericState		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetAtmosphericState( skClimatology* neutralatmosphere )
{
	if (neutralatmosphere      != NULL) neutralatmosphere->AddRef();
	if (m_backgroundatmosphere != NULL) m_backgroundatmosphere->Release();
	m_backgroundatmosphere = neutralatmosphere;
	if (m_highresopticalproperties != nullptr) m_highresopticalproperties->SetAtmosphericState( neutralatmosphere);
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetAtmosphericState		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )
{
	bool		ok;
	bool										found;
	skOpticalProperty_AdditionalStateInfoKey	index;
	iterator	iter;

	NXASSERT(( m_highresopticalproperties != NULL ) && (m_atmosphericstateinfo != NULL) );

	ok = ((m_atmosphericstateinfo != NULL) && (m_highresopticalproperties != NULL) && (m_backgroundatmosphere != nullptr));
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetLOcation, cannot SetLocation as the as the Atmospheric State Dependency object is not set");
	}
	else
	{
		ok = m_atmosphericstateinfo->KeyedIndexFromAtmosphericState( m_backgroundatmosphere, pt, &index );		// Use the atmopsheric state to produce a unique key for given conditions (typically using Temperature and maybe pressure)
		if (!ok)																							// ensure that it worked
		{																									//
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetAtmosphericState	, There were errors looking up the keyed index of the current atmopsheric state");
		}
		else
		{
			iter  = m_cachedentriestable.find(index);														// See if we have a convolved spectrum for this atmospheric state 		
			found = !(iter == m_cachedentriestable.end());													// See if we found the entry
			if (!found)																										// If it does not yet exist
			{																												// then create the entry
				value_type					dummy(index, skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable() );		// so create a new blank entry
				std::pair <iterator, bool>	pr;																				// variable to hold creation result

				pr   = m_cachedentriestable.insert( dummy );																// do the creation
				ok   = pr.second;																							// get the creation status
				iter = pr.first;																							// and iterator to the new entry
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetAtmosphericState, Error creatinga new cache entry for the given stamospheric state");
				}
			}
		}
	}
	m_entriestable = ok ? &(iter->second) : NULL;
	ok = ok && m_highresopticalproperties->SetLocation( pt, crosssectionschanged);
	if (!ok) nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::SetAtmosphericState, Error setting atmopsheric state. Make sure the highres optical properties is properly set");
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::InternalClimatology_UpdateCache		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool	ok;

	NXASSERT(( m_highresopticalproperties != NULL ));
	ok =       ( m_highresopticalproperties != NULL) ;
	ok = ok && ( m_highresopticalproperties->InternalClimatology_UpdateCache( pt) );
	if (!ok) nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::InternalClimatology_UpdateCache, Error invoking InternalClimatology_UpdateCache. Make sure the highres optical properties is properly set");
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculateCrossSections		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	skOpticalProperties_ConvolvedDiscreteWavelenEntry*	entry = nullptr;
	bool						ok;
	double						wavelennm;
	double						fwhm_nm;
	double						highresspacing;
	double						highresfwhm;
	double						adjustedfwhm;
	static std::mutex			g_mutex_xsectionslock;

	ok = (m_entriestable != NULL);																			// Do we have a valid entries table for a specific atmospheric state.
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculateCrossSections, the cached entries object is null. You probably need to successfully call SetAtmopshericState");
	}
	ok = ok && (m_entriestable->FindEntry( wavenumber, &entry ));
	if (!ok)
	{
		std::unique_lock<std::mutex>	lock(g_mutex_xsectionslock);					// The code is mutexed to make it thread safe. This is "okay" as we only call this code when we need to make mie scattering tables
																							// So we dont have to worry about parallelism and speed execution for the moment.	

		wavelennm      = 1.0E7/wavenumber;
		highresspacing = m_highresdetails->GetInstrumentPointSpacing(wavelennm);							// Get the nominal spacing of data points in the high resolution cross-section data data 
		highresfwhm    = m_highresdetails->GetInstrumentPSF_FWHM( wavelennm );								// get the FWHM spectral resolution of high res data in nm 
		ok             = GetTargetFWHM( wavelennm, &fwhm_nm);
		ok             = ok &&  (highresfwhm < fwhm_nm);													// make sure we are not convolving with a meaurement that is wider than our desired FWHM.
		if (ok)
		{
			adjustedfwhm = sqrt( fwhm_nm*fwhm_nm - highresfwhm*highresfwhm);								// adjust the desired FWHM so when convolved with desired cross-section it comes back to fwhm_nm.
			ok    =       m_entriestable->AddEntry ( wavenumber, m_highresopticalproperties, highresspacing, adjustedfwhm);
			ok    = ok && m_entriestable->FindEntry( wavenumber, &entry );
		}
	}
	if (ok)
	{
		*extxs   = entry->ConvolvedExtinction ();
		*absxs   = entry->ConvolvedAbsorption ();
		*scattxs = entry->ConvolvedScattering ();
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculateCrossSections, Error calculating convolved optical properties, Setting all cross sections to zero");
		*extxs   = 0.0;
		*absxs   = 0.0;
		*scattxs = 0.0;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculatePhaseMatrix		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	bool ok;
	NXASSERT(( m_highresopticalproperties != NULL ));
	ok =       ( m_highresopticalproperties != NULL) ;
	ok = ok && ( m_highresopticalproperties->CalculatePhaseMatrix( wavenumber, cosscatterangle, phasematrix) );
	if (!ok) nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::CalculatePhaseMatrix, Error invoking CalculatePhaseMatrix. Make sure the highres optical properties is properly set");
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsScatterer		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsScatterer() const
{
	bool	ok;
	bool	isscatter = false;

	NXASSERT(( m_highresopticalproperties != NULL ));
	ok =       ( m_highresopticalproperties != NULL) ;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsScatterer, Error invoking IsScatterer. Make sure the highres optical properties is properly set");
		isscatter = false;
	}
	else
	{
		isscatter = m_highresopticalproperties->IsScatterer();
	}
	return isscatter;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsAbsorber		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsAbsorber() const
{
	bool	ok;
	bool	isabsorber = false;

	NXASSERT(( m_highresopticalproperties != NULL ));
	ok =       ( m_highresopticalproperties != NULL) ;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsAbsorber, Error invoking IsAbsorber. Make sure the highres optical properties is properly set");
		isabsorber = false;
	}
	else
	{
		isabsorber = m_highresopticalproperties->IsAbsorber();
	}
	return isabsorber;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsDeltaFunctionForwardScatter		2012-5-4*/
/** **/
/*---------------------------------------------------------------------------*/

double skOpticalProperties_ConvolvedDiscreteWavelenCachedState::DeltaFunctionForwardScatterFraction() const
{
	bool	ok;
	double	delta;

	NXASSERT(( m_highresopticalproperties != NULL ));
	ok =     ( m_highresopticalproperties != NULL) ;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenCachedState::IsDeltaFunctionForwardScatter, Error invoking IsDeltaFunctionForwardScatter. Make sure the highres optical properties is properly set");
		delta = 0.0;
	}
	else
	{
		delta = m_highresopticalproperties->DeltaFunctionForwardScatterFraction();
	}
	return delta;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM::CreateClone		2012-8-2*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM::CreateClone( skOpticalProperties** clone) const
{
	skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM*	aclone;
	bool													ok;

	aclone = new skOpticalProperties_ConvolvedDiscreteWavelenCachedStateFixedFWHM;
	ok = (aclone != NULL);
	ok = ok && aclone->skOpticalProperties_ConvolvedDiscreteWavelenCachedState::DeepCopy( *this );
	if (ok)	aclone->m_fwhm = m_fwhm;
	return ok;
}
*/
