#include <skopticalproperties21.h>


/*
class skOpticalProperties_MultipleOverlappingSpectra : public skOpticalProperties
{
	private:
		std::list< entry_struct >						m_entries;		

	private:
		bool											DeepCopy( const skOpticalProperties_MultipleOverlappingSpectra& other );

	public:
		bool skOpticalProperties_MultipleOverlappingSpectra::SetAtmosphericState					( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt,  bool* crosssections_changed );			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		bool skOpticalProperties_MultipleOverlappingSpectra::CalculateCrossSections				( double wavenumber );			//!< Calculate cross-sections at the specified wave-number.
		bool skOpticalProperties_MultipleOverlappingSpectra::IsScatterer							() const;						//!< Returns true if this particles scatters radiation
		bool skOpticalProperties_MultipleOverlappingSpectra::IsAbsorber							() const;						//!< Returns true if this particles absorbs radiation radiation
		bool skOpticalProperties_MultipleOverlappingSpectra::CalculatePhaseMatrix				( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix );
		bool skOpticalProperties_MultipleOverlappingSpectra::CreateClone							( skOpticalProperties** clone) const;
};

*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::entry_struct::entry_struct		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MultipleOverlappingSpectra::entry_struct::entry_struct()
{
	m_lowerwavelength = 0.0;
	m_upperwavelength = 0.0;
	m_crosssection    = NULL;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::entry_struct::~entry_struct		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MultipleOverlappingSpectra::entry_struct::~entry_struct()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::entry_struct:ReleaseResources		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_MultipleOverlappingSpectra::entry_struct::ReleaseResources()
{
	if (m_crosssection != NULL)
	{
		m_crosssection->Release();
		m_crosssection = NULL;
	}
	m_lowerwavelength = 0.0;
	m_upperwavelength = 0.0;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::entry_struct::SetValues		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::entry_struct::SetValues( double lowerwavelength_nm, double upperwavelength_nm, skOpticalProperties* crosssection)
{
	bool	ok;

	m_lowerwavelength = lowerwavelength_nm;											// Copy over the lower wavenumber from the upper wavelength
	m_upperwavelength = upperwavelength_nm;											// Copy over the upper wavenumber from the lower wavelength
	if (crosssection   != NULL) crosssection->AddRef();
	if (m_crosssection != NULL) m_crosssection->Release();
	m_crosssection = crosssection;

	NXASSERT(  ( m_lowerwavelength <= m_upperwavelength) );
	ok = ( m_lowerwavelength <= m_upperwavelength);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::entry_struct::SetValues, Error Setting the values, either cloning or setting wavenumber ranges");
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::entry_struct::DeepCopy		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_MultipleOverlappingSpectra::entry_struct::DeepCopy( const entry_struct& other )
{
	bool	ok;

	ok = SetValues( other.m_lowerwavelength, other.m_upperwavelength, other.m_crosssection );
	return ok;
}
*/



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::skOpticalProperties_MultipleOverlappingSpectra		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MultipleOverlappingSpectra::skOpticalProperties_MultipleOverlappingSpectra()
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::~skOpticalProperties_MultipleOverlappingSpectra		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MultipleOverlappingSpectra::~skOpticalProperties_MultipleOverlappingSpectra()
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::DeepCopy		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_MultipleOverlappingSpectra::DeepCopy( const skOpticalProperties_MultipleOverlappingSpectra& other )
{
	const_iterator	iter;
	entry_struct	dummy;
	bool			ok = true;
	bool			ok1;

	m_entries.clear();																	// Delete all the entries
	for (iter = other.m_entries.begin(); !(iter == other.m_entries.end()); ++iter)		// for all of the entries 
	{																					// in the other list
		m_entries.push_back( dummy );													// Add a dummy entry onto the end of teh list
		ok1 = m_entries.back().DeepCopy( *iter );										// Now copy the entry over to the last value
		ok = ok && ok1;																	// and update the status
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::DeepCopy, There was an error making a DeepCopy of the spectra");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::SetAtmosphericState		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::SetAtmosphericState( skClimatology* neutralatmosphere)
{
	iterator						iter;
	skOpticalProperties*		crosssection;
	bool							ok = true;
	bool							ok1;

	for (iter = m_entries.begin(); !(iter == m_entries.end()); ++iter)								// Iterate through 
	{																								// all of the entries in the list
		crosssection = (*iter).m_crosssection;														// get the cross-section object
		if (crosssection != NULL)																	// if it has a valid pointer
		{																							// then
			ok1      = crosssection->SetAtmosphericState( neutralatmosphere);						// Change the atmospheric state
			ok       = ok && ok1;																	// make sure it worked
		}																							// otherwise ignore if cross-section object is NULL
	}																								// do all of the cross-section objects in this class
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::SetAtmosphericState, There was an error setting the atmospheric state in at least one of the underlying cross-sections");
	}
	return ok;																						// return the overall status
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::SetLocation		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::SetLocation( const GEODETIC_INSTANT& pt,  bool* crosssections_changed )
{
	iterator						iter;
	skOpticalProperties*		crosssection;
	bool							haschanged;
	bool							ischange = false;
	bool							ok = true;
	bool							ok1;

	for (iter = m_entries.begin(); !(iter == m_entries.end()); ++iter)								// Iterate through 
	{																								// all of the entries in the list
		crosssection = (*iter).m_crosssection;														// get the cross-section object
		if (crosssection != NULL)																	// if it has a valid pointer
		{																							// then
			ok1      = crosssection->SetLocation( pt,  &haschanged);								// Change the atmospheric state
			ok       = ok && ok1;																	// make sure it worked
			ischange = ischange || haschanged;														// and see if it has changed the cross-section data
		}																							// otherwise ignore if cross-section object is NULL
	}																								// do all of the cross-section objects in this class
	if (crosssections_changed != NULL) *crosssections_changed = ischange;							// set the cross-sections changed flag if requested by user
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::SetLocation, There was an error calling SetLocation in at least one of the underlying cross-sections");
	}
	return ok;																						// return the overall status
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::InternalClimatology_UpdateCache		2011-8-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt )
{
	iterator						iter;
	skOpticalProperties*		crosssection;
	bool							ok = true;
	bool							ok1;

	for (iter = m_entries.begin(); !(iter == m_entries.end()); ++iter)								// Iterate through 
	{																								// all of the entries in the list
		crosssection = (*iter).m_crosssection;														// get the cross-section object
		if (crosssection != NULL)																	// if it has a valid pointer
		{																							// then
			ok1      = crosssection->InternalClimatology_UpdateCache( pt );		// Change the atmospheric state
			ok       = ok && ok1;																	// make sure it worked
		}																							// otherwise ignore if cross-section object is NULL
	}																								// do all of the cross-section objects in this class
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::InternalClimatology_UpdateCache, There was an error updating the atmospheric state in at least one of the underlying cross-sections");
	}
	return ok;																						// return the overall status
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::IsScatterer		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::IsScatterer() const
{
	const_iterator						iter;
	const skOpticalProperties*	crosssection;
	bool								isscatterer = false;

	for (iter = m_entries.begin(); !(iter == m_entries.end()); ++iter)								// Iterate through 
	{																								// all of the entries in the list
		crosssection = (*iter).m_crosssection;														// get the cross-section object
		if (crosssection != NULL)																	// if it has a valid pointer
		{																							// then
			isscatterer = isscatterer || crosssection->IsScatterer();
		}																							// otherwise ignore if cross-section object is NULL
	}																								// do all of the cross-section objects in this class
	return isscatterer;																						// return the overall status
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::IsAbsorber		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::IsAbsorber() const
{
	const_iterator						iter;
	const skOpticalProperties*	crosssection;
	bool								isabsorber = false;

	for (iter = m_entries.begin(); !(iter == m_entries.end()); ++iter)								// Iterate through 
	{																								// all of the entries in the list
		crosssection = (*iter).m_crosssection;														// get the cross-section object
		if (crosssection != NULL)																	// if it has a valid pointer
		{																							// then
			isabsorber = isabsorber || crosssection->IsAbsorber();
		}																							// otherwise ignore if cross-section object is NULL
	}																								// do all of the cross-section objects in this class
	return isabsorber;																				// return the overall status

}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	entry_struct*		entry;
	bool				ok;
	
	ok = FindBoundingEntry(wavenumber, &entry);				// This function must be thread safe (It is)
	ok = ok && (entry->m_crosssection != NULL);
	if (ok)
	{
		ok   = entry->m_crosssection->CalculateCrossSections( wavenumber, absxs, extxs, scattxs );
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::CalculateCrossSections, No spectrum found for wavenumber %e", (double)wavenumber);
		*absxs   = 0.0;
		*extxs   = 0.0;
		*scattxs = 0.0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::CalculatePhaseMatrix		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	entry_struct*		entry;
	bool				ok;
	
	ok = FindBoundingEntry(wavenumber, &entry);
	ok = ok && (entry->m_crosssection != NULL);
	if (ok)
	{
		ok   = entry->m_crosssection->CalculatePhaseMatrix( wavenumber, cosscatterangle, phasematrix);
	}
	if (!ok)
	{
		phasematrix->SetTo(0.0);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::CreateClone		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_MultipleOverlappingSpectra::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_MultipleOverlappingSpectra*	clone;
	bool										ok;

	clone = new skOpticalProperties_MultipleOverlappingSpectra;
	ok = (clone != NULL);
	if (ok)
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::CreateClone, Error creating clone, thats probably going to be a problem");
	}
	*userclone = clone;
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::FindBoundingEntry		2009-11-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::FindBoundingEntry( double wavenumber, entry_struct** entry)
{
	double			wavelen;
	bool			ok;
	iterator		iter;

	wavelen = 1.0E7/wavenumber;
	iter    = std::find_if( m_entries.begin(), m_entries.end(), entry_locator(wavelen) );
	ok      = !(iter == m_entries.end());
	*entry  =  ok ? &(*iter) : NULL;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MultipleOverlappingSpectra::AddEntry		2009-11-23*/
/** Adds**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MultipleOverlappingSpectra::AddEntry( double startlambda_nm, double endlambda_nm, skOpticalProperties* crosssection )
{
	entry_struct	dummy;
	bool			ok;

	m_entries.push_back( dummy );												
	ok = m_entries.back().SetValues( startlambda_nm, endlambda_nm, crosssection );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MultipleOverlappingSpectra::AddEntry, There was an error adding this entry to the table of overlapping spectra");
	}
	return ok;
}
