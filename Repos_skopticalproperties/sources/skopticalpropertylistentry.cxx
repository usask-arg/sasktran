#include <skopticalproperties21.h>


/*---------------------------------------------------------------------------
 *  skOpticalProperties_ListEntry::skOpticalProperties_ListEntry  2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ListEntry::skOpticalProperties_ListEntry( double startwavenumber, double endwavenumber,  skOpticalProperties* optprop)
	: m_startwavenumber( startwavenumber),
	  m_endwavenumber  ( endwavenumber  )
{
	if (optprop != nullptr) optprop->AddRef();
	m_opticalproperty = optprop;
}

/*---------------------------------------------------------------------------
 *  skOpticalProperties_ListEntry::Copy Constructor  2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ListEntry::skOpticalProperties_ListEntry( const skOpticalProperties_ListEntry& other )
	: m_startwavenumber( other.m_startwavenumber),
	  m_endwavenumber( other.m_endwavenumber)
{
	if (other.m_opticalproperty != nullptr) other.m_opticalproperty->AddRef();
	m_opticalproperty = other.m_opticalproperty;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_ListEntry::~skOpticalProperties_ListEntry  2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ListEntry::~skOpticalProperties_ListEntry()
{
	if (m_opticalproperty != nullptr) m_opticalproperty->Release();
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_ListEntries::skOpticalProperties_ListEntries2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ListEntries::skOpticalProperties_ListEntries()
{
}

/*---------------------------------------------------------------------------
 * skOpticalProperties_ListEntries::~skOpticalProperties_ListEntries2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ListEntries::~skOpticalProperties_ListEntries()
{
}

/*---------------------------------------------------------------------------
 *           skOpticalProperties_ListEntries::FindEntry           2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties* skOpticalProperties_ListEntries::FindEntry( double wavenum )
{
	skOpticalProperties* optprop;
	bool				 ok;
	iterator			iter = std::lower_bound( m_entries.begin(), m_entries.end(), wavenum);
	
	ok       = (iter != m_entries.end());
	ok       = ok && (wavenum >= iter->m_startwavenumber);
	optprop  = ok  ? iter->m_opticalproperty : nullptr;
	return optprop;
}

/*---------------------------------------------------------------------------
 *      skOpticalProperties_ListEntries::CheckNonOverlapping      2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::CheckNonOverlapping(double wavenum)
{
	bool	ok = true;
	bool	overlap;
	for (iterator iter = m_entries.begin(); iter != m_entries.end(); ++iter)
	{
		overlap = ( wavenum  >= iter->m_startwavenumber) && ( wavenum <= iter->m_endwavenumber);
		ok = ok && !overlap;
	}
	return ok;
}
/*---------------------------------------------------------------------------
 *           skOpticalProperties_ListEntries::AddEntry            2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::AddEntry( double startwavenum, double endwavenum, skOpticalProperties* optprop)
{
	bool	ok;

	ok = CheckNonOverlapping(startwavenum) && CheckNonOverlapping( endwavenum);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ListEntries::AddEntry, Cannot add entry %e to %e as it overlaps with existing entries", (double)startwavenum, (double)endwavenum);
	}
	else
	{
		iterator iter   = std::lower_bound( m_entries.begin(), m_entries.end(), endwavenum);
		iterator status = m_entries.insert( iter, skOpticalProperties_ListEntry( startwavenum, endwavenum,optprop) );
		ok = status != m_entries.end();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ListEntries::AddEntry, Error inserting entry %e to %e.", (double)startwavenum, (double)endwavenum);
		}
	}
	return	ok;
}

/*---------------------------------------------------------------------------
 *      skOpticalProperties_ListEntries::SetAtmosphericState      2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::SetAtmosphericState( skClimatology* neutralatmosphere)
{
	bool	ok = true;
	bool	ok1;
	for (iterator iter = m_entries.begin(); iter != m_entries.end(); ++iter)
	{
		ok1 = iter->m_opticalproperty->SetAtmosphericState( neutralatmosphere);
		ok  = ok && ok1;
	}
	return ok;

}

/*---------------------------------------------------------------------------
 *          skOpticalProperties_ListEntries::SetLocation          2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )
{
	bool	ok = true;
	bool	ok1;
	bool	xschanged;

	*crosssectionschanged = false;
	for (iterator iter = m_entries.begin(); iter != m_entries.end(); ++iter)
	{
		ok1 = iter->m_opticalproperty->SetLocation( pt, &xschanged);
		*crosssectionschanged = *crosssectionschanged || xschanged;
		ok  = ok && ok1;
	}
	return ok;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_ListEntries::InternalClimatology_UpdateCache 2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool	ok = true;
	bool	ok1;

	for (iterator iter = m_entries.begin(); iter != m_entries.end(); ++iter)
	{
		ok1 = iter->m_opticalproperty->InternalClimatology_UpdateCache( pt );
		ok  = ok && ok1;
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *    skOpticalProperties_ListEntries::CalculateCrossSections     2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	bool					ok = true;
	skOpticalProperties*	optprop  = FindEntry( wavenumber);

	if (optprop == nullptr)
	{
		*absxs = 0.0;
		*extxs = 0.0;
		*scattxs = 0.0;
	}
	else
	{
		ok = optprop->CalculateCrossSections( wavenumber, absxs, extxs, scattxs);
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *     skOpticalProperties_ListEntries::CalculatePhaseMatrix      2020-01-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ListEntries::CalculatePhaseMatrix( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool					ok = true;
	skOpticalProperties*	optprop  = FindEntry( wavenumber);

	if (optprop == nullptr)
	{
		phasematrix->SetTo(0.0);
		phasematrix->At(1,1) = 1.0;
	}
	else
	{
		ok = optprop->CalculatePhaseMatrix( wavenumber, cosscatterangle, phasematrix);
	}
	return ok;
}

