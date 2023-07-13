#include <skopticalproperties21.h>
#include <limits>


/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_TableArray::AddEntry		2012-7-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skWavelengthToPSF_TableArray::AddPSFEntry( double nm, double fwhm)
{
	iterator	iter;
	bool		ok;


	iter = m_entries.find( nm );
	ok = !(iter == m_entries.end());
	if (ok)
	{
		(*iter).second = fwhm;
	}
	else
	{
		std::pair<iterator, bool>	status;
		std::pair<double, double>	value( nm, fwhm);

		status = m_entries.insert( value );
		ok    = status.second;
	}
	if (!ok) nxLog::Record(NXLOG_WARNING,"skWavelengthToPSF_TableArray::AddEntry, Error adding entry wavelength=%f, PSF_FWHM=%f to table", (double)nm , (double)fwhm);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_TableArray::GetInstrumentPSF_FWHM		2012-7-24*/
/** Lookup the spectral resolution from the table for the given wavelength.
 *	We just need to lookup the element containing this wavelength. No interpolation
 *	required.
 **/
/*---------------------------------------------------------------------------*/

double skWavelengthToPSF_TableArray::GetInstrumentPSF_FWHM( double nm ) const
{
	const_iterator					iter;
	bool							ok;
	double							fwhm = 0.0;

	ok = (m_entries.size() > 0);
	if (ok)
	{
		iter = m_entries.lower_bound( nm );
		if (!(iter == m_entries.begin())) --iter;
		fwhm = (*iter).second;
	}
	return fwhm;
}
	

/*-----------------------------------------------------------------------------
 *				skWavelengthToPSF_TableConstantWavenumber::GetInstrumentPSF_FWHM		2012-7-24*/
	
/** Calculate the spectral resolution in nanometers at wavelength nm given the
 *	spectral resolution in wavenumber.
 *		-# lambda  = 1.0E7/k 
 *		-# dlambda = 1.0E7/k^2*deltak
 *		-# dlambda = lambda^2/1.0E7*deltak
 **/
/*---------------------------------------------------------------------------*/

double skWavelengthToPSF_TableConstantWavenumber::GetInstrumentPSF_FWHM( double nm  ) const
{
	return (nm*nm*m_deltak*1.0E-7);
}


/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_TableConstantWavenumber::GetInstrumentPointSpacing		2012-7-27*/
/** Get the instrument spacing of points in wavelength given that it is at fixed 
 *	wavenumber spacing
 **/
/*---------------------------------------------------------------------------*/

double skWavelengthToPSF_TableConstantWavenumber::GetInstrumentPointSpacing( double nm  ) const
{
	return (nm*nm*m_pointspacingcm1*1.0E-7);
}

/*-----------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::sk_AbsorptionTabulatedTableEntry		2012-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

sk_AbsorptionTabulatedTableEntry::sk_AbsorptionTabulatedTableEntry	()
{
	m_t     = 0.0;
	ClearMinMaxRange();
}


/*-----------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::sk_AbsorptionTabulatedTableEntry		2012-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

sk_AbsorptionTabulatedTableEntry::sk_AbsorptionTabulatedTableEntry( double t)
{
	m_t     = t;
	ClearMinMaxRange();
}

/*-----------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::ClearTandMinMaxRange		2012-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

void sk_AbsorptionTabulatedTableEntry::ClearMinMaxRange()
{
	m_minnm = 0.0;
	m_maxnm = -1.0; 
}

/*-----------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::CheckWavelengths		2012-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_AbsorptionTabulatedTableEntry::CheckWavelengths()
{
	bool	ok;
	bool	ok1;
	size_t	i;
	double	lastnm;
	double	nm;
	
	m_minnm =  999999.0;
	m_maxnm =  -1;
	ok = (m_nm.size() == m_xsection.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_AbsorptionTabulatedTableEntry::CheckWavelengths, The wavelength and cross-section ararys are not equal in size. Thats not good");
	}
	else
	{
		lastnm = -1.0;
		for ( i = 0; i < m_nm.size(); i++ )
		{
			nm      = m_nm.At(i);
			ok1     = (nm > lastnm);
			if (!ok1)
			{
				nxLog::Record(NXLOG_WARNING, "sk_AbsorptionTabulatedTableEntry::CheckWavelengths, Cross sections wavelengths not in ascending order near %f nm", (double)nm ) ;
			}
			m_minnm = (nm < m_minnm) ? nm : m_minnm;
			m_maxnm = (nm > m_maxnm) ? nm : m_maxnm;
			ok      = ok && ok1;
			lastnm  = nm;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"sk_AbsorptionTabulatedTableEntry::CheckWavelengths, The wavelengths are not in ascending order or we have duplicate wavelengths");
		}
	}
	return ok;
}



/*---------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::Configure	 2003-11-28 */
/** Configures the tabulated entry so it uses the actual tables passed in by the user.*/
/*-------------------------------------------------------------------------*/

bool sk_AbsorptionTabulatedTableEntry::Configure( double t, double* nm, intptr_t nmstrides, double *xs, intptr_t xsstrides, size_t npts )
{
	bool	ok;
	size_t	nmstride = nmstrides;
	size_t	xsstride = xsstrides;

	m_nm.erase();
	m_xsection.erase();
	m_t = t;
	ok  =       m_nm.nxArrayLinear<double>::Attach      ( 1, &npts, nm, NULL, &nmstride);
	ok  = ok && m_xsection.nxArrayLinear<double>::Attach( 1, &npts, xs, NULL, &xsstride);
	ok =  ok && CheckWavelengths(); 
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_AbsorptionTabulatedTableEntry::Configure, error attchaing to cross-section arrays. This will create problems. I am destroying Entries");
		m_nm.Detach();
		m_xsection.Detach();
		m_t = 0.0;
		ClearMinMaxRange();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::Configure		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_AbsorptionTabulatedTableEntry::Configure( double t, nx1dArray<double>& nm, nx1dArray<double>& xsection )
{
	bool	ok;

	m_nm.erase();
	m_xsection.erase();
	m_t = t;
	ok  =       m_nm.DeepCopy( nm );
	ok  = ok && m_xsection.DeepCopy( xsection );
	ok  = ok && CheckWavelengths(); 
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_AbsorptionTabulatedTableEntry::Configure, error copying cross-section arrays. This will create problems");
		m_nm.erase();
		m_xsection.erase();
		m_t = 0.0;
		ClearMinMaxRange();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					sk_AbsorptionTabulatedTableEntry::Configure		 2015- 4- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_AbsorptionTabulatedTableEntry::Configure( double t, const std::vector<double>& nm, const std::vector<double>& xsection )
{
	bool	ok;

	m_nm.erase();
	m_xsection.erase();
	m_t = t;
	ok = (nm.size() == xsection.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_AbsorptionTabulatedTableEntry::Configure the nm and xsection arrays must be the same size");
	}
	else
	{
		m_nm.SetSize( nm.size() );
		m_xsection.SetSize( nm.size() );
		for (size_t i = 0; i < nm.size(); i++)
		{
			m_nm.At(i)       = nm.at(i);
			m_xsection.At(i) = xsection.at(i);
		}
		ok  = CheckWavelengths(); 
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"sk_AbsorptionTabulatedTableEntry::Configure, error copying cross-section arrays. This will create problems");
			m_nm.erase();
			m_xsection.erase();
			m_t = 0.0;
			ClearMinMaxRange();
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					sk_AbsorptionTabulatedTableEntry::GetCrossSection		2003-11-28
 *-------------------------------------------------------------------------*/

bool sk_AbsorptionTabulatedTableEntry::GetCrossSection( double nanometer, double* crosssection ) const
{
	nxArrayIter<double>		uppern;
	nxArrayIter<double>		lowern;
	nxArrayIter<double>		begin  = m_nm.begin();
	nxArrayIter<double>		finish = m_nm.end();
	ptrdiff_t				idx1,idx2;
	double					y1 = 0.0;
	double					y2 = 0.0;
	double					y = std::numeric_limits<double>::quiet_NaN( );
	double					x1,x2,dx;
	bool					ok;

	ok = ( nanometer >= m_minnm) && ( nanometer <= m_maxnm);		// Make sure we are asking for a value in our range
	if (ok)
	{
		uppern = std::lower_bound( begin, finish, nanometer );		// Find the first entry in our table bigger than our wavelength
		if (uppern == finish) --uppern;								// if the first entry is beyond our table  then (nanometer == m_maxnm) so popint to point before
		lowern = uppern;											// copy the pointer
		if (lowern != begin)  --lowern;								// If we are not at the start then get the value before

		idx1 = (lowern-begin);
		idx2 = (uppern-begin);

		x1 = *lowern;
		x2 = *uppern;
		y1 = m_xsection.At(idx1);
		y2 = m_xsection.At(idx2);
		dx = (x2-x1);
		if (dx == 0) y = y1;
		else
		{
			y = y1 + (y2-y1)*(nanometer - x1)/dx;
		}
	}
	*crosssection = ok ? y : std::numeric_limits<double>::quiet_NaN( );
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_UserDefinedAbsorption::Set_Temperature		2003-11-28
 *	Overides the default implementation provided in nxVRTEScatterbase
 *	The 
 *-------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::Set_Temperature( double kelvin)
{
	m_temperature =  kelvin;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::FetchNewOrExistingEntryAtTemperature		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

sk_AbsorptionTabulatedTableEntry* skOpticalProperties_UserDefinedAbsorption::FetchNewOrExistingEntryAtTemperature( double t )
{
	iterator								iter;
	sk_AbsorptionTabulatedTableEntry		dummy(t);
	sk_AbsorptionTabulatedTableEntry*		entry;
	bool									isequal;

	iter    = std::lower_bound( m_temperature_entries.begin(), m_temperature_entries.end(), dummy );		// Find the lowest entry point >= to current temperature
	isequal = !(iter == m_temperature_entries.end()) && (dummy == *iter);									// See if they are identical
	if (!isequal ) iter = m_temperature_entries.insert(iter, dummy);										// If they are not then insert a new entry
	entry = &(*iter);																						// either way return the new entry
	return entry;
}

/*---------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::AddEntry		2003-11-28 */
/**	Add this cross-section data to the list of entries. The entries are
 *	sorted into ascending temperature.  I have assumed that this 
 *	activity is carfeully done inside a constructor and have not checked to see
 *	if the entry for this temperature already exists
 **/
/*-------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::AddEntry( double t, double* nm, int nmstride, double *xs, int xsstride, int npts )
{
	sk_AbsorptionTabulatedTableEntry*		entry;

	entry = FetchNewOrExistingEntryAtTemperature( t );
	entry->Configure( t, nm, nmstride, xs, xsstride, npts );
	return true;
}


/*---------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::AddEntry		2003-11-28 */
/**	Add this cross-section data to the list of entries. The entries are
 *	sorted into ascending temperature.  I have assumed that this 
 *	activity is carfeully done inside a constructor and have not checked to see
 *	if the entry for this temperature already exists
 **/
/*-------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::AddAscendingWavenumberEntry( double t, double* usercm, int cmstrides, double *userxs, int xsstrides, int npts )
{
	sk_AbsorptionTabulatedTableEntry*		entry;
	nx1dArray<double>						cm;
	nx1dArray<double>						revxs;
	nx1dArray<double>						nm;
	nx1dArray<double>						xs;
	bool									ok;
	size_t									j;
	size_t									N        = npts;
	size_t									cmstride = cmstrides;
	size_t									xsstride = xsstrides;

	ok  =          cm.nxArrayLinear<double>::Attach( 1, &N, usercm, NULL, &cmstride);
	ok  = ok && revxs.nxArrayLinear<double>::Attach( 1, &N, userxs, NULL, &xsstride);
	ok  = ok && nm.SetSize(N);
	ok  = ok && xs.SetSize(N);
	if (ok)
	{
		j = N-1;
		for (size_t i = 0; i < N; i++, j--)
		{
			nm.at(j) = 1.0E7/cm.at(i);
			xs.at(j) = revxs.at(i);
		}
		entry = FetchNewOrExistingEntryAtTemperature( t );
		ok = entry->Configure( t, nm, xs);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::AddUserEntry		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::AddUserEntry( double kelvin, nx1dArray<double>& nm, nx1dArray<double>& xsection)
{
	sk_AbsorptionTabulatedTableEntry*		entry;

	entry = FetchNewOrExistingEntryAtTemperature( kelvin);
	entry->Configure( kelvin, nm, xsection );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::AddUserEntry		 2015- 4- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::AddUserEntry( double kelvin, const std::vector<double>& usernm, const std::vector<double>& userxsection)
{
	nx1dArray<double>						nm;
	nx1dArray<double>						xsection;
	double*									usernmptr       = (double *)&usernm.front();
	double*									userxsectionptr = (double *)&userxsection.front();
	sk_AbsorptionTabulatedTableEntry*		entry;

	nm.Attach( usernm.size(),usernmptr);
	xsection.Attach( userxsection.size(), userxsectionptr );
	entry = FetchNewOrExistingEntryAtTemperature( kelvin);
	entry->Configure( kelvin, nm, xsection );
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables		2008-3-4*/
/** Copies over the object without the tables
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables( const skOpticalProperties_UserDefinedAbsorption& other )
{
	bool	ok;

	ok = skOpticalProperties::DeepCopy(other);
	m_temperature = other.m_temperature;
	m_isabsorber  = other.m_isabsorber;
	m_refractiveindexair = other.m_refractiveindexair;
	m_wavelengthinvacuum = other.m_wavelengthinvacuum;		// True if the wavelengthtables are in vacuum
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::ShallowCopy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::CopyEntries( const skOpticalProperties_UserDefinedAbsorption& other )
{
	const_iterator		iter;

	m_quietwavelengthtruncation = other.m_quietwavelengthtruncation;
	m_temperature_entries.clear();
	for (iter = other.m_temperature_entries.begin(); !(iter == other.m_temperature_entries.end()); ++iter)
	{
		m_temperature_entries.push_back( *iter );
	}
	return (m_temperature_entries.size() == other.m_temperature_entries.size());
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::DeepCopy		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::DeepCopy( const skOpticalProperties_UserDefinedAbsorption& other )
{
	bool	ok;

	ok =       CopyEntries( other );
	ok = ok && DeepCopyWithoutTables ( other);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_UserDefinedAbsorption::CreateClone, Error copying this object over to clone");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::CreateClone		2008-11-17*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_UserDefinedAbsorption::CreateClone(skOpticalProperties ** userclone) const
{
	skOpticalProperties_UserDefinedAbsorption*	clone;
	bool				ok;

	clone = new skOpticalProperties_UserDefinedAbsorption;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_UserDefinedAbsorption::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = ok && clone->DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_UserDefinedAbsorption::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::InterpolateCrossSection		2012-7-19*/
/** This code interpolates the cross-section data in temperature. The code only uses
 *	temperature records that cover the desired wavelength. This feature was first
 *	introduced to support the Dumont, Brion, Malicet O3 tables where different
 *	temperature ranges cover different wavelength intervals.
 *
 *	The temperature interpolation scheme truncates the interpolation at the
 *	upper and lower limit bounds.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::InterpolateCrossSectionInTemperature( double t, double nm, double* crosssection) const
{

	const_iterator							iter;
	const_iterator							loweriter;
	const_iterator							upperiter;
	double									y;
	double									y1 = 0.0;
	double									t1;
	double									y2 = 0.0;
	double									t2;
	double									dt;
	double									thistemp;
	double									dlower = 999999.0;
	double									dl;
	bool									foundlower = false;
	double									dupper = 999999.0;
	double									du;
	bool									foundupper = false;
	bool									ok;
	bool									valid;

	dlower     = 999999.0;
	dupper     = 999999.0;
	foundlower = false;
	foundupper = false;
	
	if (  m_temperature_entries.size() == 1)
	{
		ok = m_temperature_entries.front().GetCrossSection( nm, crosssection );
	}
	else
	{


		for ( iter = m_temperature_entries.begin(); !(iter ==  m_temperature_entries.end()); ++iter)
		{
			valid =  (nm >= (*iter).MinimumWavelengthEntry())  && ( nm <= (*iter).MaximumWavelengthEntry());	// Does this entry support this wavelength
			if (valid)
			{
				thistemp = (*iter).GetTemperature();
				dl       = t - thistemp;								// Get the distance of this temp from the desired temp
				if ((dl >= 0.0) && ( dl < dlower))						// If we are below the desired temp and closer  than current
				{														// then
					dlower = dl;										// save the new closest entry
					loweriter = iter;									// save the index of the "best" entry
					foundlower = true;									// and flag that we have found an entry
				}

				du  = thistemp - t;										// Get the distance of this temp from the desired temp
				if ((du >= 0.0) && ( du < dupper))						// If we are below the desired temp and closer  than current
				{														// then
					dupper     = du;									// save the new closest entry
					upperiter  = iter;									// save the index of the "best" entry
					foundupper = true;									// and flag that we have found an entry
				}
			}
		}
		ok = (foundupper || foundlower);								// Di we find at least one valid entry.
		if (ok)
		{
			if (!foundupper) upperiter = loweriter;
			if (!foundlower) loweriter = upperiter;

			ok = ok && (*loweriter).GetCrossSection( nm, &y1 );
			ok = ok && (*upperiter).GetCrossSection( nm, &y2 );
			t1 = (*loweriter).GetTemperature();
			t2 = (*upperiter).GetTemperature();	   
			dt = t2-t1;
			if (dt == 0) y = y1; 
			else         y = y1 + (y2-y1)*(t-t1)/dt;
			*crosssection = y;
		}
		if (!ok)
		{
			if (m_quietwavelengthtruncation)
			{
				*crosssection = 0;
				ok = true;
			}
			else
			{
				*crosssection = std::numeric_limits<double>::quiet_NaN();
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::AirToVacuumCorrection		2012-7-26*/
/** **/
/*---------------------------------------------------------------------------*/

double skOpticalProperties_UserDefinedAbsorption::AirToVacuumCorrection( double  wavenum_inairatSTP ) const
{
	return m_refractiveindexair.AirWavenumberToVacuum( wavenum_inairatSTP );
}


/*---------------------------------------------------------------------------
 *'					skOpticalProperties_UserDefinedAbsorption::CalculateCrossSections		2003-11-28
 *	This is the method called by the base class skOpticalProperties when it
 *	decided the current cross-section data are dirty and need to be refereshed.
 *	The cross-section data are dependent upon wavenumber and temperature.
 *
 *	In the first instance find out which two temperature ranges straddle the
 *	rwquired temperature.  Calculate the cross-section at the wavenumber at
 *	both of the two straddling temperatures and linearly interpolate to the current
 *	temperature.
 *-------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs)
{
	bool									ok;
	double									y;
	double									nm;
	double									t  = Temperature();

	if (m_wavelengthinvacuum) wavenum = AirToVacuumCorrection( wavenum );
	nm = 1.0E7/wavenum;
	ok = InterpolateCrossSectionInTemperature( t, nm, &y);
	if (!ok)
	{
		y = 0;
		ok = m_quietwavelengthtruncation;
		if (!ok)
		{
			nxLog::Verbose( NXLOG_WARNING, "skOpticalProperties_UserDefinedAbsorption::CalculateCrossSections, There were no cross-section entries for any range of temperature at the requested wavelength (%f nm). Settinng cross-section to 0.0", (double)nm);
		}
	}
	if (m_isabsorber)
	{
		*absxs = y;
		*extxs = y;
		*scattxs = 0.0;
	}
	else
	{
		*absxs = 0.0;
		*extxs = y;
		*scattxs = y;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::SetAtmosphericState		2008-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::SetAtmosphericState( skClimatology* neutralatmosphere)
{
	if (neutralatmosphere      != NULL) neutralatmosphere->AddRef();
	if (m_backgroundatmosphere != NULL) m_backgroundatmosphere->Release();
	m_backgroundatmosphere = neutralatmosphere;
	return true;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_UserDefinedAbsorption::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_UserDefinedAbsorption::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )
{
	double		kelvin= 0.0;
	bool		ok;
	
	ok = (m_backgroundatmosphere != nullptr);
	ok  = ok && m_backgroundatmosphere->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, pt, &kelvin, false );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_UserDefinedAbsorption::SetLocation, There was an error fetching the temperature. Have you called SetAtmopshericState properly?");
	}
	else
	{
		Set_Temperature(kelvin);
	}
	if (crosssectionschanged != NULL) *crosssectionschanged = true;					// Return flag stating whether cross-sections have changed
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::skOpticalProperties_TabulatedExtinction_HeightWavelength		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedExtinction_HeightWavelength::skOpticalProperties_TabulatedExtinction_HeightWavelength()
{
	m_isabsorber = false;
	m_idxh1      = 0;
	m_idxh2      = 0;
	m_h1         = 0.0;
	m_h2         = 0.0;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::LookupIndicesAndWeights		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::LookupIndicesAndWeights( const nx1dArray<double>& h, double value, double* w1, size_t* idx1, double* w2, size_t* idx2 ) const
{
	nxArrayIter<double>	start;
	nxArrayIter<double>	finish;
	nxArrayIter<double>	iter2;
	bool				ok;
	double				h1;
	double				h2;
	double				f;


	ok     = (h.size() > 0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedExtinction_HeightWavelength::LookupIndicesAndWeights, Cannot lookup up height indices as the height array is empty");
		*idx1 = 0;
		*w1   = 0.0;
		*idx2 = 0;
		*w2   = 0;
	}
	else
	{
		start  = h.begin();
		finish = h.end();
		iter2 = std::upper_bound( start, finish, value);
			
		if (iter2 == start )										// If we are befor ethe first element
		{															// then truncate to the first element
			*idx1 = 0;
			*w1   = 1.0;
			*idx2 = 0;
			*w2   = 0;
		}
		else if (iter2 == finish)									// if we are at the end of the array
		{															// then truncate to the last element
			*idx1 = h.size() -1 ;
			*w1   = 1.0;
			*idx2 = 0;
			*w2   = 0.0;
		}
		else
		{
			*idx2   = (iter2 - start);
			*idx1   = *idx2 - 1;
		
			h2      = h.At(*idx2); 
			h1      = h.At(*idx1);
			f       = (value-h1)/(h2-h1);
			*w1     = (1-f);
			*w2     = f;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::SetAtmosphericState		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::SetAtmosphericState( skClimatology* /* neutralatmosphere*/ )
{
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::SetAtmosphericState		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )
{
	bool	ok;

	ok = LookupIndicesAndWeights( m_heights,pt.heightm, &m_h1, &m_idxh1, &m_h2, &m_idxh2 );
	if (crosssectionschanged != NULL) *crosssectionschanged = true;					// Return flag stating whether cross-sections have changed
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs) const
{
	double		wavelen;
	double		w1;
	double		w2;
	size_t		idx1;
	size_t		idx2;
	double		exth1;
	double		exth2;
	double		ext;
	bool		ok;

	wavelen = 1.0E7/wavenumber;
	ok      = LookupIndicesAndWeights( m_wavelennm, wavelen, &w1, &idx1, &w2, &idx2 );
	exth1   = w1*m_extinction.At(idx1, m_idxh1) + w2*m_extinction.At(idx2,m_idxh1);			// Get the extinction at required wavelength at lower altitude
	exth2   = w1*m_extinction.At(idx1, m_idxh2) + w2*m_extinction.At(idx2,m_idxh2);			// Get the extinction at required wavelength at upper altitude
	ext     = m_h1*exth1 + m_h2*exth2;
	if (m_isabsorber)
	{
		*absxs   = ext;
		*extxs   = ext;
		*scattxs = 0.0;
	}
	else
	{
		*absxs   = 0.0;
		*extxs   = ext;
		*scattxs = ext;
	}
	return ok;
}

/*----------------------------------  -------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::CalculateCrossSections		2009-5-19*/
/** [THREAD-SAFE]. Calculates the cross-section at the specified wavenumber. Code
 *	is instrinsically thread-safe as it only uses const methods 
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	return CalculateCrossSections( wavenumber, absxs, extxs, scattxs);
}


/*-----------------------------------------------------------------------------
 *					ReleaseResources		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedExtinction_HeightWavelength:: ReleaseResources()
{
	m_idxh1      = 0;
	m_idxh2      = 0;
	m_h1         = 0.0;
	m_h2         = 0.0;
	m_heights.erase();
	m_wavelennm.erase();
	m_extinction.erase();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::CreateClone		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_TabulatedExtinction_HeightWavelength::CreateClone( skOpticalProperties ** userclone ) const
{
	bool													ok;
	skOpticalProperties_TabulatedExtinction_HeightWavelength*	clone;

	clone = new skOpticalProperties_TabulatedExtinction_HeightWavelength;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedExtinction_HeightWavelength::CreateClone, Error creating allocating memory for clone");
	}
	else
	{
		clone->AddRef();
		clone->m_isdirty    = true;
		clone->m_isabsorber = m_isabsorber;
		clone->m_idxh1      = m_idxh1;
		clone->m_idxh2      = m_idxh2;
		clone->m_h1         = m_h1;
		clone->m_h2         = m_h2;
		ok =       clone->skOpticalProperties::DeepCopy( *this );
		ok = ok && clone->m_heights.DeepCopy	( m_heights);
		ok = ok && clone->m_wavelennm.DeepCopy	( m_wavelennm );
		ok = ok && clone->m_extinction.DeepCopy	( m_extinction );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedExtinction_HeightWavelength::CreateClone, Error making clone, I'm going to erase the clone to blanks");
		clone->ReleaseResources();
	}
	*userclone = clone;
	return ok;
};
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::SetExtinctionTable		2009-5-19*/
/** Sets up the table that will be used for extinction (or absorption) as a
*	function of altitude and wavelength. The user passes in a 2-D array of extinction
*	values and two 1D arrays that define the wavelength and height "axes" of the grid.
*
*	\param extinction
*	A 2D array of atmopsheric extinction (or absorption). The units are extinction per cm.
*	The array is double( numwavelens, numheights ) where numwavelens is the number of
*	wavelengths in parameter wavelens and numheights is the number of heights in parameter
*	heights_meets. The code linearly interpolates the table in both height and wavelength.
*	It truncates the interpolation at the end points of the array.
*
*	\param wavelens
*	The wavelengths at which the table is specified. The units are nm. Must be 1 or more. All elements must
*	be provided in ascending order.
*
*	\param heights_meters
*	The heights at which the table is specified. The units are in meters. Muts be 1 or more values. All elements
*	must be provided in ascending order.
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::SetExtinctionTable( const nx2dArray<double>& extinction, const nx1dArray<double>& wavelens, const nx1dArray<double>& heights_meters)
{
	const size_t*		dims;
	bool				ok;

	dims = extinction.ArrayRankSpecs()->Dims();
	ok   =	   ( dims[0] == wavelens.size()) 
			&& ( dims[1] == heights_meters.size())
			&& ( dims[0] > 0 )
			&& ( dims[1] > 0 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedExtinction_HeightWavelength::SetExtinctionTable, Error setting table as array sizes are incompatible with each other");
		ReleaseResources();
	}
	else
	{
		ok =       m_heights.DeepCopy	( heights_meters);
		ok = ok && m_wavelennm.DeepCopy	( wavelens );
		ok = ok && m_extinction.DeepCopy( extinction );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedExtinction_HeightWavelength::SetExtinctionTable, Error copying arrays over. I'm clearing array to zero");
		}
	}
	if (!ok)
	{
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedExtinction_HeightWavelength::LoadHeightWavelengthProfileFromFile		2009-6-26*/
/** Loads in the tabulated extinction from a file.
 *
 *	Line 1:					numwave	numheights			number of wavelengths and number of heights
 *	Line 2:					wavelengths(numwavelen)		expressed in nanometers.
 *	Line 3:					heights(numheights)			expressed in meters.
 *	Line 4:					extinction(numwave)			extinction for wavelengths at height (0)
 *	Line 5:					extinction(numwave)			extinction for wavelengths at height (1)
 *	..
 *	Line (numheights+3):	extinction(numwave)			extinction for wavelengths at height (numheights-1)
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedExtinction_HeightWavelength::LoadHeightWavelengthProfileFromFile( const char* filename )
{
	std::ifstream			strm;
	size_t					numh;
	size_t					numwave;
	nx1dArray<double>		wavelens;
	nx1dArray<double>		heights;
	nx2dArray<double>		extinction;
	bool					ok;

	strm.open(filename, std::ios_base::in );
	strm  >> numwave >> numh;							// Load in the number of heights and the number of wavelengths
	ok = (numh > 0) && (numwave > 0) && !strm.fail();	// Make sure it is good
	ok = ok && wavelens.SetSize(numwave);
	ok = ok && heights.SetSize (numh );
	ok = ok && extinction.SetSize( numwave, numh );
	if (ok)
	{
		strm >> wavelens; 
		strm >> heights;
		strm >> extinction;
		ok = ok && !strm.fail();
	}
	strm.close();
	ok = ok && SetExtinctionTable( extinction, wavelens, heights);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedExtinction_HeightWavelength::LoadHeightWavelengthProfileFromFile, error loading and setting the height extinction table from file <%s>", (const char*) filename);
		ReleaseResources();
	}
	return ok;
}

