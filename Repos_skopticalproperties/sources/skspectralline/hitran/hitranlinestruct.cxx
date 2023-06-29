#include <skopticalproperties21.h>
#include <cctype> 
#include <omp.h>


/*---------------------------------------------------------------------------
 *                      Class argsort_index                       2020-09-15 */
/** **/
/*---------------------------------------------------------------------------*/

class argsort_index
{
	private:
		size_t						m_index;
		const std::vector<double>*	m_values;

	public:
				argsort_index( size_t index, const std::vector<double>* values )	{ m_index = index; m_values = values;}
		double	Value		() const												{ return m_values->at(m_index); }
		size_t	Index		() const												{ return m_index;}
		bool	operator <	(const argsort_index& other) const						{ return Value() < other.Value(); }
};


/*---------------------------------------------------------------------------
 *                            argsort                             2020-09-15 */
/** Sorts an array into ascending order and returns the indices of the
 *	sorted elements in the original array.
 **/
/*---------------------------------------------------------------------------*/

static bool argsort( const std::vector<double>& uservalues, std::vector<double>* sortedarray, std::vector<size_t >* sortedindices)
{
	std::vector<argsort_index>		values;
	size_t							N = uservalues.size();

	values.reserve( N );
	sortedindices->resize( N );
	sortedarray->resize(N);

	for (size_t i = 0; i < N; i++ )
	{
		values.push_back( argsort_index( i, &uservalues) );
	}
	std::sort( values.begin(), values.end() );
	for (size_t i = 0; i < N; i++ )
	{
		sortedindices->at(i) = values.at(i).Index();
		sortedarray->at(i)   = values.at(i).Value();
	}
	return true;
};


/*-----------------------------------------------------------------------------
 *					skSpectralLine_HitranLine		2013-3-5*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine_HitranLine::skSpectralLine_HitranLine()
{
	ClearRecord();
}


/*---------------------------------------------------------------------------
 *      skSpectralLine_HitranLine::skSpectralLine_HitranLine      2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine_HitranLine::skSpectralLine_HitranLine( const HitranLineStruct& record)
	                      : m_entry(record)
{
	m_moleculeinfo		= NULL;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLine_HitranLine::~skSpectralLine_HitranLine		2013-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine_HitranLine::~skSpectralLine_HitranLine()
{
}


/*-----------------------------------------------------------------------------
 *					skSpectralLine_HitranLine::skSpectralLine_HitranLine		2014-2-28*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine_HitranLine::skSpectralLine_HitranLine( const skSpectralLine_HitranLine& other)
{
	size_t i;

	m_moleculeinfo		= other.m_moleculeinfo;
	m_entry.m_molNum			= other.m_entry.m_molNum;
	m_entry.m_isoNum			= other.m_entry.m_isoNum;
	m_entry.m_nuTrans			= other.m_entry.m_nuTrans;
	m_entry.m_intensity			= other.m_entry.m_intensity;
	m_entry.m_einsteinA			= other.m_entry.m_einsteinA;
	m_entry.m_airBroad			= other.m_entry.m_airBroad;
	m_entry.m_selfBroad			= other.m_entry.m_selfBroad;
	m_entry.m_nuLow				= other.m_entry.m_nuLow;
	m_entry.m_nAir				= other.m_entry.m_nAir;
	m_entry.m_deltaAir			= other.m_entry.m_deltaAir;
	m_entry.m_upperStatWt		= other.m_entry.m_upperStatWt;
	m_entry.m_lowerStatWt		= other.m_entry.m_lowerStatWt;
	for ( i = 0; i < 16; i++)
	{
		m_entry.m_qGlobU[i]     = other.m_entry.m_qGlobU[i];
		m_entry.m_qGlobL[i] 	= other.m_entry.m_qGlobL[i];
		m_entry.m_qLocU [i] 	= other.m_entry.m_qLocU [i];
		m_entry.m_qLocL [i]		= other.m_entry.m_qLocL [i];
	}
	m_entry.m_qGlobU[16]    = '\0';
	m_entry.m_qGlobL[16] 	= '\0';
	m_entry.m_qLocU [16] 	= '\0';
	m_entry.m_qLocL [16]		= '\0';

	for ( i = 0; i < 7;  i++) m_entry.m_errCodes [i]  = other.m_entry.m_errCodes[i];
	for ( i = 0; i < 13; i++) m_entry.m_refCodes [i]  = other.m_entry.m_refCodes[i];
	for ( i = 0; i < 2;  i++) m_entry.m_lineFlag [i]  = other.m_entry.m_lineFlag[i];
	m_entry.m_errCodes [7] = '\0';
	m_entry.m_refCodes [13] ='\0';
	m_entry.m_lineFlag [2]  = '\0';
}

/*-----------------------------------------------------------------------------
 *					skSpectralLine_HitranLine::ClearRecord		2013-3-5*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLine_HitranLine::ClearRecord()
{
	m_moleculeinfo		= NULL;
	m_entry.m_molNum			= -1;
	m_entry.m_isoNum			= -1;
	m_entry.m_nuTrans			= -9999.0;
	m_entry.m_intensity			= -9999.0;
	m_entry.m_einsteinA			= -9999.0;
	m_entry.m_airBroad			= -9999.0;
	m_entry.m_selfBroad			= -9999.0;
	m_entry.m_nuLow				= -9999.0;
	m_entry.m_nAir				= -9999.0;
	m_entry.m_deltaAir			= -9999.0;

	m_entry.m_qGlobU[0]			= '\0';
	m_entry.m_qGlobL[0]			= '\0';
	m_entry.m_qLocU [0]			= '\0';
	m_entry.m_qLocL [0]			= '\0';
	m_entry.m_errCodes[0]			= '\0';
	m_entry.m_refCodes[0]			= '\0';
	m_entry.m_lineFlag[0]			= '\0';
	m_entry.m_upperStatWt			= -9999.0;
	m_entry.m_lowerStatWt			= -9999.0;
}

/*-----------------------------------------------------------------------------
 *			HitranLinesText::ParseAndWriteBin		2010-4-5*/
/** 
 *	Reads in each record, parses the data, and writes the record out to the binary file in struct format
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLine_HitranLine::Parse160CharRecord( const char* lineString, size_t	numchar)
{
	bool	ok;

	ok = numchar == 160;
	if (ok)
	{
		m_entry.m_molNum			= IntegerValFromString  ( lineString,   0,  2 );
		m_entry.m_isoNum			= ExtendedHexValFromChar( lineString,   2     );
		m_entry.m_nuTrans			= DoubleValFromString   ( lineString,   3, 12 );
		m_entry.m_intensity			= DoubleValFromString   ( lineString,  15, 10 );
		m_entry.m_einsteinA			= DoubleValFromString   ( lineString,  25, 10 );
		m_entry.m_airBroad			= DoubleValFromString   ( lineString,  35,  5 );
		m_entry.m_selfBroad			= DoubleValFromString   ( lineString,  40,  5 );
		m_entry.m_nuLow				= DoubleValFromString   ( lineString,  45, 10 );
		m_entry.m_nAir				= DoubleValFromString   ( lineString,  55,  4 );
		m_entry.m_deltaAir			= DoubleValFromString   ( lineString,  59,  8 );

		SubstringFromString ( m_entry.m_qGlobU,   lineString,  67, 15 ) ;
		SubstringFromString ( m_entry.m_qGlobL,   lineString,  82, 15 ) ;
		SubstringFromString ( m_entry.m_qLocU,    lineString,  97, 15 ) ;
		SubstringFromString ( m_entry.m_qLocL,    lineString, 112, 15 ) ;
		SubstringFromString ( m_entry.m_errCodes, lineString, 127,  6 ) ;		
		SubstringFromString ( m_entry.m_refCodes, lineString, 133, 12 ) ;
		SubstringFromString ( m_entry.m_lineFlag, lineString, 145,  1 );

		m_entry.m_upperStatWt		= DoubleValFromString ( lineString, 146,  7 );
		m_entry.m_lowerStatWt		= DoubleValFromString ( lineString, 153,  7 );
	}
	else
	{
		ClearRecord();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *			HitranLinesText::IntegerValFromString		2010-4-5*/
/** 
 *	Must supply c-string, first index of substring. and length of substring
 **/
/*---------------------------------------------------------------------------*/

int	skSpectralLine_HitranLine::IntegerValFromString( const char* str, int firstIdx,  int subStrLen )
{
	char	substring[30];
	int		val;
	int		i;
	bool	ok = true;

	str = str + firstIdx;
	for(i=0; i <subStrLen; i++)
	{
		substring[i] = str[i];
		ok = ok && ( ((str[i] >= '0' ) && (str[i] <= '9')) || (str[i] == ' '));
	}
	substring[subStrLen] = '\0';
	if (!ok) nxLog::Record(NXLOG_WARNING, "skSpectralLine_HitranLine::IntegerValFromString, non-numeric values found in integer string [%s]", (const char*)substring);
	val = atoi(substring);
	return val;
}


/*-----------------------------------------------------------------------------
 *			HitranLinesText::IntegerValFromString		2010-4-5*/
/** 
 *	Must supply c-string, first index of substring. and length of substring
 **/
/*---------------------------------------------------------------------------*/

int	skSpectralLine_HitranLine::ExtendedHexValFromChar( const char* str, int firstIdx )
{
	char	c;
	int		val;
	bool	ok = true;

	c = str[firstIdx];
	if (( c >= '0') && ( c <= '9'))
	{
		val = c - '0';
		if (val == 0) val = 10;
	}
	else 
	{
		c = toupper(c);
		ok = (c >= 'A') && ( c <= 'Z');
		val = ok ? (c - 'A') + 11: -1;
	}
	if (!ok) nxLog::Record(NXLOG_WARNING, "skSpectralLine_HitranLine::ExtendedHexValFromChar, non-hex character values found in character [%c]", (char)c);
	return val;
}

/*-----------------------------------------------------------------------------
 *			HitranLinesText::DoubleValFromString		2010-4-5*/
/** 
 *	Must supply c-string, first index of substring. and length of substring
 **/
/*---------------------------------------------------------------------------*/

double skSpectralLine_HitranLine::DoubleValFromString( const char* str, int firstIdx,  int subStrLen )
{
	double	 val;
	char	 substring[50];
	int		 i;

	str = str + firstIdx;
	for(i=0; i <subStrLen; i++) substring[i] = str[i];
	substring[subStrLen] = '\0';
	val = atof(substring);
	return val;
}


/*-----------------------------------------------------------------------------
 *			HitranLinesText::SubstringFromString		2010-4-5*/
/** 
 *	Must supply c-string, first index of substring. and length of substring
 **/
/*---------------------------------------------------------------------------*/

void skSpectralLine_HitranLine::SubstringFromString( char* outputstr, const char* inputstr, int firstIdx,  int subStrLen )
{
	int i;

	inputstr = inputstr + firstIdx;
	for( i=0; i < subStrLen; i++) *outputstr++ = *inputstr++;
	*outputstr = '\0';
}


size_t skSpectralLineCollection_HitranIsotope::g_numinstances = 0;

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranIsotope::skSpectralLineCollection_HitranIsotope		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection_HitranIsotope::skSpectralLineCollection_HitranIsotope( const skSpectralLineCollection_HitranChemical* parentchemical, size_t molnum, size_t isotopeid, const skHitranPartitionTableEntry*	moleculeinfo)
	                                   : skSpectralLineCollection( parentchemical->MicroWindowMinWavenum(), parentchemical->MicroWindowMaxWavenum() ), 
	                                     m_partitioncache(moleculeinfo)
{
	bool ok;
	g_numinstances++;
	m_moleculeinfo = moleculeinfo;
	m_molNum       = molnum;
	m_isoNum       = isotopeid;

	m_parentchemical = parentchemical;
	ok = (m_molNum == moleculeinfo->MoleculeNumber()) && ( m_isoNum == moleculeinfo->IsotopeNumber());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranIsotope::skSpectralLineCollection_HitranIsotope, Molecule number and isotope numbers are not consistent. Thats a problem that needs debugging");
	}
}


/*---------------------------------------------------------------------------
 * skSpectralLineCollection_HitranIsotope::skSpectralLineCollection_HitranIsotope 2020-07-29 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection_HitranIsotope::skSpectralLineCollection_HitranIsotope( const skSpectralLineCollection_HitranIsotope& other)
	                                   :skSpectralLineCollection( other.MicroWindow_MinWavenum(), other.MicroWindow_MaxWavenum() ),
	                                    m_partitioncache(other.m_moleculeinfo)
{
	g_numinstances++;
	m_moleculeinfo	 = other.m_moleculeinfo;
	m_molNum		 = other.m_molNum;
	m_isoNum		 = other.m_isoNum;
	m_parentchemical = other.m_parentchemical;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranIsotope::~skSpectralLineCollection_HitranIsotope		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection_HitranIsotope::~skSpectralLineCollection_HitranIsotope()
{
	g_numinstances--;
}



/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranIsotope::QPartition		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineCollection_HitranIsotope::QPartition ( double T ) const
{
	return m_partitioncache.InternalPartition( T );
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranIsotope::MassAMU		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineCollection_HitranIsotope::MassAMU	( ) const
{
	return m_moleculeinfo != NULL ? m_moleculeinfo->m_molarmass : 0.0;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranIsotope::PartialPressure		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineCollection_HitranIsotope::PartialPressure( const GEODETIC_INSTANT& geopt, double airpressure, double airtemp) const
{
	double pp = 0.0;

	if (m_parentchemical != NULL)
	{
		pp = m_parentchemical->PartialPressure( geopt, airpressure, airtemp);
	}
	return pp;
}



/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::skSpectralLineCollection_HitranChemical		2013-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

//skSpectralLineCollection_HitranChemical::skSpectralLineCollection_HitranChemical()
//{
//	m_molnum = 0;
//	m_manager = skHitranMoleculeManager::CreateManagerInstance();
//}
//
size_t skSpectralLineCollection_HitranChemical::g_numinstances = 0;

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::skSpectralLineCollection_HitranChemical		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection_HitranChemical::skSpectralLineCollection_HitranChemical( const char* chemicalname, double lowerwavenumber, double upperwavenumber, double windowmargin_wavenumber, bool hapicompliant, int isotopefilter, const char* lowerstate_globalquanta_filter, const char* upperstate_globalquanta_filter)
{
	bool	ok;
	g_numinstances++;

	m_hapicompliant = hapicompliant;
	m_userdefined_maxlinestrength = std::numeric_limits<double>::quiet_NaN();	// A User defined maximum line strength
	m_use_userdefined_maxlinestrength = false;
	m_numberdensityguid        = SKCLIMATOLOGY_UNDEFINED;		// The number desnity guid for retrieving number density of this chmical species from the internal climatology for partial pressure
	m_selfbroadeningclimatology = NULL;
	m_maxlinestrength = 0.0;
	m_molnum          = 0;
	m_isotopeidfilter = isotopefilter;
	m_upperstate_globalquanta_filter = GlobalStateToString(upperstate_globalquanta_filter);	// Filter used to select only lines that belong to the given upper state global quanta
	m_lowerstate_globalquanta_filter = GlobalStateToString(lowerstate_globalquanta_filter);	// Filter used to select only lines that belong to the given lower state global quanta

	m_manager         = skHitranMoleculeManager::CreateManagerInstance(m_hapicompliant);
	m_lowerwavenumber = lowerwavenumber;
	m_upperwavenumber = upperwavenumber;
	m_windowmargin    = windowmargin_wavenumber;
	ok  = (lowerwavenumber < upperwavenumber);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::Constructor, lower wavenumber is bigger than upper wavenumber. No lines will be selected. Check numbers.");
	}
	ok  =       SetChemicalName(chemicalname);
	ok  = ok && LoadFile();

}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::~skSpectralLineCollection_HitranChemical		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection_HitranChemical::~skSpectralLineCollection_HitranChemical()
{
	m_isotopes.clear();
	if (m_manager != NULL) m_manager->Release();
	if (m_selfbroadeningclimatology != NULL) m_selfbroadeningclimatology->Release();
	g_numinstances--;

}



/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::ReleaseResources		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineCollection_HitranChemical::ReleaseResources()
{
	m_isotopes.clear();
	m_molnum = 0;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetAtmosphericState		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::UpdateLocation( const GEODETIC_INSTANT& geopt, skClimatology* atmosphere )
{
	iterator	iter;
	bool		ok = true;
	bool		ok1;
	double		temperature;
	double		pressure = 0.0;

	m_maxlinestrength = 0.0;			// The maximum Snm line strength in the window	

	ok =       atmosphere->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, geopt, &temperature, false );
	ok = ok && atmosphere->GetParameter( SKCLIMATOLOGY_PRESSURE_PA,   geopt, &pressure,    false );

	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1 = (*iter).second.UpdateLocation( temperature, pressure, geopt, atmosphere );
		ok = ok && ok1;
		if (ok1) m_maxlinestrength = nxmax( m_maxlinestrength, (*iter).second.MaxLineStrength() );	// Keep track of the maximum line strength loaded in.
	}
	ok = ok && SetLineLimitsfromMaxLineStrength();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetUserDefinedMaxLineStrength		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SetUserDefinedMaxLineStrength( double maxlinestrength)
{
	bool ok;

	ok =  (maxlinestrength >= 0.0);
	if (ok)
	{
		m_use_userdefined_maxlinestrength = (maxlinestrength > 0.0);
		m_userdefined_maxlinestrength     = maxlinestrength;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::SetUserDefinedMaxLineStrength, The maximum line strength (%e) must be 0.0 or positive. Nothing has been changed", (double)maxlinestrength );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetLineTolerance		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/


bool skSpectralLineCollection_HitranChemical::SetLineTolerance( double tolerance )
{
	bool		ok  = true;
	bool		ok1;
	iterator	iter;

	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1 = (*iter).second.SetLineTolerance( tolerance);
		ok = ok && ok1;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetLineLimitsfromMaxLineStrength		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SetLineLimitsfromMaxLineStrength()
{
	bool		ok  = true;
	bool		ok1;
	iterator	iter;
	double		maxstrength = m_use_userdefined_maxlinestrength ? m_userdefined_maxlinestrength : m_maxlinestrength;

	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1 = (*iter).second.SetLineLimitsfromMaxLineStrength( maxstrength );
		ok = ok && ok1;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::AbsorptionCrossSection		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

 bool skSpectralLineCollection_HitranChemical::AbsorptionCrossSection( double nu, double* absxsec) const
{
	const_iterator	iter;
	bool			ok1;
	bool			ok = true;
	double			xsec;
	
	*absxsec = 0.0;
	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1       = (*iter).second.AbsorptionCrossSectionOrEmission( nu, &xsec);
		*absxsec += xsec;
		ok        = ok && ok1;

	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::AddAbsorptionCrossSection		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

 bool skSpectralLineCollection_HitranChemical::AddAbsorptionCrossSectionArray( const std::vector<double>& userwavenum, std::vector<double>* userabsxs)
 {
	iterator	iter;
	bool			ok1;
	bool			ok = true;
	std::vector<double> sortedwavenum;
	std::vector<double> absxs;
	std::vector<size_t> indices;

	argsort( userwavenum, &sortedwavenum, &indices );													// Sort the wavenumnbers into ascending order, get the indices of the sorted elements into the original array
	absxs.resize(userabsxs->size());
	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1       = (*iter).second.AddAbsorptionCrossSectionOrEmissionArray( sortedwavenum, &absxs);	// Do the calculation with ascending wavenumber
		ok        = ok && ok1;
	}
	
	userabsxs->resize( indices.size());																	// Now copy the sorted wavenumbers backl to the original order
	for (size_t i = 0; i < indices.size(); i++)
	{
		userabsxs->at( indices.at(i) ) = absxs.at(i);
	}

	return ok;
 }


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetNumThreads		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

 bool skSpectralLineCollection_HitranChemical::SetNumThreads( size_t numthreads)
 {
	iterator	iter;
	bool			ok1;
	bool			ok = true;
	
	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1       = (*iter).second.SetNumThreads(numthreads);
		ok        = ok && ok1;
	}
	return ok;
 }


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetChemicalName		2013-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SetChemicalName( const char* chemicalname )
{
	bool	ok;
	size_t	molid;

	ReleaseResources();
	ok = m_manager->FindMoleculeId( chemicalname, &molid);
	if (ok)
	{
		m_molnum = (int)molid;
		m_chemicalname = chemicalname;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical, Error loading Hitran database for <%s>", (const char*)chemicalname);
		ReleaseResources();
	}
	return ok;

}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetLineShapeObject		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SetLineShapeObject	( skSpectralLineShape* lineshapeobject )
{
	iterator	iter;
	bool		ok = true;
	bool		ok1;

	for (iter =m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	{
		ok1 = (*iter).second.SetLineShapeObject(lineshapeobject);
		ok = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetIsotopeFilter		 2014- 10- 15*/
/** Used to load only one isotope.**/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SetIsotopeFilter( size_t filterid)
{
	bool	ok;

	ok = (m_isotopeidfilter == filterid);
	if (!ok)
	{
		ok = (m_isotopes.size() == 0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::SetIsotopeFilter, You cannot change the isotope selection filter once HITRAN records have been loaded, Call this function during intialization");
		}
		else
		{
			m_isotopeidfilter = filterid;
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *  skSpectralLineCollection_HitranChemical::MatchesGlobalQuanta  2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::MatchesGlobalQuanta( const HitranLineStruct* spectralline ) const
{
	bool matches;

	matches =  (m_upperstate_globalquanta_filter.size() == 0) && (m_lowerstate_globalquanta_filter.size() == 0 );
	if (!matches )
	{
		matches =     ( (m_upperstate_globalquanta_filter.size() == 0) || m_upperstate_globalquanta_filter == GlobalStateToString( spectralline->m_qGlobU ))
			      &&  ( (m_lowerstate_globalquanta_filter.size() == 0) || m_lowerstate_globalquanta_filter == GlobalStateToString( spectralline->m_qGlobL ));
	}
	return matches;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::InsertSpectralLineEntry		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::InsertSpectralLineEntryNoFilter( size_t isotopeid, skSpectralLine_HitranLine* spectralline)
{
	size_t									key;
	iterator								iter;
	skSpectralLineEntry*					lineentry;
	skSpectralLineCollection_HitranIsotope*	isotope = NULL;
	const skHitranPartitionTableEntry*		moleculeinfo;
	bool									ok;
	
	key = isotopeid;						// Get the key for the isotopes
	iter = m_isotopes.find(key);			// The array of isotopes associated with this chemical
	if (iter==m_isotopes.end())
	{
		m_manager->FindMoleculeEntry( m_molnum, isotopeid, &moleculeinfo );
		m_isotopes.insert( value_type( key, skSpectralLineCollection_HitranIsotope(this, m_molnum, isotopeid,moleculeinfo )) );
		iter = m_isotopes.find(key);
		if ( iter == m_isotopes.end() )
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::FindIsotopeEntry, there were erors creating a new entry for isotope %d", (int)isotopeid);
		}
		else
		{
			isotope = &(*iter).second;
		}
	}
	else
	{
		isotope = &(*iter).second;
	}
	ok =       (isotope != NULL);
	ok = ok && spectralline->SetParentMolecule( isotope );
	ok = ok && spectralline->SetMoleculeInfo  ( isotope->MoleculeInfo() );
	if (ok)
	{
		lineentry = new skSpectralLineEntry;
		ok = ok && lineentry->SetSpectralLine     ( spectralline );
		ok = ok && isotope->AddEntry(lineentry);
	}
	
	return ok;
}


/*---------------------------------------------------------------------------
 * skSpectralLineCollection_HitranChemical::InsertAllSpectralLines   2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::InsertAllSpectralLines( size_t isotopeid, const std::vector<HitranLineStruct>& spectrallines)
{
	size_t									key;
	iterator								iter;
	skSpectralLineEntry*					lineentry;
	skSpectralLineCollection_HitranIsotope*	isotope = NULL;
	const skHitranPartitionTableEntry*		moleculeinfo;
	bool									ok, ok1;
	bool									addthisisotope;
	skSpectralLine_HitranLine*				spectralline;

	
	addthisisotope = ((m_isotopeidfilter  == 0) || ( m_isotopeidfilter == isotopeid)) ;
	ok = !addthisisotope;
	if (!ok)
	{
		key = isotopeid;						// Get the key for the isotopes
		iter = m_isotopes.find(key);			// The array of isotopes associated with this chemical
		if (iter==m_isotopes.end())
		{
			m_manager->FindMoleculeEntry( m_molnum, isotopeid, &moleculeinfo );
			m_isotopes.insert( value_type( key, skSpectralLineCollection_HitranIsotope(this, m_molnum, isotopeid,moleculeinfo )) );
			iter = m_isotopes.find(key);
			if ( iter == m_isotopes.end() )
			{
				nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::InsertAllSpectralLines, there were errors creating a new isotope entry for isotope %d", (int)isotopeid);
			}
			else
			{
				isotope = &(*iter).second;
			}
		}
		else
		{
			isotope = &(*iter).second;
		}
		ok = (isotope != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::InsertAllSpectralLines, There were errors obtainging an isotope instance for molecule %d isotope %d", (int)m_molnum, (int)isotopeid);
		}
		else
		{
			isotope->ClearLines(spectrallines.size());
			for (size_t i = 0; i < spectrallines.size(); i++)
			{
				if ( MatchesGlobalQuanta(&spectrallines.at(i)))
				{

					spectralline = new skSpectralLine_HitranLine( spectrallines.at(i) );
					ok1 =        spectralline->SetParentMolecule( isotope );
					ok1 = ok1 && spectralline->SetMoleculeInfo  ( isotope->MoleculeInfo() );
					if (ok1)
					{
						lineentry = new skSpectralLineEntry; 
						ok1 = ok1 && lineentry->SetSpectralLine     ( spectralline );
						ok1 = ok1 && isotope->AddEntry(lineentry);
					}
					ok  = ok  && ok1;
				}
			}
			if (!ok) nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::InsertAllSpectralLines, There were errors inserting all the lines for molecule %d isotope %d", (int)m_molnum, (int)isotopeid);
		}
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::FindFile		2013-3-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::FindFile( nxString* filename)
{
	nxString	basename;
	nxString	moleculedir;
	bool		ok;

	if (m_hapicompliant)
	{
		basename.sprintf("%s.data", (const char*)m_chemicalname.c_str());
	}
	else
	{
		basename.sprintf("%02d_hit.par", (int)m_molnum);
	}
	ok = m_manager->FindHitranMoleculeDirectory( &moleculedir);
	if (ok)
	{
		moleculedir.EnsureLastCharIsDirectoryChar();
		nxDirectory	files;
		files.ScanDirectory(basename, false, moleculedir);

		ok = files.List().GetSize() == 1;
		if (ok)
		{
			filename->sprintf( "%s", (const char*)files.List().GetAt(0) );
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::FindFile, could not find any files matching <%s> in <%s>", (const char*)basename, (const char*)moleculedir);
		}
	}
	if (!ok)
	{
		*filename = "Not Found";
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::FindIsotopeId		 2014- 10- 15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::FindIsotopeId( size_t isotopeid, const std::vector<const skHitranPartitionTableEntry*>& isotopetable)
{
	bool ok;

	ok = (isotopeid == 0);
	if (!ok)
	{
		auto iter = isotopetable.begin();
		while (!ok && !(iter == isotopetable.end() ))
		{
			ok = ((*iter)->m_isotopeid == isotopeid);
			++iter;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::FindIsotopeId, could not find the requested isotope (%d) in the list of available isotopes", (int)isotopeid);
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *                     global_state_to_string                     2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/

std::string  skSpectralLineCollection_HitranChemical::GlobalStateToString( const char* globalstate) const
{
	std::string  answer;

	if ( globalstate != nullptr)
	{
		nxStringArray fields;
		int n = fields.Strtok(globalstate);
		for (int i=0; i < n; i++)
		{
			answer += (const char*)(fields.GetAt(i));
		}
	}
	return answer;
}
/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::LoadFile		2013-3-8*/
/** Loads in the hitran database for the moecule specified in member m_molnum.
 *	This code loads in the hitran "by Molecule" *.par file. The code was written
 *	and tested on the Hitran 2008 database and, if they keep the same file and
 *	directory structure, should work for later editions but who knows.
 *
 *	The code does do sanity checks but we dont go crazy, e.g. if you manually
 *	edit the HITRAN files into a crazy format it probably will crash.
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::LoadFile()
{
	std::ifstream										hitranTextData;
	nxString											filename;
	char												str[164];
	skSpectralLine_HitranLine*							spectralline;
	skSpectralLine_HitranLine							dummyline;
	bool												ok;
	bool												ok1;
	//bool												ok2;
	bool												ok3;
	size_t												numchar;
	size_t												isotopeorder;
	size_t												isotopeid;
	size_t												numbadlines;
	size_t												numgoodlines;
	std::vector<const skHitranPartitionTableEntry*>		isotopetable;
	HitranLineStructCache								iocache;


	numbadlines = 0;
	numgoodlines = 0;
	dummyline.SetStatic();
	ok =       m_manager->FetchAllIsotopeEntries( m_molnum, &isotopetable );				// Get the list of isotopes in ascending abundance order
	ok = ok && FindIsotopeId( m_isotopeidfilter, isotopetable);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::LoadFile, Errors loading isotope information from the HITRAN database files");
	}
	else
	{
		ok = iocache.LoadSpectralLines( this, isotopetable, m_lowerwavenumber-m_windowmargin, m_upperwavenumber+m_windowmargin );			// The upper wavenumber for selecting lines from the HITRAN/spectral line database
		if (!ok)
		{
			ok = FindFile( &filename );
			if (ok)
			{
				hitranTextData.open( filename, std::ifstream::in );
				int linecount = 0;
				ok = !hitranTextData.fail();
				if (ok)
				{
					while ( !hitranTextData.eof() ) //data storage loop; repeat to end of file
					{
						hitranTextData.getline( str, 162, '\n' );
						linecount++;
						numchar = strlen( str );
						while(str[numchar-1] == '\r' || str[numchar-1] == '\n')
						{
							// Line endings don't count for the 160 char filesize
							numchar--;
						}
						if (numchar > 0)
						{
							ok1 = (numchar == 160);
							if (!ok1)
							{
								ok1 = (numchar < 2);
							}
							else
							{
								ok1 = ok1 && dummyline.Parse160CharRecord( str, numchar );
								ok3 = dummyline.MoleculeNumber() == m_molnum; 
								if (!ok3) nxLog::Record(NXLOG_WARNING, " skSpectralLineCollection_HitranChemical::LoadFile, Error decoding line %d of [%s]. Molecule number is incorrect, expected %d but got %d", (int)linecount, (const char*)filename, (int)m_molnum, (int)dummyline.MoleculeNumber());
								ok1 = ok1 && ok3;
								if (ok1)
								{
									isotopeorder = dummyline.IsotopeOrderNumber();					// Get the isotope order from the spectral line. This is a 1-based index of the isotopes in isotopetable.
									--isotopeorder;														// We are not really interested in this 1 based index, we want the actual isotope id (eg numbers like 626)
									ok1    = (isotopeorder < isotopetable.size());						// Make sure we can index the array
									if (ok1)															// and if we can
									{																	// then
										spectralline = new skSpectralLine_HitranLine(dummyline);
										spectralline->AddRef();
										isotopeid = isotopetable.at(isotopeorder)->m_isotopeid;					// get the actual isotope id 
										ok1       = InsertSpectralLineEntryNoFilter( isotopeid, spectralline );	// and insert this spectral line into that isotope.
										spectralline->Release();
									}
									else
									{
										nxLog::Record(NXLOG_WARNING, " skSpectralLineCollection_HitranChemical::LoadFile, Error decoding line %d of [%s]", (int)linecount, (const char*)filename);
									}
								}
							}
							if (ok1) numgoodlines++;
							else     numbadlines++;
							ok = ok && ok1;
						}
					}
				}
				if (ok)
				{
					ok1 =        iocache.WriteSpectralLines(this);
					ok  = ok1 && iocache.LoadSpectralLines( this, isotopetable, m_lowerwavenumber-m_windowmargin, m_upperwavenumber+m_windowmargin );			// The upper wavenumber for selecting lines from the HITRAN/spectral line databas
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::LoadFile, Cannot create hitran line cache as there were error loading records from HITRAN HAPI database text files ");
				}
			}
		}
	}
	UpdateMaxLineStrength();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"bool skSpectralLineCollection_HitranChemical::LoadFile, There were errors loading the hitran database for molecule id %d. %d good lines, %d bad lines", (int)m_molnum, (int)numgoodlines, (int)numbadlines);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::UpdateMaxLineStrength		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::UpdateMaxLineStrength()
{
	iterator	iter;
	 m_maxlinestrength = 0.0;
	 for (iter = m_isotopes.begin(); !(iter == m_isotopes.end()); ++iter)
	 {
		 m_maxlinestrength = nxmax( m_maxlinestrength, (*iter).second.MaxLineStrength() );
	 }
	 return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::PartialPressure		2013-3-20*/
/** Calculate the partial pressure of the species of interest. This is only
 *	done if the user has previously supplied a climatology and specified the
 *	number density climatology
 **/
/*---------------------------------------------------------------------------*/

double skSpectralLineCollection_HitranChemical::PartialPressure( const GEODETIC_INSTANT& geopt, double /*airpressure*/, double airtemp) const
{
	double	partialpressure = 0.0;
	double	N = 0.0;
	bool	ok;

	ok = (m_selfbroadeningclimatology != NULL);
	if (ok)
	{
		/** THIS MAY NOT BE THREAD SAFE. The "const" on the function call applies to the numberdensityclimatology pointer, not to the object it calls e.g. GetParameter **/
		ok = ok && m_selfbroadeningclimatology->GetParameter( m_numberdensityguid, geopt, &N, false );// Get number density in molecules/CM3
		if (ok)
		{
			partialpressure = N*1.0E6*airtemp*nxSI::KBOLTZMAN;			// Get the partial pressure in Pascals. Convert number density to molecules/M3
		}
	}
	return partialpressure;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::SetInternalNumberDensityClimatology		2013-3-20*/
/** Sets the internal climatology used to calculate the number density of the
 *	HITRAN chemical. This is solely used for the calculation of partial pressure.
 *	If you are not concerned about the differences between self and air broadening
 *	then you can never call this function and partial pressure is always set to 0.0.
 *
 *	To work properly the climatology passed in must retrieve the number density
 *	(in molecules/cm3) of the chemical. The user is completely responsible for making
 *	sure the right climatology and parameter guid come in as there is no practical test 
 *	we can perform to verify integrity.
 *
 *	In this version we dont differentiate between the various isotopes of a chemical species. 
 *	All isotopes of the same chemical get the same partial pressure (ie the same self-broadening)
 *	I think this is correct but let me know if I am wrong.
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SetSelfBroadeningClimatology( const CLIMATOLOGY_HANDLE& parameterguid, skClimatology* numberdensityclimatology )
{
	bool	ok = true;

	if ( numberdensityclimatology   != NULL) numberdensityclimatology->AddRef();
	if ( m_selfbroadeningclimatology != NULL) m_selfbroadeningclimatology->Release();
	m_selfbroadeningclimatology = numberdensityclimatology;

	m_numberdensityguid	 = parameterguid;
	if (m_selfbroadeningclimatology != NULL)
	{
		ok = m_selfbroadeningclimatology->IsSupportedSpecies( parameterguid );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::SetInternalNumberDensityClimatology, The climatology object does not support the requested number density species");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection_HitranChemical::InternalClimatology_UpdateCache		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection_HitranChemical::SelfBroadeningClimatology_UpdateCache( const GEODETIC_INSTANT& geopt )
{
	bool	ok;

	ok = (m_selfbroadeningclimatology == NULL);
	if (!ok)
	{
		ok = m_selfbroadeningclimatology->UpdateCache( geopt );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralLineCollection_HitranChemical::InternalClimatology_UpdateCache, There were errors updating theinternal climatology");
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_HitranChemical2008::skOpticalProperties_HitranChemical2008 2019-05-24 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HitranChemical2008::skOpticalProperties_HitranChemical2008()
{
	SetHitran2008Compliant();
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_HitranChemical2008::skOpticalProperties_HitranChemical2008 2019-05-24 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HitranChemical2008::skOpticalProperties_HitranChemical2008(const char* chemicalname, double lowerwavenumber, double upperwavenumber)
{
	SetHitran2008Compliant();
	SetChemicalName(chemicalname);
	SetWavenumberRange(lowerwavenumber, upperwavenumber);
}

/*---------------------------------------------------------------------------
 * skOpticalProperties_HitranChemical2008::~skOpticalProperties_HitranChemical2008 2019-05-24 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HitranChemical2008::~skOpticalProperties_HitranChemical2008()
{
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::skOpticalProperties_HitranChemical		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HitranChemical::skOpticalProperties_HitranChemical()
{
	init();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::skOpticalProperties_HitranChemical		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HitranChemical::skOpticalProperties_HitranChemical	(const char* chemicalname )
{

	init();
	SetChemicalName( chemicalname);
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::init		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_HitranChemical::init()
{
	m_xs_optimizer              = nullptr;
	m_hapicompliant             = true;
	m_isdirty                   = true;
	m_lowwavenum                = 0.0;
	m_hihwavenum                = 9999999.0;
	m_manual_microwindow_isset  = false;
	m_margin_width_wavenum      = 10.0;					// By default add 10 wavenumbers to the edges of any micro-window
	m_hitranchemical            = NULL;
	m_selfbroadeningclimatology = NULL;
	m_chemicalnumberdensityguid = SKCLIMATOLOGY_UNDEFINED;
	m_lineshapeobject           = NULL;
	m_atmosphericstateclimatology      = NULL;
	m_manualmaxlinestrength     = std::numeric_limits<double>::quiet_NaN();
	m_manualtolerance           = std::numeric_limits<double>::quiet_NaN();
	m_isotopefilter             = 0;
	m_numthreads                = omp_get_num_procs();
	m_lastpoint.latitude = 0.0;
	m_lastpoint.longitude = 0.0;
	m_lastpoint.heightm   = 0.0;
	m_lastpoint.mjd       = 0.0;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::~skOpticalProperties_HitranChemical		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_HitranChemical::~skOpticalProperties_HitranChemical()
{
	if ( m_xs_optimizer        != nullptr) delete m_xs_optimizer;
	if ( m_selfbroadeningclimatology != NULL) m_selfbroadeningclimatology->Release();
	if ( m_lineshapeobject     != NULL) m_lineshapeobject->Release();
	if (m_hitranchemical       != NULL) delete m_hitranchemical;
	if ( m_atmosphericstateclimatology != nullptr) m_atmosphericstateclimatology->Release();

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetDirty		2013-3-20*/
/** Flag the obejct as dirty, ie the user has changed some configuration info.
 *	Configuration is normally done during the early stages of the object life. If
 *	for some reasona the user trie sto change the configuration after all of the 
 *	internal objects are created we issue a warning saying that that is not the
 *	right way to do things as we get large computational hits. 
**/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_HitranChemical::SetDirty()
{
	m_isdirty = true;
	if ( m_hitranchemical != NULL)
	{
		delete m_hitranchemical;
		m_hitranchemical = NULL;
		nxLog::Record(NXLOG_INFO,"skOpticalProperties_HitranChemical::SetDirty, It is very inefficient to reset parameters of an skOpticalProperties_HitranChemical object after the HITRAN object is created. It makes more sense to create a brand new instance");
	}
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::CheckDirtyAndUpdate		2013-3-20*/
/**	We use CheckDirtyandUpdate to delay creation of the Hitran object until
 *	the user actually needs to do something with the obejct. This way we avoid
 *	quite large computational overheads if the object is created but never
 *	used.
 *
 *	Typically the underlying Hitran objects are not created until the user
 *	makes a call to SetAtmopshericState or InternalClimatology_UpdateCache, both of which
 *	must be called before a call to CalculateCrossSections.
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::CheckDirtyAndUpdate()
{
	bool	ok;

	ok = !m_isdirty;
	if (!ok)
	{
		NXASSERT( m_hitranchemical == NULL );
		ok = !m_chemicalname.IsEmpty();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::CheckDirtyAndUpdate, No hitran chemical species is defined. Try calling SetChemicalName()");
		}
		else
		{
			if (m_lineshapeobject == NULL)											// IF we have no line shape object passed in by the suer
			{																		// then
				m_lineshapeobject = new skSpectralLineShape_VoigtKuntz;			// Use the tabulated voigt as the default
				m_lineshapeobject->AddRef();
			}
			m_hitranchemical = new skSpectralLineCollection_HitranChemical( m_chemicalname, m_lowwavenum, m_hihwavenum, m_margin_width_wavenum, m_hapicompliant, m_isotopefilter, m_lowerstate_global_quanta_filter.c_str(), m_upperstate_global_quanta_filter.c_str()  );
			ok               = (m_hitranchemical != NULL);
			ok               = ok && m_hitranchemical->SetLineShapeObject ( m_lineshapeobject );
			ok               = ok && m_hitranchemical->SetSelfBroadeningClimatology( m_chemicalnumberdensityguid, m_selfbroadeningclimatology);
			if ( NXFINITE((m_manualmaxlinestrength)) )
			{
				ok = ok && m_hitranchemical ->SetUserDefinedMaxLineStrength(m_manualmaxlinestrength);
			}
			double tol = NXFINITE((m_manualtolerance)) ? m_manualtolerance : 0.0;
			ok = ok && m_hitranchemical ->SetLineTolerance( tol );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::CheckDirtyAndUpdate, There were errors creating the Hitran Chemical Instance");
			}
		}
		m_isdirty = !ok;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetAtmosphericState		2013-3-20*/
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

bool skOpticalProperties_HitranChemical::SetAtmosphericState( skClimatology* atmosphericstate )
{
	if ( atmosphericstate != nullptr) atmosphericstate->AddRef();
	if ( m_atmosphericstateclimatology != nullptr) m_atmosphericstateclimatology->Release();
	m_atmosphericstateclimatology = atmosphericstate;
	return true;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged )
{
	bool	ok;

	m_lastpoint = pt;
	ok = !m_manual_microwindow_isset;
	if (!ok)
	{
		ok = CheckDirtyAndUpdate();
		if (crosssectionschanged != NULL) *crosssectionschanged = true;
		if (m_xs_optimizer != nullptr)
		{
			ok = ok && m_xs_optimizer->SetLocation(pt);
		}
		else
		{
			ok = ok && m_hitranchemical->UpdateLocation(pt, m_atmosphericstateclimatology);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::CalculateCrossSectionsInternal		2013-3-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::CalculateCrossSectionsInternal( double wavenumber, double* absxs, double* extxs, double* scattxs ) const
{
	double		absxsec = 0.0;
	bool		ok;

	ok =       (m_hitranchemical != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_HitranChemical::CalculateCrossSections, The internal hitran object is not yet loaded. Try calling SetAtmopshericState first, otherwise get the debugger going");
	}
	ok = ok && m_hitranchemical->AbsorptionCrossSection(wavenumber, &absxsec);
	if (!ok)
	{
		absxsec = std::numeric_limits<double>::quiet_NaN();
	}
	*absxs   = absxsec;
	*extxs   = absxsec;
	*scattxs = 0.0;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	bool ok = false;

	if (!m_manual_microwindow_isset )
	{
		SetWavenumberRange( wavenumber - 10, wavenumber +10 );
		SetLocation( m_lastpoint, nullptr );
	}
	if (m_xs_optimizer != nullptr) ok = m_xs_optimizer->CalculateCrossSections( wavenumber, absxs, extxs, scattxs);	// Are we using the optimized pre-cached cross-section configuration
	else                           ok = CalculateCrossSectionsInternal( wavenumber, absxs, extxs, scattxs );		// Guarantees thread safety as we call the const implementation
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::CheckWavenumberIsAscending		 2014- 10- 27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::CheckWavenumberIsAscending( const std::vector<double>&	wavenumber) const
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


/*---------------------------------------------------------------------------
 *                          min_and_max                           2020-09-15 */
/** **/
/*---------------------------------------------------------------------------*/

static void min_and_max( const double*	x, int nx, double* minval, double* maxval )
{
	*minval = x[0];
	*maxval = x[0];
	for (int i = 1; i < nx; i++)
	{
		if (x[i] < *minval) *minval = x[i];
		if (x[i] > *maxval) *maxval = x[i];
	}
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::CalculateCrossSectionsArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::CalculateCrossSectionsArray( const double*	userwavenumber,
																	  int			numwave,
																	  double*		userabsxs,
																	  double*		userextxs,
																	  double*	    userscattxs)
{
	bool	ok;
	std::vector<double>	absxs;
	std::vector<double> extxs;
	std::vector<double> scattxs;
	std::vector<double>	wavenumber;

	wavenumber.assign( userwavenumber,  userwavenumber+numwave);

	if (!m_manual_microwindow_isset )												// If the user has not explicitly called SetWavenumberRange, then do it now using the array of wavelengths passed in.
	{
		double	minwave;
		double	maxwave;

		min_and_max( userwavenumber, numwave, &minwave, &maxwave);					// Find the minimum and maximum value sin the array of wavenumbers
		SetWavenumberRange( minwave, maxwave);										// and set the wavenumber range
		SetLocation( m_lastpoint, nullptr );										// Set the location to force the update.
	}
	ok = (m_hitranchemical != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_HitranChemical::CalculateCrossSectionsArray, The internal hitran object is not yet loaded. Try calling SetAtmopshericState first, otherwise get the debugger going");
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_HitranChemical::CalculateCrossSectionsArray, Error calculating cross-sections. Returning a zeroed array");
		absxs.assign( wavenumber.size(), 0.0);
		extxs.assign( wavenumber.size(), 0.0);
		scattxs.assign( wavenumber.size(), 0.0);
	}
	else
	{

		absxs.assign( wavenumber.size(), 0.0);
		ok = m_hitranchemical->AddAbsorptionCrossSectionArray(wavenumber, &absxs);
		if (userscattxs != nullptr) scattxs.assign(wavenumber.size(), 0.0);
		extxs = absxs;
	}
	std::copy( absxs.begin(), absxs.end(), userabsxs );
	std::copy( extxs.begin(), extxs.end(), userextxs );
	if ( userscattxs != nullptr) std::copy( scattxs.begin(), scattxs.end(), userscattxs );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::InternalClimatology_UpdateCache		2013-3-20*/
/** Implements the skOpticalProperties interface function. The HITRAN optical
 *	properties do use a secondary internal climatology to calculate partial
 *	pressure of the chemical species using an ideal gas formulation. This function
 *	allows this skOpticalProperties to update its chemical number density information
 *	which is then used to calculate partial pressure.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool	ok;

	ok = CheckDirtyAndUpdate();
	ok = ok && m_hitranchemical->SelfBroadeningClimatology_UpdateCache( pt );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::InternalClimatology_UpdateCache, There were errors updating the internal climatology of the hitran chemical. Make sure it is properly defined");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetNumberDensityClimatology		2013-3-20*/
/** Sets the number density climatology associated with the selected chemical species.
 *	The current implementation uses this number density to calculate the species
 *	partial pressure. This in turn is used to work out the self broadening term used
 *	in the voigt, lorenz and other line shape objects. As such this is a fairly
 *	detailed refinement and is not required for medium precision calculations: the
 *	point being that it can be difficult for users to get profiles of various
 *	chemical species. If you choose not to use this function the partial pressure
 *	is set to 0.0 and self broadening of the line shape is ignored.
 *
 *	/param parameterguid
 *		The SKCLIMATOLOGY that provides the chemical species number density in molecules per cm3.
 *		The partial pressure is calculated using the ideal gas law p = NkT where T is the
 *		atmospheric temperature and is taken from the primary atmospheric state climatology.
 *
 *	/param numberdensityclimatology
 *		An instance of skClimatology* that provides the number density of the chemical species in
 *		molecules per cm3. The code keeps a reference count on this object and releases it when
 *		it is finished.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetSelfBroadeningClimatology( skClimatology* numberdensityclimatology  )
{
	if ( numberdensityclimatology != NULL) numberdensityclimatology->AddRef();
	if ( m_selfbroadeningclimatology    != NULL) m_selfbroadeningclimatology->Release();
	m_selfbroadeningclimatology = numberdensityclimatology;
	SetDirty();
	return true;
}


bool skOpticalProperties_HitranChemical::SetSelfBroadeningClimatologyHandle( const CLIMATOLOGY_HANDLE& parameterguid )
{
	m_chemicalnumberdensityguid =  parameterguid;
	SetDirty();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetLineShapeObject		2013-3-21*/
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

bool skOpticalProperties_HitranChemical::SetLineShapeObject( skSpectralLineShape* lineshapeobject )
{
	if (lineshapeobject != NULL) lineshapeobject->AddRef();
	if (m_lineshapeobject != NULL) m_lineshapeobject->Release();
	m_lineshapeobject = lineshapeobject;
	SetDirty();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetIsotopeFilter		 2014- 10- 15*/
/** Sets the Hitran chemical to use a single isotope **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetIsotopeFilter( size_t filterid)
{
	m_isotopefilter = (int)filterid;
	SetDirty();
	return true;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_HitranChemical::SetLowerStateGlobalQuantaFilter2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetLowerStateGlobalQuantaFilter( std::string lowerstateglobalquanta)
{
	m_lowerstate_global_quanta_filter = lowerstateglobalquanta;
	return true;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_HitranChemical::SetUpperStateGlobalQuantaFilter2020-08-25 */
/** **/
/*---------------------------------------------------------------------------*/
bool skOpticalProperties_HitranChemical::SetUpperStateGlobalQuantaFilter( std::string upperstateglobalquanta)
{
	m_upperstate_global_quanta_filter = upperstateglobalquanta;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetUserDefinedMaxLineStrength		 2016- 1- 25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetUserDefinedMaxLineStrength		( double maxlinestrength)
{
	bool	ok;

	ok = (maxlinestrength >= 0.0);
	if (ok)
	{
		m_manualmaxlinestrength = maxlinestrength;
		ok = (m_hitranchemical != NULL);
		ok = ok && m_hitranchemical->SetUserDefinedMaxLineStrength(maxlinestrength);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::SetUserDefinedMaxLineStrength, There were errors manually setting the micro-window line strength to %e", (double)maxlinestrength);
		}
	}
	else
	{
		m_manualmaxlinestrength = std::numeric_limits<double>::quiet_NaN();
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::SetUserDefinedMaxLineStrength, Negative values of max line strength are not allowed, %e", (double)maxlinestrength);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetLineTolerance		 2016- 1- 25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetLineTolerance( double tolerance )
{
	bool	ok = true;

	if (tolerance >= 0.0)
	{
		m_manualtolerance = tolerance;
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::SetLineTolerance, The line tolerance value (%20.10e) must be greater than or equal to 0.0", (double)tolerance);
	}
	if (m_hitranchemical != NULL) ok = ok && m_hitranchemical->SetLineTolerance(tolerance);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetChemicalName		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetChemicalName( const char* chemicalname )
{
	m_chemicalname = chemicalname;
	SetDirty();
	return true;
}


/*---------------------------------------------------------------------------
 *    skOpticalProperties_HitranChemical::SetMicroWindowMargin    2019-11-14 */
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetMicroWindowMargin( double margin_wavenum)
{
	bool ok;

	SetDirty();
	ok = (margin_wavenum >= 0.0) && (margin_wavenum < 1000000.0);
	if (ok)
	{
		m_margin_width_wavenum = margin_wavenum;
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::SetMicroWindowMargin. The micro-window marging must be a value between 0.0 and 1000000.0. We got [%e]", (double)margin_wavenum);
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_HitranChemical::SetWavenumberRange		2013-3-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::SetWavenumberRange( double lowwavenum, double highwavenum )
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
	m_manual_microwindow_isset = ok;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::SetWavenumberRange, The low wavenumber (%g) is bigger than the high wavenumber (%g). Thats not good.", (double)m_lowwavenum, (double)m_hihwavenum);
	}
	SetDirty();
	return ok;
}


/*---------------------------------------------------------------------------
 * skOpticalProperties_HitranChemical::EnableCachedCrossSections  2019-11-08 */
/** Enables the Hitran_CrossSection_Cache class which calculates and caches
 *	cross-sections for the given wavenumbers at every call to SetLocation
 *	after this point.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_HitranChemical::EnableCachedCrossSections( double* wavenumbers, size_t numwave )
{
	bool ok  =  (m_xs_optimizer != nullptr);

	if (ok)																				// If the cache already exists
	{																					// then
		for ( size_t i = 0; ok && (i < numwave); i++)									// see if all the wavenumbers 
		{																				// are already
			ok = ok && m_xs_optimizer->HasWavenumberAlreadyInCache(wavenumbers[i]);		// in the cache
		}																				// If every wavenumber is in the cache
	}																					// then we dont need to do anything else
	if (!ok)																			// otherwise we will do a flush
	{
		SetDirty();
		if ( m_xs_optimizer != nullptr) 
		{ 
			delete m_xs_optimizer; 
			m_xs_optimizer = nullptr;
		}
		if (numwave > 0)
		{
			nx1dArray<double>	wavenum( numwave, wavenumbers );
			double minval = 1.0E20;
			double maxval = -99999.0;
		
			for (size_t i = 0; i < numwave; i++)
			{
				minval = nxmin( minval, wavenumbers[i]);
				maxval = nxmax( maxval, wavenumbers[i]);
			}
			ok = (maxval > minval) && (minval > 0) && (maxval < 100000.0);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::EnableCachedCrossSections, The range of wavelengths from %e to %e is invalid.", (double)minval, (double)maxval);
			}
			else
			{
				ok = SetWavenumberRange( minval, maxval);		
				if (ok)
				{
					m_xs_optimizer = new Hitran_CrossSection_Cache(this);
					ok =       m_xs_optimizer != nullptr;
					ok = ok && m_xs_optimizer->SetCachedWavenumbers( wavenum );
				}
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"skOpticalProperties_HitranChemical::EnableCachedCrossSections, There were errors enabling cross-section caching. This will need debugging.");
				}
			}
		}
	}
	return ok;

}


