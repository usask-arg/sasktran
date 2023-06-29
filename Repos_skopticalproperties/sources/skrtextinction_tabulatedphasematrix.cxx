#include <skopticalproperties21.h>

//#include <sasktranv21.h>
//#include "skrtextinction_tabulatedphasematrix.h"
#include <fstream>



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry		2009-6-16*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry()
{
	m_roundedwavelen = 0.0;
	m_phasematrix.SetStatic();
	m_cosscatterangle.SetStatic();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::ConfigureEntry		2009-6-17*/
/** Configure sthis entry. Sets up the wavelength by rounding. Copies over the 
 *	phasematricx and cos(scattering angle) arrays. Checksx the cos(scattering angle) 
 *	to make sure it is in proper range and order.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::ConfigureEntry( double wavelennm, const double* phasematrix, const double* cosangle, size_t numpoints)
{
	bool		ok;
	size_t		idx;

	m_roundedwavelen = RoundTheWavelength(wavelennm);

	ok = (cosangle[0] >= -1.00001) && (cosangle[0] <= 1.00001);						// Make sure the cosine of scattering angle is in proper range. Check first element
	for (idx = 1; idx < numpoints; idx++ )											// also check that it is in ascending order (this is easy for the user to get wrong)
	{																				// as scattering anlgle 0 to 180 is opposite to cos(scatter angle)
		ok = ok && (cosangle[idx] >= -1.00001) && (cosangle[idx] <= 1.00001) && (cosangle[idx] > cosangle[idx-1]);
	}
	ok = ok && m_phasematrix.CopyGridArray    ( phasematrix, numpoints );			// Create the phase matrix grid
	ok = ok &&  m_cosscatterangle.CopyGridArray( cosangle,    numpoints );			// and the scattering angle grid. They must be the same size.

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::ConfigureEntry, The entry is invalid as the cos(scattering angles) are not in ascending order (-1 to +1) or are out of the range (-1 to +1) or there was a memory allocation error.");
		ReleaseResources();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::GetPhaseMatrix		2009-6-16*/
/**	Returns the scalar phasematrix using linear interpolation of the
 *	scattering angle grid. The 
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::GetPhaseMatrix( double cosangle, double* phasematrix ) const
{
	SKTRAN_GridIndex	lowercell;
	SKTRAN_GridIndex	uppercell;
	double				lowerweight;
	double				upperweight;
	bool				ok;
	
	*phasematrix = 0;
	ok = (cosangle >= -1.00001) && (cosangle <= 1.00001);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::GetPhaseMatrix, The cosine(scatteringing) is not in the range -1 to 1. Make sure you are using cosine(scattering angle). Do OT use radians or degree!");
	}
	else
	{
		ok = m_cosscatterangle.FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowercell, &lowerweight, &uppercell, &upperweight);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::PhaseMatrix, There was an error looking up phase matrix for cos(scatter angle)  = %g, setting value to 0.0", (double)cosangle);
		}
		else
		{
			*phasematrix = lowerweight* m_phasematrix.At(lowercell) + upperweight*m_phasematrix.At(uppercell);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::ReleaseResources		2009-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::ReleaseResources()
{
	m_phasematrix.erase();
	m_cosscatterangle.erase();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::DeepCopy		2009-6-17*/
/** Make a cloned copy of the other object.
 **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::DeepCopy( const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry& other )
{
	bool	ok;

	m_roundedwavelen = other.m_roundedwavelen;
	ok =       m_phasematrix.DeepCopy( other.m_phasematrix );
	ok = ok && m_cosscatterangle.DeepCopy( other.m_cosscatterangle );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::DeepCopy, Error copying wavelength entries, thats a problem ");
		ReleaseResources();
	}
	return ok;
}
*/



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		2009-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedPhaseMatrix_HeightEntry	()
{
	m_heightmeters = -999999.0;
	m_numwavelens = 0;
	m_maxwavelens = 0;
	m_wavelenentries = NULL;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		2009-11-26*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedPhaseMatrix_HeightEntry( double heightm )
{
	m_heightmeters = heightm;
	m_numwavelens = 0;
	m_maxwavelens = 0;
	m_wavelenentries = NULL;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::~skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		2009-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::~skOpticalProperties_TabulatedPhaseMatrix_HeightEntry()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::ReleaseResources		2009-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::ReleaseResources()
{
	if (m_wavelenentries != NULL)
	{
		delete [] m_wavelenentries;
		m_wavelenentries = NULL;
	}
	m_numwavelens = 0;
	m_maxwavelens =0;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::AllocateEntries		2009-6-17*/
/** Allocates N entries. Uses the small optimation os using the same memory if the
 *	the current emory allocation is sufficient.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::AllocateEntries	( size_t maxentries )
{
	bool	ok;

	ok = (maxentries <= m_maxwavelens);
	if (!ok)
	{
		ReleaseResources();
		m_wavelenentries = new skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry[maxentries];
		ok = (m_wavelenentries != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightEntry, Error allocating memeory for %u wavelength entries", (unsigned int)maxentries );
		}
		else
		{
			m_maxwavelens  = maxentries;
		}
	}
	if (ok) m_numwavelens = maxentries;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile		2009-6-17*/
/** Reads in the phase matrix entries for a set of wavelengths from a
 *	single text file associated with one altitude. The description of the file
 *	format can be found in the class description.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile( double heightm, const char* filename )
{
	std::ifstream		strm;
	bool				ok;
	bool				ok1;
	size_t				numwavelengths;
	size_t				numangles;
	double				wavelen_nm;
	nx1dArray<double>	phasematrix;
	nx1dArray<double>	cosangles;
	size_t				i;
	

	m_heightmeters = heightm;
	strm.open(filename, std::ios_base::in );
	ok = strm.is_open();

	strm >> numwavelengths >> numangles;						// Read in the number of wavelengths
	ok = !strm.fail();											// MAke sure the stream is ok
	ok = ok && AllocateEntries( numwavelengths );				// Allocate space for the wavelength entries
	ok = ok && phasematrix.SetSize (numangles);					// Allocate a temporary buffer for holding the phasematrix
	ok = ok && cosangles.SetSize(numangles );					// and another for holding the cos(scatter angle)
	
	strm >> cosangles; ok = ok && !strm.fail();					// read in the cos(scatter angle) from the file

	for (i = 0; i < numwavelengths; i++ )
	{
		strm >> wavelen_nm;
		strm >> phasematrix;
		ok1 = !strm.fail();
		ok1 = ok1 && m_wavelenentries[i].ConfigureEntry( wavelen_nm, phasematrix.UnsafeArrayBasePtr(), cosangles.UnsafeArrayBasePtr(), numangles );
		if (!ok1)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile, There was an error loading phase-matrix entries for the wavelength entry %u, wavelength = %g", (unsigned int)i, (double)wavelen_nm);
		}
		ok = ok && ok1;
	}
	strm.close();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile, There was an error loading the pahse matrix tables from file <%s>, I'm deleting the  whole structuire form memory");
		ReleaseResources();
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::DeepCopy		2009-6-17*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::DeepCopy( const skOpticalProperties_TabulatedPhaseMatrix_HeightEntry& other )
{
	bool	ok;

	m_heightmeters = other.Height();
	ok = AllocateEntries( other.m_numwavelens );
	for (size_t idx = 0; idx < other.m_numwavelens; idx++)
	{
		ok = ok && m_wavelenentries[idx].DeepCopy( other.m_wavelenentries[idx] );
	}
	m_numwavelens = other.m_numwavelens;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::DeepCopy, Error copying one intsance to another. Thats a problem");
		ReleaseResources();
		m_heightmeters = -9999999.0;
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					entries_equal		2009-6-18*/
/** \internal
 *	A quick and dirty comparator Predicate class used by the stl find_if in
 *	method FindWavelengthEntry
**/
/*---------------------------------------------------------------------------*/

class entries_equal
{
	private:
		double	m_wavelen;

	public:
					entries_equal	(double roundedwavelen)														{ m_wavelen = roundedwavelen;}
		bool		operator		() (const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry& entry) const { return m_wavelen == entry.WavelengthIndex();}
};


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::FindWavelengthEntry		2009-6-17*/
/** Find teh wavelength entry at this altitudde**/
/*------------------------------------------------------------------------c---*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::FindWavelengthEntry( double wavelennm, const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry** entry ) const
{
	double														roundednm;
	const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry*	start;
	const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry*	iter;
	const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry*	finish;	


	start  = m_wavelenentries;
	iter   = start;
	finish = start + m_numwavelens;

	roundednm = skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry::RoundTheWavelength( wavelennm );
	entries_equal	comparator(roundednm);

	iter = std::find_if( start, finish, comparator );
	if (iter == finish)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::FindWavelengthEntry, could not find a match for wavelength %g", (double)wavelennm);
		iter = NULL;
	}
	*entry = iter;
	return (iter != NULL);
}





/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength		2009-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength()
{
	m_crosssectionobject = NULL;
	m_lowerheight        = NULL;
	m_upperheight        = NULL;
	m_f1                 = 0.0;
	m_f2                 = 0.0;							// the interpolating factor for the upepr height entry.
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::~skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength		2009-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::~skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength()
{
	if (m_crosssectionobject != NULL) m_crosssectionobject->Release();
	m_crosssectionobject = NULL;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetCrossSectionObject		2009-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetCrossSectionObject( skOpticalProperties* crosssectionobject )
{
	if (crosssectionobject   != NULL) crosssectionobject->AddRef();
	if (m_crosssectionobject != NULL) m_crosssectionobject->Release();
	m_crosssectionobject = crosssectionobject;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CheckCrossSectionValid		2009-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CheckCrossSectionValid() const
{
	bool	ok;

	ok = (m_crosssectionobject != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CheckCrossSectionValid, You must successfully call SetCrossSectionObject with a real object before doing any sort of real calculations");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::ConfigureHeightInterpolation		2009-6-19*/
/** Finds the entries that bound the specified height and configures the linear combination.
 *	The algorithm sets phase matrices to zero for heights outside the upper and lower bounds
 *	of the set of heights. This is important if the phase matrix is only defined 
 *	for mesospheric altitudes as we do not want it being used all the way to the ground.
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::ConfigureHeightInterpolation( double heightm )
{
	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		dummy(heightm);
	iterator											iter;
	bool												ok;
	double												h1;
	double												h2;
	double												dh;

	iter = std::lower_bound( m_heightentries.begin(), m_heightentries.end(), dummy);		// Find the entry that is "bigger" than this entry
	ok   = !(iter == m_heightentries.begin()) && !(iter == m_heightentries.end());			// ok if we are inside the range of altitudes
	if (ok )
	{
		m_upperheight = &(*iter);
		--iter;
		m_lowerheight = &(*iter);
		h2            = m_upperheight->Height();
		h1            = m_lowerheight->Height();
		dh            = h2-h1;
		m_f1          = (h2-heightm)/dh;
		m_f2          = 1-m_f1;
	}
	else
	{
		m_upperheight = NULL;
		m_lowerheight = NULL;
		m_f1          = 0.0;
		m_f2          = 0.0;
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculatePhaseMatrix		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	return CalculatePhaseMatrix( wavenumber, cosscatterangle, phasematrix );
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculatePhaseMatrix		2009-6-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix ) const
{
	bool															ok;
	double															wavelennm;
	const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry*		upperentry = nullptr;
	const skOpticalProperties_TabulatedPhaseMatrix_WavelengthEntry*		lowerentry = nullptr;
	double															upperp = 0.0;
	double															lowerp = 0.0;
	

	wavelennm = 1.0E7/wavenumber;
	phasematrix->SetTo(0.0);

	ok =       (m_upperheight != NULL) && (m_lowerheight != NULL);
	ok = ok && m_upperheight->FindWavelengthEntry(wavelennm, &upperentry);
	ok = ok && m_lowerheight->FindWavelengthEntry(wavelennm, &lowerentry);
	ok = ok && upperentry->GetPhaseMatrix( cosscatterangle, &upperp );
	ok = ok && lowerentry->GetPhaseMatrix( cosscatterangle, &lowerp );
	
	phasematrix->At(1,1) = (SKRTFLOAT)(m_f1*lowerp + m_f2*upperp);
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetAtmosphericState		2009-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetAtmosphericState( skClimatology* neutralatmosphere )
{
	bool	ok;

	ok =       CheckCrossSectionValid();
	ok = ok && m_crosssectionobject->SetAtmosphericState( neutralatmosphere );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetAtmosphericState, There was an error setting the atmospheric state for the specified point. That might be a subtle problem as it disables this objects influence on any calculation");
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetAtmosphericState		2009-6-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetLocation( const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )
{
	bool	ok;

	ok =       CheckCrossSectionValid();
	ok = ok && m_crosssectionobject->SetLocation( pt, crosssectionschanged );
	ok = ok && ConfigureHeightInterpolation( pt.heightm );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::SetAtmosphericState, There was an error setting the atmospheric state for the specified point. That might be a subtle problem as it disables this objects influence on any calculation");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::InternalClimatology_UpdateCache		2011-8-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt )
{
	bool	ok;

	ok =       CheckCrossSectionValid();
	ok = ok && m_crosssectionobject->InternalClimatology_UpdateCache( pt );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::InternalClimatology_UpdateCache, There was an error updating the internal climatology. That might be a subtle problem as it disables this objects influence on any calculation");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculateCrossSections		2009-6-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool	ok;

	ok = CheckCrossSectionValid();
	ok = ok && m_crosssectionobject->CalculateCrossSections(wavenumber, absxs, extxs, scattxs);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CalculateCrossSections, Error calculating cross-sections. Thats a problem.");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::IsScatterer		2009-6-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::IsScatterer() const
{
	return ((m_crosssectionobject != NULL) && m_crosssectionobject->IsScatterer());
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::IsAbsorber		2009-6-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::IsAbsorber() const
{
	return ((m_crosssectionobject != NULL) && m_crosssectionobject->IsAbsorber());
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::ReleaseResources		2009-6-19*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::ReleaseResources()
{
	m_heightentries.clear();
	m_lowerheight = NULL;
	m_upperheight = NULL;
	m_f1 = 0.0;
	m_f2 = 0.0;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::LoadHeightWavelengthProfileFromMasterFile		2009-6-19*/
/** Loads in the set of phasematrices asa function of height, wavelength and scattering angle.
 *	This method is responsible for managing the entire load and it controls
 *	loading the height profiles. This consists of reading 
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::LoadHeightWavelengthProfileFromMasterFile( const char* masterfilename )
{
	std::ifstream										strm;
	double												h;
	double												lasth = -9999999.0;
	char												buffer[1000];
	nxString											filename;
	nxString											fullfilename;
	nxFileSpec											spec(masterfilename);
	bool												ok=false;
	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		dummy;
	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry*	entry;

	ReleaseResources();

	strm.open(masterfilename, std::ios_base::in );
	while ((!strm.eof() && !strm.fail()) )
	{
		strm >> h;										// Input the height
		if (!strm.eof() && !strm.fail())				// IF we have not hit the end of file or had an error
		{												// then
			ok = (lasth < h);								
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "LoadHeightWavelengthProfileFromMasterFile, the height profiles are not specified in ascending order, this is going to mess up future height indexing and interpolation");
			}
			else
			{
				lasth = h;
				strm.getline( buffer, 1000);
				filename = buffer;
				filename.RemoveWhiteSpace();
				fullfilename = spec.FullDirSpec() + filename;

				#if defined(NX_WINDOWS)
					char fullexpandedname[_MAX_PATH+1];
					if( _fullpath( fullexpandedname,fullfilename, _MAX_PATH ) != NULL ) fullfilename = fullexpandedname;
				#else
					NXTRACE_ONCEONLY(firsttime, ("skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::LoadHeightWavelengthProfileFromMasterFile, Got to check that filenames are properly resolved in the Lunix code\n"));
				#endif

				ok = nxDirectory::FileExists(fullfilename);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::LoadHeightWavelengthProfileFromMasterFile, Cannot load in phase matrix data for height %g as the file <%s> does not exist", (double)h, (const char*)fullfilename);
				}
				else
				{
					m_heightentries.push_back(dummy);
					entry = &m_heightentries.back();
					ok = entry->LoadWavelengthEntriesFromHeightFile( h, fullfilename );
					if (!ok) m_heightentries.pop_back();
				}
			}
		}
	}
	strm.close();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CreateClone		2009-6-19*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CreateClone( skOpticalProperties ** userclone ) const
{
	skOpticalProperties*								xsectionclone;
	skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength*	clone = NULL;
	bool													ok;
	bool													ok1;
	const_iterator											iter;
	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry			newentry;
	skOpticalProperties_TabulatedPhaseMatrix_HeightEntry*		entry;

	ok = CheckCrossSectionValid();
	if (ok)
	{
		clone = new skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength;
		ok    = (clone != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CreateClone, Memory allocation error creating clone. Thats a problem");
		}
		else
		{
			clone->AddRef();
			ok =       m_crosssectionobject->CreateClone( &xsectionclone );
			ok = ok && clone->SetCrossSectionObject( xsectionclone );
			if (xsectionclone != NULL ) xsectionclone->Release();
			for (iter = m_heightentries.begin(); !(m_heightentries.end() == iter); ++iter)
			{
				clone->m_heightentries.push_back(newentry);
				entry = &(clone->m_heightentries.back());
				ok1 = entry->DeepCopy( *iter );
				ok = ok && ok1;
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength::CreateClone, error making copy. Thats a problem.");
				clone->ReleaseResources();
			}
		}
	}
	*userclone = clone;
	return ok;
}
*/



