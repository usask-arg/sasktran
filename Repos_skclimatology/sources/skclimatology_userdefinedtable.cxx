#include <skclimatology21.h>

/*
HISTORY
2007-04-03	NDL.
   I changed the caching so default profiles are only loaded from files if the user makes a call for a species that
   has not been previously defined.  I also did some tweaking replacing some of the 'if's with 'if else'

*/



/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::skClimatology_UserDefinedTable		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserDefinedTable::skClimatology_UserDefinedTable( )
{
	m_climatetypes = NULL;
	m_numclimate = 0;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::skClimatology_UserDefinedTable		2008-3-5*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserDefinedTable::skClimatology_UserDefinedTable	( CLIMATOLOGY_HANDLE id, const char* filename )
{
	m_climatetypes = NULL;
	m_numclimate = 0;
	LoadProfileFromTextFile( &id, 1, filename );
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::~skClimatology_UserDefinedTable		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_UserDefinedTable::~skClimatology_UserDefinedTable( )
{
	ReleaseResources();
}


///*-----------------------------------------------------------------------------
// *					skClimatology_UserDefinedTable::DeepCopy		2008-3-4*/
///** **/
///*---------------------------------------------------------------------------*/
//
//bool skClimatology_UserDefinedTable::DeepCopy( const skClimatology_UserDefinedTable& other )
//{
//	bool	ok1;
//	bool	ok2;
//	bool	ok;
//	size_t	idx;
//
//	ok1 =        skClimatology::DeepCopy(other);
//	ok1 = ok1 && m_profile.DeepCopy( other.m_profile );
//	ok2 = AllocateClimatologyTable( other.m_numclimate );				// Allocate the array for the climatology table
//	if (ok2)																	// If that worked
//	{																		// then copy the tables
//		for (idx = 0; idx < m_numclimate; idx++)							// over to the new instance
//		{
//			m_climatetypes[idx] = other.m_climatetypes[idx];
//		}
//	}
//	ok = ok1 && ok2;
//	if (!ok)
//	{
//		nxLog::Record( NXLOG_WARNING,"skClimatology_UserDefinedTable::CreateClone, error initializing the clone ");
//		ReleaseResources();
//	}
//	return ok;
//
//}

/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::CreateClone		2008-3-3*/
/** **/
/*---------------------------------------------------------------------------*/

//bool skClimatology_UserDefinedTable::CreateClone( skClimatology** userclone ) const
//{
//	skClimatology_UserDefinedTable*	clone;
//	bool							ok;
//
//	clone = new skClimatology_UserDefinedTable;							// Create a new instance
//	ok = (clone != NULL);
//	if (!ok)
//	{
//		nxLog::Record(NXLOG_WARNING,"skClimatology_UserDefinedTable::CreateClone, Error allocating clone object");
//	}
//	else
//	{
//		clone->AddRef();														// add a reference that the user will release
//		ok = clone->DeepCopy(*this);
//	}
//	*userclone = clone;
//	return ok;
//}
//
/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::ReleaseResources		2008-2-29*/
/** **/
/*---------------------------------------------------------------------------*/

void skClimatology_UserDefinedTable::ReleaseResources()
{
	if (m_climatetypes != NULL ) delete [] m_climatetypes;
	m_climatetypes = NULL;
	m_numclimate = 0;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::AllocateClimatologyTable		2008-2-29*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::AllocateClimatologyTable( size_t numclimates )
{
	bool	ok;
	size_t	idx;

	ok = (numclimates == m_numclimate);
	if (!ok)
	{
		ReleaseResources();
		ok = (numclimates == 0);
		if (!ok)
		{
			m_climatetypes = new CLIMATOLOGY_HANDLE [numclimates];
			ok = (m_climatetypes != NULL);
			if (ok)
			{
				m_numclimate = numclimates;
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skClimatology_UserDefinedTable::AllocateClimatologyTable, error allocating space for %Iu climatology handles", (size_t)numclimates);
		ReleaseResources();
	}
	for (idx = 0; idx < m_numclimate; idx++) m_climatetypes[idx] = SKCLIMATOLOGY_UNDEFINED;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::LoadProfileFromTextFile		2006-7-3*/
/** Inputs a "N" column text file. The first column of the text file is
 *	height in meters. The second and sunbsequent columns the required species
 *	value. The number of species passed in must match the number of species columns
 *	in the text file.
 **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::LoadProfileFromTextFile( const CLIMATOLOGY_HANDLE* species, size_t numspecies, const char* filename )
{
	bool				ok;
	nx1dArray<double>	h;
	nx2dArray<double>	profile;

	ok = profile.InputColumnMajorText( filename, numspecies+1 );
	if (ok)
	{
		profile.XSlice( 0, &h  );
		if (nxarray::Max(h) < 250.0)
		{
			nxLog::Verbose(NXLOG_INFO, "skClimatology_UserDefinedTable::LoadProfileFromTextFile, The text file looks like it is in kilometers rather than meters, Im converting first column to meters");
			h *= 1000.0;	
		}
		ok = LoadProfileFrom2dArray( species, numspecies, profile );
	}
	else
	{
		ReleaseResources();
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefinedTable::LoadProfileFromTextFile, Error loading species profile from text file, setting this species to all zeros");
	}
	return( ok );
}



/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::LoadProfileFrom2dArray		2008-2-29*/
/** Loads the species profiles from a 2-D array. The object can load several
 *	species on one altitude grid.  The array is stored as as an array
 *	[numheights x numspecies +1]
 *
 *	The first column is the altitutude in meters.  The second and subsequent
 *	columns are profiles of various species. The number of species must match
 *	number of columns in the array-1 (the first column is always height)
 **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::LoadProfileFrom2dArray( const CLIMATOLOGY_HANDLE* species, size_t numspecies, const nx2dArray<double>& profile )
{
	bool				ok;
	nx1dArray<double>	h;
	size_t				idx;

	m_profile.DeepCopy( profile );
	ok =      (m_profile.YSize() == (numspecies+1));
	ok = ok && AllocateClimatologyTable(numspecies);
	if (ok)
	{
		for (idx = 0; idx < numspecies; idx++ ) m_climatetypes[idx] = species[idx];
	}
	else
	{
		ReleaseResources();
		m_profile.erase();
		nxLog::Record(NXLOG_WARNING, "skClimatology_UserDefinedTable::LoadProfileFromTextFile, Error loading species profile from profile, either memory allocation or column mismatch");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::InterpolateValueAtAltitude		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::InterpolateValueAtAltitude( double* value, double altitude, size_t columnidx )
{
	size_t					size;
	bool					ok;
	nx1dArray<double>		h;
	nx1dArray<double>		v;

	size =  m_profile.XSize( );							// See if we have a profile loaded for this species
	ok   = (size > 1);
	if (!ok)
	{
		*value = 0.0;
		ok     = true; 
	}
	else
	{
		m_profile.XSlice( 0,              &h );
		m_profile.XSlice( (int)columnidx, &v );
		nxArrayIter<double> hstart = h.begin();
		nxArrayIter<double> vstart = v.begin();
		*value = nxLinearInterpolate::EvaluateYatX( altitude, hstart, vstart, v.size(), nxLinearInterpolate::ENUM_TRUNCATE, MissingValue() );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::UpdateCache		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	bool ok;

	ok = (m_profile.size() > 0) && (m_numclimate > 0);
	SetCacheIsLoaded( ok );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::GetParameter		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	bool					ok;
	CLIMATOLOGY_HANDLE*		iter;
	CLIMATOLOGY_HANDLE*		finish;
	size_t					columnidx;

	ok = (m_numclimate > 0);
	if (ok)
	{
		if (updatecache) UpdateCache(placeandtime);
		finish = m_climatetypes + m_numclimate;
		iter   = std::find( m_climatetypes, finish, species );
		ok     = (iter != finish);
		if (ok)
		{
			columnidx = (iter - m_climatetypes) + 1;
			ok        = InterpolateValueAtAltitude( value, placeandtime.heightm , columnidx );
		}
	}
	if (!ok)
	{
		*value = MissingValue();
	}
	return( ok );
}


/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedTable::IsSupportedSpecies		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_UserDefinedTable::IsSupportedSpecies(const CLIMATOLOGY_HANDLE& species )
{
	bool					ok;
	CLIMATOLOGY_HANDLE*		iter;
	CLIMATOLOGY_HANDLE*		finish;

	ok = (m_numclimate > 0);
	if (ok)
	{
		finish = m_climatetypes + m_numclimate;
		iter   = std::find( m_climatetypes, finish, species );
		ok     = (iter != finish);
	}
	return ok;
}

