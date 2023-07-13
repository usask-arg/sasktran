#include <skclimatology21.h>

/*
HISTORY
2008-02-11	CZR.
   Added a zero climatology for O3, Aerosol, and NO2 (sets the corresponding species to zero density).

2007-04-03	NDL.
   I changed the caching so default profiles are only loaded from files if the user makes a call for a species that
   has not been previously defined.  I also did some tweaking replacing some of the 'if's with 'if else'

*/

//	currently supported species: Static //
//										//
//	SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3	//
//	SKCLIMATOLOGY_O3_CM3				//
//	SKCLIMATOLOGY_AEROSOL_CM3			//
//	SKCLIMATOLOGY_TEMPERATURE_K			//
//	SKCLIMATOLOGY_NO2_CM3				//
//										//
//	currently supported species: Zero   //
//										//
//	SKCLIMATOLOGY_O3_CM3				//
//	SKCLIMATOLOGY_AEROSOL_CM3			//
//	SKCLIMATOLOGY_NO2_CM3				//
//										//
//	currently supported species: OSIRIS //
//										//
//	SKCLIMATOLOGY_O3_CM3				//
//	SKCLIMATOLOGY_AEROSOL_CM3			//
//	SKCLIMATOLOGY_NO2_CM3				//
//										//
//	currently supported species: Pratmo //
//										//
//	SKCLIMATOLOGY_NO2_CM3				//




/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::skClimatology_TextFileDatabase		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_TextFileDatabase::skClimatology_TextFileDatabase( )
{
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::~skClimatology_TextFileDatabase		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_TextFileDatabase::~skClimatology_TextFileDatabase( )
{
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::SetPath		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::SetPath( const char* path )
{
	m_path  = path;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::LoadProfileFromTextFile		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::LoadProfileFromTextFile( const char* file, nx2dArray<double>* profile )
{
	bool				ok;
	nx1dArray<double>	h;
	nxString			filename;

	filename.sprintf( "%s/%s", (const char *)m_path, (const char *)file );
	ok = profile->InputColumnMajorText( filename );
	if (ok)
	{
		profile->XSlice( 0, &h  );
		if ( nxarray::Max(h) < 900.0 )
		{
			nxLog::Verbose(NXLOG_WARNING,"skClimatology_TextFileDatabase::LoadProfileFromTextFile, the heights in file <%s> appear to be in Kms. I am converting the heights to meters.", (const char*)file);
			h *= 1000.0;
		}
	}
	return( ok );
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::LoadProfileFromArray		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::LoadProfileFromArray( int numalts, double* alts, double* speciesarray, nx2dArray<double>* profile )
{
	int altctr;
	bool	ok;

	ok = profile->SetSize( numalts, 2 );
	if (ok)
	{
		for ( altctr=0; altctr<numalts; altctr++ )
		{
			profile->At(altctr,0) = alts[altctr];//*1000.0;
			profile->At(altctr,1) = speciesarray[altctr];
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::SetSpeciesProfile		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::SetSpeciesProfile( const CLIMATOLOGY_HANDLE& species, int numalts, double* alts_km, double* profile )
{
	bool ok = false;

	if      ( species == SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3 ) ok = LoadProfileFromArray( numalts, alts_km, profile, &m_neutraldensity );
	else if ( species == SKCLIMATOLOGY_O3_CM3 )               ok = LoadProfileFromArray( numalts, alts_km, profile, &m_ozonedensity   );
	else if ( species == SKCLIMATOLOGY_AEROSOL_CM3 )          ok = LoadProfileFromArray( numalts, alts_km, profile, &m_aerosoldensity );
	else if ( species == SKCLIMATOLOGY_TEMPERATURE_K )        ok = LoadProfileFromArray( numalts, alts_km, profile, &m_temperatureprofile ); 
	else if ( species == SKCLIMATOLOGY_NO2_CM3 )              ok = LoadProfileFromArray( numalts, alts_km, profile, &m_NO2density ); 
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::InterpolateValueAtAltitude		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::InterpolateValueAtAltitude( double* value, double altitude, nx2dArray<double>& profile, const char* filename_of_cache, const char* speciesname )
{
	int		size;
	bool	ok;
	nx1dArray<double>	h;
	nx1dArray<double>	v;
	nxArrayIter<double>	iter2;
	nxArrayIter<double>	iter1;
	nxArrayIter<double>	start;
	nxArrayIter<double>	finish;
	size_t				idx;
	double				h1,h2;
	double				v1,v2;
	double				f;

	size = (int)profile.XSize( );													// See if we have a profile loaded for this species
	if (size < 1)																	// if we dont
	{																				// then
		ok   = LoadProfileFromTextFile( filename_of_cache,   &profile  );			// load in from the filename cache
		size = (int)profile.XSize( );												// and see if it works
	}
	ok = (size > 1);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skClimatology_TextFileDatabase::InterpolateValueAtAltitude, I need at least 2 points in the height profile. The current profile for <%s> only has %d points", (const char*)speciesname, (int)size);
	}
	else
	{
		profile.XSlice( 0, &h );
		profile.XSlice( 1, &v );

		start  = h.begin();
		finish = h.end();
		iter2 = std::upper_bound( start, finish, altitude );
		
		if      (iter2 == start ) ++iter2;
		else if (iter2 == finish) --iter2;

		idx   = (iter2 - start);
		h2    = h.At(idx); 
		v2    = v.At(idx);
		--idx;
		h1    = h.At(idx);
		v1    = v.At(idx); 
		f     = (altitude-h1)/(h2-h1);
		ok = ( f > -0.6 ) && ( f < 1.6);			// Allow a little extraploation freedom at the ends of the array (usually fixes 1/2 km issues
		if (ok) *value = (1-f)*v1 + f*v2;
	}
	if (!ok)
	{
		*value = 0.0;							// Make sure 
		nxLog::Record(NXLOG_WARNING,"skClimatology_TextFileDatabase::InterpolateValueAtAltitude, failed to get a value at height %f meters, returning zero. You probably need to extend the table to this altitude to make the code work properly", (double)altitude);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::UpdateCache		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
/*
   This is now performed during calls to InterpolateValueAtAltitude.
   Avoids annoying warning messages being created for species the user is never actually interested in.

	bool	ok;

	m_species_ozone			= LoadProfileFromTextFile( "OzoneDensity.txt",   &m_ozonedensity       );
	m_species_neutral		= LoadProfileFromTextFile( "NeutralDensity.txt", &m_neutraldensity     ); 
	m_species_aerosol		= LoadProfileFromTextFile( "AerosolDensity.txt", &m_aerosoldensity     );
	m_species_temperature	= LoadProfileFromTextFile( "Temperature.txt",    &m_temperatureprofile );
	m_species_NO2			= LoadProfileFromTextFile( "NO2Density.txt",     &m_NO2density         );

	ok =     m_species_ozone
		  && m_species_neutral
		  && m_species_aerosol
		  && m_species_temperature
		  && m_species_NO2;
*/
//	bool ok;										// If it loads cache late rthen we must flag that update cache actually works.

//	ok =     ( m_aerosoldensity.size() > 0)
//		  && ( m_neutraldensity.size() > 0 );

	SetCacheIsLoaded( true );

	return true;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::GetParameter		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	bool ok;
	
	if      ( species == SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3 ) ok = InterpolateValueAtAltitude(  value,  placeandtime.heightm, m_neutraldensity,     "NeutralDensity.txt", "Neutrals" );
	else if ( species == SKCLIMATOLOGY_O3_CM3               ) ok = InterpolateValueAtAltitude(  value, placeandtime.heightm, m_ozonedensity,       "OzoneDensity.txt" ,  "O3"   );
	else if ( species == SKCLIMATOLOGY_AEROSOL_CM3          ) ok = InterpolateValueAtAltitude(  value, placeandtime.heightm, m_aerosoldensity,     "AerosolDensity.txt", "Aerosols"  );
	else if ( species == SKCLIMATOLOGY_TEMPERATURE_K        ) ok = InterpolateValueAtAltitude(  value, placeandtime.heightm, m_temperatureprofile, "Temperature.txt",    "Temperature"     );
	else if ( species == SKCLIMATOLOGY_NO2_CM3              ) ok = InterpolateValueAtAltitude(  value, placeandtime.heightm, m_NO2density,         "NO2Density.txt",     "NO2"      );
	else
	{
		*value = MissingValue();
		ok     = false; 
		nxLog::Record( NXLOG_WARNING, "skClimatology_TextFileDatabase::GetParameter, Not a supported species of skClimatology_TextFileDatabase. You might want to call SetPath or SetSpeciesProfile \n" );
	}

	return( ok );
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::IsSupportedSpecies		2006-7-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_TextFileDatabase::IsSupportedSpecies(const CLIMATOLOGY_HANDLE& species )
{
	bool ok;

	ok =       ( (species == SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3)) // && ( !m_neutraldensity.IsEmpty()     ))
            || ( (species == SKCLIMATOLOGY_O3_CM3) ) //              && ( !m_ozonedensity.IsEmpty()       ))
	        || ( (species == SKCLIMATOLOGY_NO2_CM3))//              && ( !m_NO2density.IsEmpty()         ))
	        || ( (species == SKCLIMATOLOGY_TEMPERATURE_K))//        && ( !m_temperatureprofile.IsEmpty() ))
	        || ( (species == SKCLIMATOLOGY_AEROSOL_CM3));//          && ( !m_aerosoldensity.IsEmpty()     ));		  

	return( ok );
}


/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::DeepCopy		2008-12-11*/
/** **/
/*---------------------------------------------------------------------------*/

//bool skClimatology_TextFileDatabase::DeepCopy( const skClimatology_TextFileDatabase& other )
//{
//	bool	ok;
//
//	m_path = other.m_path;
//	ok = skClimatology::DeepCopy(other);
//	ok = ok && m_ozonedensity.DeepCopy      ( other.m_ozonedensity );
//	ok = ok && m_neutraldensity.DeepCopy    ( other.m_neutraldensity );
//	ok = ok && m_aerosoldensity.DeepCopy    ( other.m_aerosoldensity );
//	ok = ok && m_temperatureprofile.DeepCopy( other.m_temperatureprofile );
//	ok = ok && m_NO2density.DeepCopy        ( other.m_NO2density);
//	if (!ok)
//	{
//		nxLog::Record(NXLOG_WARNING,"skClimatology_TextFileDatabase::DeepCopy, There was an error copying the original object. That is probably going to cause problems");
//	}
//	return ok;
//}

/*-----------------------------------------------------------------------------
 *					skClimatology_TextFileDatabase::CreateClone		2008-12-11*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool skClimatology_TextFileDatabase::CreateClone( skClimatology** userclone ) const
//{
//	skClimatology_TextFileDatabase*	clone;
//	bool							ok;
//
//	clone = new skClimatology_TextFileDatabase;
//	ok = (clone != NULL );
//	if (!ok)
//	{
//		nxLog::Record(NXLOG_WARNING,"skClimatology_TextFileDatabase::CreateClone, Error allocating clone object");
//	}
//	else
//	{
//		clone->AddRef();
//		ok = clone->DeepCopy( *this);
//		if (!ok)
//		{
//			nxLog::Record(NXLOG_WARNING,"skClimatology_TextFileDatabase::CreateClone, Error copying this instance to the clone");
//		}
//	}
//	*userclone = clone;
//	return ok;
//}
//
