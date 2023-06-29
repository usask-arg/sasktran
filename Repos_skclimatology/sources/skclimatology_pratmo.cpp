#include <skclimatology21.h>

//using namespace std;

static const size_t			NUMTIMESPERYEAR    = 24;						// There are two time bins per month, 24 per year
static const size_t			NUMHEIGHTS         = 101;						// There are 101 height bins
static const size_t			NUMLATITUDE        = 71;						// There are 71 latitude bins


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::skClimatology_Pratmo					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Pratmo::skClimatology_Pratmo( )
                     : m_config( "USask-ARG","/Climatology/Pratmo/",nxRegistryConfiguration::GLOBAL_INI, true)

{
	m_cachemjd       = -99999.0;
	m_cachelatitude  = -99999.0;
	m_cachelongitude = -99999.0;
	ClearCurrentIndex();
	NXTRACE(("***** skClimatology_Pratmo, This code still has lots of bugs etc NEED TO CHECK OUT CreateClone, UpdateCache, DeepCopy etc. etc.*****\n"));
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::~skClimatology_Pratmo					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_Pratmo::~skClimatology_Pratmo( )
{
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::ClearCurrentIndex		2008-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

void skClimatology_Pratmo::ClearCurrentIndex()
{

	m_timeIndex = -1;
	m_angleIndex = -1;
	m_sza        = -9999.0;
	m_isPM       = false;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::SetPath					2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

nxBOOL skClimatology_Pratmo::SetPath( const char* path )
{
	m_path  = path;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOsirisL2_FileLocator::CheckBaseDirectory		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::FetchBaseDirectory( )
{
	bool	ok;

	ok = !m_path.IsEmpty();
	if (!ok)
	{
		ok = m_config.LocateDirectoryFromKey("Storage", &m_path, true, true, "Browse for the directory containg the PRATMO database file <skClimatology_PratmoNO2.bin>.");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::TimeIndex		2008-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

int skClimatology_Pratmo::TimeIndex( const GEODETIC_INSTANT& geopt )
{
	int				year;
	int				month;
	int				day;
	nxTimeStamp		timeStamp(geopt.mjd);
	int				timeIndex;
	
	
	timeStamp.GetDate( &day, &month, &year );

	timeIndex = ((month-1)*2);
	if( day >= 8 )					// Calculate the time index (0..23)
	{
		if( day <= 22 ) timeIndex += 1;
		else
		{
			timeIndex += 2;
			if( month == 12) timeIndex = 0;
		}
	}
	return timeIndex;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::IsPM		2008-12-12*/
/** alculate if it is local am or local pm**/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::IsPM( const GEODETIC_INSTANT& geopt )
{
	bool	ispm;
	double	localtime;

	localtime  = geopt.mjd + geopt.longitude*12.0/180.0;
	localtime -= floor(localtime);

	ispm = ( localtime > 0.5 )? true : false;
	return ispm;
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::SZA		2008-12-12*/
/** Calculate the Solar Zenith Angle using the Geodetic Instant**/
/*---------------------------------------------------------------------------*/

double skClimatology_Pratmo::SZA(const GEODETIC_INSTANT& geopt)
{
	PlanetSun   sun;
	nxGeodetic  geoid;
	nxVector    sununit;
	nxVector    zenith;
	nxVector    south;
	nxVector    west;
	double		sza;

	sun.UpdateECIPosition(geopt.mjd);
	sununit = sun.Location().EquatorialToGeographic( geopt.mjd ).UnitVector();

	geoid.FromGeodetic( geopt.latitude, geopt.longitude, 0 );
	geoid.GetGeodeticWestSouthUp  (  &west, &south, &zenith );
	sza = 180.0/nxmath::Pi*acos(zenith & sununit);
	return sza;
}

/*-----------------------------------------------------------------------------
 *					size_tskClimatology_Pratmo::LatitudeIndex		2008-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

int skClimatology_Pratmo::LatitudeIndex( const GEODETIC_INSTANT& geopt)
{
	int		angleIndex;
	double	latitude;
	double	angle;

	latitude   = geopt.latitude;
	angleIndex = 0;
	angle      = -87.5;
	
	while( latitude > (angle + 1.25) )
	{
		angle += 2.5;
		++angleIndex;
	}
	if( 71 == angleIndex)  angleIndex = 70;
	return angleIndex;
}

/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::InitializeMemory		2008-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::InitializeMemory()
{
	bool	ok;
	size_t	altCtr;
	
	ok = !m_NO2density.IsEmpty();
	if (!ok)
	{
		ok = m_NO2density.SetSize( NUMHEIGHTS, 2 );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skClimatology_Pratmo::CheckMemoryAllocation(), error allocating memory");
		}
		else
		{
			for( altCtr = 0; altCtr < NUMHEIGHTS; ++altCtr )
			{
				m_NO2density.At(altCtr,0) = altCtr*1000.0;
				m_NO2density.At(altCtr,1) = MissingValue();
			}
		}
		ClearCurrentIndex();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::GoToStartOfTimeIndex		2009-2-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::GoToStartOfTimeIndex( size_t timeIndex, std::ifstream& nfile)
{
	bool	ok;

	nxDWORD					numProfilesPerTime[NUMTIMESPERYEAR];
	std::ifstream::pos_type	offset;
	size_t					idx;

	ok = (nfile.read( (char *)&numProfilesPerTime[0], NUMTIMESPERYEAR*sizeof(numProfilesPerTime[0]) )).good();	// Read in the index that tells us the total number of profiles per time bin (there are 24 of them)

	if (ok)																					// If we are good then find the start of the 
	{																						// requested record
		offset = 0;																			// so reset the offset counter
		for( idx = 0; idx < timeIndex; ++idx)												// and for each record before ours
		{																					// get
			offset += numProfilesPerTime[idx]*(2*(NUMHEIGHTS+1))*sizeof(double) + NUMLATITUDE*sizeof(nxDWORD);		// the size of the record
		}																					// and doo all of the records
		ok =  ( nfile.seekg( offset, std::ios_base::cur )).good();								// then seek this position
	}	
	return ok;																	
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::GoToStartOfLatitudeIndex		2009-2-20*/
/** Skip for the start of the current time record to the start of the **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::GoToStartOfLatitudeIndex( size_t latitudeIndex, bool isPM, std::ifstream& nfile, size_t* numszaprofiles )
{
	bool						ok;
	nxDWORD						numProfilesPerLat [NUMLATITUDE];
	std::ifstream::pos_type		offset;
	size_t						idx;


	ok = ( nfile.read( (char *)&numProfilesPerLat[0], NUMLATITUDE*sizeof(numProfilesPerLat[0]))).good();
	if (ok)
	{
		*numszaprofiles = numProfilesPerLat[latitudeIndex];
		offset = 0;										
		for( idx = 0; idx < latitudeIndex; ++idx)
		{
			offset += numProfilesPerLat[idx]*(2*(NUMHEIGHTS+1));
		}
		if (isPM) offset += numProfilesPerLat[latitudeIndex]*(NUMHEIGHTS+1);
		ok = ( nfile.seekg( offset*sizeof(double), std::ios_base::cur )).good();
	}
	if (!ok) *numszaprofiles = 0;
	return ok;
}
			

/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::InterpolateSZAprofiles		2009-2-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::InterpolateSZAprofiles( size_t numszaprofiles, double sza, std::ifstream& nfile )
{
	double*						szas;
	double*						szasend;
	double*						szaIter;
	double						profile1[NUMHEIGHTS];
	double						profile2[NUMHEIGHTS];
	bool						outofbounds;
	bool						ok;
	size_t						idx1;
	std::ifstream::pos_type		offset;
	size_t						hidx;
	double						w1,w2;

	ok =  (numszaprofiles > 0);
	if (ok)
	{
		szas    = new double[numszaprofiles];
		szasend = szas + numszaprofiles;
		ok      = ( nfile.read( (char *)szas, (std::streamsize)(numszaprofiles*sizeof(*szas)))).good();		// Read in the array of SZA angles
		if (ok)
		{
			szaIter     = std::upper_bound( szas, szasend, sza );					// Find the first entry where (entry > sza)
			idx1        = (szaIter-szas);											// Get index of first entry greater than 
			outofbounds = (idx1 >= numszaprofiles) || (idx1 == 0);					// See if we are before the first or after the last entry
			if (outofbounds)														// if we are
			{																		// then 
				if (idx1 > 0) idx1 = numszaprofiles -1;								// if we are not before the first then we are goind to use the last record
			}																		// and that is that
			else																	// otherwise we 
			{																		// are straddled by two entries
				--idx1;																// so point to the lower bounding value
			}
			offset = (std::streamsize)(idx1*NUMHEIGHTS*sizeof(double));				// so skip over the number of profiles
			ok     =       ( nfile.seekg( offset, std::ios_base::cur )).good();
			ok     = ok && ( nfile.read( (char *)&profile1[0], NUMHEIGHTS*sizeof(double))).good();
			if (!outofbounds)
			{
				ok  = ok && ( nfile.read ( (char *)&profile2[0], NUMHEIGHTS*sizeof(double)) ).good();
				if (ok)
				{
					w2  = ( sza - szas[idx1] ) / ( szas[idx1+1] - szas[idx1] );
					w1  = 1-w2;

					for( hidx = 0; hidx < NUMHEIGHTS; ++hidx )
					{
						 profile1[hidx] = profile1[hidx] * w1 + profile2[hidx] * w2;
					}
				}
			}
		}
		delete [] szas;
	}
	if (ok) for( hidx = 0; hidx < NUMHEIGHTS; ++hidx ) m_NO2density.At(hidx,1) = profile1[hidx];
	else    for( hidx = 0; hidx < NUMHEIGHTS; ++hidx ) m_NO2density.At(hidx,1) = MissingValue();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::LoadProfileFromBinaryFile			2008-2-11*/
/**  skClimatology_Pratmo data is indexed by:
 *		time (bi-monthly),				exactly 24
 *		latitude (-87.5:2.5:87.5),		exactly 71
 *		heights	 (0 to 101)				exactly 102
 *		am/pm							exactly 2
 *		and sza (varies),
 *		Load the number of profiles in each time
 *
 *	The file starts with a 24 (4 byte integer) element array that says
 *	how many records are in each timebin 
 *
 *	Each time bin record consists
 *	a 71 (4 byte integer) element array that says how many Am and PM height profiles
 *  are in each latitude bin (ie it is the value of the sza index)
**/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::LoadProfileFromBinaryFile( const GEODETIC_INSTANT& geopt )
{
	nx1dArray<double>			h;
	nxString					filename;
	bool						ok;
	std::ifstream				nfile;
	size_t 						timeIndex;
	size_t						latitudeIndex;
	double						sza;
	bool						isPM;
	size_t						numszaprofiles;

	ok = InitializeMemory();
	if (ok)
	{
		timeIndex     = TimeIndex    ( geopt );
		latitudeIndex = LatitudeIndex( geopt );
		sza           = SZA			 ( geopt );
		isPM          = IsPM		 ( geopt );

		ok = ok && ( timeIndex     == m_timeIndex)			// See if we need to actually updat ethe climatology
				&& ( latitudeIndex == m_angleIndex)		// These are the 4 indices.
				&& ( sza           == m_sza )				// If they are the same then we dont need
				&& (isPM           == m_isPM );			// to do anything
		if (!ok)
		{
			ok = FetchBaseDirectory();
			if (ok)
			{
				filename.sprintf( "%s/skClimatology_PratmoNO2.bin", (const char *)m_path );

				nfile.open( filename, std::ios::binary );
				ok =       nfile.is_open();
				ok = ok && GoToStartOfTimeIndex    ( timeIndex, nfile);
				ok = ok && GoToStartOfLatitudeIndex( latitudeIndex, isPM, nfile, &numszaprofiles ); 
				ok = ok && InterpolateSZAprofiles  ( numszaprofiles, sza, nfile );
				if (ok)
				{
					m_timeIndex  = timeIndex;
					m_angleIndex = latitudeIndex;
					m_sza        = sza;
					m_isPM       = isPM;
				}
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skClimatology_Pratmo::LoadProfileFromBinaryFile, Error loading the PRATMO climatology. Erasing the profile");
		m_NO2density.erase();
		ClearCurrentIndex();
	}
	return	ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::InterpolateValueAtAltitude		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::InterpolateValueAtAltitude( const GEODETIC_INSTANT& placeandtime, double* value )
{
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

	ok = CacheIsLoaded();
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skClimatology_Pratmo::InterpolateValueAtAltitude, Cannot interpolate as the PRATMO climatology is not loaded.");
	}
	else
	{
		m_NO2density.XSlice( 0, &h );
		m_NO2density.XSlice( 1, &v );

		start  = h.begin();
		finish = h.end();
		iter2 = std::upper_bound( start, finish, placeandtime.heightm );
		
		if      (iter2 == start ) ++iter2;
		else if (iter2 == finish) --iter2;

		idx   = (iter2 - start);
		h2    = h.At(idx); 
		v2    = v.At(idx);
		--idx;
		h1    = h.At(idx);
		v1    = v.At(idx); 
		f     = (placeandtime.heightm-h1)/(h2-h1);
		ok = ( f > -0.6 ) && ( f < 1.6);			// Allow a little extraploation freedom at the ends of the array (usually fixes 1/2 km issues
		if (ok) *value = (1-f)*v1 + f*v2;
		
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::UpdateCache		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

nxBOOL skClimatology_Pratmo::UpdateCache( const GEODETIC_INSTANT& placeandtime )
{
	bool	ok;
	
	ok =    ( m_cachemjd       == placeandtime.mjd)
		 && ( m_cachelatitude  == placeandtime.latitude)
		 && ( m_cachelongitude == placeandtime.longitude)
		 && CacheIsLoaded();

	if (!ok)
	{
		ok  = LoadProfileFromBinaryFile( placeandtime );				// load in from the binary cache
		SetCacheIsLoaded( ok );
		m_cachemjd       = placeandtime.mjd;
		m_cachelatitude  = placeandtime.latitude;
		m_cachelongitude = placeandtime.longitude;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::GetParameter			2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_Pratmo::GetParameter( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	bool ok = true;
	
	if (updatecache || !CacheIsLoaded() ) ok = UpdateCache( placeandtime );
	else             ok = CheckCache();
	if (ok)
	{
		if  ( species == SKCLIMATOLOGY_NO2_CM3  ) ok = InterpolateValueAtAltitude( placeandtime, value);
		else
		{
			*value = MissingValue();
			ok     = false; 
			nxLog::Record( NXLOG_WARNING, "skClimatology_Pratmo::GetParameter, Not a supported species of skClimatology_Pratmo." );
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::IsSupportedSpecies		2008-2-11*/
/** **/
/*---------------------------------------------------------------------------*/

nxBOOL skClimatology_Pratmo::IsSupportedSpecies(const CLIMATOLOGY_HANDLE& species )
{
	bool ok;

	ok =  (species == SKCLIMATOLOGY_NO2_CM3) == TRUE;		  

	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::DeepCopy		2008-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

//bool skClimatology_Pratmo::DeepCopy( const skClimatology_Pratmo& other)
//{
//	bool	ok;
//
//	ok           =       skClimatology::DeepCopy(other);
//	ok           = ok && m_NO2density.DeepCopy(other.m_NO2density);
//	m_path       = other.m_path;
//	m_timeIndex  = other.m_timeIndex;					// The index of the current time loaded in m_NO2Density
//	m_angleIndex = other.m_angleIndex;					// The index o fthe latitude loaded into m_NO2Density
//	m_sza        = other.m_sza;							// The solar zenith angle of the data loaded in m_NO2density.
//	m_isPM       = other.m_isPM;
//	return ok;
//}
//

/*-----------------------------------------------------------------------------
 *					skClimatology_Pratmo::CreateClone		2008-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

//bool skClimatology_Pratmo::CreateClone(skClimatology** userclone)	const
//{
//	skClimatology_Pratmo*	clone;
//	bool					ok;
//
//	clone = new skClimatology_Pratmo;
//	ok = (clone != NULL);
//	if (!ok)
//	{
//		nxLog::Record(NXLOG_WARNING,"skClimatology_Pratmo::CreateClone, Error allocating memory for pratmo climatology");
//	}
//	else
//	{
//		ok = clone->DeepCopy(*this);
//	}
//	if (ok)
//	{
//		clone->AddRef();
//		*userclone = clone;
//	}
//	else
//	{
//		*userclone = NULL;
//		nxLog::Record(NXLOG_WARNING,"skClimatology_Pratmo::CreateClone, Error creating and/or copying cloned instance of pratmo climatology. Thats a problem!");
//	}
//	return ok;
//}
