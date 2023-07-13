#include <skclimatology21.h>
#include <float.h>
#include "sknetcdf4.h"
//#include "skosirismoderadius.h"

/*-----------------------------------------------------------------------------
 *					skModeRadius_WorldMapEntry::skModeRadius_WorldMapEntry		2013-2-6*/
/** **/
/*---------------------------------------------------------------------------*/

skModeRadius_WorldMapEntry::skModeRadius_WorldMapEntry()
{
	m_missingvalue = -9999.0;
}

/*-----------------------------------------------------------------------------
 *					skModeRadius_WorldMapEntry::CheckEntry		2013-2-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_WorldMapEntry::CheckEntry	()
{
	bool	ok;

	ok =       m_worldmap.XSize() == m_latitude.size();
	ok = ok && m_worldmap.YSize() == m_longitude.size();
	ok = ok && m_worldmap.ZSize() == m_altitude.size();

	for (size_t i = 0; i < m_longitude.size(); i++)
	{
		if (m_longitude.at(i) < 0) m_longitude.at(i) += 360.0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_WorldMapEntry::ReleaseResources		2013-2-7*/
/** **/
/*---------------------------------------------------------------------------*/

void skModeRadius_WorldMapEntry::Clear()
{
	m_mjd = 0.0;
	m_latitude.clear();
	m_longitude.clear();
	m_altitude.clear();
	m_worldmap.erase();
}

/*-----------------------------------------------------------------------------
 *					skEcmwf_IO::InterpolateLatLongAllPressures		2005-6-23*/
/** Interpolates the primary value (Temperature or geopotential height) to the
 *	specified latitiude and longitude usning linear interpolation of the 4 lat/long
 *	points enclosing the required place. The interpolation is made at all
 *	of the pressure levels.
 **/
/*---------------------------------------------------------------------------*/


bool skModeRadius_WorldMapEntry::Profile_InterpolateLatLong	( double latitude, double longitude, std::vector<double>* profile)
{
	size_t						iy0, iy1; 
	size_t						ix0, ix1;
	double						y0,y1,y;
	double						x0,x1,x;
	double						vs[4];
	bool						ok;
	bool						ok1;
	size_t						pidx;
	size_t						nlevels;


	if (longitude < 0.0) longitude += 360.0;							// get Longitude in range 0-360.0;
	NXASSERT(( (longitude >= 0.0)   && (longitude <= 360.0) ));			// Check input variables are in range
	NXASSERT(( (latitude  >= -90.0) && (latitude  <= 90.0)  ));

	y   = latitude;
	x   = nxmath::inrange(longitude, 360.0);

	nlevels = m_altitude.size();
	profile->resize( nlevels );

	ok      =       nxLinearInterpolate::FindBoundingIndicesAscending      ( m_latitude,  y,          &iy0, &iy1, &y0, &y1 );		// Interpolate latitude
	ok      = ok && nxLinearInterpolate::FindBoundingIndicesAscendingCyclic( m_longitude, x,  360.0,  &ix0, &ix1, &x0, &x1 );		// Interpolate longitude
	ok      = ok && (profile->size() == nlevels );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"VariableProfile_InterpolateLatLong, There was an error interpolating in latitude and longitude. Thats rather odd");
	}
	else
	{
		for (pidx = 0; pidx < nlevels; pidx++)
		{
			vs[0] = m_worldmap.At( iy0, ix0, pidx );
			vs[1] = m_worldmap.At( iy1, ix0, pidx );
			vs[2] = m_worldmap.At( iy1, ix1, pidx );
			vs[3] = m_worldmap.At( iy0, ix1, pidx );

			ok1 =	   ( NXFINITE(vs[0]) && (vs[0] > 0.0) )
					&& ( NXFINITE(vs[1]) && (vs[1] > 0.0) )
					&& ( NXFINITE(vs[2]) && (vs[2] > 0.0) )
					&& ( NXFINITE(vs[3]) && (vs[3] > 0.0) );
			if (!ok1)
			{
				profile->at(pidx) = m_missingvalue;
			}
			else
			{
				profile->at(pidx) = nxLinearInterpolate::FromSquare( x,y,x0,x1,y0,y1, vs);
			}
		}
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"VariableProfile_InterpolateLatLong, There was an error generating the requested profile. It is set to empty");
		profile->clear();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::LoadCacheFromFile		2013-2-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_WorldMapEntry::LoadFromFile( const char* filename )
{
	sknetcdf_IONetCDF				file;
	bool							ok;
	bool							ok1;
	const nxArrayLinear<double>*	dataptr;
	size_t							indexlo[4] = { NXARRAY_STARSELECT, NXARRAY_STARSELECT, NXARRAY_STARSELECT, 0};
	size_t							indexhi[4] = { NXARRAY_STARSELECT, NXARRAY_STARSELECT, NXARRAY_STARSELECT, 0};
	nx1dArray<double>				latitude;
	nx1dArray<double>				longitude;
	nx1dArray<double>				altitude;
	nx1dArray<double>				mjd;
	

	ok = file.LoadFile( filename, "mjd;latitude;longitude;altitude;ModeRadiusV6", -9999.0);
	if (ok)
	{
		dataptr = file.Field("ModeRadiusV6")->Data();
		ok = ok && mjd      .DeepCopy( *file.Field("mjd"      )->Data() );
		ok = ok && latitude .DeepCopy( *file.Field("latitude" )->Data() );
		ok = ok && longitude.DeepCopy( *file.Field("longitude")->Data() );
		ok = ok && altitude .DeepCopy( *file.Field("altitude" )->Data() );
		
		if (ok && (mjd.size() == 1) )
		{
			SetMjd      ( mjd.At(0) );
			SetAltitude ( altitude );
			SetLatitude ( latitude );
			SetLongitude( longitude );	
			indexlo[3] = 0;
			indexhi[3] = 0;
			ok1 = dataptr->Slice( indexlo, indexhi, 4, &m_worldmap );
			ok1 = ok1 && CheckEntry();
			ok = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"skClimatology_OsirisAerosolModeRadiusV600::LoadCacheFromFile, Error loading mode radius map from file %s ", (const char*)filename);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::skClimatology_OsirisAerosolModeRadiusV600		2013-2-6*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_OsirisAerosolModeRadiusV600::skClimatology_OsirisAerosolModeRadiusV600()
{
	ResetCurrentSettings();
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::~skClimatology_OsirisAerosolModeRadiusV600		2013-2-6*/
/** **/
/*---------------------------------------------------------------------------*/

skClimatology_OsirisAerosolModeRadiusV600::~skClimatology_OsirisAerosolModeRadiusV600()
{
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::ResetCurrentSettings		2013-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

void skClimatology_OsirisAerosolModeRadiusV600::ResetCurrentSettings()
{
	m_currentlatitude  = -9999.0;
	m_currentlongitude = -9999.0;
	m_currentmjd       = -9999.0;
	m_isAscendingNode  = false;
	m_isAscendingSet   = false;
}

	
/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::LoadCacheForMjd		2013-2-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::LoadCacheForMjd( double usermjd )
{
	bool		ok;
	nxString	beforename;
	nxString	aftername;

	ok =    ( m_beforentry.Mjd() >  40000.0 )
		 && ( m_beforentry.Mjd() <= usermjd )
		 && ( m_afterentry.Mjd() >  40000.0 )
		 && ( m_afterentry.Mjd() >= usermjd );

	if (!ok)
	{
		ok = m_isAscendingSet;
		if (!ok) nxLog::Record(NXLOG_WARNING,"skClimatology_OsirisAerosolModeRadiusV600::LoadCacheForMjd, You haven't set the orbital track to load - using descending" );
		ok = ok && m_filelocator.LocateBoundingModeRadiusFiles( usermjd, &beforename, &aftername, m_isAscendingNode );
		ok = ok && m_beforentry.LoadFromFile( beforename );
		ok = ok && m_afterentry.LoadFromFile( aftername );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skClimatology_OsirisAerosolModeRadiusV600::LoadCacheForMjd, There were errors loading the mode radius climatology files straddling mjd %f", (double)usermjd );
			m_beforentry.Clear();
			m_afterentry.Clear();
		}
		ResetCurrentSettings();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::CheckCurrentModeRadiusProfile		2013-2-7*/
/** Fetches the interpolated profile of mode radius at the specified
 *	latitude, longitude and time and updates the current cached profile
 **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::CheckCurrentModeRadiusProfile( double usermjd, double latitude, double longitude )
{
	bool	ok;
	std::vector<double>	y0;
	std::vector<double>	y1;
	double				t;
	double				t0;
	double				t1;
	double				dt;
	double				f0;
	double				f1;
	size_t				i;

	ok =    ( latitude  == m_currentlatitude )
		 && ( longitude == m_currentlongitude)
		 && ( usermjd   == m_currentmjd);
	if (!ok)
	{

		ok    = LoadCacheForMjd( usermjd );
		ok    = ok && m_beforentry.Profile_InterpolateLatLong( latitude, longitude, &y0);
		ok    = ok && m_afterentry.Profile_InterpolateLatLong( latitude, longitude, &y1);
		if (ok)
		{
			m_currentprofile.resize( y0.size() );
			t   = usermjd;
			t0  = m_beforentry.Mjd();
			t1  = m_afterentry.Mjd();
			dt  = t1-t0;
			f0  = (t1-t)/dt;
			f1  = (t-t0)/dt;
			for (i = 0; i < y0.size(); i++)
			{
				m_currentprofile.at(i) = y0[i]*f0 + y1[i]*f1;
			}
			m_currentlatitude  = latitude;
			m_currentlongitude = longitude;
			m_currentmjd       = usermjd;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skClimatology_OsirisAerosolModeRadiusV600::FetchModeRadiusProfile, There were errors fetching the mode radius profile");
			ResetCurrentSettings();
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::CacheIsValid		2013-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::CacheIsValid()
{
	return (m_currentmjd > 0.0) && (m_currentprofile.size() > 1 );
}

/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::UpdateCache		2013-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::UpdateCache( const GEODETIC_INSTANT& placeandtime)
{
	bool	ok;
	ok = CheckCurrentModeRadiusProfile( placeandtime.mjd, placeandtime.latitude, placeandtime.longitude );
	SetCacheIsLoaded(ok);					// Notify the base class that a cache is loaded or not
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::GetParameter		2013-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::GetParameter( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache)
{
	bool	ok = false;

	if (species == SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH)
	{
		*value = 1.6;
		ok     = true;
	}
	else if (species == SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS)
	{
		if (updatecache) UpdateCache(placeandtime);
		ok = CacheIsValid();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skClimatology_OsirisAerosolModeRadiusV600::GetParameter, You cannot call GetParameter until you have successfully called UpdateCache");
		}
		else
		{
			*value = nxLinearInterpolate::EvaluateYatX(	placeandtime.heightm/1000.0,
														m_beforentry.Altitudes(),
														m_currentprofile,
														nxLinearInterpolate::ENUM_TRUNCATE,
														-9999.0);
			ok = (*value > 0 );
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "skClimatology_OsirisAerosolModeRadiusV600::GetParameter, The mode radius at altitude %f meters was invalid", (double)placeandtime.heightm);
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::IsSupportedSpecies		2013-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::IsSupportedSpecies( const CLIMATOLOGY_HANDLE& species )
{
	bool	ok;

	ok  =    (species == SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS)
		  || (species == SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH);
		  
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skClimatology_OsirisAerosolModeRadiusV600::CreateClone		2013-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool skClimatology_OsirisAerosolModeRadiusV600::CreateClone( skClimatology** clone)
{
	bool	ok;

	*clone = new skClimatology_OsirisAerosolModeRadiusV600;
	ok = (*clone != NULL);
	if (ok)
	{
		(*clone)->AddRef();
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skClimatology_OsirisAerosolModeRadiusV600::CreateClone, Error creating clone");
	}
	return ok;
}

void skClimatology_OsirisAerosolModeRadiusV600::SetIsAscendingNode( bool isAscendingNode)
{
	m_isAscendingNode = isAscendingNode;
	m_isAscendingSet  = true;
	//nxLog::Record(NXLOG_WARNING, "skClimatology_OsirisAerosolModeRadiusV600::SetIsAscendingNode, ascending node is %s", (m_isAscendingNode) ? "true" : "false");
}
