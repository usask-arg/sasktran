#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"
#include <boost/algorithm/string.hpp>



/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::ISKGeodetic_Stub_std		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKGeodetic_Stub_std::ISKGeodetic_Stub_std()
{
	MakeStringSetFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::~ISKGeodetic_Stub_std		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKGeodetic_Stub_std::~ISKGeodetic_Stub_std()
{
}
 

/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::FromGeodetic		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic_Stub_std::FromGeodetic( double latitude, double longitude,  double Height )
{
	m_geoid.FromGeodetic( latitude, longitude, Height );
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::FromGeocentric		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic_Stub_std::FromGeocentric( const nxVector& geocentric  )
{
	m_geoid.FromGeocentricVector( geocentric );
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::FromTangentPointLocation		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic_Stub_std::FromTangentPointLocation	( const nxVector& r, const nxVector& lookv )
{
	m_geoid.FromTangentPointLocation( r, lookv);
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::FromTangentAltitude		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic_Stub_std::FromTangentAltitude( double required_height, const nxVector& spacecraftlocation, const nxVector& boresightplane, nxVector* requiredlookvector)
{
	return m_geoid.FromTangentAltitude( required_height, spacecraftlocation, boresightplane, requiredlookvector);
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::GeodeticWest		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic_Stub_std::GeodeticWest( )
{
	nxVector	west;
	nxVector	south;
	nxVector	up;
	m_geoid.GetGeodeticWestSouthUp( &west, &south, &up);
	return west;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::GeodeticSouth		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic_Stub_std::GeodeticSouth( )
{
	nxVector	west;
	nxVector	south;
	nxVector	up;
	m_geoid.GetGeodeticWestSouthUp( &west, &south, &up);
	return south;

}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::GeodeticUp		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic_Stub_std::GeodeticUp( )
{
	nxVector	west;
	nxVector	south;
	nxVector	up;
	m_geoid.GetGeodeticWestSouthUp( &west, &south, &up);
	return up;

}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::Location		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic_Stub_std::Location( )
{
	return m_geoid.Location();
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::GeodeticLongitude		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

double ISKGeodetic_Stub_std::GeodeticLongitude( )
{
	return m_geoid.GeodeticLongitude();
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::GeodeticLatitude		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

double ISKGeodetic_Stub_std::GeodeticLatitude( )
{
	return m_geoid.GeodeticLatitude();
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::Height		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

double ISKGeodetic_Stub_std::Height( )
{
	return m_geoid.Height();
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic_Stub_std::GetShellHeightLocation		 2014- 5- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic_Stub_std::GetShellHeightLocation( double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint)
{
	return m_geoid.GetShellHeightLocation( H, observerposition, look, entrypoint, exitpoint);
}

nxVector ISKGeodetic_Stub_std::OsculatingSpheroidCenter()
{
	nxVector offset;
	double radius;
	m_geoid.GetOsculatingSpheroid(&radius, &offset);
	return offset;
}

double ISKGeodetic_Stub_std::OsculatingSpheroidRadius()
{
	nxVector offset;
	double radius;
	m_geoid.GetOsculatingSpheroid(&radius, &offset);
	return radius;
}



/*---------------------------------------------------------------------------
 *          ISKGeodetic_Stub_std::MakeStringSetFunctions          2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKGeodetic_Stub_std::MakeStringSetFunctions()
{
	AddSetStringFunction( "geoid",
		[&, this](const char *value )
		{
			bool ok = true;
			std::string name( value);
			boost::algorithm::to_upper(name);
			if      (name == "IAU1976") m_geoid.SelectGeoid( nxGeodetic::IAU1976 );
			else if (name == "GRS80")   m_geoid.SelectGeoid( nxGeodetic::GRS80 );
			else if (name == "MERIT83") m_geoid.SelectGeoid( nxGeodetic::MERIT83 );
			else if (name == "WGS84")   m_geoid.SelectGeoid( nxGeodetic::WGS84);
			else
			{
				nxLog::Record( NXLOG_WARNING,"ISKGEodetic: unsupported geoid %s", (const char*) name.c_str() );
				ok = false;
			}
			return ok;
		}
	);
	
	AddSetScalarFunction( "iau1976",
		[&, this](double )
		{
			 m_geoid.SelectGeoid( nxGeodetic::IAU1976 );
			 return true;
		}
	);

	AddSetScalarFunction( "grs80",
		[&, this](double )
		{
			 m_geoid.SelectGeoid( nxGeodetic::GRS80 );
			 return true;
		}
	);

	AddSetScalarFunction( "merit83",
		[&, this](double )
		{
			 m_geoid.SelectGeoid( nxGeodetic::MERIT83 );
			 return true;
		}
	);

	AddSetScalarFunction( "wgs84",
		[&, this](double )
		{
			 m_geoid.SelectGeoid( nxGeodetic::WGS84 );
			 return true;
		}
	);


}

