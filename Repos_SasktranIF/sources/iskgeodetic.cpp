#include "sasktranif_internals.h"

/*-----------------------------------------------------------------------------
 *					ISKGeodetic::ISKGeodetic		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKGeodetic::ISKGeodetic()
{
	SasktranIF_ClassFactoryLocator	classfactory;

	classfactory.CreateISKGeodetic( "STANDARD", &m_geoid,DllNamePtr() );
}


/*-----------------------------------------------------------------------------
 *					ISKEngine::~ISKEngine		2014-2-8*/
/** Destroys the ISKEngine object and all of the associated resources.
**/
/*---------------------------------------------------------------------------*/

ISKGeodetic::~ISKGeodetic()
{
	if (m_geoid != nullptr) m_geoid->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::SetLocationLatLonAlt		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetLocationLatLonAlt( double latitude, double longitude,  double Height )
{
	return (m_geoid != nullptr) ? m_geoid->FromGeodetic( latitude, longitude, Height) : false;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::SetLocationXYZ		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetLocationXYZ( const nxVector& geocentric )
{
	return (m_geoid != nullptr) ? m_geoid->FromGeocentric( geocentric) : false;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::SetLocationFromTangentPoint		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetLocationFromTangentPoint( const nxVector& r, const nxVector& lookv )
{
	return (m_geoid != nullptr) ? m_geoid->FromTangentPointLocation( r, lookv) : false;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::FromTangentAltitude		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetLocationFromTangentAltitude( double required_height, const nxVector& spacecraftlocation, const nxVector& boresightplane, nxVector* requiredlookvector)
{
	return (m_geoid != nullptr) ? m_geoid->FromTangentAltitude( required_height,spacecraftlocation, boresightplane, requiredlookvector) : false;
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetLocalWest		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic::GetLocalWest( )
{
	return (m_geoid != nullptr) ? m_geoid->GeodeticWest() : nxVector( 0.0, 0.0, 0.0 );
}



/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetLocalSouth		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic::GetLocalSouth( )
{
	return (m_geoid != nullptr) ? m_geoid->GeodeticSouth() : nxVector( 0.0, 0.0, 0.0 );
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetLocalUp		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic::GetLocalUp( )
{
	return (m_geoid != nullptr) ? m_geoid->GeodeticUp() : nxVector( 0.0, 0.0, 0.0 );
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetLocationXYZ		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

nxVector ISKGeodetic::GetLocationXYZ( )
{
	return (m_geoid != nullptr) ? m_geoid->Location() : nxVector( 0.0, 0.0, 0.0 );
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetLongitude		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double ISKGeodetic::GetLongitude( )
{
	return (m_geoid != nullptr) ? m_geoid->GeodeticLongitude() : std::numeric_limits<double>::quiet_NaN(); 
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetLatitude		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double ISKGeodetic::GetLatitude( )
{
	return (m_geoid != nullptr) ? m_geoid->GeodeticLatitude() : std::numeric_limits<double>::quiet_NaN(); 
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetAlt		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double ISKGeodetic::GetAlt( )
{
	return (m_geoid != nullptr) ? m_geoid->Height() : std::numeric_limits<double>::quiet_NaN(); 
}


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::GetAltitudeIntercepts		 2014- 5- 6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::GetAltitudeIntercepts( double H, const nxVector& observerposition, const nxVector& look, nxVector* entrypoint, nxVector* exitpoint)
{
	return (m_geoid != nullptr) ? m_geoid->GetShellHeightLocation( H, observerposition, look, entrypoint, exitpoint) : false; 
}

 nxVector ISKGeodetic::GetOsculatingSpheroidCenter()
 {
	 return (m_geoid != nullptr) ? m_geoid->OsculatingSpheroidCenter() : nxVector(0.0, 0.0, 0.0);
 }

 double ISKGeodetic::GetOsculatingSpheroidRadius()
 {
	 return (m_geoid != nullptr) ? m_geoid->OsculatingSpheroidRadius() : std::numeric_limits<double>::quiet_NaN();
 }


/*-----------------------------------------------------------------------------
 *					ISKGeodetic::SetPropertyScalar		 2014- 5- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetPropertyScalar( const char* propertyname, double value )
{
	return (m_geoid != nullptr) ? m_geoid->SetPropertyScalar( propertyname, value ) : false;
}



/*-----------------------------------------------------------------------------
 *					ISKGeodetic::SetPropertyArray		 2015- 11- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	nxLog::Record(NXLOG_WARNING,"ISKGeodetic does not have any array based properties");
	return false;
}

/*-----------------------------------------------------------------------------
 *					ISKGeodetic::SetPropertyArray		 2015- 11- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKGeodetic::SetPropertyObject( const char* propertyname, ISKModuleBase* object)
{
	nxLog::Record(NXLOG_WARNING,"ISKGeodetic does not have any SasktranIF object based properties");
	return false;
}

bool ISKGeodetic::SetPropertyString( const char* propertyname, const char* str ) 
{
	nxLog::Record(NXLOG_WARNING,"ISKGeodetic does not have any string properties");
	return false;
}




