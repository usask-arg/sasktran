/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_geodesy.h"
using namespace::nxmath;

/*	Planetary Positions
 *	-------------------
 *	Algorithms unashamedly copied from "Astronomy on the Personal Computer"
 *	O. Montenbruck and T. Pleger, Springer Verlag. QB51.3 E43 M6613 1991
 *	ISBN 3-540-52754-0 or ISBN 0-387-52754-0.
 *
 *	Their code is in PASCAL and I have ported that code to a C++ environment.
 *	I have checked the Sun algorithms against the Astronomical Almanac at
 *	randomly chosen times.  It was always with 1 arc second of the Almanac's
 *	value. Pretty good huh?
 *
 *	PlanetaryBody Useful Methods
 *	----------------------------
 *
 *	UpdateECI	        Updates the bodies ECI coordinates at the given epoch.
 *	Nutate			Nutate the current coords (assume current are wrt mean equinox of date).
 *	asEcliptic		Return body coordinates in Ecliptic X,Y,Z
 *	asRADEC			Return body coordinates in RA and DEC.
 *	asECI   		Return body coordinates in Earth Centred Inertial, wrt mean equinox of date.
 *
 *	PlanetaryBody Internal Methods
 *	------------------------------
 *	FRAC			Fractional part of number
 *	SN			Sine in degrees
 *	CS			Cosine in degrees
 *	TN			Tangent in degrees
 *	ASN			Arc Sine in degrees
 *	ACS			Arc Cosine in degrees.
 *	ATN2			Arc Tangent in degrees.
 *	ATN			Arc tangent in degrees
 *	ADDTHE			Add cos terms
 *	SINE			Sine in degrees
*/


//---------------------------------------------------------------------------
//			Table logging difference in seconds between
//			nxTRUE Dynamical Time and Universal Time
//			see pages K8-K9 of Astronomcal Almanac 1994
//
// HISTORY:
// 11-Feb-1999	NDL Updated table to exact values for 1995 using AA 1997 (page K9)
//
//---------------------------------------------------------------------------

const  int    TDTStartYear = 1620;						// The TDT DT table begins in 1620.
const  double TDTTable [] = 							// Table to convert UT to TDT, page K8-K9
       {												// Astronomical Almanac 1994
	 124.0, 119.0, 115.0, 110.0, 106.0, 102.0,  98.0,  95.0,  91.0,  88.0,	// 1620-1629
	  85.0,  82.0,  79.0,  77.0,  74.0,  72.0,  70.0,  67.0,  65.0,  63.0,	// 1630-1639
	  62.0,  60.0,  58.0,  57.0,  55.0,  54.0,  53.0,  51.0,  50.0,  49.0,	// 1640-1649
	  48.0,  47.0,  46.0,  45.0,  44.0,  43.0,  42.0,  41.0,  40.0,  38.0,	// 1650-1659
	  37.0,  36.0,  35.0,  34.0,  33.0,  32.0,  31.0,  30.0,  28.0,  27.0,	// 1660-1669
	  26.0,  25.0,  24.0,  23.0,  22.0,  21.0,  20.0,  19.0,  18.0,  17.0,	// 1670-1679
	  16.0,  15.0,  14.0,  14.0,  13.0,  12.0,  12.0,  11.0,  11.0,  10.0,	// 1680-1689
	  10.0,  10.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,	// 1690-1699
	   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,  10.0,  10.0,	// 1700-1709
	  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  10.0,  11.0,  11.0,  11.0,	// 1710-1719
	  11.0,  11.0,  11.0,  11.0,  11.0,  11.0,  11.0,  11.0,  11.0,  11.0,	// 1720-1729
	  11.0,  11.0,  11.0,  11.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,	// 1730-1739
	  12.0,  12.0,  12.0,  12.0,  13.0,  13.0,  13.0,  13.0,  13.0,  13.0,	// 1740-1749
	  13.0,  14.0,  14.0,  14.0,  14.0,  14.0,  14.0,  14.0,  15.0,  15.0,	// 1750-1759
	  15.0,  15.0,  15.0,  15.0,  15.0,  16.0,  16.0,  16.0,  16.0,  16.0,	// 1760-1769
	  16.0,  16.0,  16.0,  16.0,  16.0,  17.0,  17.0,  17.0,  17.0,  17.0,	// 1770-1779
	  17.0,  17.0,  17.0,  17.0,  17.0,  17.0,  17.0,  17.0,  17.0,  17.0,	// 1780-1789
	  17.0,  17.0,  16.0,  16.0,  16.0,  16.0,  15.0,  15.0,  14.0,  14.0,	// 1790-1799
	  13.7,  13.4,  13.1,  12.9,  12.7,  12.6,  12.5,  12.5,  12.5,  12.5,	// 1800-1809
	  12.5,  12.5,  12.5,  12.5,  12.5,  12.5,  12.5,  12.4,  12.3,  12.2,	// 1810-1819
	  12.0,  11.7,  11.4,  11.1,  10.6,  10.2,   9.6,   9.1,   8.6,   8.0,	// 1820-1829
	   7.5,   7.0,   6.6,   6.3,   6.0,   5.8,   5.7,   5.6,   5.6,   5.6,	// 1830-1839
	   5.7,   5.8,   5.9,   6.1,   6.2,   6.3,   6.5,   6.6,   6.8,   6.9,	// 1840-1849
	   7.1,   7.2,   7.3,   7.4,   7.5,   7.6,   7.7,   7.7,   7.8,   7.8,	// 1850-1859
	   7.88,  7.82,  7.54,  6.97,  6.40,  6.02,  5.41,  4.10,  2.92,  1.82,	// 1860-1869
	   1.61,  0.10, -1.02, -1.28, -2.69, -3.24, -3.64, -4.54, -4.71, -5.11,	// 1870-1879
	  -5.40, -5.42, -5.20, -5.46, -5.46, -5.79, -5.63, -5.64, -5.80, -5.66,	// 1880-1889
	  -5.87, -6.01, -6.19, -6.64, -6.44, -6.47, -6.09, -5.76, -4.66, -3.74,	// 1890-1899
	  -2.72, -1.54, -0.02,  1.24,  2.64,  3.86,  5.37,  6.14,  7.75,  9.13,	// 1900-1909
	  10.46, 11.53, 13.36, 14.65, 16.01, 17.20, 18.24, 19.06, 20.25, 20.95,	// 1910-1919
	  21.16, 22.25, 22.41, 23.03, 23.49, 23.62, 23.86, 24.49, 24.34, 24.08,	// 1920-1929
	  24.02, 24.00, 23.87, 23.95, 23.86, 23.93, 23.73, 23.92, 23.96, 24.02,	// 1930-1939
	  24.33, 24.83, 25.30, 25.70, 26.24, 26.77, 27.28, 27.78, 28.25, 28.71,	// 1940-1949
	  29.15, 29.57, 29.97, 30.36, 30.72, 31.07, 31.35, 31.68, 32.18, 32.68,	// 1950-1959
	  33.15, 33.59, 34.00, 34.47, 35.03, 35.73, 36.54, 37.43, 38.29, 39.20,	// 1960-1969
	  40.18, 41.17, 42.23, 43.37, 44.49, 45.48, 46.46, 47.52, 48.53, 49.59,	// 1970-1979
	  50.54, 51.38, 52.17, 52.96, 53.79, 54.34, 54.87, 55.32, 55.82, 56.30,	// 1980-1989
	  56.86, 57.57, 58.31, 59.12, 59.98, 60.78, 61.63, 62.29, 62.97, 63.47,	// 1990-1999
	  63.83, 64.09, 64.  , 65.  , 66.  , 66.  , 67. 						// 2000-2001 exact, 2002 onwards extrapolated
	};

//-------------------------------------------------------------------------------
//			PlanetaryBody::Constructor
//-------------------------------------------------------------------------------

PlanetaryBody::PlanetaryBody()
{
   m_DistantObject = 1;
   m_time  = -1.0;
}

//-------------------------------------------------------------------------------
//			PlanetaryBody::FRAC
//-------------------------------------------------------------------------------

double PlanetaryBody::FRAC( double X)
{
   return (X-floor(X));
}

//-------------------------------------------------------------------------------
//			PlanetaryBody::SN
//-------------------------------------------------------------------------------

double PlanetaryBody::SN( double X )
{
   return sin(X*ONE_DEGREE);
}

//-------------------------------------------------------------------------------
//			PlanetaryBody::CS
//-------------------------------------------------------------------------------

double PlanetaryBody::CS( double X )
{
   return cos(X*ONE_DEGREE);
}

//--------------------------------------------------------------------------------
//			PlanetaryBody::TN
//--------------------------------------------------------------------------------

double PlanetaryBody::TN( double X )
{
   return tan(X*ONE_DEGREE);
}

//--------------------------------------------------------------------------------
//			PlanetaryBody::ASN
//--------------------------------------------------------------------------------

double PlanetaryBody::ASN( double X )
{
   return asin( X )/ONE_DEGREE;
}

//---------------------------------------------------------------------------------
//			PlanetaryBody::ACS
//---------------------------------------------------------------------------------

double PlanetaryBody::ACS(double X)
{
   return acos(X)/ONE_DEGREE;
}

//--------------------------------------------------------------------------------
//			PlanetaryBody::ATN2
//--------------------------------------------------------------------------------

double PlanetaryBody::ATN2( double Y, double X )
{
   return atan2(Y,X)/ONE_DEGREE;
}

//--------------------------------------------------------------------------------
//			PlanetaryBody::ATN
//--------------------------------------------------------------------------------

double PlanetaryBody::ATN( double X )
{
   return atan(X)/ONE_DEGREE;
}

//-------------------------------------------------------------------------------
//			PlanetaryBody::ADDTHE
//-------------------------------------------------------------------------------

void PlanetaryBody::ADDTHE( double C1, double S1, double C2, double S2,
			    double &C, double &S)
{
   C = C1*C2 - S1*S2;
   S = S1*C2 + C1*S2;
}


//--------------------------------------------------------------------------------
//			PlanetaryBody::SINE
//--------------------------------------------------------------------------------


double PlanetaryBody::SINE ( double PHI)			// calculate sin(phi); phi in units of 1 revolution = 360 degrees
{
   return sin(TWOPI*FRAC(PHI));
}


//---------------------------------------------------------------------------
//						PlanetaryBody::GAST
//	Calculate the Greenwich Apparent Sideral Time (GAST) at a given UTC.
//	This requires the equation of the equinoxes and is defined in Astronomical
//	Almanac 1997, Page B6.  I ignore effects due to the ascending node of the
//	Moon < 0.003 of an arc second
//---------------------------------------------------------------------------

double PlanetaryBody::GAST( nxTimeStamp utc )
{
	double		EPS;
	double		DPSI,DEPS, equequinox;
	nxTimeStamp	tdt;
	double		gast;

	tdt = TDT( utc);							// convert the universal time to TERRESTIAL DYNAMICAL TIME
	EPS    = Ecliptic( tdt, nxFALSE );			// get the mean obliquity of the ecliptic in degrees.
	Nutation( tdt, &DPSI, &DEPS );				// get nutation in longitude and obliquity in degrees.
	equequinox =  DPSI*cosd(EPS);				// This is the equation of the equinox in "days"
	gast = utc.GMST() + (equequinox/360.0);		// GAST = GMST + equation of equinoxes.
	return gast;
}

//----------------------------------------------------------------------------
//			Nutation
//
//	Calculates Nutation effects necessary to calculate apparent sidereal
//	time.  Nutation is the short term periodic (18.6 years) movement
//	of the Earth's pole about the Earth's mean pole.  Nutation produces
//	variations of upto 17 arc-seconds in longitude and 9 arc-seconds
//	in latitude.
//
//	Nutation is used to convert from the coordinates using the <mean
//	equinox of date> to the <(nxTRUE) equinox of date>.
//
//	DPSI is returned as the total nutation in longitude in degrees.
//	DPES is returned as the total nutation in obliquity in degrees.
//
//	The last calculation is always cached as the code may make multiple
//	calls for a given time.
//
//----------------------------------------------------------------------------

void PlanetaryBody::Nutation( const nxTimeStamp &tdt, double *DPSI, double *DEPS )
{
   double T;
   double LS,D,F,N;
   static nxTimeStamp	last_tdt(-9999.0);
   static double	last_dpsi  = 0.0;
   static double	last_deps  = 0.0;

	if ( tdt == last_tdt)
	{
		*DPSI = last_dpsi;
		*DEPS = last_deps;
	}
	else
	{
		T   = tdt.JD2000Centuries();
		LS  = TWOPI*(0.993133+  99.997306*T); 		// mean anomaly Sun
		D   = TWOPI*(0.827362+1236.853087*T); 		// diff. longitude Moon-Sun
		F   = TWOPI*(0.259089+1342.227826*T); 		// mean argument of latitude
		N   = TWOPI*(0.347346-   5.372447*T); 		// longit. ascending node

		*DPSI = ( -17.200*sin(N)   - 1.319*sin(2*(F-D+N)) - 0.227*sin(2*(F+N))
				+ 0.206*sin(2*N) + 0.143*sin(LS) )/3600.0;

		*DEPS = ( + 9.203*cos(N)   + 0.574*cos(2*(F-D+N)) + 0.098*cos(2*(F+N))
				- 0.090*cos(2*N))/3600.0;

		last_dpsi = *DPSI;
		last_deps = *DEPS;
		last_tdt  = tdt;
	}
}


//------------------------------------------------------------------------------
//			PlanetaryBody::Ecliptic
//	returns the obliquity of the ecliptic in degrees at nxTRUE Dynamical Time
//	tdt.
//	The caller can either request the mean obliquity of date or the
//	nxTRUE obliquity of date (ie corrected for Nutation).
//
//------------------------------------------------------------------------------

double PlanetaryBody::Ecliptic( const nxTimeStamp &tdt, nxBOOL UsenxTRUEEclipticOfDate )
{
   double ecliptic;
   double T;
   static nxTimeStamp	last_tdt(-9999.0);
   static double	last_ecliptic  = 0.0;

	if ( tdt == last_tdt)
	{
		ecliptic = last_ecliptic;
	}
	else
	{
		T        = tdt.JD2000Centuries();									// get the time in Julian Centru ies
		ecliptic = (23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0);	// work out the mean Ecliptic of date.
		last_ecliptic = ecliptic;
		last_tdt      = tdt;
	}
	if (UsenxTRUEEclipticOfDate)							// if caller wants the nxTRUE obliquity
	{														// then
		double DeltaEcliptic, DeltaLongitude;				// we must do a
		Nutation( tdt, &DeltaLongitude, &DeltaEcliptic );	// nutation calculation, returns radians!
		ecliptic += DeltaEcliptic;							// and add this to the current estimate of the ecliptic (convert to degrees).
	}
   return ecliptic;
}


//----------------------------------------------------------------------------
//	PlanetaryBody::TDT
//
//	Returns Terrestial Dynamical Time at the specified ut.
//	time of interest.  TDT is independent of Earths rotation, unlike
//	Universal Time.
//
//	Terrestial Dynamical time is used for calculations of the Sun's orbit etc.
//	Offsets for TDT are only known after the fact, hence future predictions
//	are somewhat unreliable.
//
//	I look up TDT to the nearest YEAR.  It is only accurate to a second or so
//	as I dont interpolate between days in the year.  This might change in the
//	future if needed.
//
//	This code is not thread-safe!!!!
//
//-----------------------------------------------------------------------------

nxTimeStamp PlanetaryBody::TDT( const nxTimeStamp &utc)
{
	int		day, month, year, hour, mins, secs;
	double	ticks;
	double	dt;
	double  t;
	int		TDTEndYear;
	static	nxTimeStamp	last_utc(-9999.0);
	static	nxTimeStamp	last_tdt(-9999.0);

	if ( utc != last_utc )
	{
		TDTEndYear  =  TDTStartYear + sizeof(TDTTable)/sizeof(double) - 1;		// Get the last year in the TDT table
		utc.GetUTC( &day, &month, &year, &hour, &mins, &secs, &ticks );					// Get the year of this time.

		nxBOOL ok = (year >= TDTStartYear) && (year <= TDTEndYear  );
		if (!ok)
		{
			if (year < TDTStartYear)
			{
				dt = 0;
				ok = true;
			}
			else
			{
				ok = (year < 2051);
				if (!ok)
				{
					static bool firsttime  = true;
					if ( firsttime )
					{
						nxLog::Record( NXLOG_WARNING, "PlanetaryBody::TDT, Requested year (%d) is out of range", (int)year);
						firsttime = false;
					}
				}
				t = year-2000;
				dt = 62.92 + 0.32217 *t + 0.005589 * t*t;			// Good for the range 2005-2050
			}
		}
		else
		{
			if (year > TDTEndYear) year = TDTEndYear;
			int index;

			index = year - TDTStartYear;
			dt    = TDTTable[index];
		}
		last_tdt = utc + (dt/86400.0);
		last_utc = utc;
	}
    return last_tdt;
}


//------------------------------------------------------------------------------
// 			PlanetaryBody::NutateEquatorialCoords
//
//	Assumes the current "m_location" at instant "time" is with
//	reference to the mean equinox of date.  Adjusts the coordinates to
//	reflect the nxTRUE equinox of date, by account for the short-term
//	periodic nutation of the Earth's pole.
//	I.E. applies the equation of the equinoxes etc.
//	Follows O. Montenbruck algorithm NUTEQU.
//
//	Nutation is a short term (18.6 years) periodic variation of the
//	Earth's nxTRUE pole about the mean pole.  The long term
//	movement of the Earth's pole is precession
//
//	Also see page B20 Astronomical Almanac 1990 for details.
//	Nutation should only be applied to Equatorial coordinates.
//
//      including terms >0.1" according to IAU 1980
//
//-------------------------------------------------------------------------------

void PlanetaryBody::NutateEquatorialCoords( nxVector* location, nxTimeStamp tnow )
{
   double EPS;
   double DPSI,DEPS,C,S;
   double DX,DY,DZ;

   EPS    = Ecliptic( tnow, nxFALSE );		// get the mean obliquity of the ecliptic in degrees.
   Nutation( tnow, &DPSI, &DEPS );			// get nutation in longitude and obliquity in degrees.

   DPSI *= ONE_DEGREE;
   DEPS *= ONE_DEGREE;
   C = DPSI*cosd(EPS);
   S = DPSI*sind(EPS);

   DX = location->X() - (C*location->Y() + S*location->Z());
   DY = location->Y() + (C*location->X() - DEPS*location->Z());
   DZ = location->Z() + (S*location->X() + DEPS*location->Y());

   location->SetCoords( DX, DY, DZ );
}

//----------------------------------------------------------------------------
//			PlanetaryBody::ConvertToEclipticCoords
//
//	Returns the ecliptic (x,y,z) coordinates of the planetary body.
//	Assumes the coordinates were previously in equatorial coordinates.
//
//	x is from the geocentric earth toward the mean vernal equinox
//	z is perpendicular to ecliptic (in northward direction)
//	y is perpendicular to x and z.
//
//
// based on EQUECL function in O MontenBruck.
//
// EQUECL: Conversion of equatorial into ecliptic coordinates
//         (T: equinox in Julian centuries since J2000)
//	   Note X coordinate is unchanged.
//-----------------------------------------------------------------------------

void PlanetaryBody::ConvertToEclipticCoords( nxBOOL UseTrueEclipticOfDate)
{
   m_location.RotateAboutXaxis( Ecliptic( m_time, UseTrueEclipticOfDate) );
}


//----------------------------------------------------------------------------
//			PlanetaryBody::ConvertToEquatorialCoords
//
//	Converts the current "m_location" vector at instant "time" from Ecliptic
//	coordinates to equatorial coordinates.
//
//	There are two possible conversion
//
// 	based on O. Montenbruck ECLEQU.
// 	ECLEQU: Conversion of equatorial into ecliptic coordinates
//           (T: equinox in Julian centuries since J2000)
//	     Note X coordinate is unchanged.
//-----------------------------------------------------------------------------


void PlanetaryBody::ConvertToEquatorialCoords( nxBOOL UseTrueEclipticOfDate )
{
   m_location.RotateAboutXaxis( -Ecliptic(m_time, UseTrueEclipticOfDate));
}

//----------------------------------------------------------------------------
//			PlanetaryBody::ApparentECIPosition
//
//	Calculates the apparent ECI position of a planetary body
//	from an observer.
//
//----------------------------------------------------------------------------

nxVector PlanetaryBody::ApparentECIPosition( const nxGeodetic &Site)
{
   nxVector ApparentPosition;
   nxVector SiteCoords;
   double gst;

   if (!m_DistantObject)								// If this is not a distant object
   {													// then account for Earths parallax.
      gst = m_time.GMST()*360.0;						// angle between Geographic X axis and vernal equinox.
      SiteCoords = Site.Location();						// Copy geographic Geocentric m_location of observer.
      SiteCoords.TransformToNewPole( -gst, 90.0 );		// Get Sitecoords in Equatorial coordinate system.
      ApparentPosition = m_location - SiteCoords;		// Get apparent position of planetary body from observer.
   }													// and that is that
   else													// otherwise
   {
      ApparentPosition = m_location;
   }
   return ApparentPosition;
}


//----------------------------------------------------------------------------
//			PlanetaryBody::Topocentric
//
//	Calculates the topocentric coordinates of the planetary body
//	given the observers m_location and the current time.
//	Topcentric coordinates are represented by a right angled coordinate system
//	centred on the observer.
//
//	Coordinates of returned coordinates are:-
//
//	X = +ve due south 	meters
//	Y = +ve due east	meters
//	Z = +ve due up.		meters.
//
//----------------------------------------------------------------------------

nxVector PlanetaryBody::Topocentric( const nxGeodetic &Site )
{

   nxVector ApparentPosition;
   double lst;
   double gst;

   gst = m_time.GMST()*360.0;				// angle between Geographic X axis and vernal equinox.
   lst = gst + Site.GeodeticLongitude();			// angle between Site and Vernal equinox.
   ApparentPosition = ApparentECIPosition( Site );	// Get the apparent postion of this object.
   ApparentPosition.TransformToNewPole( lst, Site.GeodeticLatitude() );
   return ApparentPosition;
}

//----------------------------------------------------------------------------
//			PlanetaryBody::RiseSet
//
//	Calculates the time of the rise and set of an astronomical body
//	after time Tnow.  Rise and set will both occur after Tnow if the
//	body rises and sets. If the body never rises or sets then rise or set
//	will be before Tnow depending upon condition.
//	This is not nxTRUE for the SUN or the MOON which change RA and DEC
//	throught the day.
//
//	I currently use the geodetic latitude and longitude.
//
//	During Daytime  a body will set before it rises. (SET < RISE)
//	During darkness a body will rise before it sets. (RISE < SET)
//
//	If the body never rises then RISE = 00:00:00 UT,  SET  = RISE + 2;
//	If the body never sets  then SET  = 00:00:00 UT,  RISE = SET  + 2;
//
//	Returns nxTRUE  if body rises and sets.
//	Returns nxFALSE if body never rises or never sets.
//
//----------------------------------------------------------------------------


nxBOOL PlanetaryBody::RiseSet( const nxTimeStamp &Tnow, nxTimeStamp* rise, nxTimeStamp* sets,
			      const nxGeodetic  &Place, double zenith)
{
   double		ra;					// Right ascension of body.
   double		dec;				// declination of body.
   double		gmst;				// Greenwich Mean sidereal time of transit across meridian.
   double		lat;				// Latitude of observer.
   double		h;					// Hour angle of rise/set
   nxBOOL		status;
   double		theta;
   nxVector		ApparentECI;

	UpdateECIPosition( Tnow );
	ApparentECI = ApparentECIPosition( Place );						// Get the apparent coordinates from this place.
	ra    = ApparentECI.Longitude();								// get the Right Ascension
	dec   = ApparentECI.Latitude();									// and declination of this body.
	lat   = Place.GeodeticLatitude();
	theta =  (cosd(zenith) - sind(dec)*sind(lat))/( cosd(dec)*cosd(lat) );

	if (fabs(theta) >= 1.0)											// if the object never rises or never sets
	{																// then
		status = nxFALSE;											// return nxFALSE.
		if (theta <  1.0)											// if the object is always above the horizon.
		{															// then return
			*sets = m_time.ZeroUT();								// make the sun set before it
			*rise = *sets + 1.0;									// rises (which means it daytime now)
		}
		else														// if the is always below horizon
		{															// so
			*rise = m_time.ZeroUT();								//make the sun appear to rise
			*sets = *rise + 1.0;									// before it sets (which means its dark now)
		}
	}
	else															// object does rise and set
	{																// so
		status = nxTRUE;											// return nxTRUE  status
		gmst = (ra - Place.GeodeticLongitude())/360.0;				// GMST of transit across
		gmst = inrange( gmst, 1.0);									// the observers meridian.
		h    = acosd( theta )/360.0;								// hour angle of rise/set phenomena in days.
		*rise = m_time.fromGMST( gmst - h );						// Get rise time
		*sets = m_time.fromGMST( gmst + h );						// Get set  time
		while (*rise < m_time) *rise = *rise + nxTimeStamp::ONESIDEREALDAY; 		// Make sure rise time is bigger than current time.
		while (*sets < m_time) *sets = *sets + nxTimeStamp::ONESIDEREALDAY;		// Make sure set  time is bigger than current time.
	}
	return status;
}

//----------------------------------------------------------------------------
//			PlanetaryBody::BelowHorizon
//
//	Returns nxTRUE if a planetary Body is below a horizon and nxFALSE if it
//	is above.
//
//----------------------------------------------------------------------------

nxBOOL PlanetaryBody::BelowHorizon( const nxTimeStamp &Tnow, const nxGeodetic &Site, double zenith)
{
   nxTimeStamp Rise,Set;

   UpdateECIPosition( Tnow );
   RiseSet( Tnow, &Rise, &Set, Site, zenith );
   return ( Rise < Set );
}


//----------------------------------------------------------------------------
//			PlanetaryBody::InContactWithGround
//	Returns nxTRUE if the planetarybody is in contact with the specified site at the
//	specified time.
//----------------------------------------------------------------------------

nxBOOL PlanetaryBody::InContactWithGround( nxTimeStamp &Tnow, const nxGeodetic &Site, double minelevation)
{
   nxBOOL   AboveHorizon;
   nxVector locationAtSite;

   UpdateECIPosition( Tnow );
   locationAtSite = Topocentric( Site );
   AboveHorizon   = locationAtSite.Latitude() > minelevation;
   return AboveHorizon;
}

//----------------------------------------------------------------------------
//			TimeofContactChange
//	Determines the time of contact change to 0.1 seconds accuracy
//	by recursively calling itself.  This may give stack problems
//	but generally it seems to work just fine.
//
//	function returns nxTRUE if a contact change was found between times
//	TStart and Tend.
//----------------------------------------------------------------------------

nxBOOL PlanetaryBody::TimeOfContactChange( nxTimeStamp TStart, nxTimeStamp Tend, nxTimeStamp *ChangeTime,
			                 const nxGeodetic &Site, double minelevation, double stepsize )
{
   nxBOOL      InContact;
   nxBOOL      InContactNow;
   nxTimeStamp Tnow = TStart;
   nxTimeStamp Tlast= TStart;
   nxBOOL	     FoundChange = nxFALSE;

   *ChangeTime = TStart;					// place a default value into the start time
   InContact   = InContactWithGround( TStart, Site, minelevation );		// get the status of contact at the beginning
   Tnow        = Tnow + stepsize;				// and set up the next step

   while (( Tnow < Tend) && (!FoundChange))			// now scan forwards until
   {								// and look at the
      InContactNow = InContactWithGround(Tnow, Site, minelevation);		// status of contact with the ground
      FoundChange  = (InContactNow != InContact);		// and see if it has changed.
      if (!FoundChange)						// if it has not changed
      {								// then
	 Tlast = Tnow;						// the last time becomes the current time
	 Tnow  = Tnow + stepsize;				// and step forward looking for the change
      }								// and keep scanning  until
   }								// until a change is found or we have timed out

   if (FoundChange)						// if we have found a change in the coarse interval
   {								// then
      *ChangeTime = Tnow;					// the best estimate of the change time is the current time
	  if (stepsize > 0.1*nxTimeStamp::ONESECOND)				// now if we have not reached the ultimate resolution
      {								// then
	 stepsize    = stepsize/10.0;				// decrement the stepsize
	 Tnow        = Tnow + stepsize;				// add one step to the current time (to make sure it works ok)
	 FoundChange = TimeOfContactChange( Tlast, Tnow, ChangeTime, Site, minelevation, stepsize ); // and recurse to the next finer iteration, should always pass!
      }								// the Changetime is modified directly at each recursion
   }								// and after each level the status is returned
   return FoundChange;                  			// and then return status to caller.
}


//----------------------------------------------------------------------------
//			StartOfContactAtSite
//  Returns when a planetary body makes contact with a given site between times
//  Tnow and TEnd.  If the satellite is in contact at time Tnow then the
//  contact time is stepped back to the start of contact, otherwise it steps
//  forward to the next contcat time.
//
//  If no contact in the period specified then time 0.0 is returned.
//
//----------------------------------------------------------------------------

nxTimeStamp PlanetaryBody::StartOfContact( const nxGeodetic &Site, double minelevation,  nxTimeStamp Tnow, nxTimeStamp TEnd, nxTimeStamp *EndOfContact)
{
   nxTimeStamp BeginContact;
   nxTimeStamp T;
   nxBOOL      ok;

   while (InContactWithGround( Tnow, Site, minelevation))								// Now if we are already in contact with the ground
   {																					// then step back a little way until we are
	   Tnow = Tnow - nxTimeStamp::ONEMINUTE;															// out of contact with the ground
   }																					// so we never get partial passages analysed.
   if (EndOfContact != NULL) *EndOfContact = 0.0;										// default to no end of contact
   ok = TimeOfContactChange( Tnow, TEnd, &BeginContact, Site, minelevation, nxTimeStamp::ONEMINUTE);	// search for the start of contact with
   if (!ok)
   {
      BeginContact = 0.0;
   }
   else
   {
      if (EndOfContact != NULL)													// if we were in contact and user requests end of contact
      {																			// then
		  Tnow = BeginContact + nxTimeStamp::ONESECOND;										// step forward a little way
         while (InContactWithGround( Tnow, Site,minelevation ))					// and now search
         {																		// for the time when the
			 Tnow = Tnow + nxTimeStamp::ONEMINUTE;											// satellite loses contact with
         }																		// this station
         T    = Tnow - nxTimeStamp::ONEMINUTE;                                               // determine the longest time to search for end of contact
         Tnow = Tnow + (nxTimeStamp::ONEMINUTE/10.0);
         ok   = TimeOfContactChange( T, Tnow, EndOfContact, Site, minelevation, (nxTimeStamp::ONEMINUTE/10.0)); // and look for time when we lose contact
         if (!ok) *EndOfContact = 0.0;
      }
   }
   return BeginContact;
}
