/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"
#include "nxbase_math.h"
#include <time.h>
#if  !defined(NX_WINDOWS)
#include <sys/time.h>
#endif

static const double Dec30_1899      = 15018.0;			// MJD of Dec 30 1899 0:0:0 GMT (Microsoft DATE struct)
static const double Jan1_1970       = 40587.0;			// MJD of Jan 1  1970 0:0:0 GMT
static const double Jan1_1993       = 48988.0;			// MJD of Jan 1  1993 0;0:0 CMT (Used for EOS TAI time)



/*-----------------------------------------------------------------------------
 *					nxTimeStamp::nxTimeStamp		2005-8-4*/
/** Default constructor for the nxTimeStamp class
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp::nxTimeStamp(void)
{
   m_Val = 0.0;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::nxTimeStamp		2005-8-4*/
/** Constructor for the nxTimeStamp class. Sets the
 *	time to the supplied NewVal, where NewVal is the modified julian Date.
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp::nxTimeStamp(const double mjd )
{
   MJD( mjd );
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::nxTimeStamp		2005-8-4*/
/** */
/*---------------------------------------------------------------------------*/

nxTimeStamp::nxTimeStamp(const int  mjd )
{
   MJD( (double)mjd );
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::nxTimeStamp		2005-8-4*/
/** */
/*---------------------------------------------------------------------------*/

nxTimeStamp::nxTimeStamp( const nxTimeStamp &tnew )
{
   m_Val = tnew.m_Val;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::nxTimeStamp		2005-8-4*/
/** Initialized constructor for the nxTimeStamp class. Str must be in the following
 *  format:  YYYY/MM/DD hh:mm:ss.ttt  to be valid.
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp::nxTimeStamp(const char *NewStr)
{
   SetToUTC(NewStr);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::nxTimeStamp						2005-8-4*/
/** Sets the MJD from the specified instant in time.
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp::nxTimeStamp( int day, int month, int year,
		      int hour, int min,  int secs, double ticks)
{
   SetToUTC( day, month, year, hour, min, secs, ticks);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::JD		2005-8-4*/
/** Returns the Julian Date of the time.
 */
/*---------------------------------------------------------------------------*/

double nxTimeStamp::JD() const
{
   return MJD()+2400000.5;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::JD2000								2005-8-4*/
/** Returns the number of days since JD2000. ie JD 2451545.0.
*/
/*---------------------------------------------------------------------------*/

double nxTimeStamp::JD2000() const
{
   return MJD() - 51544.5;
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::JD2000Centuries		2005-8-4*/
/**	Returns number of Julian centuries since JD2000. Often used in star catlogue
 *	calculations.
 */
/*---------------------------------------------------------------------------*/

double nxTimeStamp::JD2000Centuries() const
{
   return JD2000()/36525.0;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::MJD1950		2005-8-4*/
/** Returns the number of days since JD1950.
 */
/*---------------------------------------------------------------------------*/

double nxTimeStamp::MJD1950() const
{
   return MJD() - 33282.0;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::GMST		2005-8-4*/
/**	Returns the Greenwich Mean Sidereal Time at the UT time represented internally.
 *	Returns it as a number of days. Follows the formula given in Astronomical
 *  Almanac 1990, page B6. I assume the time of this object is UT. 
 *	This is normally not TRUE as I normally have the member represent UTC. 
 *	Apart from that the conversion is accurate.
 *
 *	\par HISTORY
 *	11-FEB-1999, NDL,  Updated coefficeints using 1997 Astronomical Almanac
 */
/*---------------------------------------------------------------------------*/

double nxTimeStamp::GMST()  const
{
   double ut;
   double ut0;
   double frac;
   static double gmst    = -9999.0;
   static double lastmjd = -9999.0;
   double Tu;

   if (lastmjd == MJD()) return gmst;
   lastmjd = MJD();

   ut   = JD2000();					// get time in days since JD2000
   ut0  = floor(ut-0.5)+0.5;		// Get JD2000 time at 0UT
   frac = (ut - ut0);				// Get fractional part of day.
   Tu   = ut0/36525.0;				// Time in centuries from JD2000 to 0 hours UT

   gmst = (    24110.54841
            + (8640184.812866 
	        + (0.093104
	        - (0.0000062)*Tu)*Tu)*Tu)/86400.0
	        + frac*1.00273790935;

   gmst = gmst - floor(gmst);
 return gmst;   
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::fromGMST		2005-8-4*/
/** Sets the UT time from the Greenwich Mean Sidereal
 *	Time.  The gmst is expressed as a number of days.  This modifies
 *	the fractional part of the day number.  The integer part of the day number
 *	is unmodified.
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::fromGMST(double gmst)
{
   nxTimeStamp Today;
   double    frac;

   Today = ZeroUT();						// Get UT at Zero hours.
   frac  = gmst - Today.GMST();				// Get difference in sidereal time
   frac  = nxmath::inrange(frac, 1.0 );		// watch out for day wraparound problems. 
   frac *= ONESIDEREALDAY;					// Convert sideral difference to UT difference.
   Today = Today + frac;					// Add UT difference onto today.
   return Today;							// return cuurent Date and UT of the GMST.
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::ZeroUT		2005-8-4*/
/** Return the timestamp at zero hours UT. Its a floor function call.
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::ZeroUT() const
{
   nxTimeStamp Tnow;

   Tnow.MJD( floor(MJD()) );
   return Tnow;
}

    

//-------------------------------------------------------------------------------
//	nxTimeStamp::SetToUTC
//-------------------------------------------------------------------------------



/*-----------------------------------------------------------------------------
 *					nxTimeStamp::SetToUTC		2005-8-4*/
/** Set the current time to the UTC specified in the variables passed.
 *	Calculates the Modified Julian Date (MJD) at the given instant of time (UTC).
 *	All units are self explanatory except ticks.  Parameter ticks represents
 *	the fraction of a second for the specified UTC. It is expressed as a fraction of a second
 *	(i.e. it should be between 0 and 1).  The code properly handles the discontinuity in the
 *	Gregorian calendar in 1582.
 *
 *	\par 0 A.D to 100 A.D.
 *	The code does not handle years between 0 A.D. and 100 A.D. properly.  It assumes a year between
 *	40 and 99 is really 1940 A.D. to 1999 A.D. and years between 0 and 39 are 2000 A.D to 2039 A.D.
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::SetToUTC( int day, int month, int year,
			  int hour, int min,  int secs, double ticks)
{
   int      y;
   int      B;
   double   A;
   double   mjd;

   if ((year < 100) && (year >= 0))
   {
      if (year > 40) year +=1900; else year +=2000;
   }
   y = year;


   if ( month <= 2)
   {
      y -= 1;
      month += 12;
   }

   if (year < 1582 || (year == 1582 && (month < 10 || ((month == 10) && (day < 15)))) )
   {
      B = -2 + (y+4716)/4 - 1179;
   }
   else
   {
      B = (y/400) - (y/100)+ (y/4);
   }
   A =365.0*y - 679004.0;
   mjd =  A+B+((int)(30.6001*(month+1)))+day
	    +(hour/24.0)+min/1440.0+(secs+ticks)/86400.0;
   MJD(mjd);
}

//----------------------------------------------------------------------------
// 	nxTimeStamp::SetToUTC		"
//----------------------------------------------------------------------------


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::SetToUTC		2005-8-4*/
/** Set the current time to the UTC time in the string. The string should be in the
 *	format "YYYY-MM-DD hh:mm:ss.xxx
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::SetToUTC(const char *utstring1)
{
   int    value[5] = { 1992, 1, 1, 0, 0};
   double secval   = 0.0;
   int    secs;
   char   utstring[80];
   char   *token; 													// A local copy of the string.
   char   c;
   unsigned i;

   strncpy( utstring, utstring1, sizeof(utstring) );		// Copy passed string to local copy.
   utstring[sizeof(utstring)-1] = 0;							// ensure local copy is null terminated.
   for (i= 0; i < strlen(utstring); i++ )			// Scan through the string
   {																		// and convert any characters
      c = utstring[i];												// other than '0'-'9' and '.'
      if ( ( c <  '0'  ||  c > '9' ) && ( c != '.' ) )	// to whitespace.
      {
	 utstring[i] = ' ';
      }
   }
   token = strtok( utstring, " ");
   for (i=0; i < 5; i++)
   {
      if (token != NULL) value[i] = atoi( token );
      token = strtok( NULL, " ");
   }
   if (token !=NULL) secval = atof( token );

   secs = (int)secval;
   secval = secval-secs;
   SetToUTC( value[2], value[1], value[0], value[3], value[4], secs, secval);
}																			// return.

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::SetToDeltaTime		2005-8-4*/
/**	Set the current time to the UTC time in the string. Format should be of form
 *	day:hour:min:secs.xxx.  Not used very often.
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::SetToDeltaTime(const char *deltatime)
{
   int    value[3] = { 0, 0, 0};
   double secval   = 0.0;
   char   utstring[80];
   char   *token; 					// A local copy of the string.
   char   c;
   double dt;
   unsigned	i;

   strncpy( utstring, deltatime, sizeof(utstring) );	// Copy passed string to local copy.
   utstring[sizeof(utstring)-1] = 0;			// ensure local copy is null terminated.
   for ( i= 0; i < strlen(utstring); i++ )		// Scan through the string
   {							// and convert any characters
      c = utstring[i];					// other than '0'-'9' and '.'
      if ( ( c < '0'  ||  c > '9' ) && ( c != '.' ) )	// to whitespace.
      {
	 utstring[i] = ' ';
      }
   }
   token = strtok( utstring, " ");
   for (i=0; i < 3; i++)
   {
      if (token != NULL) value[i] = atoi( token );
      token = strtok( NULL, " ");
   }
   if (token !=NULL) secval = atof( token );

   dt = value[0] + (value[1]/24.0) + (value[2]/1440.0) + secval/86400.0;
   MJD(dt);
}								// return.


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::GetUTC		2005-8-4*/
/** Returns the internal time format as a UTC set of values constituting a date
 *	plus UTC time. Finds the civil calendar date for a given value               
 *	of the Modified Julian Date (MJD). Julian calendar is used up to
 *	1582 October 4, Gregorian calendar is used from 1582 October 15 onwards.
 *	Follows Algorithm given in "Astronomy on the Personal Computer" 
 * O. Montenbruck and T. Pfleger.
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::GetUTC( int* UDAY,  int* UMONTH, int* UYEAR,
		     int* UHOUR, int* UMIN,   int* USECS,
		     double* UTICKS) const
{
   long         B, D, F;
   double       jd, JD0, C, E, Dayfrac;
   double       mjd;
   int			DAY,MONTH,YEAR, HOUR,MIN,SECS;
   double		TICKS;

   mjd     = MJD();
   Dayfrac = mjd - floor(mjd);				// Get the fraction of UTC day.
   jd      = floor(mjd) + 2400000.5;		// Get the Julian Date 
   JD0     = (long)(jd+0.5);					// and add a half.
   if ((JD0 < 2299161.0))						// Determine the calendar
   {													// Its Julian
 //     B = 0;
      C = JD0 + 1524.0;
   }													// Otherwise its
   else												// Gregorian.
   {							
      B = (long)((JD0 - 1867216.25) / 36524.25); 
      C = JD0 + (B - (B/4)) + 1525.0; 
   } 
   D     = (long)((C - 122.1) / 365.25);
   E     = 365.0 * D + (D/4);
   F     = (long)((C - E) / 30.6001);
   DAY   = (int) ( (long)(C - E + 0.5) - (long)(30.6001 * F));
   MONTH = (int) (F - 1 - 12*(F/14));
   YEAR  =(int)( D - 4715 - ((7+MONTH)/10));

   HOUR     = (int)(Dayfrac*24.0);
   MIN      = (int)(fmod(Dayfrac*1440.0,  60.0));
   Dayfrac *= 86400.0;
   TICKS    = (fmod(Dayfrac, 60.0));
   SECS     =(int)TICKS;
   TICKS    = TICKS- SECS;

   if (UDAY  ) *UDAY   = DAY;
   if (UMONTH) *UMONTH = MONTH;
   if (UYEAR ) *UYEAR  = YEAR;
   if (UHOUR ) *UHOUR  = HOUR;
   if (UMIN  ) *UMIN   = MIN;
   if (USECS ) *USECS  = SECS;
   if (UTICKS) *UTICKS = TICKS;
} 

/*-----------------------------------------------------------------------------
 *					*nxTimeStamp::UTCStr		2005-8-4*/
/** Returns the time as a ASCII string in the format YYYY/MM/DD hh:mm:ss.ttt
 *	The parameter Str must point to a buffer that is at least 24 characters in size.
 *	Function fills the Str buffer with the time string and also returns a pointer to this
 *	buffer for convenience.  A safer implementation using nxString is also provided.
 */
/*---------------------------------------------------------------------------*/

char *nxTimeStamp::UTCStr( char *Str, int DoMillisecs ) const
{
   int Year;
   int Month;
   int Day;
   int Hour;
   int Minute;
   int Second;
   double ticks;
   int MSec;
   nxTimeStamp T0;

   T0 = *this + 0.0005*ONESECOND;	//round off the to nearest millisecond.


   T0.GetUTC( &Day, &Month, &Year, &Hour, &Minute, &Second, &ticks);
   if (DoMillisecs)
   {
      MSec = (int)(ticks*1000.0);
      sprintf(Str,"%04d-%02d-%02d %02d:%02d:%02d.%03d",
	       Year,Month,Day,Hour,Minute,Second,MSec);
   }
   else
   {
      sprintf(Str,"%04d-%02d-%02d %02d:%02d:%02d",
	       Year,Month,Day,Hour,Minute,Second);
   }
   return Str;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::UTCStr		2005-8-4*/
/** Returns the time as a ASCII string in the format YYYY/MM/DD hh:mm:ss.ttt
 */
/*---------------------------------------------------------------------------*/

nxString nxTimeStamp::UTCStr( int DoMillisecs ) const
{
   char    buffer[50];
   nxString text;
   
   text = UTCStr( buffer, DoMillisecs );
   return text;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::TimeStr		2005-8-4*/
/** Returns the time as a ASCII string in the format hh:mm:ss.ttt.
 *	The parameter Str must point to a buffer that is at least 14 characters in size.
 *	Function fills the Str buffer with the time string and also returns a pointer to this
 *	buffer for convenience.  A safer implementation using nxString is also provided.
 */
/*---------------------------------------------------------------------------*/

char* nxTimeStamp::TimeStr( char *Str, int DoMillisecs ) const
{
   int Year;
   int Month;
   int Day;
   int Hour;
   int Minute;
   int Second;
   double ticks;
   int MSec;
   nxTimeStamp T0;

   T0 = *this + 0.0005*ONESECOND;				//round off the to nearest millisecond.
   T0.GetUTC( &Day, &Month, &Year, &Hour, &Minute, &Second, &ticks);
   if (DoMillisecs)
   {
      MSec = (int)(ticks*1000.0);
      sprintf(Str,"%02d:%02d:%02d.%03d",Hour,Minute,Second,MSec);
   }
   else
   {

      sprintf(Str,"%02d:%02d:%02d",Hour,Minute,Second);
   }
   return Str;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::TimeStr		2005-8-4*/
/** Returns the time as a ASCII string in the format hh:mm:ss.ttt**/
/*---------------------------------------------------------------------------*/

nxString nxTimeStamp::TimeStr( int DoMillisecs ) const
{
   char    buffer[50];
   nxString text;
   
   text = TimeStr( buffer, DoMillisecs );
   return text;
}


//-------------------------------------------------------------------------------
//	nxTimeStamp::DateStr		Returns the time as a ASCII string
//					in the format hh:mm:ss.ttt
//-------------------------------------------------------------------------------

char *nxTimeStamp::DateStr( char *Str) const
{
   int Year;
   int Month;
   int Day;
   int Hour;
   int Minute;
   int Second;
   double ticks;

   GetUTC( &Day, &Month, &Year, &Hour, &Minute, &Second, &ticks);
   sprintf(Str,"%04d/%02d/%02d",Year, Month, Day );
   return Str;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::DateStr		2005-8-4*/
/** Returns the date as a ASCII string in the format year-month-day
 */
/*---------------------------------------------------------------------------*/

nxString nxTimeStamp::DateStr() const
{
   char    buffer[50];
   nxString text;
   
   text = DateStr( buffer );
   return text;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::DayOfWeek		2005-8-4*/
/** Returns the day of the week as a number (0-6) Sunday = 0, Saturday = 6.
 */
/*---------------------------------------------------------------------------*/

int nxTimeStamp::DayOfWeek() const
{
  double today;

  today  = floor (MJD());
  return (int)fmod( today-4, 7.0);	//MJD of day 4 is a Sunday.
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::NameOfDay		2005-8-4*/
/** Returns a pointer to the name of the current day, eg "Sunday", "Monday" etc.
 */
/*---------------------------------------------------------------------------*/

const char* nxTimeStamp::NameOfDay() const
{
  static char DayName[7][10] = {"Sunday","Monday","Tuesday","Wednesday", "Thursday","Friday","Saturday"};
  int    dow;

  dow = DayOfWeek();
  return DayName[dow];
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::SetToDate		2005-8-4*/
/** Set the time to the UTC on the day specified at 0 hours UTC
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::SetToDate( int day, int month, int year )
{
   SetToUTC( day,month,year, 0,0,0,0.00);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::GetDate		2005-8-4*/
/** Return the current calendar date of the time
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::GetDate( int* day, int* month, int* year) const
{
//   int     hour, min, secs;
// double  ticks;

   GetUTC( day, month, year, NULL, NULL, NULL, NULL );
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator>		2005-8-4*/
/** Returns true if this time is greater than time2
*/
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator> (const nxTimeStamp& Time2)  const
{
   return ( MJD() > Time2.MJD());
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator>=		2005-8-4*/
/** Returns true if this time is greater than or equal to Time2**/
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator>= ( const nxTimeStamp& Time2)  const
{
   return ( MJD() >= Time2.MJD() );
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator<		2005-8-4*/
/** Returns true if this time is less than time2
**/
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator< (const nxTimeStamp& Time2)  const
{
   return ( MJD() < Time2.MJD());
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator<=		2005-8-4*/
/** returns true if this time is less than or equal to Time2
*/
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator<= (const nxTimeStamp& Time2)  const
{
  return ( MJD() <= Time2.MJD() );
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator==		2005-8-4*/
/** Returns true if this time exactly equals Time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator== ( const nxTimeStamp& Time2)  const
{
   return (MJD() == Time2.MJD());
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator!=		2005-8-4*/
/** Returns true if this time not equal Time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator!= ( const nxTimeStamp& Time2) const
{
   return (MJD() != Time2.MJD());
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator>		2005-8-4*/
/** Returns true if this time is greater than time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator> ( const double Time2) const
{
   return ( MJD() > Time2);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator>=		2005-8-4*/
/** Returns true if this time is greater than or equal to Time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator>= ( const double Time2) const
{
   return ( MJD() >= Time2 );
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator<		2005-8-4*/
/** Returns true if this time is less than time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator< ( const double Time2) const
{
   return ( MJD() < Time2);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator<=		2005-8-4*/
/** returns true if this time is less than or equal to Time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator<= ( const double Time2) const
{
  return ( MJD() <= Time2 );
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator==		2005-8-4*/
/** Returns true if this time exactly equals Time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator== ( const double Time2) const
{
   return (MJD() == Time2);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator!=		2005-8-4*/
/**	Returns true if this time not equal Time2
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::operator!= ( const double Time2)  const
{
   return (MJD() != Time2);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator+		2005-8-4*/
/** Returns (this) + Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator+ ( const nxTimeStamp& Time2)  const
{
   nxTimeStamp  Sum( m_Val + Time2.m_Val);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator+		2005-8-4*/
/** Returns (this) + Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator+ ( const double Time2)  const
{
   nxTimeStamp  Sum(m_Val + Time2);
   return Sum;
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator-		2005-8-4*/
/** returns (this) - Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator- ( const nxTimeStamp& Time2) const 
{
   nxTimeStamp  Sum(m_Val - Time2.m_Val);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator-		2005-8-4*/
/** returns (this) - Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator- (const double Time2) const 
{
   nxTimeStamp  Sum(m_Val - Time2);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator/		2005-8-4*/
/** Returns (this) / Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator/ (const nxTimeStamp& Time2) const 
{
   nxTimeStamp  Sum(m_Val / Time2.m_Val);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator/		2005-8-4*/
/** Returns (this) / Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator/ (const double Time2) const 
{
   nxTimeStamp  Sum(m_Val / Time2);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator*		2005-8-4*/
/** returns (this) *Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator* (const nxTimeStamp& Time2) const 
{
   nxTimeStamp  Sum(m_Val*Time2.m_Val);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator*		2005-8-4*/
/** returns (this) - Time2
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp nxTimeStamp::operator* (const double Time2) const 
{
   nxTimeStamp  Sum(m_Val*Time2);
   return Sum;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::FromSystem		2005-8-4*/
/** Sets the time to UTC from the local operating system.  On Windows this has
 *	the resolution of the system clock.  On Unix it has 1 second resolution.
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::FromSystem()
{
#if defined (NX_WINDOWS)
	SYSTEMTIME	utc;

	::GetSystemTime( &utc );
	SetToUTC( utc.wDay, utc.wMonth, utc.wYear, utc.wHour, utc.wMinute, utc.wSecond, utc.wMilliseconds/1000.0 );
#else

	struct timeval tp;
	double secs;

	gettimeofday( &tp, NULL);

	//time_t  timer;
	//time( &timer );			// get number of seconds since Jan 1 1970

	secs = tp.tv_sec + tp.tv_usec*1.0E-06;
	FromUnix( secs );			// Convert this Unix time to MJD.
#endif

}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::FromUnix		2005-8-4*/
/** Converts a Unix timestamp to a MJD.  Unix timestamps are the number of
 *	seconds since Jan 1 1970 0:0:0 GMT 
*/
/*---------------------------------------------------------------------------*/

void nxTimeStamp::FromUnix( double SecsSince1970)
{
   MJD( Jan1_1970 + SecsSince1970*ONESECOND );
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::asUnix		2005-8-4*/
/**	Converts the MJD to the number of seconds since Jan 1 1970 0:0:0 GMT
 */
/*---------------------------------------------------------------------------*/

double nxTimeStamp::asUnix() const 
{
   return  (MJD() - Jan1_1970)/ONESECOND;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::FromUnix		2005-8-4*/
/** Converts a EOS TAI time to  MJD.  EOS TAI timestamps are the number of
 *	seconds since Jan 1 1993 0:0:0 GMT.  This code does not account for
 *	leap seconds.
*/
/*---------------------------------------------------------------------------*/

void nxTimeStamp::FromEosTAI( double eosTAI )
{
   MJD( Jan1_1993 + eosTAI*ONESECOND );
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::asUnix		2005-8-4*/
/**	Converts the MJD to the number of seconds since Jan 1 1993 0:0:0 GMT.
 *  There may be some discrepancy with leap seconds.
 */
/*---------------------------------------------------------------------------*/

double nxTimeStamp::asEosTAI() const 
{
   return  (MJD() - Jan1_1993)/ONESECOND;
}
//----------------------------------------------------------------------------
//			nxTimeStamp::operator=
//	Equates the time to the timestamp given.
//----------------------------------------------------------------------------


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator=		2005-8-4*/
/** Copy the other time to this time. Return this time object.
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp&  nxTimeStamp::operator=( const nxTimeStamp &Time2)
{
   MJD( Time2.MJD() );
   return *this;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::operator=		2005-8-4*/
/** Copy the other time to this time. Return this time object. 
 */
/*---------------------------------------------------------------------------*/

nxTimeStamp&  nxTimeStamp::operator=( const double  Time2)
{
   MJD( Time2 );
   return *this;
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::DayOfYear		2005-8-4*/
/** Returns the daynumber (1-366) of the current time. Day 1 is Jan 1 of the
 *	current year.
 */
/*---------------------------------------------------------------------------*/

int nxTimeStamp::DayOfYear() const
{
   int    DAY;
   int    MONTH;
   int    YEAR;
   int    HOUR;
   int    MIN;
   int    SECS;
   double TICKS;
   nxTimeStamp Day1;

   GetUTC( &DAY, &MONTH, &YEAR, &HOUR, &MIN, &SECS, &TICKS);
   Day1.SetToUTC( 1, 1, YEAR, 0, 0, 0, 0.0);
   return  (int)(MJD() - Day1.MJD()) + 1;
}

//---------------------------------------------------------------------------
//			nxTimeStamp::FromVaxDateAndTimeStr
//  Convert a VAX date/time str of the form "19-JAN-1994 23:59:59.999" to
//  a timestamp form.
//---------------------------------------------------------------------------


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::FromVaxDateAndTimeStr		2005-8-4*/
/** Convert a VAX date/time str of the form "19-JAN-1994 23:59:59.999" to
 *  a timestamp
 */
/*---------------------------------------------------------------------------*/

nxBOOL nxTimeStamp::FromVaxDateAndTimeStr( const char * VaxDateAndTimeStr)
{
   nxStringArray datevals;
   nxStringArray timevals;
   nxStringArray words;
   static const char Vaxmonth[12][4] = { "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };
   nxBOOL ok = nxFALSE;

   MJD(0.0);											// default to failure
   int nw = words.Strtok( VaxDateAndTimeStr);		// See how many words are in the string
   if (nw != 2) return ok;								// there should be two, if not then error

   int nd = datevals.Strtok( words[0],  "-" );		// get the number of words in the date
   int nt = timevals.Strtok( words[1],  ":" );		// get the number of words in the time
   if ((nd != 3) || (nt !=3 )) return ok;			// if not 3 in both then quit.
     
   datevals[1].MakeUpper();					// Get an upper copy of the month string
   nxString monthstr = datevals[1];
   int     month;

   for (month = 0; month < 12; month++)				// now find out where the month string
   {								// is the same as the
      if (monthstr == Vaxmonth[month]) break;			// predefined list and break
   }								// loop until found or no more months

   if (month < 12)						// if we have a valid month
   {								// then
        month++;						// add 1 to get correct month
	int    day  = atoi( datevals[0] );			// get the date
	int    year = atoi( datevals[2] );			// get the year
	int    hour = atoi( timevals[0] );			// get the hour
	int    mins = atoi( timevals[1] );			// get the minutes
	double secs = atof( timevals[2] );			// get the seconds
	ok =     ( day  > 0  )  && ( day <  32)			// check the range
	     &&  ( hour >= 0 )  && ( hour < 24)			// of all of the 
	     &&  ( mins >= 0 )  && ( mins < 60)			//entries
	     &&  ( secs >= 0.0) && ( secs < 60.0)
	     &&   (year > 1970) && ( year < 2200);
	if (ok)	
	{
	   double ticks = secs - floor(secs);				// get the fraction of seconds
	   int    second = (int)secs;					// integerise the seconds
	   SetToUTC( day, month, year, hour, mins, second, ticks);	// and set to this mjd
	}
   }
   return ok;
}


/*-----------------------------------------------------------------------------
 *					nxTimeStamp::FromMicrosoftDATE		2005-8-4*/
/** Convert the Microsft VT_DATE variant value to an MJD 
 */
/*---------------------------------------------------------------------------*/

void nxTimeStamp::FromMicrosoftDATE( double date )
{
	MJD(date + Dec30_1899);
}

/*-----------------------------------------------------------------------------
 *					nxTimeStamp::AsVTDate		2005-8-4*/
/** Convert the internal date as a Microsft VT_DATE variant
 */
/*---------------------------------------------------------------------------*/

VARIANT nxTimeStamp::AsVTDate()
{
	VARIANT	var;

	::VariantInit(&var);
	var.vt = VT_DATE;
	var.date = MJD() - Dec30_1899;
	return var;
}


