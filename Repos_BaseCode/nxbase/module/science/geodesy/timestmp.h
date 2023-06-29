#if !defined(NXBASE_NXTIMESTAMP_H)
#define NXBASE_NXTIMESTAMP_H 1

#include "../../system/strings/nxstring.h"

/*-----------------------------------------------------------------------------
 *					class nxTimeStamp								2004-11-23*/
/** \ingroup Geodesy
	A class used to represent instants of time in a Modified Julian Date (MJD) format.
 *	which is the Julian date - 2400000.5. The object internally stores times
 *	as a "double " which means it can represent times to an accuracy of 10 micro seconds
 *	within a few hundred years of 2000 on machines that use a 64 bit
 *	double precision number. This is sufficient for most purposes.
 *
 *	\par Operators
 *	The class is an overload of a \b double and provides most of the numerical
 *	operators applicable to doubles.  The class can also be used to provide
**/
/*---------------------------------------------------------------------------*/

class nxTimeStamp
{
	public:
	static const double ONESECOND;			//!< One Second expressed in days.
	static const double ONEMINUTE;			//!< One minute expressed in days.
	static const double ONEHOUR;			//!< One Hour expressed in days
	static const double ONESIDEREALDAY; 	//!< One sidereal day expressed in UTC days
	static const int    SUNDAY;				//!< Represents Sunday as a day of the week
	static const int    MONDAY;				//!< Represents Monday as a day of the week
	static const int    TUESDAY;			//!< Represents Tuesday as a day of the week
	static const int    WEDNESDAY;			//!< Represents Wednesday as a day of the week
	static const int    THURSDAY;			//!< Represents Thursday as a day of the week.
	static const int    FRIDAY;				//!< Represents Friday as a day of the week
	static const int    SATURDAY;			//!< Represents Saturday as aday of the week

   private:
      double    	m_Val;														// Internal format is Modified julian date in Days.

   public:
					nxTimeStamp     	( );									//!< Default constructor. Initialize time to 0
					nxTimeStamp     	( const int       MJD );				//!< Set timestamp to given modified julian date.
					nxTimeStamp     	( const double    MJD );				//!< Set timestamp to given modified julian date.
					nxTimeStamp			( const nxTimeStamp &tnew );	     	//!< copy constructor
					nxTimeStamp     	( const char      *NewStr);				//!< Set timestamp to the given UTC string. "yyyy-mm-dd hh:mm:ss"
					nxTimeStamp     	( int day, int month, int year,			//!< Set timestamp to the given UT
					  				      int hour, int min,  int secs,			// specified in Day, month year, hour min
					  				      double ticks = 0.0);					// secs and ticks.
      void			MJD           	( const double mjd) { m_Val = mjd;}			//!< Set timestamp to given modified julian date.
      double 		MJD           	( void ) const { return m_Val;}				//!< Return timestamp as modified Julian Date.
      double		JD            	( void ) const; 							//!< Return timestamp as a Julian date.
      double		JD2000	        ( void ) const;								//!< Return Number of days since JD 2000.
      double    	JD2000Centuries ( void ) const;								//!< Return Number of julian centuries since JD 2000 epoch
      double		MJD1950	        ( void ) const;								//!< Return Number of days since JD 1950 epoch
      double		GMST	      	( void ) const;								//!< Return Greenwich Mean sidereal time at this UT. (fraction of day).
      nxTimeStamp 	fromGMST      	( double gmst); 							//!< Set the UTC time using the current Greenwich Mean Sidereal Time
      nxTimeStamp 	ZeroUT        	( void ) const;								//!< Get the time of Zero UT on today.
      void			SetToUTC      	( int day, int month, int year,				//!< Set timestamp to the specified UTC
					  				  int hour, int min,  int secs,
					  				  double ticks = 0.0);
      void			SetToUTC        ( const char *NewStr);						//!< Set timestamp to the specified UTC string "year-month-day hh:mm:ss.tt"
      void			SetToDeltaTime  ( const char *deltatime);					//!< 
      void			GetUTC          ( int* DAY,  int* MONTH, int* YEAR,			//!< Return the timestamp as UTC elements, day, month, year, hour, mins, secs, ticks
					  				  int* HOUR, int* MIN,   int* SECS,
					  				  double* TICKS) const;
      char*     	UTCStr          ( char *Str, int DoMillisecs = 0 ) const;	//!< Write timestamp as a UTC string in format "yyyy-mm-dd hh:mm:ss.ttt". User buffer must be big enough to hold the string
      nxString		UTCStr			( int DoMilliSeconds = 0) const;			//!< Write timestamp as a UTC string in format "yyyy-mm-dd hh:mm:ss.ttt".
      char*     	TimeStr	        ( char *Str, int DoMillisecs = 0 ) const;	//!< Write time of day section of timestamp as a string in format "hh:mm:ss.ttt"
      nxString		TimeStr			( int DoMilliSeconds = 0 ) const;			//!< Write time of day section of timestamp as a string in format "hh:mm:ss.ttt"
      char*     	DateStr         ( char *Str ) const;						//!< Write date section of timestamp as a string in format "yyyy-mm-dd"
      nxString		DateStr			() const;									//!< Write date section of timestamp as a string in format "yyyy-mm-dd"
      int			DayOfWeek       () const;									//!< Returns Day of week (0-6) 0 = Sunday, 6 = Saturday.
      const char*   NameOfDay       () const;									//!< Return English Name of day, Sunday, Monday etc.
      void			FromSystem      ();											//!< Sets timestamp object from the current system time (in UTC).
      void			SetToDate       ( int day,  int month, int year ); 			//!< Set timestamp to 0Hrs UT on the given day.
      void			GetDate          ( int *day, int *month, int *year) const;	//!< Get the date of the current timestamp.
      void 			FromUnix        ( double SecsSince1970);					//!< Set time from a "Unix time", i.e. seconds since 1970-1-1 00:00:00 UTC
      double 		asUnix() const;												//!< convert this time to equivalent Unix time.
      void 			FromEosTAI		( double eostai);							//!< Set time from an EOS TAI time i.e. seconds since 1993-1-1 00:00:00 UTC
      double 		asEosTAI() const;											//!< convert this time to equivalent EOS TAI time.
      int       	DayOfYear() const;											//!< Return the day of the year January 1 is ??
      nxBOOL 		FromVaxDateAndTimeStr( const char * VaxDateAndTimeStr);		//!< Set the timestamp object from a VAX format date and time
	  void			FromMicrosoftDATE( double date );							//!< Set the timestamp object from the Microsoft VARIANT Date value 
	  VARIANT		AsVTDate();													//!< Return this timestamp object as a Microsoft VARIANT time object
      				operator double () const {return m_Val;}					
	  nxBOOL		operator>      ( const nxTimeStamp &Time2) const;
      nxBOOL		operator>      ( const double     Time2) const;
      nxBOOL		operator>=     ( const nxTimeStamp &Time2) const;
      nxBOOL		operator>=     ( const double     Time2) const;
      nxBOOL		operator<      ( const nxTimeStamp &Time2) const;
      nxBOOL		operator<      ( const double     Time2) const;
      nxBOOL		operator<=     ( const nxTimeStamp &Time2) const;
      nxBOOL		operator<=     ( const double     Time2) const;
      nxBOOL		operator==     ( const nxTimeStamp &Time2) const;
      nxBOOL		operator==     ( const double     Time2) const;
      nxBOOL      	operator!=     ( const nxTimeStamp &Time2) const;
      nxBOOL      	operator!=     ( const double     Time2) const;
      nxTimeStamp 	operator+      ( const nxTimeStamp &Time2) const;
      nxTimeStamp 	operator+      ( const double     Time2) const;
      nxTimeStamp 	operator-      ( const nxTimeStamp &Time2) const;
      nxTimeStamp 	operator-      ( const double     Time2) const;
      nxTimeStamp 	operator/      ( const nxTimeStamp &Time2) const;
      nxTimeStamp 	operator/      ( const double     Time2) const;
      nxTimeStamp 	operator*      ( const nxTimeStamp &Time2) const;
      nxTimeStamp 	operator*      ( const double     Time2) const;
      nxTimeStamp& 	operator=      ( const nxTimeStamp &Time2);
      nxTimeStamp&  operator=      ( const double     Time2);
   };

#endif

