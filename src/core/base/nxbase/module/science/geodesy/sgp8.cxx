/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_geodesy.h"
using namespace::nxmath;

#include "module/math/zbrent.h"		// Use the template function as its more general.




//----------------------------------------------------------------------------
//       substr
//----------------------------------------------------------------------------

static char *substr( const char *s, char *sub, int i1, int len)
{
   memcpy( sub, &s[i1], len );
   sub[len] = '\0';
   return sub;
}

//----------------------------------------------------------------------------
//       Real_Value
//----------------------------------------------------------------------------

static double Real_Value( const char *s, int i1, int len )
{
  char sub[80];

  substr( s, sub, i1, len );
  return atof(sub);
}

//-----------------------------------------------------------------------------
//       Integer_Value
//----------------------------------------------------------------------------

static long Integer_Value( const char *s, int i1, int len )
{
  char sub[80];

  substr( s, sub, i1, len );
  return atoi(sub);
}

//-------------------------------------------------------------------------------
//          FMOD2P
//-------------------------------------------------------------------------------

static double FMOD2P( double X)
{
   double x;

   x = fmod( X,TWOPI );
   if (x < 0)x += TWOPI;
   return x;
}

//-------------------------------------------------------------------------------
//       ACTAN
//-------------------------------------------------------------------------------

static double ACTAN( double SINX, double COSX)
{
   double val;

   val = atan2( SINX, COSX );
   if (val < 0.0) val += TWOPI;
   return val;
}

// ---------------------------------------------------------------------------
//       Local Define Statements.
//----------------------------------------------------------------------------

static const double   AE      = 1.0;         			// Earths radius.
static const double   TOTHRD  = 2.0/((double)3.0); 		// 2/3 rd
static const double   XKMPER  = 6378.135;          		// Earth equatorial radius - kilometers (WGS '72)
static const double   GE      = 398600.8;          		// GM, Earth gravitational constant (WGS '72) distances in km.
static const double   XJ2     = 1.0826158E-3;      		// J2 harmonic (WGS '72)
static const double   XJ3     = -2.53881E-6;       		// J3 harmonic (WGS '72)
static const double   XJ4     = -1.65597E-6;       		// J4 harmonic (WGS '72)
static const double   XMNPDA  = 1440.0;            		// No of minutes in one day.
static const double   QO      = 120.0;       			// Parameter for SGP4/SGP8 density function
static const double   S0      = 78.0;        			// Parameter for SGP4/SGP8 density function.
static const double   RHO     =  0.15696615;             // Reference value of atmospheric density (2.461E-05*XKMPER)
inline double SIN ( double X) { return sin(X);   }
inline double COS ( double X) { return cos(X);   }
inline double SQRT( double X) { return sqrt(X);  }
inline double ABS ( double X) { return fabs(X);  }
inline double POW4( double X) { return (X*X*X*X);}

static const double CK2    = 0.5*XJ2*(AE*AE);         			// Agrees with NORAD doc.
static const double CK4    = -0.375*XJ4*POW4(AE);     			// Agrees with NORAD doc.
static const double QOMS2T = POW4( ((QO-S0)*AE/XKMPER) );   	// Agrees with NORAD doc.
static const double S      = AE*( 1.0+S0/XKMPER );             	// Agrees with NORAD doc.
static const double XKE    = sqrt(3600.0*GE/pow(XKMPER,3.0));  	// Agrees with NORAD doc.

//----------------------------------------------------------------------------
//       nxSatelliteBase::Constructor
//----------------------------------------------------------------------------

nxSatelliteBase::nxSatelliteBase()
{
   	m_DistantObject = nxFALSE;     	// Flag this as not a distant object.
	m_DaysPerRev 	= 0 ;
	m_StartOrbitNumber = 0;
	m_StartOfOrbit.MJD(0.0);
	m_velocity.SetCoords(0,0,0);
	m_epoch.MJD(0.0);
}


/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::Period							2004-11-23*/
/** Return the approximate period of the satellite.  For low earth orbit spacecraft
 *	this is accurate to about 0.1 seconds. The Period is returned in days.
 **/
/*---------------------------------------------------------------------------*/

double nxSatelliteBase::Period()
{
   return m_DaysPerRev;
}


//-------------------------------------------------------------------------------
//      nxSatelliteBase::ZComponent
//-------------------------------------------------------------------------------


/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::ZComponent		2004-11-23*/
/** Return the Z component of the satellite at a time deltamjd from
 *	the current epoch.  This is often used to find the ascending node crossing
 *	of the satelitte.
 **/
/*---------------------------------------------------------------------------*/

double nxSatelliteBase::ZComponent( double deltamjd )
{
  nxTimeStamp Tnow;

  Tnow.MJD( m_epoch.MJD()+deltamjd);    	// Set up a new m_time to get satellite position.
  UpdateECIPosition( Tnow );     			// Get the satellite position atthis m_time.
  return m_location.Z();         			// and return the Z component.
}

/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::EquatorCrossing		2004-11-23*/
/**	Determines the last equator crossing before or equal to
 *  m_time Tnow.
**/
/*---------------------------------------------------------------------------*/

nxTimeStamp nxSatelliteBase::EquatorCrossing( const nxTimeStamp &Tnow)
{
   double deltat = 0.2*m_DaysPerRev;
   double z;
   nxTimeStamp positiveT;
   nxTimeStamp negativeT;
   double  crosstime;

// ----  Find when Z component is positive before or equal to now.

   positiveT = Tnow - m_epoch;
   do
   {
      z = ZComponent( positiveT.MJD() );
      if ( z < 0.0) positiveT = positiveT - deltat;
   } while ( z < 0.0);

// ---- Find when Z component is negative before positiveT.

   negativeT = positiveT - deltat;
   do
   {
      z = ZComponent( negativeT.MJD() );
      if ( z >= 0.0) negativeT = negativeT - deltat;
   } while ( z >= 0.0);

//   negativeT = negativeT;
//   positiveT = positiveT;
//   func = &(nxSatelliteBase::ZComponent);
//   crosstime = zbrent( this, func, negativeT.MJD(), positiveT.MJD(), 1.0E-12);

//#ifdef __BORLANDC__
//   SATFUNCPTR func = &(nxSatelliteBase::ZComponent);
//#else
//   SATFUNCPTR func = ZComponent;
//#endif

   crosstime = zbrent( EvalZcomponent(this), negativeT.MJD(), positiveT.MJD(), 1.0E-12);
   positiveT.MJD( crosstime + m_epoch.MJD() );
   UpdateECIPosition(positiveT);
   return positiveT;
}


//----------------------------------------------------------------------------
//       nxSatelliteBase::DepthOfEclipse
//
// Returns the vertical distance (in metres) of satellite from the
// eclipse condition.  Positive values indicate the satellite is in
// sunlight.  Negative values indicate the satellite is in eclipse.
//
// Approximate only.  Possible problems,
//
// 1) Re = Radius of earth is fro spherical earth, Also uses Equatorial
// radius XKMPER, rather than a geocentric mean radius?
//
// 2) When satellite is on sunlit side of earth returns an answer which is
// positive but is not the distance of the satellite above the the eclipse
// terminator, as this is not meaningful.  May cause problems if attempting
// to find times of eclipse using a zeroing function.
//
//----------------------------------------------------------------------------


/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::DepthOfEclipse		2004-11-23*/
/**
 * Returns the vertical distance (in metres) of satellite from the
 * eclipse condition.  Positive values indicate the satellite is in
 * sunlight.  Negative values indicate the satellite is in eclipse.
 * Approximate only.  Possible problems,
 *     - Re = Radius of earth is fro spherical earth, Also uses Equatorial
 *            radius XKMPER, rather than a geocentric mean radius?
 *     - When satellite is on sunlit side of earth returns an answer which is
 *       positive but is not the distance of the satellite above the the eclipse
 *       terminator, as this is not meaningful.  May cause problems if attempting
 *       to find times of eclipse using a zeroing function.
 **/
/*---------------------------------------------------------------------------*/

double nxSatelliteBase::DepthOfEclipse( nxTimeStamp Tnow )
{
   double       xcomp;
   PlanetSun    sun;
   nxVector       sunv;
   double       xperp;
   const double Re = XKMPER*1000.0;             	// Nominal Radius of Earth in metres
   const double Rs = 696000000.0;                  	// Solar radius - meters (IAU 76)
   double       HeightOfEclipse;             		// Height of eclipse terminator
   double   	DepthOfEclipse;                  	// Height of satellite above eclipse.

   if (m_epoch == 0.0) return 0.0;
   UpdateECIPosition ( Tnow );                  // Update satellite position.
   sun.UpdateECIPosition( Tnow );               // Update suns current m_location
   sunv  = sun.Location();
   sunv  = sunv.UnitVector();              				// Get unit vector in sunward direction.
   xcomp = -( m_location & sunv );             			// get satellite component parallel to antisunward direction.
   if ( xcomp < 0.0 )                     				// if on sunward side of earth then not in eclipse
   {                          							// so return a positive
      DepthOfEclipse = (XKMPER*1000.0) - xcomp;          // depth of eclipse.
   }                          // otherwise
   else                          // on dark side of earth
   {
      HeightOfEclipse = Re - ((Rs-Re)/!sun.Location())*xcomp;    // Get height of eclipse from antisunward direction.
      xperp = !(m_location.ComponentPerpendicularTo( sunv ));             // Get distance of satellite from antisunward direction.
      DepthOfEclipse = xperp - HeightOfEclipse;
   }
   return DepthOfEclipse;
}

/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::InEclipse		2004-11-23*/
/**  Returns nxTRUE is the satellite is in eclipse
 **/
/*---------------------------------------------------------------------------*/

nxBOOL nxSatelliteBase::InEclipse( nxTimeStamp & Tnow)
{
   return (DepthOfEclipse( Tnow ) < 0.0 );
}

//----------------------------------------------------------------------------
//       TimeofTerminator
//
//----------------------------------------------------------------------------


/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::TimeOfTerminator		2004-11-23*/
/** Finds the m_time of the terminator between a start and stop m_time.
  * Recursively calls itself until the accuracy of the solution is better
  * than stepsize.
 **/
/*---------------------------------------------------------------------------*/

nxBOOL nxSatelliteBase::TimeOfTerminator( nxTimeStamp TStart, nxTimeStamp Tend, nxTimeStamp &Terminator, double stepsize )
{
   nxBOOL      Dark;
   nxBOOL      DarkNow;
   nxTimeStamp Tnow;
   nxTimeStamp Tlast;
   nxBOOL       FoundTerminator = nxFALSE;

   Tnow       = TStart;
   Tlast      = TStart;
   Dark       = InEclipse( TStart );            		// Get the eclipse status of satellite at start m_time.

   while (( Tnow < Tend) && (!FoundTerminator))       	// while we are within the range
   {                       								// then
      Tlast   = Tnow;                					// Last m_time is now
      Tnow    = Tnow + stepsize;           				// next m_time is incremented by step size.
      if (Tnow > Tend) Tnow = Tend;						// Dont overshoot the end point.
      DarkNow = InEclipse(Tnow);          				// get the current eclipse status
      FoundTerminator = (DarkNow != Dark);         		// has it changed
   }

   if (FoundTerminator)                					// if we found a terminator within valid m_time range.
   {                       								// then
      Terminator = Tnow;               					// set the m_time to be the last m_time analysed
	  if (stepsize > 0.1*nxTimeStamp::ONESECOND)     				// if the step size is greater than 0.1 seconds
      {                       							// then
         stepsize = stepsize/10.0;          			// do another recursive iteration at a finer resolution.
         if (!TimeOfTerminator( Tlast, Tnow+stepsize, Terminator, stepsize ))
         {                								// if we did not find it then there was a major
            Terminator = 0.0;							// error in my logic so simply LOG to file
            FoundTerminator = nxFALSE;					// and return no terminator found
         }
      }
   }
   return FoundTerminator;
}

/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::OrbitNumber		2004-11-23*/
/** Calculates the orbit number at time Tnow and returns the start time
 *	of that orbit.
 **/
/*---------------------------------------------------------------------------*/

long nxSatelliteBase::OrbitNumber( const nxTimeStamp &Tnow, nxTimeStamp* Start )
{
   nxTimeStamp DeltaT;
   long      norbits;

	if (m_epoch == 0.0)
	{
	  *Start = 0.0;
	  return -1;
	}
	if (m_StartOfOrbit == 0.0)											// IF the start of the orbit is undefined
	{																	// then
		nxTimeStamp tthen = m_epoch + 6.0*nxTimeStamp::ONESECOND;
		m_StartOfOrbit = EquatorCrossing( tthen );    	// work it out
	}
	*Start  = EquatorCrossing( Tnow );
	DeltaT = ( *Start - m_StartOfOrbit);
	norbits = (long) (floor( DeltaT.MJD()/m_DaysPerRev + 0.5));
	return (norbits + m_StartOrbitNumber );
}

/*-----------------------------------------------------------------------------
 *					nxSatelliteBase::OrbitNumber					2004-11-23*/
/**  Calculates the orbit number at time Tnow **/
/*---------------------------------------------------------------------------*/

long nxSatelliteBase::OrbitNumber( const nxTimeStamp &Tnow )
{
   nxTimeStamp Start;

   return OrbitNumber( Tnow, &Start );
}

/*-----------------------------------------------------------------------------
 *					NoradSatellite::NoradSatellite					2004-11-23*/
/** default constructor for the NoradSatellite Class**/
/*---------------------------------------------------------------------------*/

NoradSatellite::NoradSatellite()
{
    m_time   = -1.0;       			// set last calculated m_time is different to m_epoch m_time.
	m_XMO    = 0.0;        			// The mean anomaly at epoch         (stored as radians).
	m_XNO    = 0.0;        			// Mean motion,                (stored as radians per minute).
	m_XNODEO = 0.0;                 // The mean longitude of ascending node at epoch (stored as radians)
	m_OMEGAO = 0.0;              	// The mean argument of perigree at epoch        (stored as radians)
	m_EO     = 0.0;                 // The mean eccentricity` at epoch.
	m_XINCL  = 0.0;                 // The mean inclination at epoch              (stored as radians).
	m_XNDT20 = 0.0;                 // First m_time derivative of the Mean motion or Ballistic Coefficient. (stored as change/minute)
	m_XNDD60 = 0.0;                 // Second m_time derivative (stored as change/minute)
	m_BSTAR  = 0.0;                 // the SGP4 type drag coefficient.
	m_ideep  = nxFALSE;
	m_elset  = 0;
	m_epoch.MJD(0.0);         		// default to epoch 0, disables program until two-line elements read in.
}

/*-----------------------------------------------------------------------------
 *					NoradSatellite::TwoLine		2004-11-23*/
/** Update the Norad nxSatelliteBase model with the given two line elements **/
/*---------------------------------------------------------------------------*/

void NoradSatellite::TwoLine( const char* line1, const char* line2 )
{
  int       iexp,ibexp;
  double    a1,ao,del1,delo,xnodp,temp;
  int       year;
  int       dayno;
  double    Ut;
  //nxTimeStamp Tnow;
//
// ------ Decode line 1 of elements -------
//
	strcpy( substr(line1,m_catnr,2,5), m_catnr);		// nxSatelliteBase catalogue number
	year    = (int)Integer_Value (line1,18,2);			// get the year of epoch.
	dayno   = (int)Integer_Value (line1,20,3);			// get the day number of epoch.
	if ((year ==0) && (dayno ==0))       				// if the year and day are zero
	{                     								// then assume screwy elements
		m_epoch.MJD(0.0);             					// and set epoch equal to zero.
		return;
	}

	Ut        =      Real_Value    (line1,23,9);    		// get UT as fraction of day.
	m_XNDT20  =      Real_Value    (line1,33,10);      		// First m_time derivative of the Mean motion or Ballistic Coefficient.
	m_XNDD60  =      Real_Value    (line1,44,6)*1.0E-5;		// mantissa Second m_time derivative of Mean motion
	iexp      = (int)Integer_Value (line1,50,2);    		// exponent of 2nd m_time derivative.
	m_BSTAR   =      Real_Value    (line1,53,6)*1.0E-5;		// mantissa of m_BSTAR drag term if GP4 general perturbation theory was used.
	ibexp     = (int)Integer_Value (line1,59,2);    		// exponent of m_BSTAR
	m_elset   =      Integer_Value (line1,65,3);       		// element number
//
// -------- decode line 2 of elements.
//
	m_XINCL    = Real_Value(line2,8,8);       				// Inclination in degrees
	m_XNODEO   = Real_Value(line2,17,8);         			// Right Ascension of the ascending node degrees
	m_EO       = Real_Value(line2,26,7)*1E-7;    			// Eccentricity.
	m_OMEGAO   = Real_Value(line2,34,8);         			// Argument of perigree, degrees.
	m_XMO      = Real_Value(line2,43,8);         			// Mean anomaly degrees.
	m_XNO      = Real_Value(line2,52,11);        			// Mean motion (revs per day)

	if (year > 57) year +=1900; else year += 2000;
	m_epoch.SetToUTC( 1,1,year, 0,0,0, 0.0);     			// Set epoch to first day of year
	m_epoch = m_epoch + (dayno + Ut - 1.0);        			// add on the appropriate amount of date.
//
// -------- Convert to proper units
//
	m_XNDD60   = m_XNDD60*pow(10.0,iexp);
	m_BSTAR    = m_BSTAR*pow(10.0,ibexp)/AE;
	m_XNODEO   = m_XNODEO*ONE_DEGREE;        				// Convert to radians
	m_OMEGAO   = m_OMEGAO*ONE_DEGREE;        				// Convert to radians.
	m_XMO      = m_XMO*ONE_DEGREE;        					// Convert to radians
	m_XINCL    = m_XINCL*ONE_DEGREE;         				// Convert to radians
	m_XNO      = m_XNO*TWOPI/XMNPDA;         				// Convert to radians per minute
	m_XNDT20   = m_XNDT20*TWOPI/(XMNPDA*XMNPDA);
	m_XNDD60   = m_XNDD60*TWOPI/pow(XMNPDA,3);

//
// ---- Extract the original mean motion and semimajor axis from the
// ---- input elements following the SGP8 algorithm.
//

	a1    = pow( XKE/m_XNO,TOTHRD );
	temp  = 1.5*CK2*(3*pow(cos(m_XINCL),2)-1)/pow(1 - m_EO*m_EO,1.5);
	del1  = temp/(a1*a1);
	ao    = a1*(1 - del1*(0.5*TOTHRD + del1*(1 + 134.0/81.0*del1)));
	delo  = temp/(ao*ao);
	xnodp = m_XNO/(1 + delo);     						// xnodp = revolution rate in radians per minute.
	m_DaysPerRev = (TWOPI/(m_XNO*XMNPDA));    			// Days per revolution
	if (TWOPI/xnodp >= 225)     						// If we need the deep space model.
	{               									// then for the moment.
		m_XMO    = 0.0;         						// place bad values into elements
		m_XNO    = 0.0;         						// so we get gross errors.
		m_XNODEO = 0.0;
		m_OMEGAO = 0.0;
		m_EO     = 0.0;
		m_XINCL  = 0.0;
		m_XNDT20 = 0.0;
		m_XNDD60 = 0.0;
		m_BSTAR  = 0.0;
		m_ideep  = 1;
		m_DaysPerRev = 0.0;
	}
	else                     							// otherwise we have near Earth model
	{                     								// so
		m_ideep = 0;                  					// call the
		Model_Init();               					// model initialisation.
	}
	m_StartOrbitNumber  = Integer_Value( line2,63,5);  	// Get the orbit number of this orbit.
	m_StartOfOrbit      = 0.0;
}

/*-----------------------------------------------------------------------------
 *					NoradSatellite::GetElementsFrom		2004-11-23*/
/** Read the satellite's two line elements from the specified file for the
 *	specified satellite
 **/
/*---------------------------------------------------------------------------*/

nxBOOL NoradSatellite::GetElementsFrom( const char *Filename, const char *SatelliteName)
{
   FILE * infile;
   char   aline[80];
   char   line1[80];
   char   line2[80];
   int    length;
   nxBOOL   found;
   nxBOOL   ok;
   int     satlen;
   nxString    satname = SatelliteName;

   //strncpy( satname, SatelliteName, sizeof(satname) );
   //satname[ sizeof(satname)-1] = '\0';
   satlen = satname.GetLength();
   //strupr(satname);
   satname.MakeUpper();

   infile = fopen( Filename, "rt" );
   if (infile == NULL) return nxFALSE;

   found = nxFALSE;
   ok    = nxTRUE;
   nxString bufl;
   while (!found && ok)
   {
      ok = fgets( aline, sizeof(aline), infile ) != NULL;
      if (ok)
      {
		  bufl = aline;
		length = bufl.Find( satname ); /* THIS MIGHT NEED DEBUGGING */
		found = (length == satlen);
      }
   }
   if (found)
   {
      fgets( line1, sizeof(aline), infile );
      fgets( line2, sizeof(aline), infile );
      line1[ strlen(line1)-1 ] = '\0';
      line2[strlen(line2) -1 ] = '\0';
      TwoLine( line1, line2 );
   }
   if (infile != NULL) fclose(infile);
   return found;
}




//----------------------------------------------------------------------------
//       SGP8::Constructor
//----------------------------------------------------------------------------

SGP8::SGP8()
{
   m_location.SetCoords( 0.0, 0.0, 0.0 );
   m_velocity.SetCoords( 0.0, 0.0, 0.0 );
}


/*-----------------------------------------------------------------------------
 *					SGP8::Model_Init		2004-11-23*/
/** Initialises the SGP8 model using the orbital elements read
 *	into the nxSatelliteBase class.  This will typically happen when
 *  the orbital elements are updated.**/
/*---------------------------------------------------------------------------*/

void SGP8::Model_Init(void)
{
	double A1,AO,AODP,ALPHA2,ALDTAL;
	double B,B1,B2,B3,BETAO2,BETAO;
	double C0DTC0,C1DTC1,C4DT,C5DT;
	double C0,C1,C4,C5,C8,C9,COS2G,COSG;
	double EDDOT,EETA,EOSQ,ETA,ETA2,ETDDT,ETDT;
	double DELO,DEL1,D1DDT;
	double D1,D2,D3,D4,D5,D6,D7,D8,D9,D10;
	double D11,D12,D13,D14,D15,D16,D17,D18,D19,D20;
	double D23,D25;
	double D1DT,D2DT,D3DT,D4DT,D5DT;
	double PARDT1,PARDT2,PARDT4;
	double PO,POM2,PSIM2,PSDTPS;
	double SIN2G,SING1;
	double TEMP,THETA4,TMNDDT,TSDDTS,TSDTTS,TSI;
	double XNDDT,XNDTN,XNTRDT;


//*      RECOVER ORIGINAL MEAN MOTION (XNODP) AND SEMIMAJOR AXIS (AODP)
//*      FROM INPUT ELEMENTS --------- CALCULATE BALLISTIC COEFFICIENT
//*      (B TERM) FROM INPUT B* DRAG TERM

	m_time   = -1.0;             // reset m_time to force calculation on next call.
	A1     = pow( (XKE/m_XNO), TOTHRD );
	COSI   = COS(m_XINCL);
	THETA2 = COSI*COSI;
	TTHMUN = 3.0*THETA2-1.0;
	EOSQ   = m_EO*m_EO;
	BETAO2 = 1.0-EOSQ;
	BETAO  = SQRT(BETAO2);
	DEL1   = 1.5*CK2*TTHMUN/(A1*A1*BETAO*BETAO2);
	AO     = A1*(1.0-DEL1*(0.5*TOTHRD+DEL1*(1.0+134.0/81.0*DEL1)));
	DELO   = 1.5*CK2*TTHMUN/(AO*AO*BETAO*BETAO2);
	AODP   = AO/(1.0-DELO);
	XNODP  = m_XNO/(1.0+DELO);
	B      = 2.0*m_BSTAR/RHO;

	//*      INITIALIZATION

	ISIMP=0;
	PO=AODP*BETAO2;
	POM2=1.0/(PO*PO);
	SINI=SIN(m_XINCL);
	SING1=SIN(m_OMEGAO);
	COSG=COS(m_OMEGAO);
	TEMP=0.5*m_XINCL;
	SINIO2=SIN(TEMP);
	COSIO2=COS(TEMP);
	THETA4=THETA2*THETA2;
	UNM5TH=1.0-5.0*THETA2;
	UNMTH2=1.0-THETA2;
	A3COF=-XJ3/CK2*pow(AE,3.0);
	PARDT1=3.0*CK2*POM2*XNODP;
	PARDT2=PARDT1*CK2*POM2;
	PARDT4=1.25*CK4*POM2*POM2*XNODP;
	XMDT1=0.5*PARDT1*BETAO*TTHMUN;
	XGDT1=-0.5*PARDT1*UNM5TH;
	XHDT1=-PARDT1*COSI;
	XLLDOT=XNODP+XMDT1+0.0625*PARDT2*BETAO*(13.0-78.0*THETA2+137.0*THETA4);
	OMGDT=XGDT1+0.0625*PARDT2*(7.0-114.0*THETA2+395.0*THETA4)
	     +PARDT4*(3.0-36.0*THETA2+49.0*THETA4);
	XNODOT=XHDT1+(0.5*PARDT2*(4.0-19.0*THETA2)+2.0*PARDT4*(3.0-7.0*THETA2))*COSI;
	TSI=1.0/(PO-S);
	ETA=m_EO*S*TSI;
	ETA2=ETA*ETA;
	PSIM2=ABS(1.0/(1.0-ETA2));
	ALPHA2=1.0+EOSQ;
	EETA=m_EO*ETA;
	COS2G=2.0*COSG*COSG-1.0;
	D5=TSI*PSIM2;
	D1=D5/PO;
	D2=12.0+ETA2*(36.0+4.5*ETA2);
	D3=ETA2*(15.0+2.5*ETA2);
	D4=ETA*(5.0+3.75*ETA2);
	B1=CK2*TTHMUN;
	B2=-CK2*UNMTH2;
	B3=A3COF*SINI;
	C0=0.5*B*RHO*QOMS2T*XNODP*AODP*pow(TSI,4.0)*pow(PSIM2,3.5)/SQRT(ALPHA2);
	C1=1.5*XNODP*ALPHA2*ALPHA2*C0;
	C4=D1*D3*B2;
	C5=D5*D4*B3;
	XNDT=C1*((2.+ETA2*(3.+34.*EOSQ)+5.*EETA*(4.+ETA2)+8.5*EOSQ)+
	   D1*D2*B1+   C4*COS2G+C5*SING1);
	XNDTN=XNDT/XNODP;

//      IF DRAG IS VERY SMALL, THE ISIMP FLAG IS SET AND THE
//      EQUATIONS ARE TRUNCATED TO LINEAR VARIATION IN MEAN
//      MOTION AND QUADRATIC VARIATION IN MEAN ANOMALY

	if( ABS(XNDTN*XMNPDA) >= 2.16E-3)
	{
		D6=ETA*(30.+22.5*ETA2);
		D7=ETA*(5.+12.5*ETA2);
		D8=1.+ETA2*(6.75+ETA2);
		C8=D1*D7*B2;
		C9=D5*D8*B3;
		EDOT=-C0*( ETA*(4.+ETA2+EOSQ*(15.5+7.*ETA2))+m_EO*(5.+15.*ETA2)
		      +D1*D6*B1 + C8*COS2G+C9*SING1);
		D20=.5*TOTHRD*XNDTN;
		ALDTAL=m_EO*EDOT/ALPHA2;
		TSDTTS=2.*AODP*TSI*(D20*BETAO2+m_EO*EDOT);
		ETDT=(EDOT+m_EO*TSDTTS)*TSI*S;
		PSDTPS=-ETA*ETDT*PSIM2;
		SIN2G=2.*SING1*COSG;
		C0DTC0=D20+4.*TSDTTS-ALDTAL-7.*PSDTPS;
		C1DTC1=XNDTN+4.*ALDTAL+C0DTC0;
		D9=ETA*(6.+68.*EOSQ)+m_EO*(20.+15.*ETA2);
		D10=5.*ETA*(4.+ETA2)+m_EO*(17.+68.*ETA2);
		D11=ETA*(72.+18.*ETA2);
		D12=ETA*(30.+10.*ETA2);
		D13=5.+11.25*ETA2;
		D14=TSDTTS-2.*PSDTPS;
		D15=2.*(D20+m_EO*EDOT/BETAO2);
		D1DT=D1*(D14+D15);
		D2DT=ETDT*D11;
		D3DT=ETDT*D12;
		D4DT=ETDT*D13;
		D5DT=D5*D14;
		C4DT=B2*(D1DT*D3+D1*D3DT);
		C5DT=B3*(D5DT*D4+D5*D4DT);

		D16= D9*ETDT+D10*EDOT + B1*(D1DT*D2+D1*D2DT) +
		C4DT*COS2G+C5DT*SING1+XGDT1*(C5*COSG-2.*C4*SIN2G);

		XNDDT=C1DTC1*XNDT+C1*D16;

		EDDOT=C0DTC0*EDOT-C0*((4.+3.*ETA2+30.*EETA+EOSQ*(15.5+21.*ETA2))*ETDT
		      +(5.+15.*ETA2+EETA*(31.+14.*ETA2))*EDOT +
		      B1*(D1DT*D6+D1*ETDT*(30.+67.5*ETA2))  +
		      B2*(D1DT*D7+D1*ETDT*(5.+37.5*ETA2))*COS2G+
		      B3*(D5DT*D8+D5*ETDT*ETA*(13.5+4.*ETA2))*SING1+XGDT1*(C9*
		      COSG-2.*C8*SIN2G));

		D25=EDOT*EDOT;
		D17=XNDDT/XNODP-XNDTN*XNDTN;
		TSDDTS=2.*TSDTTS*(TSDTTS-D20)+AODP*TSI*(TOTHRD*BETAO2*D17-4.*D20*
		m_EO*EDOT+2.*(D25+m_EO*EDDOT));

		ETDDT =(EDDOT+2.*EDOT*TSDTTS)*TSI*S+TSDDTS*ETA;
		D18=TSDDTS-(TSDTTS*TSDTTS);
		D19=-(PSDTPS*PSDTPS)/ETA2-ETA*ETDDT*PSIM2-(PSDTPS*PSDTPS);
		D23=ETDT*ETDT;

		D1DDT=D1DT*(D14+D15)+D1*(D18-2.*D19+TOTHRD*D17+2.*(ALPHA2*D25/BETAO2+m_EO*EDDOT)/BETAO2);

		XNTRDT =   XNDT*(2.*TOTHRD*D17+3.*(D25+m_EO*EDDOT)/ALPHA2-6.*(ALDTAL*ALDTAL) + 4.*D18-7.*D19 )
		         	+ C1DTC1*XNDDT+C1*(C1DTC1*D16+
		         	D9*ETDDT+D10*EDDOT+D23*(6.+30.*EETA+68.*EOSQ)+
		  			ETDT*EDOT*(40.+30.*
		  			ETA2+272.*EETA)+D25*(17.+68.*ETA2) +
		  			B1*(D1DDT*D2+2.*D1DT*D2DT+D1*(ETDDT*D11+D23*(72.+54.*ETA2))) +
		  			B2*(D1DDT*D3+2.*D1DT*D3DT+D1*(ETDDT*D12+D23*(30.+30.*ETA2))) *
		  			COS2G+
		  			B3*((D5DT*D14+D5*(D18-2.*D19)) *
		  			D4+2.*D4DT*D5DT+D5*(ETDDT*D13+22.5*ETA*D23)) *SING1+XGDT1*
		  			((7.*D20+4.*m_EO*EDOT/BETAO2)*
		  			(C5*COSG-2.*C4*SIN2G)
		  			+((2.*C5DT*COSG-4.*C4DT*SIN2G)-XGDT1*(C5*SING1+4.*
		  			C4*COS2G))));

		TMNDDT=XNDDT*1.E9;
		TEMP=(TMNDDT*TMNDDT)-XNDT*1.E18*XNTRDT;
		PP=(TEMP+(TMNDDT*TMNDDT))/TEMP;
		GAMMA=-XNTRDT/(XNDDT*(PP-2.));
		XND=XNDT/(PP*GAMMA);
		QQ=1.-EDDOT/(EDOT*GAMMA);
		ED=EDOT/(QQ*GAMMA);
		OVGPP=1./(GAMMA*(PP+1.));
	}
	else
	{
		ISIMP=1;
		EDOT=-TOTHRD*XNDTN*(1.-m_EO);
	}
}

/*-----------------------------------------------------------------------------
 *					SGP8::UpdateECIPosition		2004-11-23*/
/**  Updates the satellite position and velocity using the current orbital
 * elements.
**/
/*---------------------------------------------------------------------------*/

void SGP8::UpdateECIPosition( const nxTimeStamp &TNow )
{
   double XMAM,OMGASM,XNODES,TEMP,TEMP1,XN,EM,Z1;
   double Z7,CAPE,ZC2,AM,BETA2M;
   double SINOS,COSOS,AXNM,AYNM,PM,G1,G2,G3,BETA,G4;
   double G5,SNF,CSF,FM,SNFG,CSFG,SN2F2G,CS2F2G,ECOSF;
   double G10,RM,AOVR,G13,G14,DR,DIWC,DI,SNI2DU,XLAMB;
   double Y4,Y5,R,RDOT,RVDOT,SNLAMB,CSLAMB,UX,VX,UY,VY,UZ,VZ;
   double TSINCE;
   double SINE = 0.0;
   double COSE = 0.0;
   double ZC5  = 0.0;

   nxTimeStamp Elapsed;

   if (TNow == m_time) return;           // If last calculated m_time is same then dont do calculation.
   if (m_epoch == 0.0)                      // check for valid epoch in
   {                    // the two line elements.
      m_location.SetCoords(0.0, 0.0, 0.0);
      m_velocity.SetCoords(0.0, 0.0, 0.0);
      return;
   }

   m_time = TNow;                  // Update m_time of these elements.
   Elapsed = (TNow - m_epoch);
   TSINCE = Elapsed.MJD()*XMNPDA;

//      UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG

	XMAM=FMOD2P(m_XMO+XLLDOT*TSINCE);
	OMGASM=m_OMEGAO+OMGDT*TSINCE;
	XNODES=m_XNODEO+XNODOT*TSINCE;
	if (ISIMP != 1)
	{
		TEMP=1.-GAMMA*TSINCE;
		TEMP1= pow(TEMP,PP);
		XN=XNODP+XND*(1.-TEMP1);
		EM=m_EO+ED*(1.-pow(TEMP,QQ));
		Z1=XND*(TSINCE+OVGPP*(TEMP*TEMP1-1.));
	}
	else
	{
		XN=XNODP+XNDT*TSINCE;
		EM=m_EO+EDOT*TSINCE;
		Z1=.5*XNDT*TSINCE*TSINCE;
	}


	Z7=3.5*TOTHRD*Z1/XNODP;
	XMAM=FMOD2P(XMAM+Z1+Z7*XMDT1);
	OMGASM=OMGASM+Z7*XGDT1;
	XNODES=XNODES+Z7*XHDT1;

	//      SOLVE KEPLERS EQUATION

	ZC2=XMAM+EM*SIN(XMAM)*(1.+EM*COS(XMAM));
	for (int i=1; i<=10; i++)
	{
		SINE=SIN(ZC2);
		COSE=COS(ZC2);
		ZC5=1./(1.-EM*COSE);
		CAPE=(XMAM+EM*SINE-ZC2)*ZC5+ZC2;

		if ( fabs(CAPE-ZC2) <= 1.0E-07) break;
		ZC2=CAPE;
	}

//      SHORT PERIOD PRELIMINARY QUANTITIES

	AM= pow( (XKE/XN),TOTHRD );
	BETA2M=1.-EM*EM;
	SINOS=SIN(OMGASM);
	COSOS=COS(OMGASM);
	AXNM=EM*COSOS;
	AYNM=EM*SINOS;
	PM=AM*BETA2M;
	G1=1./PM;
	G2=.5*CK2*G1;
	G3=G2*G1;
	BETA=SQRT(BETA2M);
	G4=.25*A3COF*SINI;
	G5=.25*A3COF*G1;
	SNF=BETA*SINE*ZC5;
	CSF=(COSE-EM)*ZC5;
	FM=ACTAN(SNF,CSF);      //ACTAN Returns wrong quadrant
	SNFG=SNF*COSOS+CSF*SINOS;
	CSFG=CSF*COSOS-SNF*SINOS;
	SN2F2G=2.*SNFG*CSFG;
	CS2F2G=2.*(CSFG*CSFG)-1.0;
	ECOSF=EM*CSF;
	G10=FM-XMAM+EM*SNF;
	RM=PM/(1.+ECOSF);
	AOVR=AM/RM;
	G13=XN*AOVR;
	G14=-G13*AOVR;
	DR=G2*(UNMTH2*CS2F2G-3.*TTHMUN)-G4*SNFG;
	DIWC=3.*G3*SINI*CS2F2G-G5*AYNM;
	DI=DIWC*COSI;

	//      UPDATE FOR SHORT PERIOD PERIODICS

	SNI2DU = SINIO2*(  G3*(.5*(1.-7.*THETA2)*SN2F2G-3.*UNM5TH*G10)
	    -G5*SINI*CSFG*(2.+ECOSF))-.5*G5*THETA2*AXNM/COSIO2;

	XLAMB  = FM+OMGASM+XNODES+
	  G3*(.5*(1.+6.*COSI-7.*THETA2)*SN2F2G-3.*(UNM5TH+2.*COSI)*G10)+
	  G5*SINI*(COSI*AXNM/(1.+COSI)-(2.+ECOSF)*CSFG);

	Y4=SINIO2*SNFG+CSFG*SNI2DU+.5*SNFG*COSIO2*DI;
	Y5=SINIO2*CSFG-SNFG*SNI2DU+.5*CSFG*COSIO2*DI;
	R=RM+DR;
	RDOT=XN*AM*EM*SNF/BETA+G14*(2.*G2*UNMTH2*SN2F2G+G4*CSFG);
	RVDOT=XN*(AM*AM)*BETA/RM+G14*DR+AM*G13*SINI*DIWC;

//      ORIENTATION VECTORS

	SNLAMB=SIN(XLAMB);
	CSLAMB=COS(XLAMB);
	TEMP=2.*(Y5*SNLAMB-Y4*CSLAMB);
	UX=Y4*TEMP+CSLAMB;
	VX=Y5*TEMP-SNLAMB;
	TEMP=2.*(Y5*CSLAMB+Y4*SNLAMB);
	UY=-Y4*TEMP+SNLAMB;
	VY=-Y5*TEMP+CSLAMB;
	TEMP=2.*SQRT(1.-Y4*Y4-Y5*Y5);
	UZ=Y4*TEMP;
	VZ=Y5*TEMP;

	//      POSITION AND VELOCITY

	m_location.SetCoords( R*UX, R*UY, R*UZ );
	m_location    = m_location*(XKMPER*1000.0);

	m_velocity.SetCoords( RDOT*UX+RVDOT*VX, RDOT*UY+RVDOT*VY, RDOT*UZ+RVDOT*VZ );
	m_velocity    = m_velocity*(XKMPER*1000.0*XMNPDA/86400.0);

	NutateEquatorialCoords( &m_location, m_time);
	NutateEquatorialCoords( &m_velocity, m_time);

}

