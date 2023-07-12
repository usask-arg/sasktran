/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_geodesy.h"
using namespace::nxmath;
using namespace::nxSI;				// Use the SI Physical constants

#define ARC_SECONDS_RAD   206264.8062						//@define Converts radians to arc-seconds to 9 decimal places (MACRO)


/*	Calculate position of Sun.
 *	--------------------------
 *	Algorithms unashamedly copied from "Astronomy on the Personal Computer"
 *	O. Montenbruck and T. Pleger, Springer Verlag. QB51.3 E43 M6613 1991
 *	ISBN 3-540-52754-0 or ISBN 0-387-52754-0.
 *
 *	Their code is in PASCAL and I have ported that code to a C++ environment.
 *	I have checked the Sun algorithms against the Astronomical Almanac at
 *	randomly chosen times.  It was always with 1 arc second of the Almanac's
 *	value. Pretty good huh?
 *
 *	Note that INTEGER must be "long" otherwise it bombs out.
 *
*/

//-------------------------------------------------------------------------------
//			PlanetSun::Constructor
//-------------------------------------------------------------------------------

PlanetSun::PlanetSun()
{
   m_DistantObject = 0;		//Flag object as not very distant.
}

//-------------------------------------------------------------------------------
//			PlanetarySun::TERM
//-------------------------------------------------------------------------------

void PlanetSun::TERM( int I1, int I,  int IT,
		      double DLC, double DLS, double DRC,
	              double DRS, double DBC, double DBS)
{
   if (IT == 0)
   {
      ADDTHE( C3[I1+1], S3[I1+1], C[I+8], S[I+8], U,V);
   }
   else
   {
      U = U*T;
      V = V*T;
   }
   DL = DL + DLC*U + DLS*V;
   DR = DR + DRC*U + DRS*V;
   DB = DB + DBC*U + DBS*V;
}

//-------------------------------------------------------------------------------
//			PlanetSun::PERTVEN
//-------------------------------------------------------------------------------

void PlanetSun::PERTVEN()			  // Keplerian terms and perturbations by Venus
{
   int I;

   C[0+8]  =  1.0;
   S[0+8]  =  0.0;
   C[-1+8] =  cos(M2);
   S[-1+8] = -sin(M2);
   for (I =-1; I >= -5; I--)
   {
      ADDTHE(C[I+8],S[I+8],C[-1+8],S[-1+8],C[I-1+8],S[I-1+8]);
   }
   TERM(1, 0,0,-0.22,6892.76,-16707.37, -0.54, 0.00, 0.00);
   TERM(1, 0,1,-0.06, -17.35,    42.04, -0.15, 0.00, 0.00);
   TERM(1, 0,2,-0.01,  -0.05,     0.13, -0.02, 0.00, 0.00);
   TERM(2, 0,0, 0.00,  71.98,  -139.57,  0.00, 0.00, 0.00);
   TERM(2, 0,1, 0.00,  -0.36,     0.70,  0.00, 0.00, 0.00);
   TERM(3, 0,0, 0.00,   1.04,    -1.75,  0.00, 0.00, 0.00);
   TERM(0,-1,0, 0.03,  -0.07,    -0.16, -0.07, 0.02,-0.02);
   TERM(1,-1,0, 2.35,  -4.23,    -4.75, -2.64, 0.00, 0.00);
   TERM(1,-2,0,-0.10,   0.06,     0.12,  0.20, 0.02, 0.00);
   TERM(2,-1,0,-0.06,  -0.03,     0.20, -0.01, 0.01,-0.09);
   TERM(2,-2,0,-4.70,   2.90,     8.28, 13.42, 0.01,-0.01);
   TERM(3,-2,0, 1.80,  -1.74,    -1.44, -1.57, 0.04,-0.06);
   TERM(3,-3,0,-0.67,   0.03,     0.11,  2.43, 0.01, 0.00);
   TERM(4,-2,0, 0.03,  -0.03,     0.10,  0.09, 0.01,-0.01);
   TERM(4,-3,0, 1.51,  -0.40,    -0.88, -3.36, 0.18,-0.10);
   TERM(4,-4,0,-0.19,  -0.09,    -0.38,  0.77, 0.00, 0.00);
   TERM(5,-3,0, 0.76,  -0.68,     0.30,  0.37, 0.01, 0.00);
   TERM(5,-4,0,-0.14,  -0.04,    -0.11,  0.43,-0.03, 0.00);
   TERM(5,-5,0,-0.05,  -0.07,    -0.31,  0.21, 0.00, 0.00);
   TERM(6,-4,0, 0.15,  -0.04,    -0.06, -0.21, 0.01, 0.00);
   TERM(6,-5,0,-0.03,  -0.03,    -0.09,  0.09,-0.01, 0.00);
   TERM(6,-6,0, 0.00,  -0.04,    -0.18,  0.02, 0.00, 0.00);
   TERM(7,-5,0,-0.12,  -0.03,    -0.08,  0.31,-0.02,-0.01);
}

//-------------------------------------------------------------------------------
//			PlanetSun::PERTMAR
//-------------------------------------------------------------------------------

void PlanetSun::PERTMAR()			  // perturbations by Mars
{
   int I;

   C[-1+8]=cos(M4);
   S[-1+8]=-sin(M4);
   for (I=-1; I >=-7; I--)
   {
      ADDTHE(C[I+8],S[I+8],C[-1+8],S[-1+8],C[I-1+8],S[I-1+8]);
   }
   TERM(1,-1,0,-0.22,   0.17,    -0.21, -0.27, 0.00, 0.00);
   TERM(1,-2,0,-1.66,   0.62,     0.16,  0.28, 0.00, 0.00);
   TERM(2,-2,0, 1.96,   0.57,    -1.32,  4.55, 0.00, 0.01);
   TERM(2,-3,0, 0.40,   0.15,    -0.17,  0.46, 0.00, 0.00);
   TERM(2,-4,0, 0.53,   0.26,     0.09, -0.22, 0.00, 0.00);
   TERM(3,-3,0, 0.05,   0.12,    -0.35,  0.15, 0.00, 0.00);
   TERM(3,-4,0,-0.13,  -0.48,     1.06, -0.29, 0.01, 0.00);
   TERM(3,-5,0,-0.04,  -0.20,     0.20, -0.04, 0.00, 0.00);
   TERM(4,-4,0, 0.00,  -0.03,     0.10,  0.04, 0.00, 0.00);
   TERM(4,-5,0, 0.05,  -0.07,     0.20,  0.14, 0.00, 0.00);
   TERM(4,-6,0,-0.10,   0.11,    -0.23, -0.22, 0.00, 0.00);
   TERM(5,-7,0,-0.05,   0.00,     0.01, -0.14, 0.00, 0.00);
   TERM(5,-8,0, 0.05,   0.01,    -0.02,  0.10, 0.00, 0.00);
}

//-------------------------------------------------------------------------------
//			PlanetSun::PERTJUP
//-------------------------------------------------------------------------------

void PlanetSun::PERTJUP()  // perturbations by Jupiter
{
   int I;

   C[-1+8]=cos(M5); S[-1+8]=-sin(M5);
   for (I =-1; I >= -3; I--)
   {
      ADDTHE(C[I+8],S[I+8],C[-1+8],S[-1+8],C[I-1+8],S[I-1+8]);
   }
   TERM(1,-1,0, 0.01,   0.07,     0.18, -0.02, 0.00,-0.02);
   TERM(0,-1,0,-0.31,   2.58,     0.52,  0.34, 0.02, 0.00);
   TERM(1,-1,0,-7.21,  -0.06,     0.13,-16.27, 0.00,-0.02);
   TERM(1,-2,0,-0.54,  -1.52,     3.09, -1.12, 0.01,-0.17);
   TERM(1,-3,0,-0.03,  -0.21,     0.38, -0.06, 0.00,-0.02);
   TERM(2,-1,0,-0.16,   0.05,    -0.18, -0.31, 0.01, 0.00);
   TERM(2,-2,0, 0.14,  -2.73,     9.23,  0.48, 0.00, 0.00);
   TERM(2,-3,0, 0.07,  -0.55,     1.83,  0.25, 0.01, 0.00);
   TERM(2,-4,0, 0.02,  -0.08,     0.25,  0.06, 0.00, 0.00);
   TERM(3,-2,0, 0.01,  -0.07,     0.16,  0.04, 0.00, 0.00);
   TERM(3,-3,0,-0.16,  -0.03,     0.08, -0.64, 0.00, 0.00);
   TERM(3,-4,0,-0.04,  -0.01,     0.03, -0.17, 0.00, 0.00);
}

//-------------------------------------------------------------------------------
//			PlanetSun::PERTSAT
//-------------------------------------------------------------------------------

void PlanetSun::PERTSAT()  // perturbations by Saturn
{
   C[-1+8]=cos(M6);
   S[-1+8]=-sin(M6);
   ADDTHE(C[-1+8],S[-1+8],C[-1+8],S[-1+8],C[-2+8],S[-2+8]);
   TERM(0,-1,0, 0.00,   0.32,     0.01,  0.00, 0.00, 0.00);
   TERM(1,-1,0,-0.08,  -0.41,     0.97, -0.18, 0.00,-0.01);
   TERM(1,-2,0, 0.04,   0.10,    -0.23,  0.10, 0.00, 0.00);
   TERM(2,-2,0, 0.04,   0.10,    -0.35,  0.13, 0.00, 0.00);
}

//--------------------------------------------------------------------------------
//			PlanetSun::PERTMOO
//--------------------------------------------------------------------------------

void PlanetSun::PERTMOO()  // difference between the Earth-Moon
{
 			             // barycenter and the center of the Earth
      DL = DL +  6.45*sin(D) - 0.42*sin(D-A) + 0.18*sin(D+A)
			      + 0.17*sin(D-M3) - 0.06*sin(D+M3);
      DR = DR + 30.76*cos(D) - 3.06*cos(D-A)+ 0.85*cos(D+A)
			      - 0.58*cos(D+M3) + 0.57*cos(D-M3);
      DB = DB + 0.576*sin(UU);
}



//-----------------------------------------------------------------------
//			PlanetSun::EclipticCoords
//
//	returns the apparent ecliptic coordinates of the sun wrt mean
//	equinox of date.  The apparent position of the sun accounts for
//	the annual aberration due to the velocity of the earth around the
//	sun.
//	its position in the location vector.  It returns the position
//	with reference to the mean equinox of date and assumes the
//	timestamp is Dynamical Time rather than UT.
//
//	see program  SUN200:
//	in Astronomy on the Personal Computer, O. Montenbruck anfd T. Pfleger.
//
//	l = ecliptic longitude in radians
//	r= distance in AU then converted to kilometers
//	b= ecliptic latitude in radians.
//	
//-----------------------------------------------------------------------

void PlanetSun::EclipticCoords( const nxTimeStamp Tnow)
{
    int    I;
    double l;				// ecliptic longitude
    double r;				// distance in AU
    double b;				// ecliptic latitude.
    double cosb;


    T = Tnow.JD2000Centuries();			// convert to Julian centuries since JD2000

    DL = 0.0;
    DR = 0.0;
    DB = 0.0;
    M2 = TWOPI*FRAC(0.1387306+162.5485917*T);
    M3 = TWOPI*FRAC(0.9931266+ 99.9973604*T);
    M4 = TWOPI*FRAC(0.0543250+ 53.1666028*T);
    M5 = TWOPI*FRAC(0.0551750+  8.4293972*T);
    M6 = TWOPI*FRAC(0.8816500+  3.3938722*T);
    D  = TWOPI*FRAC(0.8274+1236.8531*T);
    A  = TWOPI*FRAC(0.3749+1325.5524*T);
    UU = TWOPI*FRAC(0.2591+1342.2278*T);
    C3[0+1]  =  1.0;
    S3[0+1]  =  0.0;
    C3[1+1]  =  cos(M3);
    S3[1+1]  =  sin(M3);
    C3[-1+1] =  C3[1+1];
    S3[-1+1] = -S3[1+1];
    for ( I=2; I<= 7; I++)
    {
       ADDTHE( C3[I-1+1], S3[I-1+1], C3[1+1], S3[1+1], C3[I+1], S3[I+1]);
    }
    PERTVEN();
    PERTMAR();
    PERTJUP();
    PERTSAT();
    PERTMOO();
    DL = DL + 6.40*sin(TWOPI*(0.6983+0.0561*T))
	    + 1.87*sin(TWOPI*(0.5764+0.4174*T))
	    + 0.27*sin(TWOPI*(0.4189+0.3306*T))
	    + 0.20*sin(TWOPI*(0.3581+2.4814*T));

    l =  TWOPI*FRAC(0.7859453+M3/TWOPI+((6191.2+1.1*T)*T+DL)/1296.0E3 );	// Ecliptic longitude in radians.
    r =  (1.0001398 - 0.0000007*T  +  DR*1.0E-6);//*ASTRONOMICALUNIT;  	// Radial distance in AU.
    l -= 20.496/(r*3600.0)*ONE_DEGREE;					// account for Earth's annual aberration
    b =  DB/ARC_SECONDS_RAD;						// Ecliptic latitude in radians.
    r *= ASTRONOMICALUNIT;

   				
    cosb = r*cos(b);
    m_location.SetCoords( cosb*cos(l), cosb*sin(l), sin(b)*r );

}

//----------------------------------------------------------------------------
//			PlanetSun::UpdateECIPosition
//----------------------------------------------------------------------------

void PlanetSun::UpdateECIPosition( const nxTimeStamp &Tnow )
{
   if (Tnow == m_time) return;
   m_time = TDT( Tnow );							// Convert  the current time to True Dynamical Time
   EclipticCoords( m_time );						// Calculate position of sun using true dynamical time wrt mean equinox of date
   ConvertToEquatorialCoords(nxFALSE);				// Convert these to equatorial coordinates using mean obliquity of ecliptic.
   NutateEquatorialCoords(&m_location, m_time);		// Get apparent equatorial coords of body.
   m_time   = Tnow;									// Setup the timestamp of these coordinates in UT (not TDT).
}

//----------------------------------------------------------------------------
//			AstronomicalDark
//----------------------------------------------------------------------------

nxBOOL PlanetSun::AstronomicalDark( const nxTimeStamp &Tnow, const nxGeodetic &Site )
{
   return BelowHorizon( Tnow, Site, 102.0 );
}


/*-----------------------------------------------------------------------------
 *					PlanetSun::AUDistance		 2015- 4- 17*/
/** **/
/*---------------------------------------------------------------------------*/

double PlanetSun::AUDistance( double mjd )
{
	nxTimeStamp		tnow( mjd);
	double			x;

	UpdateECIPosition( tnow );
	x = Location().Magnitude()/149597871000.0;
	return x;
}



