// @doc SpaceAndTime

#include "nxbase_geodesy.h"
using namespace::nxmath;
//---------------------------------------------------------------------------
//						KeplerOrbit::OrbitalPeriod
//	@mfunc
//	Calculate the orbital period of this orbit.
//
//	@rdesc
//	Returns the Orbital period in seconds
//---------------------------------------------------------------------------

double KeplerOrbit::OrbitalPeriod( 
								  double apogee)	// @parm Orbital apogee in metres from orbit focus
{
	return (TWOPI*sqrt( apogee*apogee*apogee/m_mu) );
}

//---------------------------------------------------------------------------
//						KeplerOrbit::KeplerOrbit
//	@mfunc
//	Default Constructor.  Configures this object to calculate orbits around the Earth.
//	Other planets must modify the member m_mu to account for changes in gravity
//---------------------------------------------------------------------------

KeplerOrbit::KeplerOrbit()
{
	m_mu    = 3.98601210E14;		// The gravitational parameter for Earth;
	m_N0    = 0.0;					// Mean motion
	m_M0    = 0.0;					// Mean Anomaly	at the specified epoch
	m_e     = 0.0;					// eccentricity
	m_a     = 0.0;					// semi-major axis
	m_b     = 0.0;					// semi minor axis
	m_h     = 0.0;
}


//---------------------------------------------------------------------------
//						safe_acos
//---------------------------------------------------------------------------

/*
static double safe_acos( double cosa )
{
	cosa = nxmax( cosa, -1.0);
	cosa = nxmin( cosa,  1.0);
	return acos(cosa);
}
*/

//---------------------------------------------------------------------------
//						KeplerOrbit::EccentricAnomaly
//	@mfunc
//	Calculate the Eccentric anomaly by solving Kepler's equation using a
//	Newton-Raphsom method.  Iterate to a precision specified by eps which
//	is typically 1.0E-07 to 1.0E-10
//
//---------------------------------------------------------------------------

double KeplerOrbit::EccentricAnomaly(
									  double M, 		// @parm Mean Anomaly
									  double ecc,		// @parm Eccentricy
									  double eps )		// @parm Precision, typically 1.0E-07 to 1.0E-10

{
	double E;
	double delta;

	E = M;									// first guess (which is exact for a circle)
	do										// Now iterate until we have a solution
	{										// so
		delta = E - ecc * sin(E) - M;		// estimate the correction
		E -= delta / (1 - ecc * cos(E));		// get the corrected guess
	} while (fabs(delta) >= eps);			// and see if we have converged.
	return E;
}

//---------------------------------------------------------------------------
//						KeplerOrbit::UpdateECIPosition
// @mfunc
// Calculates the position and velocity at the specified instant in time
//---------------------------------------------------------------------------

void KeplerOrbit::UpdateECIPosition( 
									const nxTimeStamp& mjd )			//@parm Update position to this instant in time
{
	if (mjd == m_time) return;
	double dt = (mjd.MJD()-m_epoch.MJD())/nxTimeStamp::ONESECOND;		// Get the difference in seconds
	double M  = m_M0 + m_N0*dt;							// Get the new mean Anomaly;
	double E  = EccentricAnomaly( M, m_e, 1.0E-10 );	// Solve Keplers equation to a precision of 1.0E-10 radians.

	double x  = m_a*(cos(E) - m_e);				// Get the x coordinate of the location
	double y  = m_b*sin(E);						// Get the y coordinate of the location
	double r  = sqrt( x*x + y*y );						// Get the radial distance 
	double k  = m_mu/m_h;
	double vx = -k*y/r;
	double vy =  k*(x/r+m_e);

	m_time     = mjd;
	m_location = (m_xunit*x  + m_yunit*y);				// return the location.
	m_velocity = (m_xunit*vx + m_yunit*vy);				// return the velocity
}

//---------------------------------------------------------------------------
//						KeplerOrbit::FromElements
//@mfunc
//	Update the object so it calculates orbits derived from the specified Keplerian elements.
//---------------------------------------------------------------------------

void KeplerOrbit::FromElements( 
							    double  mjd, 				// @parm Epoch of the elements
								int     orbitnumber,		// @parm Orbit Number at epoch
								double  I0,					// @parm Inclination in degrees (between 0 and 180 )
								double  RAAN,				// @parm Right Ascension of the Ascending Node in degrees (0 to 360.0)
								double  W0,					// @parm Argument of perigree in degree (0 to 360)
								double  E0,					// @parm Eccentricity (0 to 1)
								double  N0,					// @parm Mean motion in revs per day
								double M0 )					// @parm M0	Mean anomaly in degrees at time mjd (0-360 degrees)

{

	nxVector	ecix;
	nxVector	eciy;
	nxVector	eciz(0,0,1);
	nxVector	ytemp;
	double		e2;

	ecix.SetCoords( cosd(RAAN),     sind(RAAN), 0 );				// Get vector to ascending node
	ytemp.SetCoords( cosd(RAAN+90), sind(RAAN+90), 0 );				// Get vector perpendicular to ascending node i the equatorial plane
	eciy = ytemp*cosd(I0) + eciz*sind(I0); 							// Get vector perpendicular to ascending node but in the orbital plane
	m_zunit = (ecix^eciy).UnitVector();								// Get vector perpendicular to the orbital plane.
	m_xunit = (ecix*cosd(W0) + eciy*sind(W0)).UnitVector();			// Get the unit vector pointing at the perigree (Which is the semi major axis)
	m_yunit = m_zunit^m_xunit;

	m_N0 = (N0*TWOPI)/86400.0;										// Get mean motion in radians per second
	m_M0 = M0*ONE_DEGREE;
	m_e  = E0;

	m_a = pow( m_mu/(m_N0*m_N0), 1.0/3.0);
	e2  = 1.0-m_e*m_e;
	m_b = m_a*sqrt(e2);
	m_h = sqrt( m_mu*m_a*e2);

	m_DaysPerRev  = 1.0/N0;
	m_time = 0.0;
	m_epoch = mjd;
	UpdateECIPosition( m_epoch );
	m_StartOrbitNumber = orbitnumber;			// Orbit number at time "StartofOrbit"
	m_StartOfOrbit     = 0.0;					// Start time of orbitnumber.

}

//---------------------------------------------------------------------------
//						FromStateVector
//	@mfunc
//	Update the object so it calculates orbits derived from the specified
//	State vector
//		
//---------------------------------------------------------------------------

void KeplerOrbit::FromStateVector( 
								   double	mjd,			// @parm Epoch Of the elements
								   int		orbitnumber,	// @parm Orbit number at Epoch
								   nxVector r,				// @parm Geocentric ECI coordinates at epoch (in metres)
								   nxVector v)				// @parm Geocentric ECI velocity at epoch in m/s
{
	nxVector	h;
	nxVector	ev;								// The eccentricity vector
	nxVector	runit;							// the r unit vector
	double      e2;

	m_epoch = mjd;
	h       = r^v;								// Get the angular momentum vector (which is constant for the orbit)
	runit   = r.UnitVector();					// Get the radial unit vector
	m_h     = h.Magnitude();					// Get the magnitude of the angular velocity
	m_zunit = h/m_h;							// Get the angular momementum unit vector (aka z axis of the orbital plane)
	ev      = (v^h)/m_mu - runit;				// Get the eccentricity vector
 	m_e     = ev.Magnitude();					// Get the eccentricity of the ellipse.
	if (m_e <= 0.0)								// if the eccentricity is zero (or even less due to round off
	{											// then
		m_e = 0.0;								// set the eccentricity to 0
		ev  = runit;							// and let the perigree be at this point
	}											// and that is that
	m_xunit = ev.UnitVector();					// Get the eccentricity unitvector which points to the perigree (ie 0 degrees Mean Anomaly)
	m_yunit = m_zunit^m_xunit;					// Get the y axis of the orbital plane from z cross x
	e2      = 1.0 - m_e*m_e;					// 1 - e**2
	m_a     = (h & h)/(m_mu*e2);				// Get the semi major axis	
	m_b     = m_a*sqrt( e2 );					// Get the semi-minor axis;

	m_N0  = sqrt(m_mu/(m_a*m_a*m_a));			// Get the mean motion  (radians per second)  

	double x = r & m_xunit;						// Get the x component parallel to the semi-major axis
	double y = r & m_yunit;						// Get the y component parallel to the semi-minor axis
	double sinE = y/m_b;
	double cosE = x/m_a + m_e;
	double E = atan2( sinE, cosE );				// use y = b*sin(E) and x = a*(cos(E)-e) to get E
	m_M0 = E - m_e*sinE;						// now get the mean anomaly at the epoch.

	m_DaysPerRev  = TWOPI/(m_N0*86400.0);		// Period of the satellite in days
	m_StartOrbitNumber = orbitnumber;			// Orbit number at time "StartofOrbit"
	m_StartOfOrbit     = 0.0;


	
}
