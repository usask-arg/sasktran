
#include <sasktran2/coordinates/geodetic.h>
#include <sasktran2/math/trig.h>
#include <sasktran2/math/vector.h>
#include <float.h>
#include <limits>

static const double MACHINE_PRECISION = 100.0*DBL_EPSILON;		// approx 10-14 for Microsoft double precision
using namespace sasktran2::math;

namespace sasktran2::coordinates{


    //--------------------------------------- GEOID_SPHERE    IAU1976      GRS80                MERIT83        WGS84
    const double Geodetic::m_geoid_Re[5] = {RE_GEOCENTRIC,  6378140.0,   6378137.0,           6378137.0,     6378137.0};
    const double Geodetic::m_geoid_F [5] = {0.0,            0.00335281, 1.0/298.257222101,	1.0/298.257,   1.0/298.257223563};


    void Geodetic::SelectGeoid( GEOID_MODEL id )
    {
        SetGeoid( m_geoid_Re[(int)id], m_geoid_F[(int)id] );
    }

    void Geodetic::init()
    {
        m_ReA = m_geoid_Re[(int)WGS84];			// define the Equatorial radius
        m_Ref = m_geoid_F [(int)WGS84];			// define the flattening
        m_geodeticlongitude = 0.0;					// The Geodetic/Geocentric Longitude of observer in degrees
        m_geodeticlatitude  = 0.0;					// The Geodetic latitude of the observer
        m_height            = 0.0;					// Height of observer above sea-level in metres.
        m_useexact          = true;
    }

    Geodetic::Geodetic()
    {
        init();
    }

    Geodetic::Geodetic( GEOID_MODEL id)
    {
        init();
        SelectGeoid(id);
    }

    Geodetic::Geodetic(  double latitude, double longitude, double Height )
    {
        init();
        FromGeodetic( latitude, longitude, Height );
    }

    void Geodetic::SetGeoid( double newReA, double newRef)
    {
        m_Ref = newRef;
        m_ReA = newReA;
        FromGeodetic( m_geodeticlatitude, m_geodeticlongitude, m_height );
    }

    void Geodetic::FromGeodetic( double latitude, double longitude, double aHeight )
    {
        m_geodeticlongitude = inrange( longitude, 360.0 );		// Update observers geodetic longitude
        m_geodeticlatitude  = latitude;						// geodetic latitude
        m_height            = aHeight;							// and m_height above sea level.

        double lat = DegreesToRadians(m_geodeticlatitude);		// Convert degrees to radians.
        double lon = DegreesToRadians(m_geodeticlongitude);	// Convert degrees to radians.

        double cosphi = cos( lat );
        double sinphi = sin( lat );
        double coslam = cos( lon );
        double sinlam = sin( lon );

        double fm1  = (1.0-m_Ref);
        double fm12 =  fm1*fm1;
        double C    = 1.0/sqrt( cosphi*cosphi + fm12*sinphi*sinphi );
        double S    =  fm12*C;
        double P    = (m_ReA*C+aHeight)*cosphi;

        m_location << P*coslam, P*sinlam, (m_ReA*S+aHeight)*sinphi;
    }

    static double cuberoot( double x )
    {
        const double onethird = 1.0/3.0;

        if (x >= 0.0) return pow( x, onethird);
        return -pow(fabs(x), onethird );
    }

    void Geodetic::ExactGeocentricToGeodetic( double r, double z, double* fi,  double* h)
    {
        double	G [2];
        double	DS[2];
        double	T [4];
        int		numsolutions;
        bool	getallsolns    = false;
        double	bestt;
        double  guesst = 0.0;
        double  mindt;
        double	dt;
        double	zangle;
        double	xangle;
        double	l;
        bool	gotsolution   = false;
        bool	zisnegative	  = (z < 0);						// Keep track of whether z is negative or not.

        double	s;
        double	v;
        double	a  = m_ReA;										// approx 6378137.0 for Earth;
        double	fr = m_Ref;										// approx 1.0/298.257222101 for Earth
        double	b  = a - a*fr;									// b is guaranteed same sign as z

        // ---- see if solution is close to the z axis as this needs special attention

        z      = fabs(z);										// only deal inth angles in the first quadrant (zisnegative tracks the other stuff)
        l      = sqrt( z*z + r*r );
        if (l < 1000.0*std::numeric_limits<double>::epsilon() )
        {
            *fi = 0.0;
            *h = -m_ReA;
        }
        else
        {
            zangle = (r/l);
            xangle = (z/l);
            if (zangle < 0.01)										// If we are within 1/100 of a radian
            {														// then
                getallsolns    = true;							// get all of the solutions
                guesst         = 0.0;								// t = 0 is the +ve z axis
                if (zangle  < MACHINE_PRECISION )					// If we are really close to the Z axis
                {													// within numerical precision
                    gotsolution = true;							// then we can go straight to the solution
                    *fi = Pi/2.0;									// get the positive solution
                    *h  = z-b;										// at the pole
                }
            }

            // ---- see if solution is close to the x axis as this needs special attention

            if (xangle < 0.01)										// if we are within 1/100 of a radian
            {														// then
                getallsolns    = true;							// get all of the solutions
                guesst         = 1.0;								// t=1 is the +ve x axis.
                if (xangle < MACHINE_PRECISION )					// similarly if really close to the
                {													// equator
                    gotsolution = true;
                    *fi = 0.0;											// then set the latitude
                    *h  = r-a;											// and get the height
                }
            }

            if (!gotsolution)
            {										  				//  Find solution to: t**4 + 2*E*t**3 + 2*F*t - 1 = 0
                double sqrtd;
                double E = ((z + b)*b/a - a)/r;
                double F = ((z - b)*b/a + a)/r;
                double P = (E*F + 1.0)*(4.0/3.0);						// Calculate the intermediate terms
                double Q = (E*E - F*F)*2.0;								// as specified in the Borkowski paper.
                double D = P*P*P + Q*Q;

                if(D >= 0.0)											// if (D > 0) as it usually is above ~47km from center of Earth
                {														// then we will normally want
                    sqrtd = sqrt(D);									// the cubic resolvent technique
                    s = cuberoot( sqrtd+Q );							// discussed in Borkowski paper to
                    v = P/s - s;										// to evaluate the difference of the two cube roots
                    v = -(Q + Q + v*v*v)/(3.0*P); 						// This is the formula given in the paper (equation 20)
                    if (v*v > fabs(P))									// However if accuracy is not that good
                    {													// which happens about 47km from center of Earth
                        v = cuberoot(sqrtd-Q) - cuberoot(sqrtd+Q);		// then use the direct method (equation 14a)
                    }													// which considerably improves the accuracy
                }														// otherwise if D < 0
                else													// which happens close to the
                {														// center of the Earth
                    double psqrt = sqrt(-P);							// then
                    v = 2.0*psqrt*cos(acos(Q/(P*psqrt))/3.0);			// use this alternate formula given by author (equation 14b)
                }


                double esqrt = sqrt(E*E+v);
                double g     = 0.5*(E + esqrt);							// Now get the default solution for G (equation 13 in paper)
                double ds    = g*g + (F - v*g)/(g + g - E);				// Get the contents of square root brackets in equation 12
                double ssqrt = sqrt(ds);								// and get the square root
                double t     = ssqrt - g;								// solve equation 12 for the default value of t.

                if (getallsolns)										// If we are close to the axes then we wantthe solution
                {														// that is also close to the axis, so we have to pick from solutions
                    G [0] = g;											// Get the 1st solution for G.
                    DS[0] = ds;											// Get the square root term corresponding to +ve sqrt of G
                    G [1] = 0.5*(E-esqrt);								// Get the 2nd solution for G
                    DS[1] = G[1]*G[1] + (F-v*G[1])/( G[1]+G[1] - E );	// Get the square root term correpsondig to -ve sqrt of G

                    T[0] = t;											// The first solution is the one we have already calculated
                    T[1] = -ssqrt - g;									// the second is with negative of sqrt of first value of ds

                    if (DS[1] >= 0)										// See if we have the other two real solutions or if they are complex.
                    {													// we do have them
                        T[2] =  sqrt(DS[1]) - G[1];						// so get them
                        T[3] = -sqrt(DS[1]) - G[1];						// and
                        numsolutions = 4;								// flag that we have all four real solutions
                    }													// otherwise
                    else												// we only have
                    {													// two real solutions
                        numsolutions = 2;								// and that is that
                    }
                    bestt = T[0];										// Now hunt for the solution closes to our "desired" guess estimate
                    mindt = fabs(T[0]-guesst);							// initialize with the first default solution (its often best)
                    for (int i = 1; i < numsolutions; i++)				// and then hunt through the list
                    {													// of other solutions
                        dt = fabs( T[i] - guesst );						// to see if we have any that are even better
                        if ( dt < mindt )								// if we do
                        {												// then
                            bestt = T[i];								// copy it over
                            mindt = dt;									// and set up the new minimum distance
                        }												// that is that
                    }													// check all of the solutions
                    t = bestt;											// copy over the best solution into t
                }

                double divisor   = 2.0*b*t;									// now use the best value of parameter t
                double numerator = (1.0 - t*t)*a;							// to get
                if (fabs(divisor) != 0) *fi = atan( numerator/divisor );	// the geodetic latitude
                else                    *fi = Pi/2.0;						// check for 90 degree stuff
                *h = (r - a*t)*cos(*fi) + (z - b)*sin(*fi);	// and get the geodetic altitude.
            }																// and that is that
            if (zisnegative) *fi = -(*fi);									// reverse the geodetic latitude if z was originally negative.
        }
    }

    void Geodetic::IterateGeocentricToGeodetic( double r, double z )
    {
        double e2;
        double newphi;
        double phi;
        double C;
        double sinphi;

        e2      = 2.0*m_Ref - sqr(m_Ref);						// given in the
        newphi  = atan2d( z, r );								// nautical alamanac
        if (newphi >=270.0) newphi -= 360.0;					// make sure we in -90 to +90 latitudes

        do														// now start the iteration
        {														// so
            phi    = newphi;										// start with the initial guess
            sinphi = sind(phi);									// get the sine
            C      = ( 1.0 - e2*sqr(sinphi));						// calculate C
            if ( C > 0.0)
            {
                C  = 1.0/ sqrt(C);									// calculate new C
                newphi = atan2d( (z + m_ReA*C*e2*sinphi), r );		// get improved estimate on latitude
            }
            else
            {
                newphi = 90.0*sign(sinphi);
            }

            if (newphi >= 270.0) newphi -= 360.0;					// make sure we are in -90 to +90
        } while (fabs(phi-newphi) > 0.00025);					// repeat until converged satisfactorily.

        m_geodeticlatitude = newphi;							// save the latest estimate of the Geodetic latitude

        double cp = cosd(m_geodeticlatitude);

        if (cp < MACHINE_PRECISION)
        {
            double b  = m_ReA - m_ReA*m_Ref;					// Get the semi-minor axis
            if (m_geodeticlatitude >= 0) m_height = z-b;		// then in above equator
            else 						 m_height = -b-z;
        }
        else
        {
            m_height = r/cp - m_ReA*C;							// get the m_height of the observation.
        }
    }

    void Geodetic::FromGeocentricVector( const Eigen::Vector3d &Place )
    {
        double x, y, z;
        double r;

        m_location = Place;
        x          = m_location.x();
        y          = m_location.y();
        z          = m_location.z();

        m_geodeticlongitude = inrange( atan2d( y,x ), 360.0 );				// get the geodetic longitude
        if (IsPureSphere())													// If we are using a spehrical Earth
        {																	// Then
            r = sqrt( sqr(x) + sqr(y) + sqr(z) );							// simply
            m_geodeticlatitude = asind( z/r );								// Get the latitude from the geocentric coordinates
            m_height = r- m_ReA;
        }
        else
        {
            r = sqrt( sqr(x) + sqr(y) );									// now follow the recipe
            if (m_useexact)
            {
                ExactGeocentricToGeodetic( r, z, &m_geodeticlatitude,  &m_height);
                m_geodeticlatitude *= ONE_RADIAN;
            }
            else
            {
                IterateGeocentricToGeodetic( r, z );
            }
        }
    }

    void Geodetic::FromGeocentric( double latitude, double longitude, double Height )
    {
        Eigen::Vector3d Place;
        double P;
        double H;

        P        = RE_GEOCENTRIC+Height;			// radius of m_location.
        H        = P*cosd( latitude  );			// horizontal component.
        Place << H*cosd( longitude ), H*sind( longitude ), P*sind( latitude  );
        FromGeocentricVector( Place );
    }

    void Geodetic::GetGeodeticWestSouthUp( Eigen::Vector3d *West, Eigen::Vector3d *South, Eigen::Vector3d *Up) const
    {
        Eigen::Vector3d  Horiz    ( m_location.x(), m_location.y(), 0.0 );							// vector in equatorial plane and observers meridian
        Eigen::Vector3d  Vertical ( 0.0, 0.0, 1.0 );												// vector perp to equatorial plane and in observer's meridian.
        Horiz /= Horiz.norm();														// unit vector in equatorial plane and observers meridian

        *Up    =  (Horiz*cosd( m_geodeticlatitude )) + (Vertical*sind(m_geodeticlatitude) );	// get vertical  unit vector at observer's m_location
        *South =  (Horiz*sind( m_geodeticlatitude )) - (Vertical*cosd(m_geodeticlatitude) );	// get due South unit vector at observer's m_location
        *West  =  (*South).cross((*Up));                                                        		// get due West  unit vector at observer's m_location
    }

    bool Geodetic::FromTangentAltitude( double h, const Eigen::Vector3d& spacecraftlocation, const Eigen::Vector3d& boresightplane, Eigen::Vector3d* requiredlookvector )
    {
        Eigen::Vector3d	xunit;					// x unit vector in local vertical (upward) direction
        Eigen::Vector3d	yunit;					// y unit vector paralle to velocity vector
//	nxVector	lookunit;
        double		radius;
        double		re_geocentric;
        Eigen::Vector3d	offset;
        bool		ok;
        int			numtries;

        xunit   = spacecraftlocation.normalized();									// Get unit vector to spacecraft
        yunit   = (component_perpindicular_to(boresightplane, xunit)).normalized();	// Get unit vector perpendicular to r (to define look plane)
        FromGeocentricVector( spacecraftlocation);
        GetOsculatingSpheroid( &re_geocentric, & offset );
        radius  = (spacecraftlocation-offset).norm();

        double   theta = (re_geocentric+h)/radius;									// Get angle to tangent altitude on spherical Earth
        if (theta > 1.0)
        {
            BOOST_LOG_TRIVIAL(warning) << "Geodetic::FromTangentAltitude, Your requested height is too high for this code. Tangent point is above or equal to observers altitude";
            theta = 1.0;
        }
        theta = acos(theta);
        Eigen::Vector3d lookv = yunit*cos(theta) - xunit*sin(theta);						// and initialize the look vector
        numtries = 0;
        do
        {
            FromTangentPointLocation( spacecraftlocation, lookv );					// calculate the tangent point
            double newh = Height();													// get the geodetic height o fthis point
            double dh = (newh-h);													// get the difference in height from guess to geodetic
            ok = (fabs(dh) < 0.1);													// see if we have converged
            if (!ok)																// we haven't
            {																		// so
                double dtheta = 0.8*dh/(radius*sin(theta));							// adjust the estimate of theta
                theta += dtheta;													// get a new theta
                lookv = yunit*cos(theta) - xunit*sin(theta);						// and unpdate the unit vectors
            }																		// and that is that
            numtries++;
        } while (!ok && numtries < 100);											// repeat until converged;
        *requiredlookvector = lookv;												// update the satellite look direction.
        return ok;
    }

    void Geodetic::SetTrueSphere( double radius )
    {
        m_ReA = radius;
        m_Ref = 0.0;
    }
}
