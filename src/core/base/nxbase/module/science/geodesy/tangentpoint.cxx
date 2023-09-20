#include "nxbase_geodesy.h"

using namespace::nxmath;
//
///*-----------------------------------------------------------------------------
// *					nxGeodetic::FromTangentPointLocation		2004-11-23*/
///**	Calculates the geocentric coordinates of the tangent point of a look vector
// *	"lookv" emanating from a point "r".  The earth is modelled as an oblate
// *	spheroid.  All calculations assume straight ray paths.
// *
// *	\par Updates
// *	2013-07-10, ndl303, The code checks to see if the observer is looking away from the
// *	tangent point. If he is then it places the "efefctive" tangent point at the center
// *	of the Earth.
// *
// *	\par Algorithm Details
// *	Basic idea is to transform from the standard (primed) geocentric coordinate system
// *	to another system such that the look vector defines the X direction, the
// *	spacecraft position is a a second vector that defines the plane that includes
// *	the straight line ray and the center of the earth.  The new Z axis is
// *	perpendicular to the plane.	Convenient because Z is exactly zero in the new
// *	coordinate frame.  Also convenient is the fact that the straight line ray
// *	is given by the equation y = yr  where yr is the (constant) Y coordinate
// *	of the spacecraft in the new reference frame.
// *
// *	First find the equation of the oblate spheroid projected onto this plane.
// *	Second define the tangent point as the place on the oblate spheroid that
// *	is parallel to the ray direction (which is aligned in the X direction)
// *
// *	Transform the tangent point back to the normal coordinate system.
// *	Input coordinate can be geographic or ECI.
// *	Return coordinate is in same system as input coordinate.
// *	
// *	\par Limitations
// *	lookv and r must not be parallel.  If they are then there is a whole family
// *	of solutions (probably the great circle ).
// **/
///*---------------------------------------------------------------------------*/
//
//nxVector nxGeodetic::FromTangentPointLocation( nxVector r, nxVector lookv )
//{
//	nxVector loc; 
//
//	if (IsPureSphere())								// IF we have a pure sphere
//	{												// then the solutionis simple
//		nxVector lunit = lookv.UnitVector();		// Get the look vector unit vector
//		loc  = r  - (r & lunit)*lunit;				// cos(theta) = -runit.lookv
//	}												// Location = spacecraft
//	else
//	{
//		double a = m_ReA;							// semi major axis (for earths X and Y axis)
//		double c = a*(1.0-m_Ref);					// semi minor axis (for earth's Z axis) 
//		nxVector xunit;		
//		nxVector yunit;
//		nxVector zunit;
//
//		// ----- ndl303, 2013-07-10. Modified the code so it handles locations that are looking away from the tangent point
//
//		nxVector	west;
//		nxVector	south;
//		nxVector	up;
//		double		lookzenith;
//		bool		ok;
//
//		FromGeocentric( r );								// Put in the location of the observer
//		GetGeodeticWestSouthUp( &west, &south, &up );		// get the local "UP" vector at the observer
//		lookzenith = up.AngleTo( lookv );					// Get the zenith angle of the look vector
//		ok = (lookzenith >= 90.0);							// Tangent points only make sense when we are looking "down"
//		if (!ok)											// So if we are looking up
//		{													// Then
//			loc.SetCoords(0,0,0);							// return tangent point at center of Earth. (NOT very good solution but I dont want to change all existing code)
//		}
//		else
//		{
//
//			xunit = lookv.UnitVector();					// define X unit vector as parallel to the look vector
//			zunit = (lookv^r).UnitVector();				// define Z unit vector as the ray tracing plane as perpendicular to the plane defined by the ray and the position vector.
//			yunit = (zunit^xunit).UnitVector();			// define Y unit vector so it forms a right angled system.
//
//			double w11 = xunit.X();						// get the transformation numbers from the unit vectors
//			double w12 = yunit.X();						// we have to do the transpose to get the proper array (primes are Earth coords, unprimed are plane coords)
//			double w21 = xunit.Y();						// x' = w11.x + w12.y + w13.z  
//			double w22 = yunit.Y();						// y' = w21.x + w22.y + w23.z
//			double w31 = xunit.Z();						// z' = w31.x + w32.y + w33.z
//			double w32 = yunit.Z();						// for our calculations z is exactly zero.
//
//			double	a2 = a*a;
//			double  c2 = c*c;
//
//			double A = (sqr(w11) + sqr(w21))/a2 + sqr(w31)/c2;				// Calculate coefficients of the equation of oblate spehroid
//			double B = (sqr(w12) + sqr(w22))/a2 + sqr(w32)/c2;		// projected onto the straight line ray
//			double C = 2.0*( (w11*w12 + w21*w22)/a2 + (w31*w32)/c2 );				// plane
//
//			double factor = 4.0*A*B - C*C;
//			double Px;
//			double Sy     = r & yunit;						// y coordinate of spacecraft 
//
//			if (factor <= 0.0) // Nick's notes: NB after equation 10
//			{
//		//		if (factor < 0) nxLog::Record(NXLOG_WARNING, "nxGeodetic::TangentPointLocation, 4AC-B*B <= 0 (actually = %g), that cannot be!!!!", (double)factor);
//				Px = 0.0;		// This corresponds to a circular Earth
//		//		Sy = 0.0;
//			}
//			else
//			{
//				Px = -C/A*sqrt( A/factor);
//			}
//			double Tx = Px;
//			double Ty = Sy;
//			loc = (Tx*xunit)  + (Ty*yunit);		// Now get Geocentric location in original coordinate system
//		}
//	}
//	FromGeocentric( loc );				// return answer to caller.
//	return loc; 
//}
//
// 
 
/*-----------------------------------------------------------------------------
 *					nxGeodetic::FromTangentPointLocation		2014-05-28*/
/**	Calculates the geocentric coordinates of the tangent point of a look vector
 *	"lookv" emanating from a point "r".  The earth is modelled as an oblate
 *	spheroid.  All calculations assume straight ray paths.
 *
 *
 *	\par Algorithm Details
 *	Diagonal matrices are affine transformations -- the ratio of distances
 *  is preserved under the transformation. Then it is sufficient to transform
 *  the ellipsoid into a sphere, find the point along the line of sight with
 *  minimum radius, and transform the point back into the original coorinate
 *  system. 
 **/
/*---------------------------------------------------------------------------*/
nxVector nxGeodetic::FromTangentPointLocation( const nxVector& obs, const nxVector& elluser )
{
    const int numiter = 5;

	nxVector	loc; 
	nxVector	ell;
	// Stretch coordinates along the z-direction

    double f = m_Ref;

	ell = elluser.UnitVector();

    for(int i = 0; i < numiter; ++i) {
        double transformFactor = 1.0 / (1.0 - f);
        nxVector obsTransf(obs.X(), obs.Y(), obs.Z()*transformFactor);
        nxVector ellTransf(ell.X(), ell.Y(), ell.Z()*transformFactor);

        // The radius is minimized at distance dot(look,obs) along the look direction
        double s = -ellTransf.Dot(obsTransf) / ellTransf.Dot(ellTransf);

        // Set point, transforming back to ellipsoid
        loc = obs + s*ell;

        // Set as geodetic location, return to caller
        FromGeocentricVector( loc );

        f = m_Ref * (m_ReA / (m_ReA + m_height));
    }

    return loc;
}


/*---------------------------------------------------------------------------
 *					nxGeodetic::GetOsculatingSpheroid              2003-4-24*/
/**	Calculates the spheroid that best fits the ellipsoid surface (height =0)
 *	at the current location. This algorithm best fits the spheroid in a 
 *	latitudinal direction (ie North South) 
 *
 *	The radius of curvature is given by:-
 *
 *		       [       (dy)^2 ] ^3/2
 *		       [   1 + (--)   ]
 *		       [       (dx)   ]
 *		R =  -----------------------
 *                d2y
 *				  ---
 *		          dx2
 *
 *	For an ellipse:
 *	R = 1/ab [ a^2y^2/b^2 + b^2x^2/a^2]^3/2
 *
 *	\par Returns
 *	The code returns the radius o ftyhe best fit osculating spheroid and the
 *	offset of the origin of the osculating sphere from the origin of the 
 *	true geoid.
 *
*/
/*-------------------------------------------------------------------------*/

bool nxGeodetic::GetOsculatingSpheroid( double* Radius, nxVector* offset )
{
	double		oldheight = m_height;
	double		y0;
	double		x0;
	double		a;
	double		b;
	double		r;
	double		theta;
	double		dx;
	double		dy;
	double		a2y0;
	double		b2x0;
	double		a2;
	double		b2;
	nxVector	xunit( cosd(m_geodeticlongitude), sind(m_geodeticlongitude), 0);
	nxVector	yunit( 0, 0, 1 );

	FromGeodetic( m_geodeticlatitude, m_geodeticlongitude, 0.0 );		// Get the location at the surface of the geoid.

	y0    = m_location.Z();												// Get the vertical component of the point on the surface of the geoid
	x0    = sqrt( sqr(m_location.X()) + sqr(m_location.Y()) );			// Get the horizontal component of the point on the surface of the geoid
	a     = m_ReA;
	b     = m_ReA*(1-m_Ref);
	a2    = a*a;
	b2    = b*b;
	a2y0  = a2*y0;
	b2x0  = b2*x0;
	r     = 1.0/(a*b)* pow( (a2y0*y0/b2 + b2x0*x0/a2), 1.5 );			// Get the radius of curvature at this location
	theta = atan2( a2y0, b2x0 );										// Get the angle of the gradient of the vertical at the surface of the geoid
	dx    = r*cos(theta);												// Get the X offset of the center of curvature from the surface point.
	dy    = r*sin(theta);												// Get the Y offset of the center of curvature from the surface point.
	*offset = m_location - dy*yunit - dx*xunit;
	*Radius = r;
	FromGeodetic( m_geodeticlatitude, m_geodeticlongitude, oldheight );		// restore the original location
	return nxTRUE;
}


