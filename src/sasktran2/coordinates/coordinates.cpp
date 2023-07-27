#include <sasktran2/coordinates/coordinates.h>
#include <sasktran2/coordinates/geodetic.h>
#include <sasktran2/math/vector.h>

using namespace sasktran2::math;

namespace sasktran2::coordinates{

    HeliodeticUnitVector& HeliodeticUnitVector::from_vector	(const HeliodeticVector &v)
    {
        double	f = 1.0/v.Magnitude();;

        set_coords( v.X()*f, v.Y()*f, v.Z()*f );
        return *this;
    }

    void HeliodeticVector::SetCoords	( double x, double y, double z )
    {
        m_data(0) = x;
        m_data(1) = y;
        m_data(2) = z;
    }

    void HeliodeticVector::SetCoords( const HeliodeticUnitVector& unit, double magnitude )
    {
        m_data(0) = unit.X()*magnitude;
        m_data(1) = unit.Y()*magnitude;
        m_data(2) = unit.Z()*magnitude;
    }

    void HeliodeticVector::SetCoords	( const HeliodeticUnitVector& unit)
    {
        m_data(0) = unit.X();
        m_data(1) = unit.Y();
        m_data(2) = unit.Z();
    }

    HeliodeticUnitVector HeliodeticVector::UnitVector() const
    {
        HeliodeticUnitVector	unit;
        double					f;

        f = Magnitude();
        if (f > 0.0) f = 1.0/f;
        unit.set_coords( m_data[0]*f, m_data[1]*f, m_data[2]*f );
        return unit;
    }

    void HeliodeticPoint::Initialize( const HeliodeticUnitVector& unit, double magnitude, const CoordinateTransform* coords )
    {
        m_direction = unit;
        m_radius    = magnitude;
        m_heightm   = (coords != nullptr) ? coords->RadiusToAltitude(m_radius) : std::numeric_limits<double>::quiet_NaN();
    }


    void HeliodeticPoint::FromVector( const HeliodeticVector&	v, const CoordinateTransform* coords )
    {
        HeliodeticUnitVector	unit;
        double					r;
        double					f;

        r = v.Magnitude();
        f = (1.0/r);
        unit.set_coords( v.X()*f, v.Y()*f, v.Z()*f );
        Initialize( unit, r, coords );
    }

    void HeliodeticPoint::Clear()
    {
        m_direction.set_zero();
        m_radius =-9999999.0;
        m_heightm = -9999999.0;
    }

    double HeliodeticPoint::CosZenithAngle( const HeliodeticUnitVector& look ) const
    {
        double	coszen;

        coszen = look & m_direction;
        if (coszen > 1.0)
        {
            assert(( coszen < 1.000001 ));
            coszen = 1.0;
        }
        if (coszen < -1.0)
        {
            assert(( coszen > -1.000001 ));
            coszen = -1.0;
        }
        return coszen;
    }

    bool HeliodeticPoint::LocalUnitVectors( HeliodeticUnitVector* localunitarray, size_t numv)  const
    {
        double					coszen = m_direction.Z();
        double					sinzen = sqrt( 1.0 - coszen*coszen );
        double					sinzenrecip;
        HeliodeticUnitVector*	north = localunitarray;
        HeliodeticUnitVector*	west  = localunitarray + 1;
        HeliodeticUnitVector*	up    = localunitarray + 2;
        bool					ok;

        ok = (numv == 3);
        if (!ok)
        {
            BOOST_LOG_TRIVIAL(warning) << "HeliodeticPoint::LocalUnitVectors, You must pass in exactly 3 unit vetor array. You passed in " << (int)numv;
        }
        else
        {
            if (sinzen < 0.005 )				// Sun and location are with 0.28 degrees of each other
            {
                north->set_coords( -1,  0, 0 );		// North
                west-> set_coords(  0, -1, 0 );		// West
                up->   set_coords(  0,  0, 1 );		// Up
                if( coszen < 0.0 ){
                    north->negate();
                    west-> negate();
                    up->   negate();
                }
            }
            else
            {
                sinzenrecip = 1.0/sinzen;
                west->set_coords(     m_direction.Y()*sinzenrecip,		// west = direction ^ sun
                                     -m_direction.X()*sinzenrecip,		// But sun is along Z axis = (0,0,1)
                                     0 );								// Therefore cross product is really quite simple.

                north->set_coords(	west->Y()*m_direction.Z(),			//north = west ^ up
                                     -west->X()*m_direction.Z(),
                                     west->X()*m_direction.Y() - west->Y()*m_direction.X() );

                *up = m_direction;																	// Up is same as this unit vector
            }
        }
        return ok;
    }

    HeliodeticUnitVector HeliodeticPoint::TransformToLocalZenithCoords( const HeliodeticUnitVector& v, const HeliodeticUnitVector* localunitarray ) const
    {
        const HeliodeticUnitVector&	north = localunitarray[0];
        const HeliodeticUnitVector&	west  = localunitarray[1];
        const HeliodeticUnitVector&	up    = localunitarray[2];
        HeliodeticUnitVector			newv;

        double	x,y,z;

        x = v & north;
        y = v & west;
        z = v & up;
        newv.set_coords(x,y,z);
        return newv;
    }

    HeliodeticVector HeliodeticPoint::Vector() const
    {
        HeliodeticVector	v;

        v.SetCoords( m_direction, m_radius );
        return v;
    }

    HeliodeticUnitVector HeliodeticPoint::UnitVector() const
    {
        return m_direction;
    }

    CoordinateTransform::CoordinateTransform()
    {
        double nan = std::numeric_limits<double>::quiet_NaN();

        m_earthRadius = nan;
        m_latitude    = nan;
        m_longitude   = nan;
        m_mjd         = nan;
        m_groundaltitude_meters = 0.0;				// Use sensible defaults for the moment
        m_toaaltitude_meters    = 100000.0;			// Use NaN in the future
    }

    CoordinateTransform::~CoordinateTransform()
    {

    }

    bool CoordinateTransform::ConfigureCoordinates( double latitude, double longitude, double mjd, const Eigen::Vector3d& sun )
    {
        double		sza;
        Eigen::Vector3d	location;
        bool		ok;


        do																	//	Keep looping untl we are well away from the Sun
        {																	// for each loop
            m_truegeoid.FromGeodetic( latitude, longitude );				// Get the true location of this latitude and longitude
            location = m_truegeoid.Location();								// Get the geocentric location
            sza      = angle_degrees_between(location, sun);		// Get the solar zenith angle
            ok       = (sza >= 0.04);										// make sure its ok
            if (!ok) latitude += 0.05;										// if not step in small steps until it is
        } while (!ok);														// Repeat until good.

        m_truegeoid.GetOsculatingSpheroid( &m_earthRadius, &m_centre );		// Get the center and radius of the fitted osculating sphere
        m_oscgeoid.SetTrueSphere(m_earthRadius);							// Define the osculating spheroid
        m_oscgeoid.FromGeodetic( latitude,longitude, 0.0 );					// Set the location of the point
        m_referencepoint_unit = m_oscgeoid.Location().normalized();			    // Cache the unit vector towardslocation of the reference point in osculating frame
        m_latitude       = latitude;											// Cache the latitude
        m_longitude      = longitude;											// Cache the longitude
        m_mjd            = mjd;
        assert(( (sun.norm() < 1.0000000001) && (sun.norm() > 0.9999999999)));								// Throw exception if user throws in a position vector for Sun rather than a directional unit vector
        m_sun         = sun;												// Copy the sun over, its a directional unit vector so dont do any translation to osculating sphere center

        ok = ok && ConfigureGlobalTransform();
        return true;
    }

    bool CoordinateTransform::ManuallySetOsculatingSphereRadius( double radius )
    {
        double deltar = radius - m_earthRadius;						// Get the change in radius
        m_oscgeoid.SetTrueSphere(radius);							// Define the osculating spheroid to a true spher of the new radius
        m_earthRadius = radius;										// Save the new radius
        m_centre      = m_centre - deltar*m_referencepoint_unit;	// adjust offset so surface of new sphere is still at the surface of Earth
        return true;
    }

    bool CoordinateTransform::ConfigureCoordinates( const Eigen::Vector3d& observer, const Eigen::Vector3d& look, double mjd, const Eigen::Vector3d& sun )
    {
        Geodetic		geoid;
        double			latitude;
        double			longitude;
        Eigen::Vector3d		location;
        double			sza;

        geoid.FromTangentPointLocation( observer, look );
        location = geoid.Location();								// Get the tangent point location
        sza      = angle_degrees_between(location, sun);			// get the solar zenith angle
        while ( sza < 0.05 )										// While we are pretty close to the sub-solar point
        {															// we are going to move the reference point
            location = location + 5000.0*look;						// by adding 5 km in the look direction
            sza      = angle_degrees_between(location, sun);		// update the solar zenith angle
        }															// and try again.
        latitude  = geoid.GeodeticLatitude();
        longitude = geoid.GeodeticLongitude();
        return ConfigureCoordinates( latitude, longitude, mjd, sun / sun.norm() );
    }

    bool CoordinateTransform::ConfigureGlobalTransform()
    {
        bool		ok;
        Eigen::Vector3d	sun;
        Eigen::Vector3d	refpt;
        Eigen::Vector3d	zprime;
        Eigen::Vector3d	xprime;
        Eigen::Vector3d	yprime;

        zprime    = GetSunUnitVector();								// The Z axis is the direction to the sun, expressed in geographic coords
        refpt     = m_referencepoint_unit;							// Use the reference point in the spatial grid as the definition of 0 degrees longitude.
        xprime    = component_perpindicular_to(refpt, zprime);		// Get the component perpendicular to Z
        ok        = xprime.norm() > 1.0E-4;					// Are observer and sun parallel
        if (!ok)													// Yes,
        {															// Add code later to support this option
            BOOST_LOG_TRIVIAL(warning) << "SKTRAN_TableDiffusePoints_V21::ConfigureTransform, observer and sun are parallel, This is not properly supported and may cause errors";
        }
        xprime /= xprime.norm();								     // Normalize the x prime direction
        m_heliounitvector[0] = xprime;								// The Helio X axis is perpendicular to sun in plane of reference point
        m_heliounitvector[1] = zprime.cross(xprime);			// The Y axis is a third axis of right hand system
        m_heliounitvector[2] = zprime;								// The Helio Z axis is towards the sun
        return ok;
    }


/*-----------------------------------------------------------------------------
 *					CoordinateTransform::GeographicToHelio		2009-1-26*/
/** **/
/*---------------------------------------------------------------------------*/

    HeliodeticVector CoordinateTransform::GeographicToHelio( const Eigen::Vector3d& geographic ) const
    {
        double					x, y, z;
        Eigen::Vector3d				geo;
        double					r;
        HeliodeticVector		helio;

        r = geographic.norm();
        geo = geographic / r;

        x = m_heliounitvector[0].dot(geo);
        y = m_heliounitvector[1].dot(geo);
        z = m_heliounitvector[2].dot(geo);

        helio.SetCoords( x*r, y*r, z*r );
        return helio;
    }

    HeliodeticUnitVector CoordinateTransform::GeographicToHelioUnitVector( const Eigen::Vector3d& geo ) const
    {
        double					x, y, z;
        double					m;
        HeliodeticUnitVector	helio;

        m = geo.norm();
        m = (m > 0.0) ? 1.0/m : 0.0;
        x = m_heliounitvector[0].dot(geo);
        y = m_heliounitvector[1].dot(geo);
        z = m_heliounitvector[2].dot(geo);

        helio.set_coords( m*x, m*y, m*z );
        return helio;
    }

    Eigen::Vector3d  CoordinateTransform::HelioVectorToGeographic( const HeliodeticVector& helio) const
    {
        return ( m_heliounitvector[0]*helio.X()
                 + m_heliounitvector[1]*helio.Y()
                 + m_heliounitvector[2]*helio.Z());
    }

    Eigen::Vector3d CoordinateTransform::HelioUnitVectorToGeographic( const HeliodeticUnitVector& helio) const
    {
        return ( m_heliounitvector[0]*helio.X()
                 + m_heliounitvector[1]*helio.Y()
                 + m_heliounitvector[2]*helio.Z());
    }

    bool CoordinateTransform::HelioVectorToHelioPoint( const HeliodeticVector& geo, HeliodeticPoint* point) const
    {
        double				r;

        r               = geo.Magnitude();
        point->Initialize( geo.UnitVector(), r, this );
        return true;
    }

    HeliodeticPoint CoordinateTransform::ReferencePoint( double altitude_meters ) const
    {
        Eigen::Vector3d				refgeo;
        HeliodeticUnitVector	refv;
        HeliodeticPoint		location;
        double					r;

        refgeo = ReferencePointUnitVector();
        refv   = GeographicToHelio( refgeo ).UnitVector();
        r      = AltitudeToRadius(altitude_meters);
        location.Initialize( refv, r, this );
        return location;
    }

    bool CoordinateTransform::SetAtmosphereAltitudeBounds( double groundalt_meters, double toaalt_meters)
    {
        bool	ok;

        m_groundaltitude_meters = groundalt_meters;		// Altitude of the ground in meters
        m_toaaltitude_meters    = toaalt_meters;		// Altitude of the top of the atmosphere in meters.

        ok = (m_toaaltitude_meters > m_groundaltitude_meters) && (m_toaaltitude_meters > 5000.0);
        if (!ok)
        {
            BOOST_LOG_TRIVIAL(warning) << "CoordinateTransform::SetAtmosphereAltitudeBounds, the Top of the Atmosphere altitude" << m_toaaltitude_meters << "is less than the ground altitude or it is too small (did you use kms by accident) or is undefined";
        }
        return ok;
    }


}
