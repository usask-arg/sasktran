#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/coordinates/geodetic.h>


namespace sasktran2::coordinates {
    class HeliodeticVector;

    class HeliodeticPoint;

    class CoordinateTransform;

    class HeliodeticUnitVector {
    private:
        Eigen::Vector3d m_data;
    public:
        HeliodeticUnitVector() { m_data.setConstant(-9999.0); }

        void set_coords(double x, double y, double z) {
            m_data(0) = x;
            m_data(1) = y;
            m_data(2) = z;
        }

        void set_zero() { m_data.setZero(); }

        HeliodeticUnitVector &from_vector(const HeliodeticVector &v);

        double X() const { return m_data(0); }

        double Y() const { return m_data(1); }

        double Z() const { return m_data(3); }

        double operator&(const HeliodeticUnitVector &v2) const {
            return m_data(0) * v2.m_data(0) + m_data(1) * v2.m_data(1) + m_data(2) * v2.m_data(2);
        }

        HeliodeticUnitVector &negate() {
            m_data *= -1;
            return *this;
        }
    };

    class HeliodeticVector {
    private:
        Eigen::Vector3d m_data;

    public:
        HeliodeticVector() {
            m_data(0) = -99999;
            m_data(1) = -99999;
            m_data(2) = -99999.0;
        }

        HeliodeticVector(const HeliodeticUnitVector &unit, const double &magnitude) {
            this->SetCoords(unit, magnitude);
        }

        void SetCoords(double x, double y, double z);

        void SetCoords(const HeliodeticUnitVector &unit, double magnitude);

        void SetCoords(const HeliodeticUnitVector &unit);

        double X() const { return m_data(0); }

        double Y() const { return m_data(1); }

        double Z() const { return m_data(2); }

        HeliodeticVector &operator+=(const HeliodeticVector &v2) {
            m_data += v2.m_data;
            return *this;
        }

        HeliodeticVector &operator*=(double f) {
            m_data *= f;
            return *this;
        }

        HeliodeticVector &operator-=(const HeliodeticVector &v2) {
            m_data -= v2.m_data;
            return *this;
        }

        double Magnitude() const { return m_data.norm(); }

        HeliodeticVector operator*(double f) const {
            HeliodeticVector v;
            v.SetCoords(m_data[0] * f, m_data[1] * f, m_data[2] * f);
            return v;
        }

        HeliodeticVector operator+(const HeliodeticVector &v2) const {
            HeliodeticVector v;
            v.SetCoords(m_data[0] + v2.m_data[0], m_data[1] + v2.m_data[1], m_data[2] + v2.m_data[2]);
            return v;
        }

        HeliodeticVector operator-(const HeliodeticVector &v2) const {
            HeliodeticVector v;
            v.SetCoords(m_data[0] - v2.m_data[0], m_data[1] - v2.m_data[1], m_data[2] - v2.m_data[2]);
            return v;
        }

        double operator&(const HeliodeticUnitVector &v2) const {
            return (X() * v2.X()) + (Y() * v2.Y()) + (Z() * v2.Z());
        }

        HeliodeticUnitVector UnitVector() const;
    };


    class HeliodeticPoint {
    private:
        HeliodeticUnitVector m_direction;            // Unit vector, Sun is along Z axis
        double m_radius;                // Radius from centre of Osculating sphere in meters
        double m_heightm;                // Height in meters from surface of osculating sphere

    public:
        HeliodeticPoint() {
            m_radius = -99999;
            m_heightm = -99999;
        }

        void Clear();

        void Initialize(const HeliodeticUnitVector &v, double magnitude, const CoordinateTransform *coords);

        void InitializeFromRaw(const HeliodeticUnitVector &v, double magnitude, double heightm) {
            m_direction = v;
            m_radius = magnitude;
            m_heightm = heightm;
        }

        double CosSZA() const { return m_direction.Z(); }

        double LongitudeX() const { return m_direction.X(); }

        double LongitudeY() const { return m_direction.Y(); }

        double Altitude() const { return m_heightm; }

        double Radius() const { return m_radius; }

        const HeliodeticUnitVector &LocalZenith() const { return m_direction; }

        double CosZenithAngle(const HeliodeticUnitVector &look) const;

        bool LocalUnitVectors(HeliodeticUnitVector *localunitvectors, size_t numv) const;

        HeliodeticUnitVector
        TransformToLocalZenithCoords(const HeliodeticUnitVector &v, const HeliodeticUnitVector *localunitvectors) const;

        HeliodeticVector Vector() const;

        HeliodeticUnitVector UnitVector() const;

        void FromVector(const HeliodeticVector &v, const CoordinateTransform *coords);
    };


    class CoordinateTransform {
    private:
        Geodetic m_truegeoid;            //!< The oblate spheroid used for the real world
        Geodetic m_oscgeoid;                //!< The osculating spheroid used for the spherical approximation
        Eigen::Vector3d m_centre;                //!< The offset of the center of the osculating earth in "true earth" X,Y,Z meters
        double m_earthRadius;            //!< The radius of the osculating earth in meters
        Eigen::Vector3d m_sun;                    //!< The unit vector towards the Sun in geographic coordinates
        Eigen::Vector3d m_heliounitvector[3];    //!< local unit vectors that express "Solar Coordinate" axes x', y' and z' in geographic coordinates;
        Eigen::Vector3d m_referencepoint_unit;    //!< The unit vector towards the reference point in osculating sphere coordinates
        double m_latitude;                //!< The reference geodetic latitude for this grid definition (used to define the osculating sphere)
        double m_longitude;            //!< The reference geodetic longitude for this grid definition (used to define the osculating sphere)
        double m_mjd;
        double m_groundaltitude_meters;    //!< Altitude of the ground in meters
        double m_toaaltitude_meters;        //!< Altitude of the top of the atmosphere in meters.

    private:
        bool
        ConfigureGlobalTransform();        //!< Configure the transform to convert global look directions into local zenith azimuth


    public:
        // ---- Constructors
        CoordinateTransform();

        virtual                       ~CoordinateTransform();

        bool ConfigureCoordinates(double latitude, double longitude, double mjd, const Eigen::Vector3d &sun);

        bool ConfigureCoordinates(const Eigen::Vector3d &observer, const Eigen::Vector3d &look, double mjd,
                                  const Eigen::Vector3d &sun);

        const Eigen::Vector3d &GetSunUnitVector() const { return m_sun; };

        double AltitudeToRadius(double alt_meters) const { return (alt_meters + m_earthRadius); }

        double RadiusToAltitude(double radius_meters) const {
            return floor((radius_meters - m_earthRadius) * 1000.0 + 0.5) / 1000.0;
        } // Round to the nearest micron
        const Eigen::Vector3d &SunUnit() const { return m_sun; }

        double ReferencePtLatitude() const { return m_latitude; }

        double ReferencePtLongitude() const { return m_longitude; }

        const Eigen::Vector3d &ReferencePointUnitVector() const { return m_referencepoint_unit; }

        double ReferencePointMJD() const {
            assert((m_mjd > 10000));
            return m_mjd;
        }

        HeliodeticPoint ReferencePoint(double altitude_meters) const;

        bool ManuallySetOsculatingSphereRadius(double radius);

        Eigen::Vector3d TranslateGeoidToOsculatingSphere(const Eigen::Vector3d &truecoord) const {
            return (truecoord - m_centre);
        };

        Eigen::Vector3d TranslateOsculatingSphereToGeoid(const Eigen::Vector3d &osccoord) const {
            return (osccoord + m_centre);
        };

        bool HelioVectorToHelioPoint(const HeliodeticVector &vector, HeliodeticPoint *point) const;

        HeliodeticVector GeographicToHelio(const Eigen::Vector3d &geographic) const;

        HeliodeticUnitVector GeographicToHelioUnitVector(const Eigen::Vector3d &geo) const;

        Eigen::Vector3d HelioVectorToGeographic(const HeliodeticVector &helio) const;

        Eigen::Vector3d HelioUnitVectorToGeographic(const HeliodeticUnitVector &helio) const;

        const Geodetic &
        TrueGeoid() const { return m_truegeoid; }            //!< The oblate spheroid used for the real world
        const Geodetic &
        OsculatingGeoid() const { return m_oscgeoid; }                //!< The osculating spheroid used for the spherical approximation
        bool SetAtmosphereAltitudeBounds(double groundalt_meters, double toaalt_meters);

        double GroundAltitude() const { return m_groundaltitude_meters; }
    };


}

