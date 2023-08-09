#pragma once
#include <sasktran2/internal_common.h>
#include <sasktran2/grids.h>

namespace sasktran2 {
    /** Enum defining the type of geometry used in the calculation.  Possible options are plane parralel, spherical,
     *  and ellipsoidal.  Note that currently only spherical is implemented for most classes, this is only included
     *  at the moment for future extensibility.
     */
	enum geometrytype {
		planeparallel,
		spherical,
		ellipsoidal
	};

    /** Struct defining a single location inside the atmosphere.  The core of this is a x,y,z position coordinate,
     *  and helper functions for things like radius.  The location can also hold optional information to aid in
     *  interpolation.
     */
    struct Location {
        Eigen::Vector3d position; /**< x,y,z location coordinates */

        bool on_exact_altitude; /**< OPTIONAL.  True if the location is on an exact grid point of the global altitude grid */
        int lower_alt_index;    /**< OPTIONAL.  the lower altitude index of this point in the global altitude grid.  Set to -1 if not used. */

        Location() : on_exact_altitude(false), lower_alt_index(-1) {};

        /**
         *
         * @return The norm of the location.
         */
        double radius() const {
            return position.norm();
        }

        /**
         *
         * @param other Local vector at the Location
         * @return cosine of zenith angle from the location to the other vector
         */
        double cos_zenith_angle(const Eigen::Vector3d& other) const {
            return position.dot(other) / (position.norm() * other.norm());
        }
    };

    /** Defines the coordinate system used for the radiative transfer calculation.  We define Local Coordinates
     *  as a cartesian coordinate system defined by x=reference_plane unit vector, z=reference_point unit vector,
     *  and y = z cross x.  Local cartesian coordinates are used when using plane parallel mode.
     *
     *  Local angular coordinates are defined as, for a positional vector v, |r| = |v| - reference point radius,
     *  theta = angle of the point projected into the plane formed by x and z, phi = angle of the point projected
     *  into the plane formed by y and z
     *
     *  All radiative transfer calculations internally happen on this local coordinate system, which includes
     *  specification of the atmosphere.
     *
     *  For most radiative transfer calculations we need a singular sun position, the Coordinates object is also
     *  responsible for creating this sun unit vector.
     *
     *  The preferred method of constructing the Coordinates is to only specify the solar zenith angle and solar azimuth
     *  angle at the
    */
    class Coordinates {
    private:
        Eigen::Vector3d m_x_unit;    /**< Reference plane unit vector */
        Eigen::Vector3d m_y_unit;    /**< z cross x */
        Eigen::Vector3d m_z_unit;    /**< Reference point unit vector */

        Eigen::Vector3d m_sun_unit;  /**< Unit vector towards the sun */

        double m_earth_radius;       /**< Earth radius in meters */
        geometrytype m_geotype;      /**< Enum indicating the type of geometry (spherical or plane parallel) */

    public:
        /** Constructs the Coordinates from a set of solar angles at the reference point, this is the preferred and
         *  default method of constructing the geometry.  Here solar azimuth angle is defined relative to the
         *  reference plane direction.  For 1d calculations this is typically ambiguous, and it may be preferred
         *  to simply set the saa to 0, then the relative azimuth angle between any line of sight and the sun
         *  is determined from the line of sight azimuth angles.  If all lines of sight are expected to have the
         *  same relative azimuth angle you may want to set the solar azimuth here to a non-zero value and leave
         *  the line of sight viewing azimuth at 0.
         *
         *  For 2d/3d calculations it is necessary to properly set the solar azimuth angle relative to
         *  the viewing plane, since here the viewing plane also influences the definition of the 2d/3d atmosphere.
         *  For example in a 2d calculation where all lines of sight are within the direction the atmosphere varies,
         *  you would specify the line of sight azimuth of 0, and set the solar azimuth to be the relative azimuth.
         *
         *
         * @param cos_sza Cosine of solar zenith angle
         * @param saa Solar azimuth angle in [radians]
         * @param earth_radius Earth radius in [m]
         * @param geotype Type of geometry (spherical or plane parallel) defaults to spherical
         */
        Coordinates(double cos_sza,
                    double saa,
                    double earth_radius,
                    geometrytype geotype=geometrytype::spherical
                 );

        /** Constructs the coordinates by manually specifying the unit vectors for the reference point and planes.
         *
         * @param ref_point_unit Unit vector to the reference point, becomes m_z
         * @param ref_plane_unit Unit vector along the reference plane, becomes m_x
         * @param sun_unit Unit vector towards the sun
         * @param earth_radius Earth radius in [m]
         * @param geotype Type of geometry (spherical or plane parallel) defaults to spherical
         */
        Coordinates(Eigen::Vector3d ref_point_unit,
                    Eigen::Vector3d ref_plane_unit,
                    Eigen::Vector3d sun_unit,
                    double earth_radius,
                    geometrytype geotype=geometrytype::spherical
                    );

        /**
         * @return The geometry type (plane parallel, spherical, ellipsoidal) for the coordinate system
         */
        geometrytype geometry_type() const { return m_geotype; }

        /** Calculates a unit vector in the local coordinate system from the angles (theta, phi).
         *
         * @param theta Angle in radians in the along reference plane dimension
         * @param phi Angle in radians in the across reference plane dimension
         * @return Unit vector that has coordinates theta, phi
         */
        Eigen::Vector3d unit_vector_from_angles(double theta, double phi) const;

        /**
         *
         * @param theta Angle in radians in the along reference plane dimension
         * @param phi Angle in radians in the across reference plane dimension
         * @return
         */

        std::pair<Eigen::Vector3d, Eigen::Vector3d> local_x_y_from_angles(double theta, double phi) const;

        /** Transforms a unit vector in the local coordinate system to the
         *
         * @param uv Unit vector in the local coordinate systems
         * @return Angular pair (theta, phi) both in radians
         */
        std::pair<double, double> angles_from_unit_vector(const Eigen::Vector3d& uv) const;

        /** Returns the earth radius.
         *
         * @return The radius of the earth in [m]
         */
        double earth_radius() const { return m_earth_radius; }

        /** Returns the sun unit vector in the local coordinate system
         *
         * @return The sun unit vector
         */
        const Eigen::Vector3d& sun_unit() const { return m_sun_unit; }

        /** Returns a vector pointing to a location parameterized by solar angles
         *
         *  @param cos_sza Cosine of solar zenith angle
         *  @param saa Solar azimuth angle in radians
         *  @param altitude Altitude above earth_radius in [m]
         *  @return Vector pointing to the specified location
         */
        Eigen::Vector3d solar_coordinate_vector(double cos_sza, double saa, double altitude) const;

        /** Calculates the solar angles corresponding to a geographic location
         *
         *  @param location Geographic location
         *  @return Pair (cos_sza, saa) with saa in radians
         */
        std::pair<double, double> solar_angles_at_location(const Eigen::Vector3d& location) const;

        /** Calculates the solar zenith angle at the reference point
         *
         * @return cosine of solar zenith angle at the reference point
         */
        double cos_sza_at_reference() const { return solar_angles_at_location(m_z_unit).first; }

        /** Constructs a look vector at a given location with a specified relative azimuth angle to the sun
         *
         *  @param location Geographic location
         *  @param saa Solar azmiuth angle in radians
         *  @param cos_viewing Cosine of viewing angle, positive is upwards?
         *  @return Look vector with a given relative azimuth angle to the sun
         */
         Eigen::Vector3d look_vector_from_azimuth(const Eigen::Vector3d& location, double saa, double cos_viewing) const;
    };

    /** Base class that defines the Geometry for the calculation.  The Geometry object contains a reference to the
     *  coordinate transform, as well as defines the geometry grids for the calculation.  For example, a 1d geometry
     *  object will contain an altitude grid.
     */
    class Geometry {
    private:
        const Coordinates m_coordinates;

    public:
        Geometry(Coordinates&& coordinates) : m_coordinates(coordinates) {};

        /**
         *
         * @return Reference to the coordinates object
         */
        const Coordinates& coordinates() const { return m_coordinates; }

        /**
         *
         * @return The number of atmosphere dimensions
         */
        virtual int num_atmosphere_dimensions() const = 0;

        /**
         *
         * @return Total number of grid points in the atmosphere
         */
        virtual int size() const = 0;

        /** Constructs interpolation weights on the geometry grid.
         *
         * @param loc Location to interpolate to
         * @param index_weights Pairs of index/weight objects
         */
        virtual void assign_interpolation_weights(const Location& loc, std::vector<std::pair<int, double>>& index_weights) const = 0;
    };

    /** Derived Geometry object which contains a 1D altitude grid.  Usually serves as the Base geometry object in
     *  most places since almost all calculations require at least an altitude grid to be implemented.
     */
    class Geometry1D : public Geometry {
    private:
        const grids::AltitudeGrid m_alt_grid;
    public:
        Geometry1D(Coordinates&& coordinates, grids::AltitudeGrid&& alt_grid) : Geometry(std::forward<Coordinates&&>(coordinates)), m_alt_grid(alt_grid) {}

        const grids::AltitudeGrid& altitude_grid() const { return m_alt_grid; }

        int num_atmosphere_dimensions() const { return 1; }
        int size() const { return (int)m_alt_grid.grid().size(); }

        virtual void assign_interpolation_weights(const Location& loc, std::vector<std::pair<int, double>>& index_weights) const override;
    };

    /** Not currently implemented.  Planned to be Altitude/Angle along reference plane.
     *
     */
    class Geometry2D : public Geometry1D {
        int num_atmosphere_dimensions() const { return 2; }
    };

    /** Not currently implemented.  Planned to be Altitude/Angle along reference plane/Angle across reference plane.
     *
     */
    class Geometry3D : public Geometry2D {
        int num_atmosphere_dimensions() const { return 3; }
    };

}