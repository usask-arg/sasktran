#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/geometry.h>

namespace sasktran2::viewinggeometry {

    /** struct which contains all of the necessary information to define a viewing ray internally in the model. Typically
     *  this is just the observer location and look vector.
     *
     */
    struct ViewingRay {
        sasktran2::Location observer; /**< Observer location */
        Eigen::Vector3d look_away; /**< Look vector away from the observer */

        /**
         *
         * @return Cosine of the local viewing angle
         */
        double cos_viewing() const {
            return observer.cos_zenith_angle(look_away);
        }

    };

    /** Generally an observing viewing geometry is defined as an observer position, and a local look direction.
     *  For most problems though it is not convenient for the user to specify these quantities directly as they are
     *  coordinate system dependent, while the radiative transfer calculation is mostly coordinate system independent.
     *  Therefore we define a base class ViewingGeometryBase that provides the internal functionality to uniquely
     *  specify a line of sight for the radiative transfer calculation.  Derived classes can implement different
     *  methods that are convenient for the user to uniquely specify the line of sight.
     *
    */
    class ViewingGeometryBase {

    public:
        /** Interface function where Derived classes are given information on the internal geometry/coordinate
         * system and then are responsible for constructing an observer location/viewing vector pair.
         *
         * @param geometry The internal SASKTRAN geometry
         * @return ViewingRay
         */
        virtual ViewingRay construct_ray(const sasktran2::Coordinates& geometry) = 0;

        virtual ~ViewingGeometryBase() {};

    };

    /** A singular line of sight that is defined from parameters at the tangent altitude of the measurement.
    *  This class should only be used for limb viewing lines of sight. And is only valid when operating
    *  in Spherical geometry.
    */
    class TangentAltitude : public ViewingGeometryBase {
    private:
        double m_tangentaltitude; /**< The unrefracted tangent altitude of the line of sight in [m]   */
        double m_observeraltitude; /**< The altitude of the observer in [m].  If None then the observer is assumed to be outside the atmosphere */
        double m_relative_azimuth_angle; /**< Relative azimuth angle in [radians] */

        double m_theta; /**< Angle in [radians] along the reference_plane direction */
        double m_phi; /**< Angle in [radians] along the cross reference_plane direction */

    public:
        /**
        * @param tangentaltitude The unrefracted tangent altitude of the line of sight in [m]
        * @param relative_azimuth_angle The relative azimuth angle for the line of sight in [radians]
        * @param observeraltitude The altitude of the observer in [m].
        * @param theta Angle in [radians] in the along reference plane direction
        * @param phi Angle in [radians] in the across reference plane direction
        */
        TangentAltitude(double tangentaltitude, double relative_azimuth_angle, double observeraltitude,
                        double theta=0,
                        double phi=0
                        );

        /** Constructs the ray from the user provided angles and altitudes
         *
         * @param geometry Internal sasktran Coordinates
         * @return ViewingRay
         */
        ViewingRay construct_ray(const sasktran2::Coordinates& geometry) override final;
    };

    /** A singular line of sight that is defined from solar parameters at the tangent altitude of the measurement.
     *  This class should only be used for limb viewing lines of sight. And is only valid when operating
     *  in Spherical geometry with a 1D atmosphere.
     */
    class TangentAltitudeSolar : public ViewingGeometryBase {
    private:
        double m_tangentaltitude; /**< The unrefracted tangent altitude of the line of sight in [m]   */
        double m_observeraltitude; /**< The altitude of the observer in [m].  If None then the observer is assumed to be outside the atmosphere */
        double m_relative_azimuth_angle; /**< Relative azimuth angle in [radians] */
        double m_cos_sza; /**< Cosine of solar zenith angle */

    public:
        /**
        * @param tangentaltitude The unrefracted tangent altitude of the line of sight in [m]
        * @param relative_azimuth_angle The relative azimuth angle for the line of sight in [radians]
        * @param observeraltitude The altitude of the observer in [m].
        * @param cos_sza Cosine of solar zenith angle
        */
        TangentAltitudeSolar(double tangentaltitude, double relative_azimuth_angle, double observeraltitude,
                             double cos_sza);

        /** Constructs the ray from the user provided angles and altitudes
         *
         * @param geometry Internal sasktran Coordinates
         * @return ViewingRay
         */
        ViewingRay construct_ray(const sasktran2::Coordinates& geometry) override final;
    };

    /** A singular line of sight that is defined based upon the solar angles at the ground point.
     */
    class GroundViewingSolar : public ViewingGeometryBase {
    private:
        double m_cos_sza;
        double m_relative_azimuth_angle;
        double m_observer_altitude;
        double m_cos_viewing_zenith;

    public:
        GroundViewingSolar(double cos_sza, double relative_azimuth_angle, double cos_viewing_zenith, double observer_altitude);

        /** Constructs the ray from the user provided angles and altitudes
         *
         * @param geometry Internal sasktran Coordinates
         * @return ViewingRay
         */
        ViewingRay construct_ray(const sasktran2::Coordinates& geometry) override final;
    };


    /** A container that defines all of the viewing rays
     *
     */
    class ViewingGeometryContainer {
    private:
        std::vector<std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryBase>> m_observer_rays;
    public:
        /**
         *
         * @return The viewing rays
         */
        std::vector<std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryBase>>& observer_rays() { return m_observer_rays; }

        /**
         *
         * @return The viewing rays
         */
        const std::vector<std::unique_ptr<sasktran2::viewinggeometry::ViewingGeometryBase>>& observer_rays() const { return m_observer_rays; }

    };

}