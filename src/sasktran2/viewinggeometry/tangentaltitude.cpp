#include <sasktran2.h>
#include "sasktran2/viewinggeometry.h"


namespace sasktran2::viewinggeometry {
    TangentAltitude::TangentAltitude(double tangentaltitude, double relative_azimuth_angle, double observeraltitude, double theta, double phi) :
    m_observeraltitude(observeraltitude),
    m_tangentaltitude(tangentaltitude),
    m_relative_azimuth_angle(relative_azimuth_angle),
    m_theta(theta),
    m_phi(phi)
    {
    }

    ViewingRay TangentAltitude::construct_ray(const sasktran2::Coordinates &geometry) {
        if(geometry.geometry_type() != sasktran2::geometrytype::spherical) {
            auto msg = "Error constructing ray in TangentAltitude::construct_ray, TangentAltitude ray construction can only be used in spherical geometry mode.";
            BOOST_LOG_TRIVIAL(error) << msg;
            throw std::invalid_argument(msg);
        }
        ViewingRay ray;

        // Get the unit vector pointing to the tangent altitude
        Eigen::Vector3d uv = geometry.unit_vector_from_angles(m_theta, m_phi);
        Eigen::Vector3d tangent_point = uv * (geometry.earth_radius() + m_tangentaltitude);

        // And the local coordinate system at the tangent point
        std::pair<Eigen::Vector3d, Eigen::Vector3d> x_y = geometry.local_x_y_from_angles(m_theta, m_phi);

        // Calculate the local look vector
        ray.look_away = cos(m_relative_azimuth_angle) * x_y.first +
                sin(m_relative_azimuth_angle) * x_y.second;

        // Now we need to back calculate the observer position based upon altitude
        double s = sqrt(math::sqr(geometry.earth_radius() + m_observeraltitude) -
                math::sqr(geometry.earth_radius() + m_tangentaltitude));

        ray.observer.position = tangent_point - s*ray.look_away;

        return ray;
    }
}