#include <sasktran2.h>
#include "sasktran2/viewinggeometry.h"


namespace sasktran2::viewinggeometry {

    GroundViewingSolar::GroundViewingSolar(double cos_sza, double relative_azimuth_angle, double cos_viewing_zenith, double observer_altitude) :
    m_cos_sza(cos_sza),
    m_relative_azimuth_angle(relative_azimuth_angle),
    m_observer_altitude(observer_altitude),
    m_cos_viewing_zenith(cos_viewing_zenith)
    {

    }

    ViewingRay GroundViewingSolar::construct_ray(const sasktran2::Coordinates &geometry) {
        ViewingRay result;

        // Coordinate of ground point that has the correct angles
        Eigen::Vector3d ground_vector = geometry.solar_coordinate_vector(m_cos_sza, m_relative_azimuth_angle, 0.0);

        result.look_away = -1.0*geometry.look_vector_from_azimuth(ground_vector, m_relative_azimuth_angle, m_cos_viewing_zenith);

        double distance_from_ground = 0.0;

        if(geometry.geometry_type() == sasktran2::geometrytype::planeparallel) {
            // Distance from ground is a simple scaling of the observer altitude
            distance_from_ground = m_observer_altitude / m_cos_viewing_zenith;

        } else if (geometry.geometry_type() == sasktran2::geometrytype::spherical) {
            // Law of cosines + quadratic formula
            double b = 2.0*geometry.earth_radius() * m_cos_viewing_zenith;
            double c = -(2.0*geometry.earth_radius() * m_observer_altitude + m_observer_altitude * m_observer_altitude);

            // want the positive solution
            distance_from_ground = (-b + sqrt(b*b - 4*c)) / 2;

        } else {
            BOOST_LOG_TRIVIAL(error) << "GroundViewingSolar does not support the given geometry type";
        }

        result.observer.position = ground_vector - result.look_away * distance_from_ground;

        return result;
    }


}
