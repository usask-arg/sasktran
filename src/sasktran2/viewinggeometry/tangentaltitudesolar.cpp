#include <sasktran2.h>
#include "sasktran2/viewinggeometry.h"


namespace sasktran2::viewinggeometry {
    TangentAltitudeSolar::TangentAltitudeSolar(double tangentaltitude, double relative_azimuth_angle,
                                               double observeraltitude, double cos_sza) :
                                               m_tangentaltitude(tangentaltitude),
                                               m_relative_azimuth_angle(relative_azimuth_angle),
                                               m_observeraltitude(observeraltitude),
                                               m_cos_sza(cos_sza)
                                               {
        #ifdef SASKTRAN_DEBUG_ASSERTS
            if(std::isnan(tangentaltitude) || tangentaltitude < 0) {
                BOOST_LOG_TRIVIAL(error) << "Invalid Tangent Altitude in TangentAltitudeSolar: " << tangentaltitude;
            }

            if(std::isnan(relative_azimuth_angle)) {
                BOOST_LOG_TRIVIAL(error) << "NaN Relative Azimuth Angle in TangentAltitudeSolar";
            }

            if(std::isnan(observeraltitude)) {
                BOOST_LOG_TRIVIAL(error) << "NaN Observer Altitude in TangentAltitudeSolar";
            }

            if(std::isnan(cos_sza) || std::abs(cos_sza) > 1) {
                BOOST_LOG_TRIVIAL(error) << "Invalid cos_sza in TangentAltitudeSolar: " << cos_sza;
            }
        #endif

    }

    ViewingRay TangentAltitudeSolar::construct_ray(const sasktran2::Coordinates &geometry) {
        if(geometry.geometry_type() != sasktran2::geometrytype::spherical) {
            auto msg = "Error constructing ray in TangentAltitude::construct_ray, TangentAltitude ray construction can only be used in spherical geometry mode.";
            BOOST_LOG_TRIVIAL(error) << msg;
            throw std::invalid_argument(msg);
        }
        ViewingRay ray;

        // Get the unit vector pointing to the tangent altitude
        Eigen::Vector3d tangent_point = geometry.solar_coordinate_vector(m_cos_sza, 0.0, m_tangentaltitude);

        // Calculate the local look vector
        ray.look_away = geometry.look_vector_from_azimuth(tangent_point, m_relative_azimuth_angle, 0);

        // Now we need to back calculate the observer position based upon altitude
        double s = sqrt(math::sqr(geometry.earth_radius() + m_observeraltitude) -
                        math::sqr(geometry.earth_radius() + m_tangentaltitude));

        ray.observer.position = tangent_point - s*ray.look_away;

        #ifdef SASKTRAN_DEBUG_ASSERTS
            if(!tangent_point.allFinite()) {
                BOOST_LOG_TRIVIAL(error) << "Error calculating tangent point from cos_sza: " << m_cos_sza << " tangent_alt: " << m_tangentaltitude;
            }

            if(!ray.look_away.allFinite()) {
                BOOST_LOG_TRIVIAL(error) << "Error calculating the look vector from tangent point: " << tangent_point << " rel_az: " << m_relative_azimuth_angle;
            }

            if(!ray.observer.position.allFinite()) {
                BOOST_LOG_TRIVIAL(error) << "Error estimating the observer position from: ";
                BOOST_LOG_TRIVIAL(error) << "tangent_point: " << tangent_point;
                BOOST_LOG_TRIVIAL(error) << "look_away: " << ray.look_away;
                BOOST_LOG_TRIVIAL(error) << "obs_alt: " << m_observeraltitude;
                BOOST_LOG_TRIVIAL(error) << "tangent_alt: " << m_tangentaltitude;
                BOOST_LOG_TRIVIAL(error) << "s: " << s;
            }
        #endif

        return ray;

    }

}