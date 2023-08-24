#include <sasktran2/math/unitsphere.h>

namespace sasktran2::math {
    UnitSphereGround::UnitSphereGround(std::unique_ptr<const UnitSphere>&& sphere, const Eigen::Vector3d location) :
    m_full_sphere(std::move(sphere)), m_location(location)
    {
        Eigen::Vector3d quad;

        m_contributing_map.reserve(m_full_sphere->num_points()/2);
        m_is_full_sphere_looking_up.resize(m_full_sphere->num_points());

        // We need to find the points that are looking up
        for(int i = 0; i < m_full_sphere->num_points(); ++i) {
            quad = m_full_sphere->get_quad_position(i);

            if(quad.dot(location) > 0) {
                // Looking up

                m_contributing_map.push_back(i);
                m_is_full_sphere_looking_up[i] = true;
            } else {
                m_is_full_sphere_looking_up[i] = false;
            }
        }
    }

    int UnitSphereGround::num_points() const {
        return m_contributing_map.size();
    }

    Eigen::Vector3d UnitSphereGround::get_quad_position(int index) const {
        return m_full_sphere->get_quad_position(m_contributing_map[index]);
    }

    double UnitSphereGround::quadrature_weight(int i) const {
        // TODO: Should we renormalize the quadrature weights to guarantee they integrate to the half sphere?
        return m_full_sphere->quadrature_weight(m_contributing_map[i]);
    }

    void UnitSphereGround::interpolate(const Eigen::Vector3d& direction,
                                       std::vector<std::pair<int, double>>& index_weights,
                                       int& num_interp
    ) const {
        // Start by interpolating over the full sphere
        m_full_sphere->interpolate(direction, index_weights, num_interp);

        // Now go through every interpolation index and see if they are outside the boundaries
        double total_weight;
        for(int i = 0; i < index_weights.size(); ++i) {
            if(!m_is_full_sphere_looking_up[i]) {
                // Have to remove this index from the interpolator
                // Rather than removing it we just set the index to 0 and weight to 0
                index_weights[i].first = 0;
                index_weights[i].second = 0;
            }
            total_weight += index_weights[i].second;
        }

        // Renormalize the weights
        for(int i = 0; i < index_weights.size(); ++i) {
            index_weights[i].second /= total_weight;
        }
    }

}