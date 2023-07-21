#include <sasktran2/grids.h>
#include <sasktran2/geometry.h>

namespace sasktran2::grids {
    AltitudeSZASourceLocationInterpolator::AltitudeSZASourceLocationInterpolator(AltitudeGrid &&altitude_grid,
                                                                                 Grid &&sza_grid) :
            SourceLocationInterpolator(std::move(altitude_grid)),
            m_cos_sza_grid(sza_grid)
                                                                                 {

    }

    int AltitudeSZASourceLocationInterpolator::num_ground_points() const {
        return (int)m_cos_sza_grid.grid().size();
    }

    int AltitudeSZASourceLocationInterpolator::num_interior_points() const {
        return (int)(num_ground_points() * m_altitude_grid.grid().size());
    }

    Eigen::Vector3d AltitudeSZASourceLocationInterpolator::grid_location(const sasktran2::Coordinates &coords,
                                                                         int location_index) const {
        int alt_index = location_index % m_altitude_grid.grid().size();
        int csz_index = (int)(location_index / m_altitude_grid.grid().size());

        return coords.solar_coordinate_vector(m_cos_sza_grid.grid()(csz_index), EIGEN_PI, m_altitude_grid.grid()(alt_index));
    }

    int AltitudeSZASourceLocationInterpolator::interior_linear_index(int alt_index, int sza_index) {
        return (int)(alt_index + sza_index * m_altitude_grid.grid().size());
    }

    void AltitudeSZASourceLocationInterpolator::interior_interpolation_weights(const sasktran2::Coordinates& coords,
                                                                               const sasktran2::Location &location,
                                                                               std::vector<std::pair<int, double>> &weights,
                                                                               int &num_interp) {
        std::array<int, 2> alt_index, sza_index;
        std::array<double, 2> alt_weight, sza_weight;
        int num_alt_contrib, num_sza_contrib;

        double altitude = location.radius() - coords.earth_radius();

        double cos_sza = coords.solar_angles_at_location(location.position).first;

        m_cos_sza_grid.calculate_interpolation_weights(cos_sza, sza_index, sza_weight, num_sza_contrib);
        m_altitude_grid.calculate_interpolation_weights(altitude, alt_index, alt_weight, num_alt_contrib);

        num_interp = num_sza_contrib * num_alt_contrib;

        if( weights.size() < num_interp) {
            weights.resize(num_interp);
        }

        for(int i = 0; i < num_alt_contrib; ++i) {
            for(int j = 0; j < num_sza_contrib; ++j) {
                weights[i + j*num_alt_contrib].first = interior_linear_index(alt_index[i], sza_index[j]);
                weights[i + j*num_alt_contrib].second = sza_weight[j] * alt_weight[i];
            }
        }
    }
}