#include <sasktran2/geometry.h>


namespace sasktran2 {

    void Geometry1D::assign_interpolation_weights(const Location &loc,
                                                  std::vector<std::pair<int, double>> &index_weights) const {
        std::array<double, 2> weight;
        std::array<int, 2> index;

        int num_contrib;

        m_alt_grid.calculate_interpolation_weights(loc.radius() - coordinates().earth_radius(),
                index,
                weight,
                num_contrib);

        index_weights.resize(num_contrib);

        for(int i = 0; i < num_contrib; ++i) {
            index_weights[i].first = index[i];
            index_weights[i].second = weight[i];
        }
    }
}