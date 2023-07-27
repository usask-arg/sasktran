#include <sasktran2/grids.h>
#include <sasktran2.h>

namespace sasktran2::grids {
    Grid::Grid(Eigen::VectorXd&& grid_values,
    gridspacing spacing,
    outofbounds out_of_bounds_mode,
    interpolation interp
    ) :
    m_grid_spacing(spacing),
    m_out_of_bounds_mode(out_of_bounds_mode),
    m_interp_method(interp),
    m_grid_values(grid_values)
    {
        if(m_interp_method != interpolation::linear) {
            BOOST_LOG_TRIVIAL(warning) << "Requested interpolation mode is not implemented, falling back to linear interpolation";
        }

        if(m_grid_spacing == gridspacing::constant) {
            if(grid_values.size() > 1 ) {
                m_x0 = m_grid_values(0);
                m_dx = (m_grid_values(1) - m_grid_values(0));
            } else {
                // Unused values
                m_x0 = std::numeric_limits<double>::quiet_NaN();
                m_dx = std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            // Unused values
            m_x0 = std::numeric_limits<double>::quiet_NaN();
            m_dx = std::numeric_limits<double>::quiet_NaN();
        }
    }


    void Grid::interpolate_constant_spacing(double x, std::array<int, 2> &index, std::array<double, 2> &weight,
                                            int &num_contributing) const {
        if(x < m_x0) {
            // out of bounds on the lower side
            if(m_out_of_bounds_mode == outofbounds::setzero) {
                num_contributing = 0;
                // For safety set the index and weights to 0
                index[0] = 0;
                index[1] = 0;
                weight[0] = 0;
                weight[1] = 0;
                return;
            } else {
                // Set to the first value
                index[0] = 0;
                index[1] = 0;
                weight[0] = 1;
                weight[1] = 0;
                num_contributing = 1;
                return;
            }
        } else {
            // Greater than the first value
            int i = int(floor((x - m_x0) / m_dx));

            if( i >= m_grid_values.size() - 1) {
                // out of bounds on the right side
                if(m_out_of_bounds_mode == outofbounds::setzero) {
                    num_contributing = 0;
                    // For safety set the index and weights to 0
                    index[0] = 0;
                    index[1] = 0;
                    weight[0] = 0;
                    weight[1] = 0;
                    return;
                } else {
                    // Set to the last value
                    index[0] = (int)m_grid_values.size() - 1;
                    index[1] = 0;
                    weight[0] = 1;
                    weight[1] = 0;
                    num_contributing = 1;
                    return;
                }
            } else {
                // Perform linear interpolation
                index[0] = i;
                index[1] = i+1;

                weight[1] = (x - m_grid_values(i)) / m_dx;
                weight[0] = 1 - weight[1];

                num_contributing = 2;
                return;
            }
        }
    }

    void Grid::interpolate_varying_spacing(double x, std::array<int, 2> &index, std::array<double, 2> &weight,
                                           int &num_contributing) const {
        // Start by checking the out of bounds parameters
        if(x < m_grid_values(0)) {
            // out of bounds on the lower side
            if (m_out_of_bounds_mode == outofbounds::setzero) {
                num_contributing = 0;
                // For safety set the index and weights to 0
                index[0] = 0;
                index[1] = 0;
                weight[0] = 0;
                weight[1] = 0;
                return;
            } else {
                // Set to the first value
                index[0] = 0;
                index[1] = 0;
                weight[0] = 1;
                weight[1] = 0;
                num_contributing = 1;
                return;
            }
        }

        if(x > m_grid_values(Eigen::last)) {
            // out of bounds on the lower side
            if (m_out_of_bounds_mode == outofbounds::setzero) {
                num_contributing = 0;
                // For safety set the index and weights to 0
                index[0] = 0;
                index[1] = 0;
                weight[0] = 0;
                weight[1] = 0;
                return;
            } else {
                // Set to the last
                index[0] = (int)m_grid_values.size() - 1;
                index[1] = 0;
                weight[0] = 1;
                weight[1] = 0;
                num_contributing = 1;
                return;
            }
        }

        // Do a binary search to find the index
        int i = (int)(std::lower_bound(m_grid_values.begin(), m_grid_values.end(), x) - m_grid_values.begin());

        if( i == 0) {
            // Equal to lowest layer
            i += 1;
        }

        // This points to the index with the first one bigger than x, so our two indexes are i and i-1
        index[0] = i-1;
        index[1] = i;
        // Perform linear interpolation

        weight[1] = (x - m_grid_values(i-1)) / (m_grid_values(i) - m_grid_values(i-1));
        weight[0] = 1 - weight[1];

        num_contributing = 2;
    }

    void Grid::calculate_interpolation_weights(double x, std::array<int, 2> &index, std::array<double, 2> &weight,
                                               int &num_contributing) const {
        if(m_grid_values.size() == 1) {
            // Special case for constant grid
            index[0] = 0;
            index[1] = 0;

            weight[0] = 1;
            weight[1] = 0;

            num_contributing = 1;
            return;
        }

        if (m_grid_spacing == gridspacing::constant) {
            return interpolate_constant_spacing(x, index, weight, num_contributing);
        } else {
            return interpolate_varying_spacing(x, index, weight, num_contributing);
        }
    }
}