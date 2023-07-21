#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_geometrylayerarray.h"

template <int NSTOKES, int CNSTR>
sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>::GeometryLayerArray(
        const PersistentConfiguration<NSTOKES, CNSTR> &config, const OpticalState<NSTOKES> *opticalstate) : m_config(config),
                                               GeometryLayerArrayROP<NSTOKES>(config),
                                               m_opticalstate(opticalstate) {

    // Load in the altitude grid from the user config
    std::vector<double> altitude_grid = m_config.userSpec()->getAltitudeGrid();
    // Reverse since we construct layers from the TOA downwards
    std::reverse(altitude_grid.begin(), altitude_grid.end());

    // Initialize the chapman factor storage
    m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
    m_chapman_factors.setZero();

    // Calculate the ceiling heights and floor heights of each layer
    m_ceiling_h.resize(this->M_NLYR);
    m_floor_h.resize(this->M_NLYR);

    if (config.layer_construction_method() == SKTRAN_DO_UserSpec::LayerConstructionMethod::uniform_pressure) {
        // Evenly spaced layers in pressure space
        double pressure_difference = (opticalstate->pressure()(Eigen::last) - opticalstate->pressure()(0)) / this->M_NLYR;  // negative

        for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
            double layer_p_ceil = opticalstate->pressure()(Eigen::last) - p * pressure_difference;
            double layer_p_floor = std::min(opticalstate->pressure()(Eigen::last) - (p + 1) * pressure_difference, opticalstate->pressure()(0));

            m_ceiling_h[p] = opticalstate->altitude_at_pressure(layer_p_ceil);
            m_floor_h[p] = opticalstate->altitude_at_pressure(layer_p_floor);
        }
    }
    else if (config.layer_construction_method() == SKTRAN_DO_UserSpec::LayerConstructionMethod::uniform_height) {
        // Evenly spaced layers in height
        double altitude_diff = this->m_userspec->getTopAltitude() - this->m_userspec->getBottomAltitude();
        double h_difference = altitude_diff / this->M_NLYR;
        for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
            m_ceiling_h[p] = altitude_diff - p * h_difference;
            m_floor_h[p] = std::max(altitude_diff - (p + 1)*h_difference, 0.0);
        }
    }
    else if (config.layer_construction_method() == SKTRAN_DO_UserSpec::LayerConstructionMethod::match_altitudegrid) {
        if (this->M_NLYR != altitude_grid.size() - 1) {
            throw("Number of layers must equal the length of the altitude grid - 1 when using layerconstructionmethod::match_altitudegrid");
        }
        for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
            m_ceiling_h[p] = altitude_grid[p];
            m_floor_h[p] = altitude_grid[p + 1];
        }
    }
    else if (config.layer_construction_method() == SKTRAN_DO_UserSpec::LayerConstructionMethod::manual) {
        const auto& alts = m_config.userSpec()->getManualLayerAltitudes();

        if (this->M_NLYR != alts.size() - 1) {
            throw("Number of layers must equal the length of the manual altitude grid - 1");
        }
        for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
            m_ceiling_h[p] = alts[alts.size() - p - 1];
            m_floor_h[p] = alts[alts.size() - p - 2];
        }
    }

    // Now we have the layer locations, we can calculate the chapman factors
    calculate_chapman_factors(m_config.coords()->AltitudeToRadius(0.0));

    // Next we can write layer OD = A * optical_state.extinction, and calculate the A matrix using pure geometry factors
    calculate_optical_table_matrices();

}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>::calculate_optical_table_matrices() {
    const auto& alt_grid = m_config.userSpec()->getAltitudeGrid();

    m_optical_interpolator.resize(this->M_NLYR, alt_grid.size());
    m_optical_interpolator.setZero();

    for(LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double ceiling_h = m_ceiling_h[p];
        double floor_h = m_floor_h[p];

        int ceil_idx, floor_idx;

        if (ceiling_h >= alt_grid[alt_grid.size() - 1]) {
            ceil_idx = alt_grid.size() - 1;
        }
        else {
            ceil_idx = std::distance(
                alt_grid.cbegin(),
                std::upper_bound(alt_grid.cbegin(), alt_grid.cend(), ceiling_h)
            );
        }

        if (floor_h <= alt_grid[0]) {
            floor_idx = 0;
        }
        else {
            floor_idx = std::distance(
                alt_grid.cbegin(),
                --std::lower_bound(alt_grid.cbegin(), alt_grid.cend(), floor_h)
            );
        }

        auto num_interior_partitions = ceil_idx - floor_idx;

        for (decltype(num_interior_partitions) i = 0; i < num_interior_partitions; ++i) {
            double low_alt, high_alt;
            if(i == 0) {
                low_alt = m_floor_h(p);
                if( num_interior_partitions == 1 ) {
                    high_alt = m_ceiling_h(p);
                } else {
                    high_alt = alt_grid[floor_idx + i + 1];
                }
            } else {
                low_alt = alt_grid[floor_idx + i];
                if( i == num_interior_partitions - 1) {
                    high_alt = m_ceiling_h[p];
                } else {
                    high_alt = alt_grid[floor_idx + i + 1];
                }
            }

            m_opticalstate->fill_interpolation_matrix(high_alt, low_alt, p, m_optical_interpolator);
        }
    }
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>::calculate_chapman_factors(double earth_rad) {
    if (m_config.use_pseudo_spherical()) {
        // Calculate the chapman factors for the layer
        double sinthetasq = 1 - this->M_CSZ * this->M_CSZ;

        for(sktran_do_detail::LayerIndex p = 0; p < this->M_NLYR; ++p) {
            double rp = earth_rad + m_floor_h(p);

            for (sktran_do_detail::LayerIndex q = 0; q <= p; ++q) {
                double rfloor = earth_rad + m_floor_h(q);
                double rceil = earth_rad + m_ceiling_h(q);

                m_chapman_factors(p, q) = (sqrt(rceil * rceil - rp * rp*sinthetasq) - sqrt(rfloor * rfloor - rp * rp*sinthetasq)) / (rceil - rfloor);
            }
        }
    }
    else {
        m_chapman_factors.setConstant(1 / this->M_CSZ);
    }
}

INSTANTIATE_TEMPLATE(sktran_do_detail::GeometryLayerArray);