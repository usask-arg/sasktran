#include <sasktran2/solartransmission.h>
#include <sasktran2/dual.h>
#include <sasktran2/config.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/math/trig.h>
#include <boost/math/quadrature/gauss.hpp>

namespace sasktran2::solartransmission {
    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere)  {
        // Store the atmosphere for later
        m_atmosphere  = &atmosphere;

        // Calculate the phase interpolators for the LOS rays
        int max_index = 0;
        for(const auto& map : m_phase_index_map) {
            for(const auto& ele : map) {
                if(ele > max_index) {
                    max_index = ele;
                }
            }
        }

        int num_phase = max_index + 1;

        m_phase_interp.resize(num_phase);
        #pragma omp parallel for
        for(int i = 0; i < (*m_los_rays).size(); ++i) {
            const auto& ray = (*m_los_rays)[i];
            for(int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];
                m_phase_interp[m_phase_index_map[i][j]].load_scattering_angle(m_atmosphere->storage().phase[0], m_geometry.coordinates().sun_unit(), layer.average_look_away);
            }
        }

        // Initialize some local memory storage
        for(int i = 0; i < m_start_source_cache.size(); ++i) {
            m_start_source_cache[i].resize(NSTOKES, atmosphere.num_deriv(), false);
            m_end_source_cache[i].resize(NSTOKES, atmosphere.num_deriv(), false);
        }

        // Precompute some parts of the source calculation
        /*
        for(int losidx = 0; losidx < m_los_rays->size(); ++losidx) {
            for(int layeridx = 0; layeridx < (*m_los_rays)[losidx].layers.size(); ++layeridx) {
                // Integrates assuming the source is constant in the layer and determined by the average of the
                // layer boundaries
                int exit_index = m_index_map[losidx][layeridx];
                int entrance_index = m_index_map[losidx][layeridx] + 1;

                auto& phase_interp = m_phase_interp[m_phase_index_map[losidx][layeridx]];

                auto& start_phase = m_start_source_cache[threadidx];
                auto& end_phase = m_end_source_cache[threadidx];

                double entrance_factor = exp(-solar_od_entrance) / (sasktran2::math::Pi * 4);
                double exit_factor = exp(-solar_od_exit) / (sasktran2::math::Pi * 4);

                start_phase.value.setZero();
                end_phase.value.setZero();

                // Calculate SSA, don't do derivatives until later
                for(auto& ele : layer.entrance.interpolation_weights) {
                    ssa_start += m_atmosphere->storage().ssa(ele.first, wavelidx) * ele.second;
                }
                for(auto& ele : layer.exit.interpolation_weights) {
                    ssa_end += m_atmosphere->storage().ssa(ele.first, wavelidx) * ele.second;
                }

                // Incident solar beam is unpolarized
                start_phase.value(0) = ssa_start * entrance_factor;
                end_phase.value(0) = ssa_end * exit_factor;

                phase_interp.scatter(m_atmosphere->storage().phase[wavelidx], layer.entrance.interpolation_weights, start_phase);
                phase_interp.scatter(m_atmosphere->storage().phase[wavelidx], layer.exit.interpolation_weights, end_phase);

                if(calculate_derivative) {
                    if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                        if (m_config->wf_precision() != sasktran2::Config::WeightingFunctionPrecision::limited) {
                            // Have to apply the solar transmission derivative factors
                            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(m_geometry_sparse,
                                                                                                entrance_index); it; ++it) {
                                start_phase.deriv(Eigen::all, it.index()) -= it.value() * start_phase.value;
                            }

                            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(m_geometry_sparse,
                                                                                                exit_index); it; ++it) {
                                end_phase.deriv(Eigen::all, it.index()) -= it.value() * end_phase.value;
                            }
                        }
                    }

                    // And SSA derivative factors
                    for(auto& ele : layer.entrance.interpolation_weights) {
                        start_phase.deriv(Eigen::all, m_atmosphere->ssa_deriv_start_index() + ele.first) += ele.second * start_phase.value / ssa_start;
                    }
                    for(auto& ele : layer.exit.interpolation_weights) {
                        end_phase.deriv(Eigen::all, m_atmosphere->ssa_deriv_start_index() + ele.first) += ele.second * end_phase.value / ssa_end;
                    }
                }


                double source_factor1 = (1 - exp(-shell_od.od));
                // Note dsource_factor = d_od * (1 - source_factor)

                // Get the phase matrix and add on the sources
                // The source factor term will only have extinction derivatives, the phase term will have
                // local SSA/scattering derivatives and is ~dense in a 1D atmosphere

                source.value.array() += source_factor1 * start_phase.value.array() * layer.od_quad_start_fraction;
                source.value.array() += source_factor1 * end_phase.value.array() * layer.od_quad_end_fraction;

                if(calculate_derivative) {
                    // Now for the derivatives, start with dsource_factor which is sparse
                    for(auto it = shell_od.deriv_iter; it; ++it) {
                        source.deriv(Eigen::all, it.index()).array() += it.value() * (1 - source_factor1) * start_phase.value.array() * layer.od_quad_start_fraction;
                        source.deriv(Eigen::all, it.index()).array() += it.value() * (1 - source_factor1) * end_phase.value.array() * layer.od_quad_end_fraction;
                    }
                    // And add on d_phase
                    source.deriv.array() += source_factor1 * start_phase.deriv.array() * layer.od_quad_start_fraction + source_factor1 * end_phase.deriv.array() * layer.od_quad_end_fraction;
                }

            }
        }
         */



    };

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_config(const sasktran2::Config &config) {
        m_config = &config;

        this->m_solar_transmission.initialize_config(config);

        // Set up storage for each thread
        //m_solar_trans.resize(config.num_threads());
        m_thread_index_cache_one.resize(config.num_threads());
        m_thread_index_cache_two.resize(config.num_threads());

        m_solar_trans.resize(config.num_threads());

        m_start_source_cache.resize(config.num_threads());
        m_end_source_cache.resize(config.num_threads());

        m_phase_interp.resize(config.num_threads());
    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::calculate(int wavelidx, int threadidx) {
        // Don't have to do anything here

        // Calculate the solar transmission at each cell
        if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
            // Faster to use the dense matrix if most of the elements are nonzero
            // TODO: Is 0.25 a good number?
            if(double(m_geometry_sparse.nonZeros()) / double(m_geometry_matrix.size()) > 1) {
                m_solar_trans[threadidx] = m_geometry_matrix * m_atmosphere->storage().total_extinction(Eigen::all, wavelidx);
            } else {
                m_solar_trans[threadidx] = m_geometry_sparse * m_atmosphere->storage().total_extinction(Eigen::all, wavelidx);
            }
        }

        if constexpr (std::is_same_v<S, SolarTransmissionTable>) {
            m_solar_trans[threadidx] = m_geometry_sparse * (m_solar_transmission.geometry_matrix() * m_atmosphere->storage().total_extinction(Eigen::all, wavelidx));
        }

        m_solar_trans[threadidx] = exp(-m_solar_trans[threadidx].array());
        for(int i = 0; i < m_ground_hit_flag.size(); ++i) {
            if(m_ground_hit_flag[i]) {
                m_solar_trans[threadidx][i] = 0;
            }
        }

    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        // TODO: BRDF

        if(m_los_rays->at(losidx).ground_is_hit) {
            double albedo = m_atmosphere->surface().albedo()[wavelidx];

            // Single scatter ground source is solar_trans * cos(th) * albedo / pi
            int exit_index = m_index_map[losidx][0];

            double solar_trans = m_solar_trans[threadidx](exit_index);

            double cos_theta = m_los_rays->at(losidx).layers[0].exit.cos_zenith_angle(m_geometry.coordinates().sun_unit());

            // TODO: Linearization
            double source_value = solar_trans * albedo * cos_theta / EIGEN_PI;

            source.value(0) += source_value;
        }
    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {
        this->m_solar_transmission.initialize_geometry(los_rays);

        if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
            // Generates the geometry matrix so that matrix * extinction = solar od at grid points
            this->m_solar_transmission.generate_geometry_matrix(los_rays, m_geometry_matrix, m_ground_hit_flag);

            // Usually faster to calculate the matrix densely and then convert to sparse
            m_geometry_sparse = m_geometry_matrix.sparseView();
        }
        if constexpr (std::is_same_v<S, SolarTransmissionTable>) {
            this->m_solar_transmission.generate_interpolation_matrix(los_rays, m_geometry_sparse, m_ground_hit_flag);
        }


        // We need some mapping between the layers inside each ray to our calculated solar transmission
        m_index_map.resize(los_rays.size());
        m_phase_index_map.resize(los_rays.size());
        m_num_cells = 0;
        int c = 0;
        int cp = 0;
        for(int i = 0; i < los_rays.size(); ++i) {
            m_index_map[i].resize(los_rays[i].layers.size());
            m_phase_index_map[i].resize(los_rays[i].layers.size());

            for(int j = 0; j < m_index_map[i].size(); ++j) {
                m_index_map[i][j] = c;
                ++c;
            }
            // Final exit layer
            ++c;

            if(los_rays[i].is_straight) {
                for(int j = 0; j < m_index_map[i].size(); ++j) {
                    m_phase_index_map[i][j] = cp;
                }
                ++cp;
            } else {
                for(int j = 0; j < m_index_map[i].size(); ++j) {
                    m_phase_index_map[i][j] = cp;
                    ++cp;
                }
            }

            m_num_cells += (int)los_rays[i].layers.size();
        }

        // Store the rays for later
        m_los_rays = &los_rays;
    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source_quadrature(int wavelidx, int losidx, int layeridx,
                                                                       int threadidx,
                                                                       const sasktran2::raytracing::SphericalLayer &layer,
                                                                       const sasktran2::SparseODDualView& shell_od,
                                                                       sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {

        // TODO: Retry quadrature calculation?
    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source_constant(int wavelidx, int losidx, int layeridx,
                                                                           int threadidx,
                                                                           const sasktran2::raytracing::SphericalLayer &layer,
                                                                           const sasktran2::SparseODDualView& shell_od,
                                                                           sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        bool calculate_derivative = source.derivative_size() > 0;

        // Integrates assuming the source is constant in the layer and determined by the average of the
        // layer boundaries
        int exit_index = m_index_map[losidx][layeridx];
        int entrance_index = m_index_map[losidx][layeridx] + 1;

        double solar_trans_exit = m_solar_trans[threadidx](exit_index);
        double solar_trans_entrance = m_solar_trans[threadidx](entrance_index);

        auto& phase_interp = m_phase_interp[m_phase_index_map[losidx][layeridx]];

        auto& start_phase = m_start_source_cache[threadidx];
        auto& end_phase = m_end_source_cache[threadidx];


        if (calculate_derivative) {
            start_phase.deriv.setZero();
            end_phase.deriv.setZero();
        }

        double ssa_start = 0;
        double ssa_end = 0;

        double entrance_factor = solar_trans_entrance / (sasktran2::math::Pi * 4);
        double exit_factor = solar_trans_exit / (sasktran2::math::Pi * 4);

        start_phase.value.setZero();
        end_phase.value.setZero();

        // Calculate SSA, don't do derivatives until later
        for(auto& ele : layer.entrance.interpolation_weights) {
            ssa_start += m_atmosphere->storage().ssa(ele.first, wavelidx) * ele.second;
        }
        for(auto& ele : layer.exit.interpolation_weights) {
            ssa_end += m_atmosphere->storage().ssa(ele.first, wavelidx) * ele.second;
        }

        // Incident solar beam is unpolarized
        start_phase.value(0) = ssa_start * entrance_factor;
        end_phase.value(0) = ssa_end * exit_factor;

        phase_interp.scatter(m_atmosphere->storage().phase[wavelidx], layer.entrance.interpolation_weights, start_phase);
        phase_interp.scatter(m_atmosphere->storage().phase[wavelidx], layer.exit.interpolation_weights, end_phase);

        if(calculate_derivative) {
            if constexpr (std::is_same_v<S, SolarTransmissionExact>) {
                if (m_config->wf_precision() != sasktran2::Config::WeightingFunctionPrecision::limited) {
                    // Have to apply the solar transmission derivative factors
                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(m_geometry_sparse,
                                                                                        entrance_index); it; ++it) {
                        start_phase.deriv(Eigen::all, it.index()) -= it.value() * start_phase.value;
                    }

                    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(m_geometry_sparse,
                                                                                        exit_index); it; ++it) {
                        end_phase.deriv(Eigen::all, it.index()) -= it.value() * end_phase.value;
                    }
                }
            }

            // And SSA derivative factors
            for(auto& ele : layer.entrance.interpolation_weights) {
                start_phase.deriv(Eigen::all, m_atmosphere->ssa_deriv_start_index() + ele.first) += ele.second * start_phase.value / ssa_start;
            }
            for(auto& ele : layer.exit.interpolation_weights) {
                end_phase.deriv(Eigen::all, m_atmosphere->ssa_deriv_start_index() + ele.first) += ele.second * end_phase.value / ssa_end;
            }
        }


        double source_factor1 = (1 - shell_od.exp_minus_od);
        // Note dsource_factor = d_od * (1 - source_factor)

        // Get the phase matrix and add on the sources
        // The source factor term will only have extinction derivatives, the phase term will have
        // local SSA/scattering derivatives and is ~dense in a 1D atmosphere

        source.value.array() += source_factor1 * (start_phase.value.array() * layer.od_quad_start_fraction + end_phase.value.array() * layer.od_quad_end_fraction);

        if(calculate_derivative) {
            // Now for the derivatives, start with dsource_factor which is sparse
            for(auto it = shell_od.deriv_iter; it; ++it) {
                source.deriv(Eigen::all, it.index()).array() += it.value() * (1 - source_factor1) * (start_phase.value.array() * layer.od_quad_start_fraction + end_phase.value.array() * layer.od_quad_end_fraction);
            }
            // And add on d_phase
            source.deriv.array() += source_factor1 * start_phase.deriv.array() * layer.od_quad_start_fraction + source_factor1 * end_phase.deriv.array() * layer.od_quad_end_fraction;
        }

        #ifdef SASKTRAN_DEBUG_ASSERTS
        if(source.value.hasNaN()) {
            static bool message = false;
            if(!message) {
                BOOST_LOG_TRIVIAL(error) << "SS Source Nan" << source_factor1 << " " << start_phase.value << " " << layer.od_quad_start_fraction;
                BOOST_LOG_TRIVIAL(error) << "SS Source Nan" << source_factor1 << " " << end_phase.value << " " << layer.od_quad_end_fraction;

                message = true;
            }
        }
        #endif
    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source_linear(int wavelidx, int losidx, int layeridx,
                                                                   int threadidx,
                                                                   const sasktran2::raytracing::SphericalLayer &layer,
                                                                   const sasktran2::SparseODDualView& shell_od,
                                                                   sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
    // TODO: Go back to this?
    }

    template<typename S, int NSTOKES>
    void SingleScatterSource<S, NSTOKES>::integrated_source(int wavelidx, int losidx, int layeridx, int threadidx,
                                                   const sasktran2::raytracing::SphericalLayer &layer,
                                                   const sasktran2::SparseODDualView&  shell_od,
                                                   sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        if(layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do anything
            return;
        }

        integrated_source_constant(wavelidx, losidx, layeridx, threadidx,
                                         layer, shell_od, source);

    }


    template class SingleScatterSource<SolarTransmissionExact, 1>;
    template class SingleScatterSource<SolarTransmissionExact, 3>;

    template class SingleScatterSource<SolarTransmissionTable, 1>;
    template class SingleScatterSource<SolarTransmissionTable, 3>;
}