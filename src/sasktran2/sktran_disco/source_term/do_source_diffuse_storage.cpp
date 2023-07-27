#include "sasktran2/do_source.h"

namespace sasktran2 {
    template<int NSTOKES, int CNSTR>
    DOSourceDiffuseStorage<NSTOKES, CNSTR>::DOSourceDiffuseStorage(
            const sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR> &layer_geometry,
            const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR> &do_config,
            const sasktran2::grids::Grid& sza_grid,
            const Config &config,
            const sasktran2::Geometry1D &geometry) : m_config(config), m_geometry(geometry), m_sza_grid(sza_grid) {
        // Create an altitude grid using the layer geometry
        Eigen::VectorXd altitude_mid = (layer_geometry.layer_ceiling() + layer_geometry.layer_floor()) / 2.0;

        // Layers are defined from TOA downards, so we have to reverse it
        altitude_mid.reverseInPlace();

        m_altitude_grid = std::make_unique<sasktran2::grids::AltitudeGrid>(std::move(altitude_mid),
                                                                           sasktran2::grids::gridspacing::variable,
                                                                           sasktran2::grids::outofbounds::extend,
                                                                           sasktran2::grids::interpolation::linear
        );

        Eigen::VectorXd cos_angles = Eigen::ArrayXd::LinSpaced(40, -1, 1);

        m_cos_angle_grid = std::make_unique<sasktran2::grids::Grid>(std::move(cos_angles),
                                                                    sasktran2::grids::gridspacing::variable,
                                                                    sasktran2::grids::outofbounds::extend,
                                                                    sasktran2::grids::interpolation::linear
        );

        m_num_azi = do_config.nstr();

        m_storage.resize(config.num_threads());

        for(auto& storage: m_storage) {
            // TODO: Derivative
            for(int s = 0; s < NSTOKES; ++s) {
                storage.source_terms[s].resize((int)m_altitude_grid->grid().size() * (int)m_cos_angle_grid->grid().size() * m_num_azi * (int)m_sza_grid.grid().size(), 0, true);
                storage.singlescatter_source_terms[s].resize((int)m_altitude_grid->grid().size() * (int)m_cos_angle_grid->grid().size() * m_num_azi * (int)m_sza_grid.grid().size(), 0, true);
            }
            storage.source_terms_linear.resize((int)m_altitude_grid->grid().size() * (int)m_cos_angle_grid->grid().size() * m_num_azi * (int)m_sza_grid.grid().size()*NSTOKES, 0, true);

            storage.phase_storage.resize(m_cos_angle_grid->grid().size(), do_config.nstr());
            for(int i = 0; i < storage.phase_storage.size(); ++i) {
                storage.phase_storage[i].fill(m_cos_angle_grid->grid()(i));
            }
            storage.phase_container.resize(config.num_do_streams());
        }
    }

    template <int NSTOKES, int CNSTR>
    std::unique_ptr<sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>> DOSourceDiffuseStorage<NSTOKES, CNSTR>::geometry_interpolator(
            const std::vector<sasktran2::raytracing::TracedRay> &rays,
            bool include_azimuth_weights
    ) {
        auto result = std::make_unique<sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>>(rays.size());

        std::array<int, 2> alt_index, angle_index, sza_index;
        std::array<double, 2> alt_weight, angle_weight, sza_weight;
        int num_alt_contrib, num_angle_contrib, num_sza_contrib;

        double earth_radius = m_geometry.coordinates().earth_radius();

        for(int i = 0; i < rays.size(); ++i) {
            auto& ray = rays[i];

            (*result)[i].resize(ray.layers.size());

            for(int j = 0; j < ray.layers.size(); ++j) {
                auto& layer = ray.layers[j];
                std::array<Eigen::SparseVector<double>, NSTOKES>& sparsevec = (*result)[i][j];

                for(int s = 0; s < NSTOKES; ++s) {
                    if(include_azimuth_weights) {
                        sparsevec[s].resize(m_storage[0].source_terms[s].value.size());
                    } else {
                        sparsevec[s].resize(m_storage[0].source_terms[s].value.size() / m_num_azi);
                    }
                }

                double altitude = (layer.entrance.radius() + layer.exit.radius()) / 2.0 - earth_radius;
                double cos_angle = (layer.entrance.cos_zenith_angle(layer.average_look_away) + layer.exit.cos_zenith_angle(layer.average_look_away)) / 2.0;
                double azi = (layer.saz_entrance + layer.saz_exit) / 2.0;
                double cos_sza = (layer.cos_sza_entrance + layer.cos_sza_exit) / 2.0;

                m_altitude_grid->calculate_interpolation_weights(altitude, alt_index, alt_weight, num_alt_contrib);
                m_cos_angle_grid->calculate_interpolation_weights(cos_angle, angle_index, angle_weight, num_angle_contrib);
                m_sza_grid.calculate_interpolation_weights(cos_sza, sza_index, sza_weight, num_sza_contrib);

                for(int szaidx = 0; szaidx < num_sza_contrib; ++szaidx) {
                    for (int altidx = 0; altidx < num_alt_contrib; ++altidx) {
                        for (int angleidx = 0; angleidx < num_angle_contrib; ++angleidx) {
                            double weight = alt_weight[altidx] * angle_weight[angleidx] * sza_weight[szaidx];

                            if(include_azimuth_weights) {
                                for (int k = 0; k < m_num_azi; ++k) {
                                    double azi_factor = cos(k * azi);
                                    int index = linear_storage_index(angle_index[angleidx], alt_index[altidx], sza_index[szaidx], k);

                                    if constexpr (NSTOKES == 1) {
                                        sparsevec[0].coeffRef(index) = azi_factor * weight;
                                    } else if constexpr (NSTOKES == 3) {
                                        double sin_azi_factor = sin(k * (EIGEN_PI - azi));

                                        sparsevec[0].coeffRef(index) = azi_factor * weight;
                                        sparsevec[1].coeffRef(index) = azi_factor * weight;
                                        sparsevec[2].coeffRef(index) = sin_azi_factor * weight;
                                    }
                                }
                            } else {
                                int index = linear_storage_index(angle_index[angleidx], alt_index[altidx], sza_index[szaidx], 0);

                                if constexpr (NSTOKES == 1) {
                                    sparsevec[0].coeffRef(index) = weight;
                                } else if constexpr (NSTOKES == 3) {

                                    sparsevec[0].coeffRef(index) = weight;
                                    sparsevec[1].coeffRef(index) = weight;
                                    sparsevec[2].coeffRef(index) = weight;
                                }
                            }

                        }
                    }
                }
            }
        }
        return result;
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceDiffuseStorage<NSTOKES, CNSTR>::create_location_source_interpolator(
            const std::vector<Eigen::Vector3d> &locations, const std::vector<Eigen::Vector3d> &directions,
            Eigen::SparseMatrix<double, Eigen::RowMajor> &interpolator) const {

        interpolator.resize(locations.size()*NSTOKES, m_storage[0].source_terms_linear.value_size());

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        std::array<int, 2> alt_index, angle_index, sza_index;
        std::array<double, 2> alt_weight, angle_weight, sza_weight;
        int num_alt_contrib, num_angle_contrib, num_sza_contrib;
        double csz, saa;

        sasktran2::Location temp;

        double earth_radius = m_geometry.coordinates().earth_radius();

        for(int i = 0; i < locations.size(); ++i) {
            const auto& location = locations[i];
            const auto& direction = directions[i];

            temp.position = location;

            sasktran2::raytracing::calculate_csz_saz(m_geometry.coordinates().sun_unit(), temp, direction, csz, saa);

            double cos_angle = temp.cos_zenith_angle(-1*direction);
            double altitude = temp.radius() - earth_radius;

            m_altitude_grid->calculate_interpolation_weights(altitude, alt_index, alt_weight, num_alt_contrib);
            m_cos_angle_grid->calculate_interpolation_weights(cos_angle, angle_index, angle_weight, num_angle_contrib);
            m_sza_grid.calculate_interpolation_weights(csz, sza_index, sza_weight, num_sza_contrib);

            for(int szaidx = 0; szaidx < num_sza_contrib; ++szaidx) {
                for (int altidx = 0; altidx < num_alt_contrib; ++altidx) {
                    for (int angleidx = 0; angleidx < num_angle_contrib; ++angleidx) {
                        double weight = alt_weight[altidx] * angle_weight[angleidx] * sza_weight[szaidx];

                        for (int k = 0; k < m_num_azi; ++k) {
                            double azi_factor = cos(k * saa);
                            int index = linear_storage_index(angle_index[angleidx], alt_index[altidx], sza_index[szaidx], k);

                            if constexpr (NSTOKES == 1) {
                                tripletList.emplace_back(T(i, index, azi_factor * weight));
                            } else if constexpr (NSTOKES == 3) {
                                double sin_azi_factor = sin(k * (EIGEN_PI - saa));

                                tripletList.emplace_back(T(i*NSTOKES, index*NSTOKES, azi_factor * weight));
                                tripletList.emplace_back(T(i*NSTOKES + 1, index*NSTOKES + 1, azi_factor * weight));
                                tripletList.emplace_back(T(i*NSTOKES + 2, index*NSTOKES + 2, sin_azi_factor * weight));
                            }
                        }
                    }
                }
            }


        }
        interpolator.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceDiffuseStorage<NSTOKES, CNSTR>::generate_scattering_matrices(const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>& do_config) {
        int num_interp_angles = (int)m_cos_angle_grid->grid().size();
        int num_stream_angles = m_config.num_do_streams() / 2;

        m_scattering_matrix_interpolation_angles.resize(m_num_azi, Eigen::MatrixXd(NSTOKES * num_interp_angles, NSTOKES));
        m_scattering_matrix_stream_angles.resize(m_num_azi, Eigen::MatrixXd(NSTOKES * num_stream_angles, NSTOKES));

        // Fill the scattering matrix
        for(int m = 0; m < m_num_azi; ++m) {
            for(int i = 0; i < num_interp_angles; ++i) {
                if constexpr (NSTOKES == 1) {
                    m_scattering_matrix_interpolation_angles[m](i, 0) = m_storage[0].phase_storage[i].storage(0, m);
                }

                if constexpr (NSTOKES == 3) {
                    int start = i*NSTOKES;
                    // Fill an NSTOKES, NSTOKES block of the matrix
                    m_scattering_matrix_interpolation_angles[m](start, 0) = m_storage[0].phase_storage[i].storage(0, m);
                    m_scattering_matrix_interpolation_angles[m](start+1, 1) = m_storage[0].phase_storage[i].storage(1, m);
                    m_scattering_matrix_interpolation_angles[m](start+2, 2) = m_storage[0].phase_storage[i].storage(1, m);

                    m_scattering_matrix_interpolation_angles[m](start+1, 2) = -m_storage[0].phase_storage[i].storage(2, m);
                    m_scattering_matrix_interpolation_angles[m](start+2, 1) = -m_storage[0].phase_storage[i].storage(2, m);
                }
            }
        }

        LegendrePhaseStorage<NSTOKES, CNSTR> temp_storage(m_config.num_do_streams());
        for(int i = 0; i < num_stream_angles; ++i) {
            temp_storage.fill(do_config.quadrature_cos_angle()->at(i));

            for(int m = 0; m < m_num_azi; ++m) {
                if constexpr (NSTOKES == 1) {
                    m_scattering_matrix_stream_angles[m](i, 0) = temp_storage.storage(0, m);
                }

                if constexpr (NSTOKES == 3) {
                    int start = i*NSTOKES;
                    // Fill an NSTOKES, NSTOKES block of the matrix
                    m_scattering_matrix_stream_angles[m](start, 0) = temp_storage.storage(0, m);
                    m_scattering_matrix_stream_angles[m](start+1, 1) = temp_storage.storage(1, m);
                    m_scattering_matrix_stream_angles[m](start+2, 2) = temp_storage.storage(1, m);

                    m_scattering_matrix_stream_angles[m](start+1, 2) = -temp_storage.storage(2, m);
                    m_scattering_matrix_stream_angles[m](start+2, 1) = -temp_storage.storage(2, m);
                }
            }
        }
    }

    template <int NSTOKES, int CNSTR>
    int DOSourceDiffuseStorage<NSTOKES, CNSTR>::linear_storage_index(int angleidx, int layeridx, int szaidx,
                                                                     int aziidx) const {
        return angleidx + (int)m_cos_angle_grid->grid().size() * layeridx + (int)m_cos_angle_grid->grid().size() * (int)m_altitude_grid->grid().size() * szaidx + (int)m_cos_angle_grid->grid().size() * (int)m_altitude_grid->grid().size() * (int)m_sza_grid.grid().size() * aziidx;
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceDiffuseStorage<NSTOKES, CNSTR>::accumulate_sources(
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> &optical_layer,
            sasktran_disco::AEOrder m,
            sasktran2::DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
            int szaidx,
            int thread_idx
    ) {

        auto& storage = m_storage[thread_idx];

        for(int lidx = 0; lidx < m_altitude_grid->grid().size(); ++lidx) {
            auto* layer = &optical_layer.layer((int)m_altitude_grid->grid().size() - lidx - 1);
            sasktran_disco::LayerIndex p = layer->index();
            double altitude = m_altitude_grid->grid()(lidx);

            double layer_fraction = (layer->altitude(sasktran_disco::Location::CEILING) - altitude) / (layer->altitude(sasktran_disco::Location::CEILING) - layer->altitude(sasktran_disco::Location::FLOOR));
            double x = layer_fraction * layer->dual_thickness().value;
            const auto& transmission = layer->dual_beamTransmittance(sasktran_disco::Location::CEILING, optical_layer.inputDerivatives());

            const auto& average_secant = layer->dual_average_secant();
            const auto& ssa = layer->ssa();
            const sasktran_disco::LayerDual<double>& dual_thickness = layer->dual_thickness();


            auto& cache = thread_storage.postprocessing_cache[p];
            sasktran_disco::uint numderiv = (sasktran_disco::uint)optical_layer.inputDerivatives().numDerivativeLayer(p);
            sasktran_disco::uint numtotalderiv = (sasktran_disco::uint)optical_layer.inputDerivatives().numDerivative();
            int layerStart = (int)optical_layer.inputDerivatives().layerStartIndex(p);

            auto& solution = layer->solution(m);

            cache.resize(m_config.num_do_streams(), p, numderiv, layerStart, numtotalderiv);

            const auto& dual_Aplus = solution.value.dual_green_A_plus();
            const auto& dual_Aminus = solution.value.dual_green_A_minus();

            auto& hp = cache.hp;
            auto& hm = cache.hm;

            auto& Y_plus = cache.Y_plus;
            auto& Y_minus = cache.Y_minus;

            auto& Dm = cache.Dm;
            auto& Dp = cache.Dp;

            auto& Q = cache.Q;
            auto& Qtemp = cache.Qtemp;
            auto& temp = cache.temp;

            const auto& dual_L = solution.boundary.L_coeffs;
            const auto& dual_M = solution.boundary.M_coeffs;

            using MatrixView = Eigen::Map<Eigen::MatrixXd>;
            using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
            MatrixView Y_plus_matrix(Y_plus.value.data(), NSTOKES, m_config.num_do_streams()/2 *NSTOKES);
            MatrixView Y_minus_matrix(Y_minus.value.data(), NSTOKES, m_config.num_do_streams()/2 * NSTOKES);

            ConstMatrixView homog_plus_matrix(solution.value.dual_homog_plus().value.data(), m_config.num_do_streams()/2 * NSTOKES, m_config.num_do_streams()/2 *NSTOKES);
            ConstMatrixView homog_minus_matrix(solution.value.dual_homog_minus().value.data(), m_config.num_do_streams()/2 * NSTOKES, m_config.num_do_streams()/2 *NSTOKES);

            for(int aidx = 0; aidx < m_cos_angle_grid->grid().size(); ++aidx) {
                int sourceidx = linear_storage_index(aidx, lidx, szaidx, m);

                for(int s = 0; s < NSTOKES; ++s) {
                    storage.source_terms[s].value(sourceidx) = 0.0;
                    storage.source_terms_linear.value(sourceidx*NSTOKES + s) = 0.0;
                }

                storage.phase_storage[aidx].set_phase_container(storage.phase_container, m);

                // Calculate legendre sums multiplied by stream weights
                sasktran_disco::VectorLayerDual<double>& dual_lpsum_plus = cache.dual_lpsum_plus;
                sasktran_disco::VectorLayerDual<double>& dual_lpsum_minus = cache.dual_lpsum_minus;

                layer->vectordual_scatPhaseF(m, storage.phase_container, optical_layer.inputDerivatives(), dual_lpsum_minus, dual_lpsum_plus);

                layer->singleScatST(m, storage.phase_container, Qtemp, temp);


                if constexpr (NSTOKES == 1 ) {
                    storage.singlescatter_source_terms[0].value(sourceidx) = Qtemp.value;
                } else {
                    for(int s = 0; s < NSTOKES; ++s) {
                        storage.singlescatter_source_terms[s].value(sourceidx) = Qtemp.value(s);
                    }
                }

                // dual_lpsum_plus and minus are now linearly storage of N * NSTOKES * NSTOKES, have to multiply the blocks by
                // the weights
                for(int i = 0; i < this->m_config.num_do_streams()/2; ++i) {
                    for(int j = 0; j < NSTOKES*NSTOKES; ++j) {
                        int linearindex = i*NSTOKES*NSTOKES + j;
                        dual_lpsum_minus.value(linearindex) *= thread_storage.sza_calculators[0].persistent_config->quadrature_weights()->at(i);
                        dual_lpsum_plus.value(linearindex) *= thread_storage.sza_calculators[0].persistent_config->quadrature_weights()->at(i);
                    }
                }

                ConstMatrixView lpsum_plus_matrix(dual_lpsum_plus.value.data(), NSTOKES, m_config.num_do_streams()/2 * NSTOKES);
                ConstMatrixView lpsum_minus_matrix(dual_lpsum_minus.value.data(), NSTOKES, m_config.num_do_streams()/2 * NSTOKES);

                Y_plus_matrix.noalias() = lpsum_plus_matrix * homog_plus_matrix + lpsum_minus_matrix * homog_minus_matrix;
                Y_minus_matrix.noalias() = lpsum_plus_matrix * homog_minus_matrix + lpsum_minus_matrix * homog_plus_matrix;

                for (int i = 0; i < m_config.num_do_streams() / 2 * NSTOKES; ++i) {
                    const auto& eigval = solution.value.dual_eigval();

                    hp.value = std::exp(-1.0*eigval.value(i) * layer->dual_thickness().value * layer_fraction);
                    hm.value = std::exp(-eigval.value(i) * (layer->dual_thickness().value - layer->dual_thickness().value * layer_fraction));

                    for(int s = 0; s < NSTOKES; ++s) {
                        storage.source_terms[s].value(sourceidx) += Y_plus_matrix(s, i) * hp.value * dual_L.value(i);
                        storage.source_terms[s].value(sourceidx) += Y_minus_matrix(s, i) * hm.value * dual_M.value(i);

                        storage.source_terms_linear.value(sourceidx*NSTOKES + s) += Y_plus_matrix(s, i) * hp.value * dual_L.value(i) / ssa;
                        storage.source_terms_linear.value(sourceidx*NSTOKES + s) += Y_minus_matrix(s, i) * hm.value * dual_M.value(i) / ssa;

                    }

                    Dp.value = transmission.value / (eigval.value(i) + average_secant.value) * (exp(-x*average_secant.value) - exp(-dual_thickness.value * average_secant.value) * exp(-(dual_thickness.value - x) * eigval.value(i)));
                    Dm.value = transmission.value / (average_secant.value - eigval.value(i)) * (exp(-x*eigval.value(i)) - exp(-x*average_secant.value));

                    for(int s = 0; s < NSTOKES; ++s) {
                        storage.source_terms[s].value(sourceidx) += dual_Aplus.value(i) * Y_plus_matrix(s, i) * Dm.value + dual_Aminus.value(i) * Y_minus_matrix(s, i) * Dp.value;

                        storage.source_terms_linear.value(sourceidx*NSTOKES + s) += (dual_Aplus.value(i) * Y_plus_matrix(s, i) * Dm.value + dual_Aminus.value(i) * Y_minus_matrix(s, i) * Dp.value) / ssa;
                    }
                }
            }
        }
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceDiffuseStorage);
}