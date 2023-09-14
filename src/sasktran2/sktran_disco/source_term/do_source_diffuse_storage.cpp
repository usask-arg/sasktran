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

        m_ground_start = (int)m_altitude_grid->grid().size() * (int)m_cos_angle_grid->grid().size() * m_num_azi * (int)m_sza_grid.grid().size();

        int num_ground_points = (int)m_cos_angle_grid->grid().size() * m_num_azi * (int)m_sza_grid.grid().size()*NSTOKES;

        int num_source_points = m_ground_start*NSTOKES + num_ground_points;

        for(auto& storage: m_storage) {
            storage.source_terms_linear.resize(num_source_points, 0, true);

            storage.phase_storage.resize(m_cos_angle_grid->grid().size(), do_config.nstr());
            for(int i = 0; i < storage.phase_storage.size(); ++i) {
                storage.phase_storage[i].fill(m_cos_angle_grid->grid()(i));
            }
            storage.phase_container.resize(config.num_do_streams());
        }

        m_need_to_calculate_map.resize(num_source_points);
        m_converged_map.resize(num_source_points);
        m_need_to_calculate_map.setConstant(false);
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceDiffuseStorage<NSTOKES, CNSTR>::initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmo) {
        // TODO: This seems like the wrong place to do this, but I'm not sure where else it should go

        int numderiv = atmo.num_deriv();

        for(auto& storage : m_storage) {
            storage.source_terms_linear.deriv.resize(storage.source_terms_linear.value_size(), numderiv);
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
                        sparsevec[s].resize(m_storage[0].source_terms_linear.value.size());
                    } else {
                        sparsevec[s].resize(m_storage[0].source_terms_linear.value.size() / m_num_azi);
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
                                    double azi_factor = cos(k * (EIGEN_PI - azi));
                                    int index = linear_storage_index(angle_index[angleidx], alt_index[altidx], sza_index[szaidx], k);

                                    m_need_to_calculate_map[index] = true;

                                    if constexpr (NSTOKES == 1) {
                                        sparsevec[0].coeffRef(index) = azi_factor * weight;
                                    } else if constexpr (NSTOKES == 3) {
                                        double sin_azi_factor = sin(k * (EIGEN_PI - azi));

                                        sparsevec[0].coeffRef(index*NSTOKES) = azi_factor * weight;
                                        sparsevec[1].coeffRef(index*NSTOKES + 1) = azi_factor * weight;
                                        sparsevec[2].coeffRef(index*NSTOKES + 2) = sin_azi_factor * weight;
                                    }
                                }
                            } else {
                                int index = linear_storage_index(angle_index[angleidx], alt_index[altidx], sza_index[szaidx], 0);

                                m_need_to_calculate_map[index] = true;

                                if constexpr (NSTOKES == 1) {
                                    sparsevec[0].coeffRef(index) = weight;
                                } else if constexpr (NSTOKES == 3) {
                                    sparsevec[0].coeffRef(index*NSTOKES) = weight;
                                    sparsevec[1].coeffRef(index*NSTOKES + 1) = weight;
                                    sparsevec[2].coeffRef(index*NSTOKES + 2) = weight;
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
            const std::vector<bool>& ground_hit_flag,
            Eigen::SparseMatrix<double, Eigen::RowMajor> &interpolator) {

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


            m_sza_grid.calculate_interpolation_weights(csz, sza_index, sza_weight, num_sza_contrib);
            if(!ground_hit_flag[i]) {
                // Interior point, interpolate in angle and altitude
                m_altitude_grid->calculate_interpolation_weights(altitude, alt_index, alt_weight, num_alt_contrib);
                m_cos_angle_grid->calculate_interpolation_weights(cos_angle, angle_index, angle_weight, num_angle_contrib);
                for(int szaidx = 0; szaidx < num_sza_contrib; ++szaidx) {
                    for (int altidx = 0; altidx < num_alt_contrib; ++altidx) {
                        for (int angleidx = 0; angleidx < num_angle_contrib; ++angleidx) {
                            double weight = alt_weight[altidx] * angle_weight[angleidx] * sza_weight[szaidx];

                            for (int k = 0; k < m_num_azi; ++k) {
                                double azi_factor = cos(k * (EIGEN_PI - saa));
                                int index = linear_storage_index(angle_index[angleidx], alt_index[altidx], sza_index[szaidx], k);

                                m_need_to_calculate_map[index] = true;

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
            } else {
                // Ground point, just have to interpolate in SZA and only use m=0
                // TODO: CHange when move away from non-lambertian
                for(int szaidx = 0; szaidx < num_sza_contrib; ++szaidx) {
                    double weight = sza_weight[szaidx];
                    int index = ground_storage_index(0, sza_index[szaidx], 0);

                    if constexpr (NSTOKES == 1) {
                        tripletList.emplace_back(T(i, index, weight));
                    } else if constexpr (NSTOKES == 3) {
                        tripletList.emplace_back(T(i*NSTOKES, index*NSTOKES, weight));
                        tripletList.emplace_back(T(i*NSTOKES + 1, index*NSTOKES + 1, weight));
                        tripletList.emplace_back(T(i*NSTOKES + 2, index*NSTOKES + 2, weight));
                    }
                }
            }
        }
        interpolator.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    template <int NSTOKES, int CNSTR>
    int DOSourceDiffuseStorage<NSTOKES, CNSTR>::linear_storage_index(int angleidx, int layeridx, int szaidx,
                                                                     int aziidx) const {
        return angleidx + (int)m_cos_angle_grid->grid().size() * layeridx + (int)m_cos_angle_grid->grid().size() * (int)m_altitude_grid->grid().size() * szaidx + (int)m_cos_angle_grid->grid().size() * (int)m_altitude_grid->grid().size() * (int)m_sza_grid.grid().size() * aziidx;
    }

    template<int NSTOKES, int CNSTR>
    int DOSourceDiffuseStorage<NSTOKES, CNSTR>::ground_storage_index(int angleidx, int szaidx, int aziidx) const {
        return angleidx  + (int)m_cos_angle_grid->grid().size()  * szaidx + (int)m_cos_angle_grid->grid().size() * (int)m_sza_grid.grid().size() * aziidx + m_ground_start;
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceDiffuseStorage<NSTOKES, CNSTR>::accumulate_ground_sources(
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> &optical_layer, sasktran_disco::AEOrder m,
            sasktran2::DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage, int szaidx, int thread_idx) {
        // Here we evaluate the upwelling term at the surface.  Note that the sign convention switches compared to
        // the standard sign convention.  +<->- on W and Z

        // TODO: This entire function needs to change when switching to non-lambertian BRDF
        // Have to loop over outgoing angles, apply azimuthal factors, and get expansion coefficients
        if(m > 0) {
            return;
        }
        const sasktran_disco::uint N = m_config.num_do_streams() / 2;

        const auto& layer = optical_layer.bottom();
        const auto& solution = layer.solution(m);
        const auto& input_deriv = optical_layer.inputDerivatives();
        int layerStart = (int)input_deriv.layerStartIndex(layer.index());
        int numLayerDeriv = (int)input_deriv.numDerivativeLayer(layer.index());

        sasktran_disco::Radiance<NSTOKES> diffuse_contrib((int)input_deriv.numDerivative()), direct_contrib((int)input_deriv.numDerivative());
        diffuse_contrib.setzero();
        diffuse_contrib.setzero();

        direct_contrib.setzero();
        direct_contrib.setzero();

        // TODO: This is wrong but for lambertian albedo it doesn't matter
        // We should really be using some other function
        auto& rho = optical_layer.albedo(m).streamBDRFromStreams(0);

        // Construct the Dual quantity for the layer transmittance
        sasktran_disco::Dual<double> beam_transmittance = layer.dual_beamTransmittance(sasktran_disco::Location::FLOOR, input_deriv);
        sasktran_disco::Dual<double> stream_transmittance;

        for(sasktran_disco::StreamIndex i = 0; i < N * NSTOKES; ++i) {
            sasktran_disco::LayerDual<double> dual_rho((sasktran_disco::uint)input_deriv.numDerivativeLayer(layer.index()), layer.index(), (sasktran_disco::uint)input_deriv.layerStartIndex(layer.index()));
            // TODO: something is definitely wrong here, but it might be fine for scalar surface reflection...
            int s1 = i % NSTOKES;
            if( i % NSTOKES != 0) {
                dual_rho.value = 0.0;
                dual_rho.deriv.setZero();
            } else {
                dual_rho.value = rho[i/NSTOKES + N];

                for (sasktran_disco::uint j = 0; j < input_deriv.numDerivativeLayer(layer.index()); ++j) {
                    dual_rho.deriv(j) = input_deriv.layerDerivatives()[input_deriv.layerStartIndex(layer.index()) + j].d_albedo * sasktran_disco::kronDelta(m, 0);
                }
            }

            sasktran_disco::Radiance<NSTOKES> stream_contrib((int)input_deriv.numDerivative());
            if constexpr(NSTOKES == 1) {
                stream_contrib.value = solution.value.dual_Gplus_bottom().value(i);
            } else {
                stream_contrib.value(s1) = solution.value.dual_Gplus_bottom().value(i);
            }

            stream_contrib.deriv(Eigen::all, s1) = solution.value.dual_Gplus_bottom().deriv(Eigen::all, i);

            // Positive homogeneous solutions
            for(sasktran_disco::uint j = 0; j < N * NSTOKES; ++j) {
                stream_transmittance = layer.dual_streamTransmittance(sasktran_disco::Location::INSIDE, m, j, input_deriv);
                sasktran_disco::uint homogIndex = j*N*NSTOKES + i;

                if constexpr(NSTOKES == 1) {
                    stream_contrib.value += solution.boundary.L_coeffs.value(j) * solution.value.dual_homog_plus().value(homogIndex) * stream_transmittance.value;
                } else {
                    stream_contrib.value(s1) += solution.boundary.L_coeffs.value(j) * solution.value.dual_homog_plus().value(homogIndex) * stream_transmittance.value;
                }

                // LCoeffs have full derivatives, for some reason stream transmittance is a full dual even though it should be layer?
                for(sasktran_disco::uint k = 0; k < input_deriv.numDerivative(); ++k) {
                    stream_contrib.deriv(k, s1) += solution.boundary.L_coeffs.deriv(k, j) * solution.value.dual_homog_plus().value(homogIndex) * stream_transmittance.value;
                    stream_contrib.deriv(k, s1) += solution.boundary.L_coeffs.value(j) * solution.value.dual_homog_plus().value(homogIndex) * stream_transmittance.deriv(k);
                }
                // Homog only have layer derivs
                for(int k = 0; k < numLayerDeriv; ++k) {
                    stream_contrib.deriv(k + layerStart, s1) += solution.boundary.L_coeffs.value(j) * solution.value.dual_homog_plus().deriv(k, homogIndex) * stream_transmittance.value;
                }

                if constexpr(NSTOKES == 1) {
                    stream_contrib.value += solution.boundary.M_coeffs.value(j) * solution.value.dual_homog_minus().value(homogIndex);
                } else {
                    stream_contrib.value(s1) += solution.boundary.M_coeffs.value(j) * solution.value.dual_homog_minus().value(homogIndex);
                }

                // MCoeffs have full derivatives
                for(sasktran_disco::uint k = 0; k < input_deriv.numDerivative(); ++k) {
                    stream_contrib.deriv(k, s1) += solution.boundary.M_coeffs.deriv(k, j) * solution.value.dual_homog_minus().value(homogIndex);
                }
                // Homog only have layer derivs
                for(int k = 0; k < numLayerDeriv; ++k) {
                    stream_contrib.deriv(k + layerStart, s1) += solution.boundary.M_coeffs.value(j) * solution.value.dual_homog_minus().deriv(k, homogIndex);
                }
            }
            // Add stream i contribution
            double factor = (1.0 + sasktran_disco::kronDelta(m, 0)) * (*thread_storage.sza_calculators[0].persistent_config->quadrature_cos_angle())[i/NSTOKES] * (*thread_storage.sza_calculators[0].persistent_config->quadrature_weights())[i/NSTOKES];

            if constexpr(NSTOKES == 1) {
                diffuse_contrib.value += factor * stream_contrib.value * dual_rho.value;
            } else {
                diffuse_contrib.value(s1) += factor * stream_contrib.value(s1) * dual_rho.value;
            }

            // stream_contrib has full derivatives
            for(sasktran_disco::uint k = 0; k < input_deriv.numDerivative(); ++k) {
                diffuse_contrib.deriv(k, s1) += factor * stream_contrib.deriv(k, s1) * dual_rho.value;
            }
            // rho only have layer derivs
            for(int k = 0; k < numLayerDeriv; ++k) {
                if constexpr(NSTOKES == 1) {
                    diffuse_contrib.deriv(k + layerStart) += factor * stream_contrib.value * dual_rho.deriv(k);
                } else {
                    diffuse_contrib.deriv(k + layerStart, s1) += factor * stream_contrib.value(s1) * dual_rho.deriv(k);
                }
            }

        }

        // Always use 0 angle index to store the lambertian scattering result
        // TODO: Derivative propagation
        int source_index = ground_storage_index(0, szaidx, m);
        auto& storage = m_storage[thread_idx];

        if constexpr (NSTOKES == 3) {
            for(int s = 0; s < NSTOKES; ++s) {
                storage.source_terms_linear.value(source_index*NSTOKES + s) = diffuse_contrib.value(s);
            }
        } else {
            storage.source_terms_linear.value(source_index) = diffuse_contrib.value;
        }
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

        if(m == 0 && szaidx == 0) {
            // Have to reset some things
            storage.source_terms_linear.value.setZero();
            storage.source_terms_linear.deriv.setZero();
            m_converged_map.setConstant(false);
        }

        accumulate_ground_sources(optical_layer, m, thread_storage, szaidx, thread_idx);


        // Storage for derivatives of Y_plus and Y_minus,
        // TODO: Move these to cache
        std::vector<Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>> Y_plus_deriv;
        std::vector<Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>> Y_minus_deriv;

        const auto& input_derivatives = optical_layer.inputDerivatives();

        // TODO: Move to cache
        Eigen::Matrix<double, -1, NSTOKES> temp_deriv(input_derivatives.numDerivative(), NSTOKES);

        for(int lidx = 0; lidx < m_altitude_grid->grid().size(); ++lidx) {
            Y_minus_deriv.clear();
            Y_plus_deriv.clear();

            auto* layer = &optical_layer.layer((int)m_altitude_grid->grid().size() - lidx - 1);
            sasktran_disco::LayerIndex p = layer->index();
            double altitude = m_altitude_grid->grid()(lidx);

            double layer_fraction = (layer->altitude(sasktran_disco::Location::CEILING) - altitude) / (layer->altitude(sasktran_disco::Location::CEILING) - layer->altitude(sasktran_disco::Location::FLOOR));
            double x = layer_fraction * layer->dual_thickness().value;
            const auto& transmission = layer->dual_beamTransmittance(sasktran_disco::Location::CEILING, optical_layer.inputDerivatives());

            const auto& average_secant = layer->dual_average_secant();
            const auto& ssa = layer->dual_ssa();
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

            const auto& dual_L = solution.boundary.L_coeffs;
            const auto& dual_M = solution.boundary.M_coeffs;

            using MatrixView = Eigen::Map<Eigen::MatrixXd>;
            using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
            MatrixView Y_plus_matrix(Y_plus.value.data(), NSTOKES, m_config.num_do_streams()/2 *NSTOKES);
            MatrixView Y_minus_matrix(Y_minus.value.data(), NSTOKES, m_config.num_do_streams()/2 * NSTOKES);

            // Allocate memory for Y_plus/Y_minus derivative maps
            Y_plus_deriv.reserve(numderiv);
            Y_minus_deriv.reserve(numderiv);

            for(sasktran_disco::uint k = 0; k < numderiv; ++k) {
                Y_plus_deriv.emplace_back(&Y_plus.deriv(k, 0), NSTOKES, m_config.num_do_streams() / 2 * NSTOKES,
                                          Eigen::InnerStride<>(numderiv));
                Y_minus_deriv.emplace_back(&Y_minus.deriv(k, 0), NSTOKES, m_config.num_do_streams() / 2 * NSTOKES,
                                           Eigen::InnerStride<>(numderiv));
            }

            ConstMatrixView homog_plus_matrix(solution.value.dual_homog_plus().value.data(), m_config.num_do_streams()/2 * NSTOKES, m_config.num_do_streams()/2 *NSTOKES);
            ConstMatrixView homog_minus_matrix(solution.value.dual_homog_minus().value.data(), m_config.num_do_streams()/2 * NSTOKES, m_config.num_do_streams()/2 *NSTOKES);

            // Precalculate the layer multipliers since they are the same for every angle
            for (int i = 0; i < m_config.num_do_streams() / 2 * NSTOKES; ++i) {
                const auto &eigval = solution.value.dual_eigval();

                sasktran_disco::postprocessing::h_plus_sampled(layer->dual_thickness(),
                                                               eigval,
                                                               i,
                                                               layer_fraction,
                                                               hp[i]);

                sasktran_disco::postprocessing::h_minus_sampled(layer->dual_thickness(),
                                                                eigval,
                                                                i,
                                                                layer_fraction,
                                                                hm[i]);

                sasktran_disco::postprocessing::d_minus_sampled(layer->dual_thickness(),
                                                                eigval,
                                                                i,
                                                                layer_fraction,
                                                                transmission,
                                                                average_secant,
                                                                layerStart,
                                                                Dm[i]);

                sasktran_disco::postprocessing::d_plus_sampled(layer->dual_thickness(),
                                                               eigval,
                                                               i,
                                                               layer_fraction,
                                                               transmission,
                                                               average_secant,
                                                               layerStart,
                                                               Dp[i]);
            }

            for(int aidx = 0; aidx < m_cos_angle_grid->grid().size(); ++aidx) {
                int sourceidx = linear_storage_index(aidx, lidx, szaidx, m);

                if(!m_need_to_calculate_map[sourceidx] || m_converged_map[sourceidx]) {
                    continue;
                }

                temp_deriv.setZero();

                storage.phase_storage[aidx].set_phase_container(storage.phase_container, m);

                // Calculate legendre sums multiplied by stream weights
                sasktran_disco::VectorLayerDual<double>& dual_lpsum_plus = cache.dual_lpsum_plus;
                sasktran_disco::VectorLayerDual<double>& dual_lpsum_minus = cache.dual_lpsum_minus;

                layer->vectordual_scatPhaseF(m, storage.phase_container, optical_layer.inputDerivatives(), dual_lpsum_minus, dual_lpsum_plus);

                // dual_lpsum_plus and minus are now linearly storage of N * NSTOKES * NSTOKES, have to multiply the blocks by
                // the weights
                for(int i = 0; i < this->m_config.num_do_streams()/2; ++i) {
                    for(int j = 0; j < NSTOKES*NSTOKES; ++j) {
                        int linearindex = i*NSTOKES*NSTOKES + j;
                        dual_lpsum_minus.value(linearindex) *= thread_storage.sza_calculators[0].persistent_config->quadrature_weights()->at(i);
                        dual_lpsum_plus.value(linearindex) *= thread_storage.sza_calculators[0].persistent_config->quadrature_weights()->at(i);

                        dual_lpsum_minus.deriv(Eigen::all, linearindex) *= thread_storage.sza_calculators[0].persistent_config->quadrature_weights()->at(i);
                        dual_lpsum_plus.deriv(Eigen::all, linearindex) *= thread_storage.sza_calculators[0].persistent_config->quadrature_weights()->at(i);
                    }
                }

                ConstMatrixView lpsum_plus_matrix(dual_lpsum_plus.value.data(), NSTOKES, m_config.num_do_streams()/2 * NSTOKES);
                ConstMatrixView lpsum_minus_matrix(dual_lpsum_minus.value.data(), NSTOKES, m_config.num_do_streams()/2 * NSTOKES);

                Y_plus_matrix.noalias() = lpsum_plus_matrix * homog_plus_matrix + lpsum_minus_matrix * homog_minus_matrix;
                Y_minus_matrix.noalias() = lpsum_plus_matrix * homog_minus_matrix + lpsum_minus_matrix * homog_plus_matrix;

                // Calculate derivatives of Y_plus/Y_minus
                for(sasktran_disco::uint k = 0; k < numderiv; ++k) {
                    Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> lpsum_plus_deriv(
                            &dual_lpsum_plus.deriv(k, 0), NSTOKES, m_config.num_do_streams()/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

                    Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> lpsum_minus_deriv(
                            &dual_lpsum_minus.deriv(k, 0), NSTOKES, m_config.num_do_streams()/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

                    Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> homog_minus_deriv(
                            &solution.value.dual_homog_minus().deriv(k, 0), m_config.num_do_streams()/2 * NSTOKES, m_config.num_do_streams()/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

                    Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> homog_plus_deriv(
                            &solution.value.dual_homog_plus().deriv(k, 0), m_config.num_do_streams()/2 * NSTOKES, m_config.num_do_streams()/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

                    Y_plus_deriv[k].noalias() = lpsum_plus_deriv * homog_plus_matrix + lpsum_plus_matrix * homog_plus_deriv +
                                                lpsum_minus_deriv * homog_minus_matrix + lpsum_minus_matrix * homog_minus_deriv;

                    Y_minus_deriv[k].noalias() = lpsum_plus_deriv * homog_minus_matrix + lpsum_plus_matrix * homog_minus_deriv +
                                                 lpsum_minus_deriv * homog_plus_matrix + lpsum_minus_matrix * homog_plus_deriv;
                }

                for (int i = 0; i < m_config.num_do_streams() / 2 * NSTOKES; ++i) {
                    const auto& hpi = hp[i];
                    const auto& hmi = hm[i];
                    const auto& Dmi = Dm[i];
                    const auto& Dpi = Dp[i];

                    for(int s = 0; s < NSTOKES; ++s) {
                        storage.source_terms_linear.value(sourceidx*NSTOKES + s) += Y_plus_matrix(s, i) * hpi.value * dual_L.value(i);
                        storage.source_terms_linear.value(sourceidx*NSTOKES + s) += Y_minus_matrix(s, i) * hmi.value * dual_M.value(i);

                        // Y_plus/Y_minus and hp/hm only have layer derivatives, but dual_L/dual_M are dense derivatives

                        // Start with the dense ones
                        temp_deriv(Eigen::all, s) += dual_L.deriv(Eigen::all, i) * Y_plus_matrix(s, i) * hpi.value;
                        temp_deriv(Eigen::all, s) += dual_M.deriv(Eigen::all, i) * Y_minus_matrix(s, i) * hmi.value;

                        // And add in the layer derivatives
                        for(int k = 0; k < numderiv; ++k) {
                            // For Yplus/Yminus
                            temp_deriv(k + layerStart, s) += Y_plus_deriv[k](s, i) * hpi.value * dual_L.value(i);
                            temp_deriv(k + layerStart, s) += Y_minus_deriv[k](s, i) * hmi.value * dual_M.value(i);

                            // for Hp, Hm
                            temp_deriv(k + layerStart, s) += hpi.deriv(k) * Y_plus_matrix(s, i) * dual_L.value(i);
                            temp_deriv(k + layerStart, s) += hmi.deriv(k) * Y_minus_matrix(s, i) * dual_M.value(i);
                        }

                    }

                    for(int s = 0; s < NSTOKES; ++s) {
                        storage.source_terms_linear.value(sourceidx*NSTOKES + s) += (dual_Aplus.value(i) * Y_plus_matrix(s, i) * Dmi.value + dual_Aminus.value(i) * Y_minus_matrix(s, i) * Dpi.value);

                        // Dp/Dm has dense derivatives, the others are small

                        // Start with the dense ones
                        temp_deriv(Eigen::all, s) += Dmi.deriv * dual_Aplus.value(i) * Y_plus_matrix(s, i);
                        temp_deriv(Eigen::all, s) += Dpi.deriv * dual_Aminus.value(i) * Y_minus_matrix(s, i);

                        // And add in the layer derivatives
                        for(int k = 0; k < numderiv; ++k) {
                            // For Yplus/Yminus
                            temp_deriv(k + layerStart, s) += Y_plus_deriv[k](s, i) * Dmi.value * dual_Aplus.value(i);
                            temp_deriv(k + layerStart, s) += Y_minus_deriv[k](s, i) * Dpi.value * dual_Aminus.value(i);

                            // for Aplus/Aminu
                            temp_deriv(k + layerStart, s) += dual_Aplus.deriv(k, i) * Y_plus_matrix(s, i) * Dmi.value;
                            temp_deriv(k + layerStart, s) += dual_Aminus.deriv(k, i) * Y_minus_matrix(s, i) * Dpi.value;
                        }

                    }
                }

                // Now we have to divide out SSA from the source since we opt to multiply by actual SSA later on
                for(int s = 0; s < NSTOKES; ++s) {
                    int index = sourceidx*NSTOKES + s;

                    storage.source_terms_linear.value(index) /= ssa.value;

                    // Check for convergence
                    if(m >= 2 && s == 0) {
                        int l1_index = linear_storage_index(aidx, lidx, szaidx, m-1);
                        int l2_index = linear_storage_index(aidx, lidx, szaidx, m-2);

                        if(abs(storage.source_terms_linear.value(index) / storage.source_terms_linear.value(l1_index*NSTOKES)) < 1e-4) {
                            if (abs(storage.source_terms_linear.value(index) /
                                    storage.source_terms_linear.value(l2_index * NSTOKES)) < 1e-4) {
                                for(int azi = m+1; azi < m_config.num_do_streams(); ++azi) {
                                    int converged_index = linear_storage_index(aidx, lidx, szaidx, azi);
                                    m_converged_map[converged_index] = true;
                                }
                            }
                        }
                    }

                    for(int k = 0; k < numderiv; ++k) {
                        temp_deriv(k + layerStart, s) -= ssa.deriv(k) * storage.source_terms_linear.value(index);
                    }
                    temp_deriv(Eigen::all, s) /= ssa.value;

                    // And we also have to translate the temporary layer DO derivatives to atmosphere derivatives
                    if(numtotalderiv > 0) {
                        storage.source_terms_linear.deriv(index, Eigen::all).setZero();
                    }
                    for(int k = 0; k < numtotalderiv; ++k) {
                        for(int l = 0; l < input_derivatives.layerDerivatives()[k].group_and_triangle_fraction.size(); ++l) {
                            const std::pair<sasktran_disco::uint, double>& group_fraction = input_derivatives.layerDerivatives()[k].group_and_triangle_fraction[l];
                            const auto& extinction = input_derivatives.layerDerivatives()[k].extinctions[l];

                            storage.source_terms_linear.deriv(index, group_fraction.first) += group_fraction.second * temp_deriv(k, s) * extinction;
                        }
                    }
                }
            }
        }
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceDiffuseStorage);
}