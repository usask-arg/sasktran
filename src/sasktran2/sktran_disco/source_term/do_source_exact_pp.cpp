#include "sasktran2/do_source.h"


namespace sasktran2 {
    template<int NSTOKES, int CNSTR>
    DOSourceExactPostProcessing<NSTOKES, CNSTR>::DOSourceExactPostProcessing(const sasktran2::Geometry1D &geometry,
                                                                             const sasktran2::raytracing::RayTracerBase &raytracer) : DOSource<NSTOKES, CNSTR>(geometry, raytracer) {

    }

    template <int NSTOKES, int CNSTR>
    void DOSourceExactPostProcessing<NSTOKES, CNSTR>::accumulate_layer_boundary_sources(const sasktran2::raytracing::SphericalLayer& layer,
                                                                                        const sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& layer_array,
                                                                                        int ray_idx, int layer_idx, DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage, sasktran_disco::AEOrder m) {
        accumulate_layer_source(layer_array,
                                thread_storage.legendre_phase_container,
                                m_lp_segments[ray_idx],
                                layer.entrance,
                                layer.exit,
                                layer.saz_entrance,
                                layer.saz_exit,
                                layer_idx,
                                0, 1,
                                thread_storage.layer_sources[ray_idx],
                                m,
                                thread_storage);
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceExactPostProcessing<NSTOKES, CNSTR>::accumulate_layer_source(
            const sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> &layer_array,
            std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>> &phase_storage,
            const std::vector<LegendrePhaseStorage<NSTOKES, CNSTR>> &lp_segments,
            const sasktran2::Location &start_location, const sasktran2::Location &end_location, double start_saa,
            double end_saa, int source_idx, int sza_index, int sza_weight,
            sasktran2::Dual<double, sasktran2::dualstorage::dense> &source, sasktran_disco::AEOrder m,
            DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage) {
        double altitude = (start_location.radius() + end_location.radius())/2.0 - this->m_geometry.coordinates().earth_radius();
        const auto* layer = layer_array.layerAtAltitude(altitude);
        bool is_upwelling = start_location.radius() < end_location.radius();

        for(int start_end = 0; start_end < 2; start_end++)
        {
            lp_segments[source_idx + start_end].set_phase_container(phase_storage, m);
            const auto& location = start_end == 0 ? end_location : start_location;
            double saa = start_end == 0 ? end_saa : start_saa;

            double altitude = location.radius() - this->m_geometry.coordinates().earth_radius();

            double layer_fraction = (layer->altitude(sasktran_disco::Location::CEILING) - altitude) / (layer->altitude(sasktran_disco::Location::CEILING) - layer->altitude(sasktran_disco::Location::FLOOR));
            double x = layer_fraction * layer->dual_thickness().value;
            const auto& transmission = layer->dual_beamTransmittance(sasktran_disco::Location::CEILING, layer_array.inputDerivatives());

            const auto& average_secant = layer->dual_average_secant();
            const sasktran_disco::LayerDual<double>& dual_thickness = layer->dual_thickness();


            const auto& solution = layer->solution(m);

            sasktran_disco::LayerIndex p = layer->index();

            auto& cache = thread_storage.postprocessing_cache[p];
            sasktran_disco::uint numderiv = (sasktran_disco::uint)layer_array.inputDerivatives().numDerivativeLayer(p);
            sasktran_disco::uint numtotalderiv = (sasktran_disco::uint)layer_array.inputDerivatives().numDerivative();
            int layerStart = (int)layer_array.inputDerivatives().layerStartIndex(p);

            cache.resize(this->m_config->num_do_streams(), p, numderiv, layerStart, numtotalderiv);

            // Calculate legendre sums multiplied by stream weights
            sasktran_disco::VectorLayerDual<double>& dual_lpsum_plus = cache.dual_lpsum_plus;
            sasktran_disco::VectorLayerDual<double>& dual_lpsum_minus = cache.dual_lpsum_minus;

            layer->vectordual_scatPhaseF(m, phase_storage, layer_array.inputDerivatives(), dual_lpsum_minus, dual_lpsum_plus);

            // dual_lpsum_plus and minus are now linearly storage of N * NSTOKES * NSTOKES, have to multiply the blocks by
            // the weights
            for(int i = 0; i < this->m_config->num_do_streams()/2; ++i) {
                for(int j = 0; j < NSTOKES*NSTOKES; ++j) {
                    int linearindex = i*NSTOKES*NSTOKES + j;
                    dual_lpsum_minus.value(linearindex) *= thread_storage.sza_calculators[sza_index].persistent_config->quadrature_weights()->at(i);
                    dual_lpsum_plus.value(linearindex) *= thread_storage.sza_calculators[sza_index].persistent_config->quadrature_weights()->at(i);

                    dual_lpsum_minus.deriv(Eigen::all, linearindex) *= thread_storage.sza_calculators[sza_index].persistent_config->quadrature_weights()->at(i);
                    dual_lpsum_plus.deriv(Eigen::all, linearindex) *= thread_storage.sza_calculators[sza_index].persistent_config->quadrature_weights()->at(i);
                }
            }

            const auto& dual_particular_plus = solution.value.dual_particular_plus();
            const auto& dual_particular_minus = solution.value.dual_particular_minus();

            const auto& dual_Aplus = solution.value.dual_green_A_plus();
            const auto& dual_Aminus = solution.value.dual_green_A_minus();

            auto& hp = cache.hp;
            auto& hm = cache.hm;
            auto& J = cache.J;

            auto& Y_plus = cache.Y_plus;
            auto& Y_minus = cache.Y_minus;

            auto& Dm = cache.Dm;
            auto& Dp = cache.Dp;

            const auto& dual_L = solution.boundary.L_coeffs;
            const auto& dual_M = solution.boundary.M_coeffs;

            using MatrixView = Eigen::Map<Eigen::MatrixXd>;
            using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
            MatrixView Y_plus_matrix(Y_plus.value.data(), NSTOKES, this->m_config->num_do_streams()/2 *NSTOKES);
            MatrixView Y_minus_matrix(Y_minus.value.data(), NSTOKES, this->m_config->num_do_streams()/2 * NSTOKES);

            ConstMatrixView lpsum_plus_matrix(dual_lpsum_plus.value.data(), NSTOKES, this->m_config->num_do_streams()/2 * NSTOKES);
            ConstMatrixView lpsum_minus_matrix(dual_lpsum_minus.value.data(), NSTOKES, this->m_config->num_do_streams()/2 * NSTOKES);

            ConstMatrixView homog_plus_matrix(solution.value.dual_homog_plus().value.data(), this->m_config->num_do_streams()/2 * NSTOKES, this->m_config->num_do_streams()/2 *NSTOKES);
            ConstMatrixView homog_minus_matrix(solution.value.dual_homog_minus().value.data(), this->m_config->num_do_streams()/2 * NSTOKES, this->m_config->num_do_streams()/2 *NSTOKES);

            Y_plus_matrix.noalias() = lpsum_plus_matrix * homog_plus_matrix + lpsum_minus_matrix * homog_minus_matrix;
            Y_minus_matrix.noalias() = lpsum_plus_matrix * homog_minus_matrix + lpsum_minus_matrix * homog_plus_matrix;

            J.setzero();

            for (int i = 0; i < this->m_config->num_do_streams() / 2 * NSTOKES; ++i) {
                const auto& eigval = solution.value.dual_eigval();

                hp.value = std::exp(-1.0*eigval.value(i) * layer->dual_thickness().value * layer_fraction);
                hm.value = std::exp(-eigval.value(i) * (layer->dual_thickness().value - layer->dual_thickness().value * layer_fraction));

                if constexpr(NSTOKES==1) {
                    J.value += Y_plus_matrix(0, i) * hp.value * dual_L.value(i);
                    J.value += Y_minus_matrix(0, i) * hm.value * dual_M.value(i);
                } else {
                    J.value += Y_plus_matrix(Eigen::all, i) * hp.value * dual_L.value(i);
                    J.value += Y_minus_matrix(Eigen::all, i) * hm.value * dual_M.value(i);
                }

                Dp.value = transmission.value / (eigval.value(i) + average_secant.value) * (exp(-x*average_secant.value) - exp(-dual_thickness.value * average_secant.value) * exp(-(dual_thickness.value - x) * eigval.value(i)));
                Dm.value = transmission.value / (average_secant.value - eigval.value(i)) * (exp(-x*eigval.value(i)) - exp(-x*average_secant.value));

                if constexpr(NSTOKES==1) {
                    J.value += dual_Aplus.value(i) * Y_plus_matrix(0, i) * Dm.value + dual_Aminus.value(i) * Y_minus_matrix(0, i) * Dp.value;
                } else {
                    J.value += dual_Aplus.value(i) * Y_plus_matrix(Eigen::all, i) * Dm.value + dual_Aminus.value(i) * Y_minus_matrix(Eigen::all, i) * Dp.value;
                }
            }

            int source_start_idx = source_idx*NSTOKES;
            if constexpr(NSTOKES==1) {
                double azimuthal_factor = cos(m*saa);
                source.value(source_start_idx) += J.value * azimuthal_factor * 0.5;

                if(isnan(J.value)) {
                    BOOST_LOG_TRIVIAL(error) << "NaN In DO Source";
                }

            } else if constexpr(NSTOKES == 3) {
                double cos_daz = cos(m*saa);

                // Azimuth in expansion differs from our definition
                double sin_daz = sin(m*(EIGEN_PI -saa));

                source.value(source_start_idx) += J.value(0) * cos_daz * 0.5;
                source.value(source_start_idx+1) += J.value(1) * cos_daz * 0.5;
                source.value(source_start_idx+2) += J.value(2) * sin_daz * 0.5;
            }
        }
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceExactPostProcessing<NSTOKES, CNSTR>::accumulate_solved_azimuth(
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> &optical_layer,
            DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage, int szaidx, sasktran_disco::AEOrder m,
            int threadidx) {

        for(int i = 0; i < this->m_los_rays->size(); ++i) {
            const auto& ray = this->m_los_rays->at(i);
            for(int j = 0; j < ray.layers.size(); ++j) {
                accumulate_layer_boundary_sources(ray.layers[j], optical_layer, i, j, thread_storage, m);
            }
        }
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceExactPostProcessing<NSTOKES, CNSTR>::initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {
        DOSource<NSTOKES, CNSTR>::initialize_geometry(los_rays);

        // Calculate legendre polynomials for every cell and initialize the boundary sources
        m_lp_segments.resize(los_rays.size());

        for(int thidx = 0; thidx < this->m_thread_storage.size(); ++thidx) {
            this->m_thread_storage[thidx].boundary_sources.resize(los_rays.size());
            this->m_thread_storage[thidx].layer_sources.resize(los_rays.size());
        }

        for(int i = 0; i < m_lp_segments.size(); ++i) {
            auto& segment = m_lp_segments[i];
            const auto& ray = los_rays[i];

            // 1 more boundary than cell
            segment.resize(ray.layers.size() + 1, this->m_config->num_do_streams());

            for(int thidx = 0; thidx < this->m_thread_storage.size(); ++thidx) {
                this->m_thread_storage[thidx].boundary_sources[i].resize(NSTOKES * (int)(ray.layers.size() + 1), 0, true);
                this->m_thread_storage[thidx].layer_sources[i].resize(NSTOKES * (int)(ray.layers.size()), 0, true);
            }

            for(int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];
                double cos_viewing_entrance = layer.entrance.cos_zenith_angle(ray.observer_and_look.look_away);
                double cos_viewing_exit = layer.exit.cos_zenith_angle(ray.observer_and_look.look_away);
                if(j == 0) {
                    m_lp_segments[i][0].fill(cos_viewing_exit);
                }
                m_lp_segments[i][j+1].fill(cos_viewing_entrance);
            }
        }

    }

    template <int NSTOKES, int CNSTR>
    void DOSourceExactPostProcessing<NSTOKES, CNSTR>::integrated_source(int wavelidx, int losidx, int layeridx, int threadidx,
                                                                               const sasktran2::raytracing::SphericalLayer &layer,
                                                                               const sasktran2::SparseODDualView &shell_od,
                                                                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        if(abs(layer.entrance.radius() - layer.exit.radius()) < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do anything
            return;
        }
        double source_factor = 1-exp(-shell_od.od);

        // TODO: Linearization

        // Exact calculation?
        const auto& sources = this->m_thread_storage[threadidx].layer_sources[losidx];

        Eigen::Vector<double, NSTOKES> layer_source = sources.value(Eigen::seq((layeridx)*NSTOKES, (layeridx)*NSTOKES + NSTOKES - 1));

        source.value += source_factor * layer_source;
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceExactPostProcessing)
}