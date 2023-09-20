#include "sasktran2/do_source.h"


namespace sasktran2 {
    template<int NSTOKES, int CNSTR>
    DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::DOSourceInterpolatedPostProcessing(const sasktran2::Geometry1D &geometry,
                                                                                           const sasktran2::raytracing::RayTracerBase &raytracer,
                                                                                           bool will_integrate_sources
                                                                                           )
            : DOSource<NSTOKES, CNSTR>(geometry, raytracer), m_will_integrate_sources(will_integrate_sources) {
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::accumulate_solved_azimuth(
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> &optical_layer,
            DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage, int szaidx, sasktran_disco::AEOrder m,
            int threadidx) {
        m_diffuse_storage->accumulate_sources(optical_layer, m, thread_storage, szaidx, threadidx);
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::calculate(int wavelidx, int threadidx) {
        DOSource<NSTOKES, CNSTR>::calculate(wavelidx, threadidx);

    }

    template<int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {
        DOSource<NSTOKES, CNSTR>::initialize_geometry(los_rays);

        m_diffuse_storage = std::make_unique<sasktran2::DOSourceDiffuseStorage<NSTOKES, CNSTR>>(
                *(this->m_thread_storage)[0].sza_calculators[0].geometry_layers.get(),
                *(this->m_thread_storage)[0].sza_calculators[0].persistent_config.get(),
                *this->m_sza_grid,
                *this->m_config,
                this->m_geometry
        );

        if(m_will_integrate_sources) {
            m_los_source_interpolator = m_diffuse_storage->geometry_interpolator(los_rays);
            m_source_interpolator_view = m_los_source_interpolator.get();
        }

    }

    template<int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere) {
        m_atmosphere = &atmosphere;

        DOSource<NSTOKES, CNSTR>::initialize_atmosphere(atmosphere);

        m_diffuse_storage->initialize_atmosphere(atmosphere);
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::initialize_config(const sasktran2::Config &config) {
        DOSource<NSTOKES, CNSTR>::initialize_config(config);
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::integrated_source(int wavelidx, int losidx, int layeridx, int threadidx,
                                                     const sasktran2::raytracing::SphericalLayer &layer,
                                                     const sasktran2::SparseODDualView &shell_od,
                                                     sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        if(layer.layer_distance < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do anything
            return;
        }

        auto& interpolator = this->m_los_source_weights[losidx][layeridx];

        // Start by calculating ssa at the source point
        double omega = 0;
        for(int i = 0; i < interpolator.size(); ++i) {
            auto& index_weight = interpolator[i];
            omega += this->m_atmosphere->storage().ssa(index_weight.first, wavelidx) * index_weight.second;
        }

        double source_factor = (1-shell_od.exp_minus_od);

        for(int s = 0; s < NSTOKES; ++s) {
            // Need some temporaries
            double source_value = (*m_source_interpolator_view)[losidx][layeridx][s].dot(m_diffuse_storage->linear_source(threadidx).value);

            source.value(s) += omega * source_factor * source_value;

            if(m_atmosphere->num_deriv() > 0) {
                // Now we need dJ/dthickness
                for(auto it = shell_od.deriv_iter; it; ++it) {
                    source.deriv(s, it.index()) += it.value() * (1 - source_factor) * source_value * omega;
                }

                // And dJ/dssa
                for(auto& ele : interpolator) {
                    source.deriv(s, m_atmosphere->ssa_deriv_start_index() + ele.first) += ele.second * source_factor * source_value;
                }

                if(this->m_config->wf_precision() == sasktran2::Config::WeightingFunctionPrecision::full) {
                    // The derivatives are omega * source_factor * (interpolator @ source_deriv_storage)
                    // This seems to be the fastest way to do the calculation
                    for(auto it = Eigen::SparseVector<double>::InnerIterator((*m_source_interpolator_view)[losidx][layeridx][s]); it; ++it) {
                        source.deriv(s, Eigen::all) += omega * source_factor * it.value() * m_diffuse_storage->linear_source(threadidx).deriv(it.index(), Eigen::all);
                    }
                }
            }
        }
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceInterpolatedPostProcessing);
}