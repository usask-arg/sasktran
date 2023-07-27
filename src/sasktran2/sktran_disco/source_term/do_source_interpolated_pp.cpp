#include "sasktran2/do_source.h"


namespace sasktran2 {
    template<int NSTOKES, int CNSTR>
    DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::DOSourceInterpolatedPostProcessing(const sasktran2::Geometry1D &geometry,
                                                                                           const sasktran2::raytracing::RayTracerBase &raytracer)
            : DOSource<NSTOKES, CNSTR>(geometry, raytracer) {
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

        m_los_source_interpolator = m_diffuse_storage->geometry_interpolator(los_rays);
        m_source_interpolator_view = m_los_source_interpolator.get();

    }

    template<int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere) {
        DOSource<NSTOKES, CNSTR>::initialize_atmosphere(atmosphere);
    }

    template <int NSTOKES, int CNSTR>
    void DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::integrated_source(int wavelidx, int losidx, int layeridx, int threadidx,
                                                     const sasktran2::raytracing::SphericalLayer &layer,
                                                     const sasktran2::SparseODDualView &shell_od,
                                                     sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        if(abs(layer.entrance.radius() - layer.exit.radius()) < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do anything
            return;
        }
        double source_factor = 1-exp(-shell_od.od);

        // TODO: Linearization

        for(int s = 0; s < NSTOKES; ++s) {
            source.value(s) += source_factor * ((*m_source_interpolator_view)[losidx][layeridx][s].dot(m_diffuse_storage->source(threadidx)[s].value));
        }
    }


    template<int NSTOKES, int CNSTR>
    DOSourceInterpolatedPostProcessingView<NSTOKES, CNSTR>::DOSourceInterpolatedPostProcessingView(
            const sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>> &interpolator,
            const DOSourceDiffuseStorage<NSTOKES, CNSTR> &source_storage) : m_source_interpolator(interpolator), m_diffuse_storage(source_storage) {

    }

    template<int NSTOKES, int CNSTR>
    void
    DOSourceInterpolatedPostProcessingView<NSTOKES, CNSTR>::integrated_source(int wavelidx, int losidx, int layeridx,
                                                                              int threadidx,
                                                                              const sasktran2::raytracing::SphericalLayer &layer,
                                                                              const sasktran2::SparseODDualView &shell_od,
                                                                              sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        int m = 0; // azi order

        auto azi_seq = Eigen::seq(m*m_diffuse_storage.num_azi_terms(), (m+1)*m_diffuse_storage.num_azi_terms()-1);

        if(abs(layer.entrance.radius() - layer.exit.radius()) < MINIMUM_SHELL_SIZE_M) {
            // Essentially an empty shell from rounding, don't have to do anything
            return;
        }
        double source_factor = 1-exp(-shell_od.od);

        // TODO: Linearization?

        for(int s = 0; s < NSTOKES; ++s) {
            source.value(s) += source_factor * ((m_source_interpolator)[losidx][layeridx][s].dot(m_diffuse_storage.source(threadidx)[s].value(azi_seq)));
            source.value(s) += source_factor * ((m_source_interpolator)[losidx][layeridx][s].dot(m_diffuse_storage.singlescatter_source(threadidx)[s].value(azi_seq)));
        }
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceInterpolatedPostProcessing);

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceInterpolatedPostProcessingView);
}