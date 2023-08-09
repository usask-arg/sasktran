#include "sasktran2/do_source.h"

namespace sasktran2 {
    template <int NSTOKES, int CNSTR>
    DOSource<NSTOKES, CNSTR>::DOSource(const sasktran2::Geometry1D& geometry, const sasktran2::raytracing::RayTracerBase& raytracer) : m_geometry(geometry), m_raytracer(raytracer) {
        // Initialize pointers
        m_atmosphere = nullptr;
        m_los_rays = nullptr;
    }

    template <int NSTOKES, int CNSTR>
    void DOSource<NSTOKES, CNSTR>::calculate(int wavelidx, int threadidx)  {
        for(int i = 0; i < m_thread_storage[threadidx].boundary_sources.size(); ++i) {
            auto& source = m_thread_storage[threadidx].boundary_sources[i];

            source.value.setZero();
        }

        for(int i = 0; i < m_thread_storage[threadidx].layer_sources.size(); ++i) {
            auto& source = m_thread_storage[threadidx].layer_sources[i];

            source.value.setZero();
        }

        for(int szaidx = 0; szaidx < m_thread_storage[threadidx].sza_calculators.size(); ++szaidx) {
            auto& solver = m_thread_storage[threadidx].sza_calculators[szaidx];
            std::unique_ptr<sasktran_disco::BRDF_Base> brdf;

            // TODO: Make non-lambertian BRDF
            brdf = std::make_unique<sasktran_disco::TestBRDF>(m_atmosphere->surface().albedo()(wavelidx));

            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> optical_layer(
                    *solver.persistent_config,
                    wavelidx,
                    m_do_los,
                    std::move(brdf),
                    *solver.geometry_layers,
                    *m_atmosphere
                    );

            sasktran_disco::RTESolver<NSTOKES, CNSTR> rte(*solver.persistent_config, optical_layer);

            int num_azi = m_config->num_do_streams();
            for(int m = 0; m < num_azi; ++m) {
                rte.solve(m);
                accumulate_solved_azimuth(optical_layer, m_thread_storage[threadidx], szaidx, m, threadidx);
            }
        }
    }

    template<int NSTOKES, int CNSTR>
    void DOSource<NSTOKES, CNSTR>::generate_sza_grid() {
        // Go through all of the los rays and find the maximum and minimum SZA
        double min_cos_sza = 1;
        double max_cos_sza = -1;

        for(const auto& ray : *m_los_rays) {
            for(const auto& layer : ray.layers) {
                min_cos_sza = std::min(layer.cos_sza_entrance, min_cos_sza);
                min_cos_sza = std::min(layer.cos_sza_exit, min_cos_sza);

                max_cos_sza = std::max(layer.cos_sza_entrance, max_cos_sza);
                max_cos_sza = std::max(layer.cos_sza_exit, max_cos_sza);
            }
        }

        int num_sza = m_config->num_do_sza();
        double ref_cos_sza = m_geometry.coordinates().cos_sza_at_reference();
        Eigen::VectorXd sza_grid;

        if(num_sza == 1) {
            // Special case where we just want the SZA at the tangent point

            sza_grid.resize(1);
            sza_grid(0) = ref_cos_sza;

            m_sza_grid = std::make_unique<sasktran2::grids::Grid>(std::move(sza_grid), sasktran2::grids::gridspacing::constant,
                                                                  sasktran2::grids::outofbounds::extend,
                                                                  sasktran2::grids::interpolation::linear);
        } else {
            sza_grid = Eigen::ArrayXd::LinSpaced(num_sza, min_cos_sza, max_cos_sza);

            m_sza_grid = std::make_unique<sasktran2::grids::Grid>(std::move(sza_grid), sasktran2::grids::gridspacing::constant,
                                                                  sasktran2::grids::outofbounds::extend,
                                                                  sasktran2::grids::interpolation::linear);
        }

    }


    template <int NSTOKES, int CNSTR>
    void DOSource<NSTOKES, CNSTR>::initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {
        m_los_rays = &los_rays;

        // Create the SZA grid
        generate_sza_grid();

        for(int thidx = 0; thidx < m_thread_storage.size(); ++thidx) {
            auto& sza_calculators = m_thread_storage[thidx].sza_calculators;
            for(int i = 0; i < sza_calculators.size(); ++i) {
                auto& sza_calculator = sza_calculators[i];
                double cos_sza = m_sza_grid->grid()(i);
                sza_calculator.persistent_config->configure(
                        sza_calculator.userspec,
                        *m_config,
                        cos_sza,
                        (int)m_geometry.altitude_grid().grid().size() - 1,
                        los_rays
                        );


                sza_calculator.geometry_layers = std::make_unique<sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>>(*sza_calculator.persistent_config, m_geometry);
            }
            // These quantities aren't needed for each SZA
            m_thread_storage[thidx].postprocessing_cache.resize(m_geometry.altitude_grid().grid().size() - 1);
            m_thread_storage[thidx].legendre_phase_container.resize(m_config->num_do_streams());
        }
    }

    template <int NSTOKES, int CNSTR>
    void DOSource<NSTOKES, CNSTR>::initialize_config(const sasktran2::Config &config) {
        m_config = &config;
        // Create the thread storage
        m_thread_storage.resize(config.num_threads());

        m_nstr = config.num_do_streams();

        int num_sza = config.num_do_sza();
        for(int i = 0; i < m_thread_storage.size(); ++i) {
            auto& sza_calculators = m_thread_storage[i].sza_calculators;
            sza_calculators.resize(num_sza);
            for( auto& calculator : sza_calculators ) {
                calculator.persistent_config = std::make_unique<sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>>();
            }
        }
    }

    template <int NSTOKES, int CNSTR>
    void DOSource<NSTOKES, CNSTR>::initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere) {
        m_atmosphere = &atmosphere;
    }



    template <int NSTOKES, int CNSTR>
    void DOSource<NSTOKES, CNSTR>::end_of_ray_source(int wavelidx, int losidx, int threadidx,
                                                     sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        // Postprocessing for ground source
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSource)
}