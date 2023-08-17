#include <sasktran2.h>
#include <omp.h>

template <int NSTOKES>
void Sasktran2<NSTOKES>::construct_raytracer() {
    // TODO: Use config to determine what type of raytracer we are going to use
    m_raytracer = std::make_unique<sasktran2::raytracing::SphericalShellRayTracer>(*m_geometry);
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::construct_integrator() {
    m_source_integrator = std::make_unique<sasktran2::SourceIntegrator<NSTOKES>>(true);
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::construct_source_terms() {

    if(m_config.single_scatter_source() == sasktran2::Config::SingleScatterSource::exact) {
        m_source_terms.emplace_back(
                std::make_unique<sasktran2::solartransmission::SingleScatterSource<sasktran2::solartransmission::SolarTransmissionExact, NSTOKES>>(*m_geometry, *m_raytracer)
        );

        m_los_source_terms.push_back(m_source_terms[m_source_terms.size() - 1].get());
    }

    if(m_config.single_scatter_source() == sasktran2::Config::SingleScatterSource::solartable) {
        m_source_terms.emplace_back(
                std::make_unique<sasktran2::solartransmission::SingleScatterSource<sasktran2::solartransmission::SolarTransmissionTable, NSTOKES>>(*m_geometry, *m_raytracer)
        );

        m_los_source_terms.push_back(m_source_terms[m_source_terms.size() - 1].get());
    }

    if(m_config.occultation_source() == sasktran2::Config::OccultationSource::standard) {
        m_source_terms.emplace_back(
                std::make_unique<sasktran2::solartransmission::OccultationSource<NSTOKES>>()
        );

        m_los_source_terms.push_back(m_source_terms[m_source_terms.size() - 1].get());
    }

    if(m_config.multiple_scatter_source() == sasktran2::Config::MultipleScatterSource::discrete_ordinates) {

        #ifdef SASKTRAN_DISCO_FULL_COMPILE
        if constexpr (NSTOKES == 1) {
            if(m_config.num_do_streams() == 2) {
                m_source_terms.emplace_back(
                        std::make_unique<sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, 2>>(*m_geometry, *m_raytracer)
                );
            } else {
                m_source_terms.emplace_back(
                        std::make_unique<sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(*m_geometry, *m_raytracer)
                );
            }
        } else {
            m_source_terms.emplace_back(
                    std::make_unique<sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(*m_geometry, *m_raytracer)
            );
        }
        #else
        m_source_terms.emplace_back(
                std::make_unique<sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(*m_geometry, *m_raytracer)
        );

        #endif


        m_los_source_terms.push_back(m_source_terms[m_source_terms.size() - 1].get());
    } else if (m_config.multiple_scatter_source() == sasktran2::Config::MultipleScatterSource::hr) {
        m_source_terms.emplace_back(
                std::make_unique<sasktran2::hr::DiffuseTable<NSTOKES>>(*m_raytracer, *m_geometry)
        );
        m_los_source_terms.push_back(m_source_terms[m_source_terms.size() - 1].get());


    }

    for(auto& source : m_source_terms) {
        source->initialize_config(m_config);
    }

}


template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_geometry() {
    // Trace every ray that we are given
    m_traced_rays.clear();
    m_traced_rays.resize(m_viewing_geometry.observer_rays().size());

    for(int i = 0; i < m_viewing_geometry.observer_rays().size(); ++i) {
        const auto& viewing_ray = m_viewing_geometry.observer_rays()[i];
        sasktran2::viewinggeometry::ViewingRay ray = viewing_ray->construct_ray(m_geometry->coordinates());

        m_raytracer->trace_ray(ray, m_traced_rays[i]);
    }

    // Initialize the integrator
    m_source_integrator->initialize_geometry(m_traced_rays, *m_geometry);

    for(auto& source : m_source_terms) {
        source->initialize_geometry(m_traced_rays);
    }

}

template<int NSTOKES>
void Sasktran2<NSTOKES>::validate_input_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) const {
    // Check that we have the required legendre information for the number of NSTOKES requested

    // Check that the atmosphere parameters have the correct dimensions
    
    // Check that the atmosphere geometry matches the global geometry
    if (atmosphere.num_wavel() != atmosphere.surface().albedo().size()) {
        BOOST_LOG_TRIVIAL(error) << "Atmosphere albedo does not have the correct dimensions";
    }
}

template <int NSTOKES>
void Sasktran2<NSTOKES>::calculate_radiance(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere, sasktran2::Output<NSTOKES>& output) const {
    omp_set_num_threads(m_config.num_threads());

    validate_input_atmosphere(atmosphere);

    // Use this method for observer geometries and make a different method for interior fluxes?

    // Initialize each source term with the atmosphere
    for(auto& source : m_source_terms) {
        source->initialize_atmosphere(atmosphere);
    }

    m_source_integrator->initialize_atmosphere(atmosphere);

    // Allocate memory, should be moved to thread storage?
    std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>> radiance(m_config.num_threads(),{NSTOKES, atmosphere.num_deriv(), true});

    output.resize((int)m_traced_rays.size(), atmosphere.num_wavel(), atmosphere.num_deriv());

    #pragma omp parallel for
    for(int w = 0; w < atmosphere.num_wavel(); ++w) {
        int thread_idx = omp_get_thread_num();

        // Trigger source term generation for this wavelength
        for(auto& source : m_source_terms) {
            source->calculate(w, thread_idx);
        }

        for(int i = 0; i < m_traced_rays.size(); ++i) {
            // Set the radiance thread storage to 0
            radiance[thread_idx].value.setZero();
            radiance[thread_idx].deriv.setZero();

            // Integrate all of the sources for the ray
            m_source_integrator->integrate(radiance[thread_idx], m_los_source_terms, w, i, thread_idx);

            // And assign it to the output
            output.assign(radiance[thread_idx], i, w);
        }

        // TODO: Is this where we should generate fluxes or other quantities that aren't through the integrator?
    }
}

template class Sasktran2<1>;
template class Sasktran2<3>;