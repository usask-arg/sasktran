#include "include/sktran_me_internals.h"

namespace sktran_me {
    template<int NSTOKES>
    void AtmosphereConstructor<NSTOKES>::construct_atmospheric_state(const sasktran2::Geometry1D& geometry,
                                                                     const sasktran2::Config& config,
                                                                     const sasktran2::viewinggeometry::ViewingGeometryContainer& viewing_geo,
                                                                     const GEODETIC_INSTANT& refpt,
                                                                     const std::vector<double>& wavelengths,
                                                                     std::unique_ptr<sasktran2::atmosphere::Atmosphere<NSTOKES>>& atmosphere) {
        // Make the output
        m_nwavel = wavelengths.size();
        m_nlos = viewing_geo.observer_rays().size();

        auto* atmospheric_state = m_atmo_interface.atmospheric_state();
        auto& all_species = m_atmo_interface.species();

        int numorder = config.num_singlescatter_moments();

        if(config.apply_delta_scaling() && numorder <= config.num_do_streams()) {
            // Have to add 1 legendre moment to get the scaling factor
            ++numorder;
        }

        m_species_quantities.resize(all_species.size(), sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>(wavelengths.size(), geometry.size(), numorder));

        m_species_nd.resize(all_species.size());
        m_species_guid.resize(all_species.size());
        for(auto& nd : m_species_nd) {
            nd.resize(geometry.altitude_grid().grid().size());
        }

        GEODETIC_INSTANT loc(refpt.latitude, refpt.longitude, refpt.heightm, refpt.mjd);
        atmospheric_state->UpdateCache(loc);
        bool haschanged;

        // Storage for each thread, Probably not even used anymore
        std::vector<int> thread_numlegendre(omp_get_num_threads());

        // Temporary storage
        std::vector<double> temp(numorder);

        // Temporary wavenumber grid
        std::vector<double> wavenumber(wavelengths.size());
        for(size_t idx = 0; idx < wavelengths.size(); ++idx) {
            wavenumber[idx] = 1e7 / wavelengths[idx];
        }

        // Temporary storage for the output cross sections
        std::vector<double> absxs, scatxs, extxs;
        absxs.resize(wavelengths.size());
        scatxs.resize(wavelengths.size());
        extxs.resize(wavelengths.size());

        // First fill the species specific tables
        for (size_t i = 0; i < all_species.size(); ++i) {
            auto& species = all_species[i];
            auto* clim = species.GetClimatology();
            auto* optprop = species.ParticleOpticalProps();
            auto& species_entry = m_species_quantities[i];
            m_species_guid[i] = species.GetSpecies();

            optprop->SetAtmosphericState(atmospheric_state);

            for (int j = 0; j < geometry.altitude_grid().grid().size(); ++j) {
                loc.heightm = geometry.altitude_grid().grid()(j);

                atmospheric_state->UpdateCache(loc);
                clim->UpdateCache(loc);

                clim->GetParameter(species.GetSpecies(), loc, &m_species_nd[i](j), false);

                if (optprop->IsHeightDependent() || j == 0) {
                    // If the optical property is height dependent we have to fill the table at every height
                    optprop->SetLocation(loc, &haschanged);
                    optprop->CalculateCrossSectionsArray(&wavenumber[0], wavenumber.size(), &absxs[0], &extxs[0], &scatxs[0]);

                    // I think this was causing some problems, CalculatePhaseMatrix is not threadsafe
                    //#pragma omp parallel for schedule(dynamic, 1)
                    for(int wavidx = 0; wavidx < wavelengths.size(); ++wavidx) {
                        size_t thread_id = omp_get_thread_num();
                        auto& numlegendre = thread_numlegendre[thread_id];

                        species_entry.ssa(j, wavidx) = scatxs[wavidx] / extxs[wavidx];
                        species_entry.total_extinction(j, wavidx) = extxs[wavidx];

                        if (optprop->IsScatterer())
                        {
                            if constexpr (NSTOKES == 1) {
                                optprop->LegendreCoefficientsP11(1e7 / wavelengths[wavidx],
                                                                 &species_entry.phase[wavidx].storage()(0, j),
                                                                 numorder,
                                                                 numlegendre);
                            } else if constexpr (NSTOKES == 3) {
                                optprop->LegendreCoefficientsPolarized(1e7 / wavelengths[wavidx],
                                                                       &species_entry.phase[wavidx].storage()(0, j),
                                                                       &species_entry.phase[wavidx].storage()(numorder, j),
                                                                       &species_entry.phase[wavidx].storage()(2*numorder, j),
                                                                       &temp[0],
                                                                       &species_entry.phase[wavidx].storage()(3*numorder, j),
                                                                       &temp[0],
                                                                       numorder,
                                                                       numlegendre
                                                                       );
                            }
                        }

                    }
                }
                else {
                    for(int wavidx = 0; wavidx < wavelengths.size(); ++wavidx) {
                        species_entry.phase[wavidx].storage()(Eigen::all, j) = species_entry.phase[wavidx].storage()(Eigen::all, 0);
                    }
                    species_entry.ssa(j, Eigen::all) = species_entry.ssa(0, Eigen::all);
                    species_entry.total_extinction(j, Eigen::all) = species_entry.total_extinction(0, Eigen::all);
                }
            }
        }

        // Next construct the combined quantities
        sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES> storage(wavelengths.size(), geometry.size(), numorder);
        sasktran2::atmosphere::Surface surface;

        surface.albedo().resize(wavelengths.size());

        // TODO: Non-lambertian BRDFS
        if(!m_atmo_interface.brdf()->IsLambertian()) {
            BOOST_LOG_TRIVIAL(error) << "sasktran_CO does not support non-lambertian BRDFs";
        }
        for(int i = 0; i < wavelengths.size(); ++i) {
            m_atmo_interface.brdf()->BRDF(wavelengths[i], refpt, 0,0,0, &surface.albedo()(i));
        }

        surface.albedo() *= M_PI;

        atmosphere = std::make_unique<sasktran2::atmosphere::Atmosphere<NSTOKES>>(std::move(storage), std::move(surface), m_wf_handler.calculating_wf());

        atmosphere->storage().total_extinction.setZero();
        atmosphere->storage().ssa.setZero();

        m_output = std::make_unique<sktran_me::OutputSKIF<NSTOKES>>(m_wf_handler, m_species_quantities, atmosphere->storage(), m_species_guid);


        // TODO: Species scattering derivatives

        // For each species add in the extinctions
        // Temporarily store scattering extinction in SSA
        for(int i = 0; i < m_species_quantities.size(); ++i) {
            atmosphere->storage().total_extinction.array() += m_species_quantities[i].total_extinction.array().colwise() * m_species_nd[i].array();
            atmosphere->storage().ssa.array() += m_species_quantities[i].ssa.array() * (m_species_quantities[i].total_extinction.array().colwise() * m_species_nd[i].array());

            // Now add in the phase information
            if(m_atmo_interface.species_type()[i] != AtmosphereInterface::PurelyAbsorbing) {
                // We have a scatterer
                for(int j = 0; j < wavelengths.size(); ++j) {
                    atmosphere->storage().phase[j].storage().array() += m_species_quantities[i].phase[j].storage().array().rowwise() * (m_species_quantities[i].ssa(Eigen::all, j).array() * m_species_quantities[i].total_extinction(Eigen::all, j).array() * m_species_nd[i].array()).transpose();
                }
            }
        }

        // Normalize phase by total scattering extinction
        for(int i = 0; i < wavelengths.size(); ++i) {
            atmosphere->storage().phase[i].storage().array().rowwise() /= atmosphere->storage().ssa(Eigen::all, i).array().transpose();
        }

        // Normalize SSA by total extinction
        atmosphere->storage().ssa.array() /= atmosphere->storage().total_extinction.array();

        // sktran_co wants extinction in /m
        atmosphere->storage().total_extinction.array() *= 100;

        if(config.apply_delta_scaling()) {
            atmosphere->apply_delta_m_scaling(config.num_do_streams());
        }

    }

    template<int NSTOKES>
    void AtmosphereConstructor<NSTOKES>::assign_output_radiance(const double **radiance, int *numwavelens,
                                                                int *numlinesofsight) {
        *radiance = &m_output->radiance()(0);

        *numwavelens = m_nwavel;
        *numlinesofsight = m_nlos;
    }

    template<int NSTOKES>
    void AtmosphereConstructor<NSTOKES>::assign_wf_buffer(const double **wf) {
        auto* output = dynamic_cast<OutputSKIF<NSTOKES>*>(m_output.get());

        *wf = &output->wf()(0);
    }

    template class AtmosphereConstructor<1>;
    template class AtmosphereConstructor<3>;

}