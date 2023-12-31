#include "include/sktran_me_internals.h"

namespace sktran_me {
    template<int NSTOKES>
    void AtmosphereConstructor<NSTOKES>::construct_atmospheric_state(const sasktran2::Geometry1D& geometry,
                                                                     const sasktran2::Config& config,
                                                                     const sasktran2::viewinggeometry::ViewingGeometryContainer& viewing_geo,
                                                                     const GEODETIC_INSTANT& refpt,
                                                                     const std::vector<double>& wavelengths) {
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

        // Temporary storages
        std::vector<double> temp(numorder);
        std::vector<double> a1(numorder);
        std::vector<double> a2(numorder);
        std::vector<double> a3(numorder);
        std::vector<double> b1(numorder);


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

                        if(extxs[wavidx] == 0) {
                            // Doesn't matter anyways
                            species_entry.ssa(j, wavidx) = 0;
                        } else {
                            species_entry.ssa(j, wavidx) = scatxs[wavidx] / extxs[wavidx];
                        }
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
                                                                       &a1[0],
                                                                       &a2[0],
                                                                       &a3[0],
                                                                       &temp[0],
                                                                       &b1[0],
                                                                       &temp[0],
                                                                       numorder,
                                                                       numlegendre
                                );

                                for(int l = 0; l < numorder; ++l) {
                                    species_entry.phase[wavidx].storage()(l*4, j) = a1[l];
                                    species_entry.phase[wavidx].storage()(l*4 + 1, j) = a2[l];
                                    species_entry.phase[wavidx].storage()(l*4 + 2, j) = a3[l];
                                    species_entry.phase[wavidx].storage()(l*4 + 3, j) = b1[l];
                                }
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

        m_atmosphere = std::make_unique<sasktran2::atmosphere::Atmosphere<NSTOKES>>(std::move(storage), std::move(surface), m_wf_handler.calculating_wf());

        m_atmosphere->storage().total_extinction.setZero();
        m_atmosphere->storage().ssa.setZero();

        m_output = std::make_unique<sktran_me::OutputSKIF<NSTOKES>>(m_wf_handler, m_species_quantities, m_atmosphere->storage(), m_species_guid);

        // For each species add in the extinctions
        // Temporarily store scattering extinction in SSA
        for(int i = 0; i < m_species_quantities.size(); ++i) {
            m_atmosphere->storage().total_extinction.array() += m_species_quantities[i].total_extinction.array().colwise() * m_species_nd[i].array();
            m_atmosphere->storage().ssa.array() += m_species_quantities[i].ssa.array() * (m_species_quantities[i].total_extinction.array().colwise() * m_species_nd[i].array());

            // Now add in the phase information
            if(m_atmo_interface.species_type()[i] != AtmosphereInterface::PurelyAbsorbing) {
                // We have a scatterer
                for(int j = 0; j < wavelengths.size(); ++j) {
                    m_atmosphere->storage().phase[j].storage().array() += m_species_quantities[i].phase[j].storage().array().rowwise() * (m_species_quantities[i].ssa(Eigen::all, j).array() * m_species_quantities[i].total_extinction(Eigen::all, j).array() * m_species_nd[i].array()).transpose();
                }
            }
        }

        // Normalize phase by total scattering extinction
        for(int i = 0; i < wavelengths.size(); ++i) {
            m_atmosphere->storage().phase[i].storage().array().rowwise() /= m_atmosphere->storage().ssa(Eigen::all, i).array().transpose();
        }

        // Normalize SSA by total extinction
        m_atmosphere->storage().ssa.array() /= m_atmosphere->storage().total_extinction.array();

        // sktran_co wants extinction in /m
        m_atmosphere->storage().total_extinction.array() *= 100;

        // Set the species scattering derivatives
        m_wf_handler.set_scattering_information(m_atmo_interface, m_species_guid);

        for(int w = 0; w < wavelengths.size(); ++w) {
            m_atmosphere->storage().phase[w].resize_derivative(geometry.size(), numorder, m_wf_handler.num_scattering_derivatives(), m_atmosphere->scat_deriv_start_index());

            int scat_deriv_index = 0;
            for(int i = 0; i < m_wf_handler.num_wf_types(); ++i) {
                if(m_wf_handler.scattering_index(i) >= 0) {
                    // Have a scattering derivative
                    m_atmosphere->storage().phase[w].deriv_storage(m_wf_handler.scattering_index(i)) = m_species_quantities[m_wf_handler.species_index(i)].phase[w].storage() - m_atmosphere->storage().phase[w].storage();
                }
            }

        }

        if(config.apply_delta_scaling()) {
            m_atmosphere->apply_delta_m_scaling(config.num_do_streams());
        } else {
            m_atmosphere->storage().f.setZero();
        }

    }

    template<int NSTOKES>
    void AtmosphereConstructor<NSTOKES>::assign_output_radiance(const double **radiance, int *numwavelens,
                                                                int *numlinesofsight) {
        *radiance = &m_output->radiance()(0);

        *numwavelens = m_nwavel;
        *numlinesofsight = m_nlos * NSTOKES;
    }

    template<int NSTOKES>
    void AtmosphereConstructor<NSTOKES>::assign_wf_buffer(const double **wf) {
        auto* output = dynamic_cast<OutputSKIF<NSTOKES>*>(m_output.get());

        *wf = &output->wf()(0);
    }

    template class AtmosphereConstructor<1>;
    template class AtmosphereConstructor<3>;

}