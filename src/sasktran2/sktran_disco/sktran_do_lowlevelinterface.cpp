//
// Created by Daniel Zawada on 2022-06-29.
//

#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_lowlevelinterface.h"
#include "sktran_disco/sktran_do_specs.h"

namespace sasktran_disco_lowlevel {

    void calculate(const Atmosphere *atmosphere, const Config *config,
                   const WeightingFunctions *weightingfunctions,
                   const ViewingGeometry *geometry, const Output *output) {

        // Validate input parameters
        if(!config) {

        }

#ifdef SKTRAN_USE_MKL
        mkl_set_num_threads(1);
#endif

        if(config->nstokes == 1) {
        #ifdef SASKTRAN_DISCO_FULL_COMPILE
            if (config->nstr == 2) {
                internal_calculate<1, 2>(atmosphere, config, weightingfunctions, geometry, output);
            }
            else if (config->nstr == 4) {
                internal_calculate<1, 4>(atmosphere, config, weightingfunctions, geometry, output);
            }
            else if (config->nstr == 16) {
                internal_calculate<1, 16>(atmosphere, config, weightingfunctions, geometry, output);
            }
            else {
                internal_calculate<1, -1>(atmosphere, config, weightingfunctions, geometry, output);
            }
        #else
            internal_calculate<1, -1>(atmosphere, config, weightingfunctions, geometry, output);
        #endif

        } else if (config->nstokes == 3) {
            internal_calculate<3, -1>(atmosphere, config, weightingfunctions, geometry, output);
        }
    }

    template <int NSTOKES,  int CNSTR>
    void internal_calculate(const Atmosphere* atmosphere, const Config* config,
                            const WeightingFunctions* weightingfunctions,
                            const ViewingGeometry* geometry, const Output* output) {
        sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR> persistent_config;
        sasktran_disco::SKTRAN_DO_UserSpec userspec;
        persistent_config.configureLowLevel(userspec, *config, *geometry);

        sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR> geometry_layer(persistent_config, *atmosphere);

        std::vector<sasktran_disco::LineOfSight> los;

        sasktran_disco::VectorDim3<sasktran_disco::LegendrePhaseContainer<NSTOKES>> lp_coszen(geometry->nlos);		//< Legendre polynomials evaluated at LOS coszen
        los.resize(geometry->nlos);

        for(int i = 0; i < geometry->nlos; ++i) {
            los[i].azimuth = geometry->saa[i];
            los[i].coszenith = geometry->cos_vza[i];

            for (int m = 0; m < config->nstr; ++m) {
                lp_coszen[i].push_back(sasktran_disco::VectorDim1<sasktran_disco::LegendrePhaseContainer<NSTOKES>>(config->nstr));
                for (int l = 0; l < config->nstr; ++l) {
                    lp_coszen[i][m][l].fill(m, l, los[i].coszenith);
                }
            }
        }

        int numderiv = 0;
        if(weightingfunctions) {
            numderiv = weightingfunctions->numderiv;
        }

        int nthreads;
        if (config->nthreads > 0) {
            nthreads = config->nthreads;
        }
        else {
            nthreads = 1;
        }

        // Temporaries
        std::vector<sasktran_disco::Radiance<NSTOKES>> integral(nthreads, numderiv);
        std::vector<sasktran_disco::Radiance<NSTOKES>> m_component(nthreads, numderiv);

        std::vector<sasktran_disco::Radiance<NSTOKES>> exact_ss(nthreads, numderiv);
        std::vector<sasktran_disco::Radiance<NSTOKES>*> exact_ss_ptr(nthreads);

        int numazimuth = config->numazimuthexpansion > 0 ? config->numazimuthexpansion : config->nstr;
        if(numazimuth > config->nstr) {
            numazimuth = config->nstr;
        }

        std::vector<std::unique_ptr<sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>>> optical_layers;

        for(int i = 0; i < nthreads; ++i) {
            optical_layers.emplace_back(std::make_unique<sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>>(persistent_config, los, geometry_layer, persistent_config.pool().thread_data(i)));
        }

        #pragma omp parallel for schedule(guided) num_threads(nthreads)
        for(int i = 0; i < config->nwavel; ++i) {
            int thread_id = omp_get_thread_num();
            const auto& thread_data = persistent_config.pool().thread_data(thread_id);

            int wavoffset = NSTOKES * geometry->nlos;

            std::unique_ptr<sasktran_disco::BRDF_Base> brdf;

            brdf = std::make_unique<sasktran_disco::TestBRDF>(atmosphere->albedo[i]);

            auto& optical_layer = *optical_layers[thread_id];
            optical_layer.set_optical(i, std::move(brdf), *atmosphere, weightingfunctions);


            /*
            // TODO: Replace reconstruction of layerarray/RTE
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> optical_layer(persistent_config,
                                                                              i,
                                                                              los,
                                                                              std::move(brdf),
                                                                              *atmosphere,
                                                                              weightingfunctions,
                                                                              geometry_layer,
                                                                              thread_data);
            */


            sasktran_disco::RTESolver<NSTOKES, CNSTR> rte(persistent_config, optical_layer);

            for(int m = 0; m < numazimuth; ++m) {
                rte.solve(m);

                for (int j = 0; j < los.size(); ++j) {
                    double observeraltitude = geometry->viewingaltitude[los[j].unsorted_index];
                    double observer_opticaldepth = 0.0;

                    if(observeraltitude >= 0) {
                        // We might be inside the atmosphere, so we have to adjust the LOS optical depth
                        observer_opticaldepth = optical_layer.opticalDepthAt(observeraltitude);
                    }

                    // Start with assigning the radiance reflected from the ground
                    auto& ground_radiance = optical_layer.reflectedIntensity(m, los[j]);

                    m_component[thread_id].value = ground_radiance.value;
                    m_component[thread_id].deriv = ground_radiance.deriv;

                    // Calculate radiance recursively upwards through atmosphere
                    for (auto layer = optical_layer.template iteratorAcross<sasktran_disco::Propagating::UP>(); layer.isValid(); ++layer) {
                        if (layer.entryOpticalDepth() < observer_opticaldepth) {
                            // Layer doesn't contribute
                            continue;
                        }
                        double layerfraction = 1;
                        if(layer.exitOpticalDepth() < observer_opticaldepth) {
                            layerfraction = (layer.entryOpticalDepth() - observer_opticaldepth) / (layer.entryOpticalDepth() - layer.exitOpticalDepth());
                        }

                        // Attenuation through the layer is exp(-opticaldepth * layerfraction)
                        auto& dual_optical_depth = layer.layer().dual_thickness();

                        m_component[thread_id].value *= exp(-1.0*dual_optical_depth.value * layerfraction / los[j].coszenith);

                        m_component[thread_id].deriv *= exp(-1.0*dual_optical_depth.value * layerfraction / los[j].coszenith);
                        auto seq = Eigen::seq(dual_optical_depth.layer_start, dual_optical_depth.layer_start + dual_optical_depth.deriv.size() - 1);
                        if constexpr(NSTOKES == 1) {
                            m_component[thread_id].deriv(seq) += -layerfraction  / los[j].coszenith * m_component[thread_id].value * dual_optical_depth.deriv;
                        } else {
                            for (int k = 0; k < NSTOKES; ++k) {
                                m_component[thread_id].deriv(seq, k) +=
                                        -layerfraction  / los[j].coszenith * m_component[thread_id].value(k) * dual_optical_depth.deriv;
                            }
                        }


                        if(m == 0 && config->useexactsinglescatter) {
                            exact_ss[thread_id].setzero();

                            double lssa = atmosphere->ssa[layer.layer().index() + i*config->nlyr];
                            double lf = atmosphere->f[layer.layer().index() + i*config->nlyr];
                            double* lphase = &atmosphere->ss_phase[los[j].unsorted_index*NSTOKES + geometry->nlos * layer.layer().index() * NSTOKES + i*(geometry->nlos * config->nlyr * NSTOKES)];

                            if constexpr(NSTOKES == 1) {
                                exact_ss[thread_id].value = *lphase *
                                                 lssa / (1 - lf * lssa);
                            } else {
                                for(int k = 0; k < NSTOKES; k++) {
                                    exact_ss[thread_id].value(k) = lphase[k] *
                                                     lssa / (1 - lf * lssa);
                                }
                            }


                            // only have layer derivatives but we have a full container for some reason
                            // TODO: should be profiled and changed if necessary
                            int numlayerderiv = (int)optical_layer.inputDerivatives().numDerivativeLayer(layer.layer().index());
                            int layerstart = (int)optical_layer.inputDerivatives().layerStartIndex(layer.layer().index());
                            for(int k = 0; k < numlayerderiv; ++k) {
                                int derivindex = optical_layer.inputDerivatives().layerDerivatives()[k + layerstart].group_and_triangle_fraction[0].first;

                                double* dlphase = &weightingfunctions->d_ss_phase[derivindex * NSTOKES + i*(geometry->nlos * numderiv * NSTOKES) + los[j].unsorted_index * numderiv * NSTOKES];
                                double dlf = weightingfunctions->d_f[derivindex + i*numderiv];
                                double dlssa = weightingfunctions->d_ssa[derivindex + i*numderiv];

                                if constexpr(NSTOKES == 1) {
                                    exact_ss[thread_id].deriv(k + layerstart) += *dlphase * lssa / (1 - lf * lssa);
                                    exact_ss[thread_id].deriv(k + layerstart) += dlssa * *lphase / (1 - lf * lssa);
                                    exact_ss[thread_id].deriv(k + layerstart) += dlssa * lf * exact_ss[thread_id].value / (1 - lf * lssa);
                                    exact_ss[thread_id].deriv(k + layerstart) += dlf * lssa * exact_ss[thread_id].value / (1 - lf * lssa);
                                } else {
                                    for(int l = 0 ; l < NSTOKES; ++l) {
                                        exact_ss[thread_id].deriv(k + layerstart, l) += dlphase[l] * lssa / (1 - lf * lssa);
                                        exact_ss[thread_id].deriv(k + layerstart, l) += dlssa * lphase[l] / (1 - lf * lssa);
                                        exact_ss[thread_id].deriv(k + layerstart, l) += dlssa * lf * exact_ss[thread_id].value(l) / (1 - lf * lssa);
                                        exact_ss[thread_id].deriv(k + layerstart, l) += dlf * lssa * exact_ss[thread_id].value(l) / (1 - lf * lssa);
                                    }
                                }
                            }

                            exact_ss_ptr[thread_id] = &exact_ss[thread_id];
                        } else {
                            exact_ss_ptr[thread_id] = nullptr;
                        }

                        // Calculate and add in new source
                        layer.ptr()->integrate_source(m, los[j].coszenith, observer_opticaldepth, lp_coszen[j][m], integral[thread_id], exact_ss_ptr[thread_id], !config->useexactsinglescatter);

                        m_component[thread_id].value += integral[thread_id].value;
                        m_component[thread_id].deriv += integral[thread_id].deriv;
                    }

                    if constexpr (NSTOKES == 1) {
                        output->radiance[j + i * wavoffset] += m_component[thread_id].value * cos(m * los[j].azimuth);

                        for(int k = 0; k < numderiv; ++k) {
                            int derivgroup = optical_layer.inputDerivatives().layerDerivatives()[k].group_and_triangle_fraction[0].first;
                            output->d_radiance[(j + i * wavoffset)*numderiv + derivgroup] += m_component[thread_id].deriv(k) * cos(m * los[j].azimuth);
                        }
                    }
                    else {
                        for (int k = 0; k < NSTOKES; ++k) {
                            double azimuthal_factor;
                            if(k < 2) {
                                azimuthal_factor = cos(m * los[j].azimuth);
                            } else {
                                azimuthal_factor = sin(m * los[j].azimuth);
                            }
                            output->radiance[j * NSTOKES + k + i * wavoffset] += m_component[thread_id].value(k) * azimuthal_factor;

                            for(int l = 0; l < numderiv; ++l) {
                                int derivgroup = optical_layer.inputDerivatives().layerDerivatives()[l].group_and_triangle_fraction[0].first;

                                output->d_radiance[derivgroup + k*numderiv + j*numderiv*NSTOKES + i*numderiv*wavoffset] += m_component[thread_id].deriv(l, k) * azimuthal_factor;
                            }

                        }
                    }
                }
            }
        }
    }
}