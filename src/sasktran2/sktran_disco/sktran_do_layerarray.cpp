#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_layerarray.h"

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::computeReflectedIntensities(AEOrder m, const sasktran_disco::LineOfSight& los) {
    // Here we evaluate the upwelling term at the surface.  Note that the sign convention switches compared to
    // the standard sign convention.  +<->- on W and Z

    m_reflection_computed[m][los.unsorted_index] = true;
    if(m_albedo[m].isLambertian() && m > 0) {
        return;
    }
    const uint N = this->M_NSTR / 2;

    const auto& layer = bottom();
    const auto& solution = layer.solution(m);
    const auto& input_deriv = m_input_derivatives;
    int layerStart = (int)input_deriv.layerStartIndex(layer.index());
    int numLayerDeriv = (int)input_deriv.numDerivativeLayer(layer.index());

    Radiance<NSTOKES> diffuse_contrib((int)input_deriv.numDerivative()), direct_contrib((int)input_deriv.numDerivative());
    diffuse_contrib.setzero();
    diffuse_contrib.setzero();

    direct_contrib.setzero();
    direct_contrib.setzero();

    auto& rho = m_albedo[m].losBDRFromStreams(los.unsorted_index);

    // Construct the Dual quantity for the layer transmittance
    Dual<double> beam_transmittance = layer.dual_beamTransmittance(Location::FLOOR, input_deriv);
    Dual<double> stream_transmittance;

    for(StreamIndex i = 0; i < N * NSTOKES; ++i) {
        LayerDual<double> dual_rho((uint)input_deriv.numDerivativeLayer(layer.index()), layer.index(), (uint)input_deriv.layerStartIndex(layer.index()));
        // TODO: something is definitely wrong here, but it might be fine for scalar surface reflection...
        int s1 = i % NSTOKES;
        if( i % NSTOKES != 0) {
            dual_rho.value = 0.0;
            dual_rho.deriv.setZero();
        } else {
            dual_rho.value = rho[i/NSTOKES + this->M_NSTR / 2];

            for (uint j = 0; j < input_deriv.numDerivativeLayer(layer.index()); ++j) {
                dual_rho.deriv(j) = input_deriv.layerDerivatives()[input_deriv.layerStartIndex(layer.index()) + j].d_albedo * kronDelta(m, 0);
            }
        }

        Radiance<NSTOKES> stream_contrib((int)input_deriv.numDerivative());
        if(solution.value.use_green_function()) {
            if constexpr(NSTOKES == 1) {
                stream_contrib.value = solution.value.dual_Gplus_bottom().value(i);
            } else {
                stream_contrib.value(s1) = solution.value.dual_Gplus_bottom().value(i);
            }

            stream_contrib.deriv(Eigen::all, s1) = solution.value.dual_Gplus_bottom().deriv(Eigen::all, i);
        } else {
            if constexpr(NSTOKES == 1) {
                stream_contrib.value = solution.value.dual_particular_plus().value(i) * beam_transmittance.value;
            } else {
                stream_contrib.value(s1) = solution.value.dual_particular_plus().value(i) * beam_transmittance.value;
            }
            stream_contrib.deriv(Eigen::all, s1) = solution.value.dual_particular_plus().value(i) * beam_transmittance.deriv +
                                   beam_transmittance.value * solution.value.dual_particular_plus().deriv(Eigen::all, i);
        }
        // Positive homogeneous solutions
        for(uint j = 0; j < N * NSTOKES; ++j) {
            stream_transmittance = layer.dual_streamTransmittance(Location::INSIDE, m, j, input_deriv);
            uint homogIndex = j*N*NSTOKES + i;

            if constexpr(NSTOKES == 1) {
                stream_contrib.value += solution.boundary.L_coeffs.value(j) * solution.value.dual_homog_plus().value(homogIndex) * stream_transmittance.value;
            } else {
                stream_contrib.value(s1) += solution.boundary.L_coeffs.value(j) * solution.value.dual_homog_plus().value(homogIndex) * stream_transmittance.value;
            }

            // LCoeffs have full derivatives, for some reason stream transmittance is a full dual even though it should be layer?
            for(uint k = 0; k < input_deriv.numDerivative(); ++k) {
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
            for(uint k = 0; k < input_deriv.numDerivative(); ++k) {
                stream_contrib.deriv(k, s1) += solution.boundary.M_coeffs.deriv(k, j) * solution.value.dual_homog_minus().value(homogIndex);
            }
            // Homog only have layer derivs
            for(int k = 0; k < numLayerDeriv; ++k) {
                stream_contrib.deriv(k + layerStart, s1) += solution.boundary.M_coeffs.value(j) * solution.value.dual_homog_minus().deriv(k, homogIndex);
            }
        }
        // Add stream i contribution
        double factor = (1.0 + kronDelta(m, 0)) * (*this->M_MU)[i/NSTOKES] * (*this->M_WT)[i/NSTOKES];

        if constexpr(NSTOKES == 1) {
            diffuse_contrib.value += factor * stream_contrib.value * dual_rho.value;
        } else {
            diffuse_contrib.value(s1) += factor * stream_contrib.value(s1) * dual_rho.value;
        }

        // stream_contrib has full derivatives
        for(uint k = 0; k < input_deriv.numDerivative(); ++k) {
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
    LayerDual<double> dual_albedo_sun((uint)input_deriv.numDerivativeLayer(layer.index()), layer.index(), (uint)input_deriv.layerStartIndex(layer.index()));
    dual_albedo_sun.value = m_albedo[m].losBDRFromSun(los.unsorted_index);
    for (uint j = 0; j < input_deriv.numDerivativeLayer(layer.index()); ++j) {
        dual_albedo_sun.deriv(j) = input_deriv.layerDerivatives()[input_deriv.layerStartIndex(layer.index()) + j].d_albedo * kronDelta(m, 0);
    }

    if (m_include_direct_bounce) {
        if constexpr(NSTOKES == 1) {
            direct_contrib.value = directIntensityTOA() * this->M_CSZ / PI * beam_transmittance.value * dual_albedo_sun.value;
        } else {
            direct_contrib.value(0) = directIntensityTOA() * this->M_CSZ / PI * beam_transmittance.value * dual_albedo_sun.value;
        }

        // beam has full deriv
        for(uint k = 0; k < input_deriv.numDerivative(); ++k) {
            direct_contrib.deriv(k, 0) += directIntensityTOA() * this->M_CSZ / PI * beam_transmittance.deriv(k) * dual_albedo_sun.value;
        }
        // albedo has layer deriv
        for(int k = 0; k < numLayerDeriv; ++k) {
            direct_contrib.deriv(k + layerStart, 0) += directIntensityTOA() * this->M_CSZ / PI * beam_transmittance.value * dual_albedo_sun.deriv(k);
        }

        if (m_config.ss_only()) {
            m_ground_reflection[m][los.unsorted_index].value = direct_contrib.value;
            m_ground_reflection[m][los.unsorted_index].deriv = direct_contrib.deriv;
        }
        else {
            m_ground_reflection[m][los.unsorted_index].value = direct_contrib.value + diffuse_contrib.value;
            m_ground_reflection[m][los.unsorted_index].deriv = direct_contrib.deriv + diffuse_contrib.deriv;
        }
    }
    else {
        m_ground_reflection[m][los.unsorted_index].value = diffuse_contrib.value;
        m_ground_reflection[m][los.unsorted_index].deriv = diffuse_contrib.deriv;
    }

    if (m == 0) {
        if constexpr(NSTOKES == 1) {
            m_ground_reflection[m][los.unsorted_index].value += m_surfaceemission;
        } else {
            m_ground_reflection[m][los.unsorted_index].value(0) += m_surfaceemission;
        }
    }
}

template <int NSTOKES, int CNSTR>
sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::OpticalLayerArray(
        const PersistentConfiguration<NSTOKES, CNSTR> &config, int wavelidx,
        const std::vector<LineOfSight>& los,
        std::unique_ptr<BRDF_Base> brdf,
        const GeometryLayerArray<NSTOKES, CNSTR> &geometry_layers,
        const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere,
        const sasktran2::Config& sk_config
        ) : OpticalLayerArrayROP<NSTOKES>(config),
        m_direct_toa(this->M_SOLAR_DIRECT_INTENSITY),
        m_config(config),
        m_include_direct_bounce(false),
        m_surfaceemission(0.0),
        m_num_los(los.size()),
        m_input_derivatives(config.pool().thread_data().input_derivatives()),
        m_albedo(los, *this->M_MU, this->M_CSZ, std::move(brdf), config.userSpec()? config.userSpec()->getNumBRDFQuadratureTerms() : 64)
{
    m_wavel_index = wavelidx;
    // Allocations
    m_layers.reserve(this->M_NLYR);
    m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
    m_chapman_factors.setZero();

    // Accumulation quantities for the layers
    double ceiling_depth = 0;
    double floor_depth = 0;

    for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double ceil_h = geometry_layers.layer_ceiling()(p);
        double floor_h = geometry_layers.layer_floor()(p);

        int top_atmosphere_idx = (int)atmosphere.storage().total_extinction.rows() - p - 1;
        int bot_atmosphere_idx = top_atmosphere_idx - 1;

        double kbot = atmosphere.storage().total_extinction(bot_atmosphere_idx, wavelidx);
        double ktop = atmosphere.storage().total_extinction(top_atmosphere_idx, wavelidx);

        double layer_dh = ceil_h - floor_h;
        double ssa_bot = atmosphere.storage().ssa(bot_atmosphere_idx, wavelidx);
        double ssa_top = atmosphere.storage().ssa(top_atmosphere_idx, wavelidx);

        double f_bot = atmosphere.storage().f(bot_atmosphere_idx, wavelidx);
        double f_top = atmosphere.storage().f(top_atmosphere_idx, wavelidx);

        double od = (kbot + ktop) / 2 * layer_dh;
        double ssa = (ssa_bot * kbot + ssa_top * ktop) / (kbot + ktop);

        double f = (ssa_bot * kbot*f_bot + ssa_top * ktop*f_top) / (kbot*ssa_bot + ktop*ssa_top);

        // Copy the legendre coefficients to a new vector
        std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>> lephasef(new VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);

        const auto& phase = atmosphere.storage().phase[wavelidx].storage();

        for (uint k = 0; k < this->M_NSTR; ++k) {
            auto& temp = (*lephasef)[k];

            if constexpr (NSTOKES == 1) {
                double avg_p = (kbot*ssa_bot*(phase(k, bot_atmosphere_idx) - f_bot * (2*k+1) / (1-f_bot)) +
                        ktop*ssa_top*(phase(k, top_atmosphere_idx) - f_top * (2*k+1) / (1-f_top))) / (kbot*ssa_bot + ktop*ssa_top);

                temp.a1 = avg_p;
            } else if constexpr (NSTOKES == 3) {
                auto stokes_seq = Eigen::seq(k*4, (k+1)*4 - 1);
                Eigen::Vector<sasktran2::types::leg_coeff, 4> avg_p = (kbot*ssa_bot*phase(stokes_seq, bot_atmosphere_idx) + ktop*ssa_top*phase(stokes_seq, top_atmosphere_idx)) / (kbot*ssa_bot + ktop*ssa_top);
                temp.a1 = avg_p(0) - f * (2*k+1) / (1-f);
                temp.a2 = avg_p(1) - f * (2*k+1) / (1-f);
                temp.a3 = avg_p(2) - f * (2*k+1) / (1-f);
                temp.b1 = -(avg_p(3) - f * (2*k+1) / (1-f)); // TODO: Check
            } else {

            }
        }

        floor_depth += od;

        double total_ext = od / (ceil_h - floor_h);
        double scat_ext = total_ext * ssa;
        scat_ext = std::max(scat_ext, total_ext * this->m_userspec->getSSAEqual1Dither());

        m_layers.push_back(
                std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(new OpticalLayer<NSTOKES, CNSTR>(config, p, scat_ext, total_ext, std::move(lephasef), ceiling_depth, floor_depth, ceil_h, floor_h, m_input_derivatives))
        );

        ceiling_depth = floor_depth;
    }
    m_chapman_factors = geometry_layers.chapman_factors();

    if (atmosphere.num_deriv() > 0 && sk_config.wf_precision() == sasktran2::Config::WeightingFunctionPrecision::full) {
        int numderiv = atmosphere.num_deriv();
        int num_atmo_grid = (int)atmosphere.storage().total_extinction.rows();
        int num_scattering_groups = atmosphere.num_scattering_deriv_groups();

        if (!m_input_derivatives.is_geometry_configured()) {
            // Construct the matrix that maps atmosphere extinction to layer quantities
            // TODO: In the future we will construct this matrix elsewhere, probably in geometry layers, and then use it
            // to construct both the atmosphere and the derivatives
            Eigen::MatrixXd atmosphere_mapping(this->M_NLYR, num_atmo_grid);
            atmosphere_mapping.setZero();

            for(LayerIndex p = 0; p < this->M_NLYR; ++p) {
                int top_atmosphere_idx = (int)atmosphere.storage().total_extinction.rows() - p - 1;
                int bot_atmosphere_idx = top_atmosphere_idx - 1;

                atmosphere_mapping(p, top_atmosphere_idx) = 0.5;
                atmosphere_mapping(p, bot_atmosphere_idx) = 0.5;
            }

            // Now avg_k in layers = atmosphere_mapping @ atmosphere_extinction
            for(LayerIndex p = 0; p < this->M_NLYR; ++p) {
                auto& ptrb_layer = layer(p);

                // Go through the atmosphere mapping and add non-zero elements for SSA/ext
                for(int scatidx = 0; scatidx < num_scattering_groups; ++scatidx) {
                    LayerInputDerivative<NSTOKES>& deriv_scat = m_input_derivatives.addDerivative(this->M_NSTR, p);

                    for(int i = 0; i < num_atmo_grid; ++i) {
                        if(atmosphere_mapping(p, i) > 0) {
                            // How d atmosphere legendre influences layer legendre
                            deriv_scat.group_and_triangle_fraction.emplace_back(atmosphere.scat_deriv_start_index() + num_atmo_grid*scatidx + i, atmosphere_mapping(p, i));
                            deriv_scat.extinctions.emplace_back(1);

                            // How d atmosphere ssa influences layer legendre?
                            //deriv_scat.group_and_triangle_fraction.emplace_back(atmosphere.ssa_deriv_start_index() + i, atmosphere_mapping(p, i));
                            //deriv_scat.extinctions.emplace_back(1);

                            // How d atmosphere k influences layer legendre?
                            //deriv_scat.group_and_triangle_fraction.emplace_back(i, atmosphere_mapping(p, i));
                            //deriv_scat.extinctions.emplace_back(1);
                        }
                    }
                }

                LayerInputDerivative<NSTOKES>& deriv_ext = m_input_derivatives.addDerivative(this->M_NSTR, p);

                deriv_ext.d_optical_depth = 1;
                // Go through the atmosphere mapping and add non-zero elements for SSA/ext
                for(int i = 0; i < num_atmo_grid; ++i) {
                    if(atmosphere_mapping(p, i) > 0) {

                        // How d atmosphere k influences layer od
                        deriv_ext.group_and_triangle_fraction.emplace_back(i, atmosphere_mapping(p, i) * (ptrb_layer.altitude(Location::CEILING) - ptrb_layer.altitude(Location::FLOOR)));
                        deriv_ext.extinctions.emplace_back(1);
                    }
                }

                LayerInputDerivative<NSTOKES>& deriv_ssa = m_input_derivatives.addDerivative(this->M_NSTR, p);
                deriv_ssa.d_SSA = 1;
                // Go through the atmosphere mapping and add non-zero elements for SSA/ext
                for(int i = 0; i < num_atmo_grid; ++i) {
                    if(atmosphere_mapping(p, i) > 0) {
                        // How d atmosphere ssa influences layer ssa
                        deriv_ssa.group_and_triangle_fraction.emplace_back(i + num_atmo_grid, atmosphere_mapping(p, i));
                        deriv_ssa.extinctions.emplace_back(1);

                        // How d atmosphere k influences layer ssa
                        deriv_ssa.group_and_triangle_fraction.emplace_back(i, atmosphere_mapping(p, i));
                        deriv_ssa.extinctions.emplace_back(1);
                    }
                }

            }
            // Add one final derivative for the surface
            LayerInputDerivative<NSTOKES>& deriv_albedo = m_input_derivatives.addDerivative(this->M_NSTR, this->M_NLYR - 1);
            deriv_albedo.d_albedo = 1;
            deriv_albedo.group_and_triangle_fraction.emplace_back(atmosphere.surface_deriv_start_index(), 1);
            deriv_albedo.extinctions.emplace_back(1);

            m_input_derivatives.set_geometry_configured();
            m_input_derivatives.sort(this->M_NLYR);
        }

        // Go through the derivatives and reassign a few things
        for(int k = 0; k < m_input_derivatives.numDerivative(); ++k) {
            LayerInputDerivative<NSTOKES>& deriv = m_input_derivatives.layerDerivatives()[k];
            const OpticalLayer<NSTOKES, CNSTR>* layer = m_layers[deriv.layer_index].get();

            // Have to check what kind of derivative this is
            if(deriv.d_optical_depth > 0) {
                // OD derivative, don't have to do anything here
            } else if (deriv.d_SSA > 0) {
                // SSA derivative

                for(int l = 0; l < deriv.group_and_triangle_fraction.size(); ++l) {
                    auto& group = deriv.group_and_triangle_fraction[l];
                    auto& extinction = deriv.extinctions[l];

                    if(group.first >= num_atmo_grid) {
                        // Layer SSA contribution to dI/d atmosphere SSA, these are equal to the relative fraction of
                        // extinction contributions

                        int atmo_index = group.first - num_atmo_grid;

                        extinction = atmosphere.storage().total_extinction(atmo_index, wavelidx) / layer->totalExt();
                    } else {
                        // This is layer SSA contribution to dI/d atmosphere k,
                        int atmo_index = group.first;

                        extinction = (atmosphere.storage().ssa(atmo_index, wavelidx) - layer->ssa()) / layer->totalExt();
                    }
                }

            } else if (deriv.d_albedo == 0) {
                // Scattering derivative
                for(int l = 0; l < deriv.group_and_triangle_fraction.size(); ++l) {
                    auto& group = deriv.group_and_triangle_fraction[l];
                    auto& extinction = deriv.extinctions[l];

                    if(group.first >= 2*num_atmo_grid) {
                        // Layer legendre contribution to dI / d atmosphere legendre
                        int group_index = (group.first - 2*num_atmo_grid) / int(num_atmo_grid);
                        int atmo_index = group.first % num_atmo_grid;

                        double f = atmosphere.storage().f(atmo_index, wavelidx);

                        for(int l = 0; l < (int)this->M_NSTR; ++l) {
                            deriv.d_legendre_coeff[l].a1 = atmosphere.storage().phase[wavelidx].deriv_storage(
                                    group_index)(l, atmo_index);

                            if(atmosphere.storage().applied_f_order > 0) {
                                deriv.d_legendre_coeff[l].a1 += -(2*l+1) / (1-f) / (1-f) * atmosphere.storage().phase[wavelidx].f_deriv_storage(group_index)(atmo_index);
                            }
                        }

                        extinction = atmosphere.storage().ssa(atmo_index, wavelidx) * atmosphere.storage().total_extinction(atmo_index, wavelidx) / layer->scatExt();
                    } else if(group.first >= num_atmo_grid) {
                        // Layer legendre contribution to dI / dssa ?
                    } else {
                        // Layer legendre contribution to dI / dk ?
                    }
                }
            }
        }
    }

    // Post configure the layers, PS beam transmittances and derivative calculations
    for (auto& layer : m_layers) {
        // Start by telling each layer to calculate the derivative of thickness
        layer->configureDerivative();
    }

    configureTransmission();

    // Configure azimuthal dependencies
    for (auto& layer : m_layers) {
        registerAzimuthDependency(*layer);
    }
    registerAzimuthDependency(m_albedo);

    // resize ground reflection vector
    m_ground_reflection.resize(this->M_NSTR, std::vector<Radiance<NSTOKES>>(
            static_cast<uint>(los.size()),
            Radiance<NSTOKES>((int)m_input_derivatives.numDerivative())
    ));
    m_reflection_computed.resize(this->M_NSTR, std::vector<bool>(los.size(), false));

    // TODO: Figure out surface emission in low level interface?

    m_surfaceemission = 0.0;
}

template <int NSTOKES, int CNSTR>
sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::OpticalLayerArray(const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>& config,
                                                                     const std::vector<LineOfSight>& los,
                                                                     const sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
                                                                     const ThreadData<NSTOKES, CNSTR>& thread_data
) : OpticalLayerArrayROP<NSTOKES>(config),
    m_direct_toa(this->M_SOLAR_DIRECT_INTENSITY),
    m_config(config),
    m_include_direct_bounce(true),
    m_surfaceemission(0.0),
    m_num_los(los.size()),
    m_input_derivatives(thread_data.input_derivatives()),
    m_albedo(los, *this->M_MU, this->M_CSZ, nullptr, config.userSpec()? config.userSpec()->getNumBRDFQuadratureTerms() : 64){

    // Allocations
    m_layers.reserve(this->M_NLYR);

    for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double ceil_h = geometry_layers.layer_ceiling()(p);
        double floor_h = geometry_layers.layer_floor()(p);

        m_layers.push_back(
                std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(new OpticalLayer<NSTOKES, CNSTR>(config, p, ceil_h, floor_h, m_input_derivatives, thread_data))
        );
    }
    m_chapman_factors = geometry_layers.chapman_factors();

    // Configure azimuthal dependencies
    for (auto& layer : m_layers) {
        registerAzimuthDependency(*layer);
    }
    registerAzimuthDependency(m_albedo);

    // TODO: Figure out surface emission in low level interface?

    m_surfaceemission = 0.0;
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::set_optical(int wavelidx,
                                                                    std::unique_ptr<BRDF_Base> brdf,
                                                                    const sasktran_disco_lowlevel::Atmosphere& atmosphere,
                                                                    const sasktran_disco_lowlevel::WeightingFunctions* weightingfunctions) {
    m_albedo.injectTestingBRDF(std::move(brdf));
    m_albedo.reset();

    // Map the low level atmosphere memory to eigen matrices
    int stream_independent_offset = this->M_NLYR * wavelidx;
    int stream_dependent_offset = this->M_NLYR * this->M_NSTR * wavelidx;

    Eigen::Map<Eigen::VectorXd> layer_od(atmosphere.od + stream_independent_offset, this->M_NLYR);
    Eigen::Map<Eigen::VectorXd> layer_ssa(atmosphere.ssa + stream_independent_offset, this->M_NLYR);
    Eigen::Map<Eigen::VectorXd> layer_f(atmosphere.f + stream_independent_offset, this->M_NLYR);

    Eigen::Map<Eigen::MatrixXd> layer_a1(atmosphere.a1 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);

    // Can't if constexpr these because it breaks scope
    Eigen::Map<Eigen::MatrixXd> layer_a2(atmosphere.a2 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);
    Eigen::Map<Eigen::MatrixXd> layer_a3(atmosphere.a3 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);
    Eigen::Map<Eigen::MatrixXd> layer_b1(atmosphere.b1 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);

    Eigen::Map<Eigen::MatrixXd> layer_a4(atmosphere.a4 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);
    Eigen::Map<Eigen::MatrixXd> layer_b2(atmosphere.b2 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);

    // Accumulation quantities for the layers
    double ceiling_depth = 0;
    double floor_depth = 0;

    for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double od = layer_od(p);
        double ssa = layer_ssa(p);
        double f = layer_f(p);

        od *= 1 - ssa * f;
        ssa *= (1 - f) / (1 - ssa*f);

        // Copy the legendre coefficients to a new vector
        std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>> lephasef(new VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);

        for (uint k = 0; k < this->M_NSTR; ++k) {
            auto& temp = (*lephasef)[k];
            temp.a1 = (layer_a1(k, p) - f * (2*k+1)) / (1-f);

            if constexpr (NSTOKES > 2) {
                temp.a2 = layer_a2(k, p);
                temp.a3 = layer_a3(k, p);
                temp.b1 = layer_b1(k, p);
            }
            if constexpr (NSTOKES > 3) {
                temp.a4 = layer_a4(k, p);
                temp.b2 = layer_b2(k, p);
            }
        }

        floor_depth += od;

        double total_ext = od / (m_layers[p]->altitude(Location::CEILING) - m_layers[p]->altitude(Location::FLOOR));
        double scat_ext = total_ext * ssa;
        scat_ext = std::max(scat_ext, total_ext * this->m_userspec->getSSAEqual1Dither());

        m_layers[p]->set_optical(scat_ext, total_ext, *lephasef, ceiling_depth, floor_depth);

        ceiling_depth = floor_depth;
    }

    if (weightingfunctions != nullptr) {
        int numderiv = weightingfunctions->numderiv;
        if (!m_input_derivatives.is_geometry_configured()) {
            for (int k = 0; k < weightingfunctions->numderiv; ++k) {
                LayerInputDerivative<NSTOKES>& deriv = m_input_derivatives.addDerivative(this->M_NSTR,
                                                                                         weightingfunctions->d_layerindex[k]);
                deriv.group_and_triangle_fraction.emplace_back(std::pair<uint, double>(k, 1));
                deriv.extinctions.emplace_back(1);
            }
            m_input_derivatives.sort(this->M_NLYR);
            m_input_derivatives.set_geometry_configured();
        }

        for (int k = 0; k < numderiv; ++k) {
            LayerInputDerivative<NSTOKES>& deriv = m_input_derivatives.layerDerivatives()[k];

            int derivindex = deriv.group_and_triangle_fraction[0].first;
            int index = derivindex + wavelidx * numderiv;

            int layerindex = weightingfunctions->d_layerindex[derivindex];

            deriv.d_albedo = weightingfunctions->d_albedo[index];
            deriv.d_optical_depth = weightingfunctions->d_od[index] * (1 - layer_ssa(layerindex) * layer_f(layerindex)) - layer_od(layerindex) * layer_f(layerindex) * weightingfunctions->d_ssa[index] - layer_ssa(layerindex) * layer_od(layerindex) * weightingfunctions->d_f[index];
            deriv.d_SSA = weightingfunctions->d_ssa[index] * (1 - layer_f(layerindex)) / (1 - layer_ssa(layerindex) * layer_f(layerindex)) -
                          weightingfunctions->d_f[index] * layer_ssa(layerindex) / (1 - layer_ssa(layerindex) * layer_f(layerindex)) +
                          (weightingfunctions->d_ssa[index] * layer_f(layerindex) + weightingfunctions->d_f[index] * layer_ssa(layerindex)) *
                          layer_ssa(layerindex) * (1 - layer_f(layerindex)) / ((1 - layer_ssa(layerindex) * layer_f(layerindex)) * (1 - layer_ssa(layerindex) * layer_f(layerindex)));

            for(int l = 0; l < (int)this->M_NSTR; ++l) {
                int coeffindex = derivindex + wavelidx * numderiv * this->M_NSTR + l*numderiv;
                deriv.d_legendre_coeff[l].a1 = weightingfunctions->d_a1[coeffindex] / (1 - layer_f(layerindex)) -
                                               weightingfunctions->d_f[index] * (2*l+1) / (1 - layer_f(layerindex)) +
                                               weightingfunctions->d_f[index] * (layer_a1(l, layerindex) - layer_f(layerindex) * (2*l+1)) / (1 - layer_f(layerindex)) / (1 - layer_f(layerindex));

                if constexpr (NSTOKES > 2) {
                    deriv.d_legendre_coeff[l].a2 = weightingfunctions->d_a2[coeffindex];
                    deriv.d_legendre_coeff[l].a3 = weightingfunctions->d_a3[coeffindex];
                    deriv.d_legendre_coeff[l].b1 = weightingfunctions->d_b1[coeffindex];
                }
                if constexpr (NSTOKES > 3) {
                    deriv.d_legendre_coeff[l].a4 = weightingfunctions->d_a4[coeffindex];
                    deriv.d_legendre_coeff[l].b2 = weightingfunctions->d_b2[coeffindex];
                }

            }

        }
    }

    // Post configure the layers, PS beam transmittances and derivative calculations
    for (auto& layer : m_layers) {
        // Start by telling each layer to calculate the derivative of thickness
        layer->configureDerivative();
    }

    configureTransmission();

    if(m_ground_reflection.size() == 0) {
        // resize ground reflection vector
        m_ground_reflection.resize(this->M_NSTR, std::vector<Radiance<NSTOKES>>(
                static_cast<uint>(m_num_los),
                Radiance<NSTOKES>((int)m_input_derivatives.numDerivative())
        ));
        m_reflection_computed.resize(this->M_NSTR, std::vector<bool>(m_num_los, false));
    } else {
        for(auto& ele : m_reflection_computed) {
            std::fill(ele.begin(), ele.end(), false);
        }
    }
}

template <int NSTOKES, int CNSTR>
sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::OpticalLayerArray(const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>& config,
                                                                     int wavelidx,
                                                                     const std::vector<LineOfSight>& los,
                                                                     std::unique_ptr<BRDF_Base> brdf,
                                                                     const sasktran_disco_lowlevel::Atmosphere& atmosphere,
                                                                     const sasktran_disco_lowlevel::WeightingFunctions* weightingfunctions,
                                                                     const sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
                                                                     const ThreadData<NSTOKES, CNSTR>& thread_data
) : OpticalLayerArray<NSTOKES, CNSTR>(config, los, geometry_layers, thread_data) {

    set_optical(wavelidx, std::move(brdf), atmosphere, weightingfunctions);
    return;

    // Allocations
    m_layers.reserve(this->M_NLYR);
    m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
    m_chapman_factors.setZero();

    // Map the low level atmosphere memory to eigen matrices
    int stream_independent_offset = this->M_NLYR * wavelidx;
    int stream_dependent_offset = this->M_NLYR * this->M_NSTR * wavelidx;

    Eigen::Map<Eigen::VectorXd> layer_od(atmosphere.od + stream_independent_offset, this->M_NLYR);
    Eigen::Map<Eigen::VectorXd> layer_ssa(atmosphere.ssa + stream_independent_offset, this->M_NLYR);
    Eigen::Map<Eigen::VectorXd> layer_f(atmosphere.f + stream_independent_offset, this->M_NLYR);

    Eigen::Map<Eigen::MatrixXd> layer_a1(atmosphere.a1 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);

    // Can't if constexpr these because it breaks scope
    Eigen::Map<Eigen::MatrixXd> layer_a2(atmosphere.a2 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);
    Eigen::Map<Eigen::MatrixXd> layer_a3(atmosphere.a3 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);
    Eigen::Map<Eigen::MatrixXd> layer_b1(atmosphere.b1 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);

    Eigen::Map<Eigen::MatrixXd> layer_a4(atmosphere.a4 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);
    Eigen::Map<Eigen::MatrixXd> layer_b2(atmosphere.b2 + stream_dependent_offset, this->M_NSTR, this->M_NLYR);

    // Accumulation quantities for the layers
    double ceiling_depth = 0;
    double floor_depth = 0;

    for (LayerIndex p = 0; p < this->M_NLYR; ++p) {
        double ceil_h = geometry_layers.layer_ceiling()(p);
        double floor_h = geometry_layers.layer_floor()(p);

        double od = layer_od(p);
        double ssa = layer_ssa(p);
        double f = layer_f(p);

        od *= 1 - ssa * f;
        ssa *= (1 - f) / (1 - ssa*f);

        // Copy the legendre coefficients to a new vector
        std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>> lephasef(new VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);

        for (uint k = 0; k < this->M_NSTR; ++k) {
            auto& temp = (*lephasef)[k];
            temp.a1 = (layer_a1(k, p) - f * (2*k+1)) / (1-f);

            if constexpr (NSTOKES > 2) {
                temp.a2 = layer_a2(k, p);
                temp.a3 = layer_a3(k, p);
                temp.b1 = layer_b1(k, p);
            }
            if constexpr (NSTOKES > 3) {
                temp.a4 = layer_a4(k, p);
                temp.b2 = layer_b2(k, p);
            }
        }

        floor_depth += od;

        double total_ext = od / (ceil_h - floor_h);
        double scat_ext = total_ext * ssa;
        scat_ext = std::max(scat_ext, total_ext * this->m_userspec->getSSAEqual1Dither());

        //m_layers.push_back(
        //    std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(new OpticalLayer<NSTOKES, CNSTR>(config, p, scat_ext, total_ext, std::move(lephasef), ceiling_depth, floor_depth, ceil_h, floor_h, m_input_derivatives))
        //);

        m_layers.push_back(
                std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(new OpticalLayer<NSTOKES, CNSTR>(config, p, ceil_h, floor_h, m_input_derivatives, thread_data))
        );

        m_layers[m_layers.size() - 1]->set_optical(scat_ext, total_ext, *lephasef, ceiling_depth, floor_depth);

        ceiling_depth = floor_depth;
    }
    m_chapman_factors = geometry_layers.chapman_factors();

    if (weightingfunctions != nullptr) {
        int numderiv = weightingfunctions->numderiv;
        if (!m_input_derivatives.is_geometry_configured()) {
            for (int k = 0; k < weightingfunctions->numderiv; ++k) {
                LayerInputDerivative<NSTOKES>& deriv = m_input_derivatives.addDerivative(this->M_NSTR,
                    weightingfunctions->d_layerindex[k]);
                deriv.group_and_triangle_fraction.emplace_back(std::pair<uint, double>(k, 1));
                deriv.extinctions.emplace_back(1);
            }
            m_input_derivatives.sort(this->M_NLYR);
            m_input_derivatives.set_geometry_configured();
        }

        for (int k = 0; k < numderiv; ++k) {
            LayerInputDerivative<NSTOKES>& deriv = m_input_derivatives.layerDerivatives()[k];

            int derivindex = deriv.group_and_triangle_fraction[0].first;
            int index = derivindex + wavelidx * numderiv;

            int layerindex = weightingfunctions->d_layerindex[derivindex];

            deriv.d_albedo = weightingfunctions->d_albedo[index];
            deriv.d_optical_depth = weightingfunctions->d_od[index] * (1 - layer_ssa(layerindex) * layer_f(layerindex)) - layer_od(layerindex) * layer_f(layerindex) * weightingfunctions->d_ssa[index] - layer_ssa(layerindex) * layer_od(layerindex) * weightingfunctions->d_f[index];
            deriv.d_SSA = weightingfunctions->d_ssa[index] * (1 - layer_f(layerindex)) / (1 - layer_ssa(layerindex) * layer_f(layerindex)) -
                    weightingfunctions->d_f[index] * layer_ssa(layerindex) / (1 - layer_ssa(layerindex) * layer_f(layerindex)) +
                    (weightingfunctions->d_ssa[index] * layer_f(layerindex) + weightingfunctions->d_f[index] * layer_ssa(layerindex)) *
                    layer_ssa(layerindex) * (1 - layer_f(layerindex)) / ((1 - layer_ssa(layerindex) * layer_f(layerindex)) * (1 - layer_ssa(layerindex) * layer_f(layerindex)));

            for(int l = 0; l < (int)this->M_NSTR; ++l) {
                int coeffindex = derivindex + wavelidx * numderiv * this->M_NSTR + l*numderiv;
                deriv.d_legendre_coeff[l].a1 = weightingfunctions->d_a1[coeffindex] / (1 - layer_f(layerindex)) -
                        weightingfunctions->d_f[index] * (2*l+1) / (1 - layer_f(layerindex)) +
                        weightingfunctions->d_f[index] * (layer_a1(l, layerindex) - layer_f(layerindex) * (2*l+1)) / (1 - layer_f(layerindex)) / (1 - layer_f(layerindex));

                if constexpr (NSTOKES > 2) {
                    deriv.d_legendre_coeff[l].a2 = weightingfunctions->d_a2[coeffindex];
                    deriv.d_legendre_coeff[l].a3 = weightingfunctions->d_a3[coeffindex];
                    deriv.d_legendre_coeff[l].b1 = weightingfunctions->d_b1[coeffindex];
                }
                if constexpr (NSTOKES > 3) {
                    deriv.d_legendre_coeff[l].a4 = weightingfunctions->d_a4[coeffindex];
                    deriv.d_legendre_coeff[l].b2 = weightingfunctions->d_b2[coeffindex];
                }

            }

        }
    }

    // Post configure the layers, PS beam transmittances and derivative calculations
    for (auto& layer : m_layers) {
        // Start by telling each layer to calculate the derivative of thickness
        layer->configureDerivative();
    }
    
    configureTransmission();

    // Configure azimuthal dependencies
    for (auto& layer : m_layers) {
        registerAzimuthDependency(*layer);
    }
    registerAzimuthDependency(m_albedo);

    // resize ground reflection vector
    m_ground_reflection.resize(this->M_NSTR, std::vector<Radiance<NSTOKES>>(
        static_cast<uint>(los.size()),
        Radiance<NSTOKES>((int)m_input_derivatives.numDerivative())
        ));
    m_reflection_computed.resize(this->M_NSTR, std::vector<bool>(los.size(), false));

    // TODO: Figure out surface emission in low level interface?

    m_surfaceemission = 0.0;
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::copyLegendre(VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>& container, const Eigen::Matrix<double, Eigen::Dynamic, 6>& legendre) {
    for (uint k = 0; k < this->M_NSTR; ++k) {
        const auto& temp = container[k];
        container[k].a1 = legendre(k, 0);

        if constexpr(NSTOKES > 2) {
            container[k].a1 = legendre(k, 0);
            container[k].a2 = legendre(k, 1);
            container[k].a3 = legendre(k, 2);
            container[k].b1 = legendre(k, 4);
        }
        if constexpr(NSTOKES > 3) {
            container[k].a4 = legendre(k, 3);
            container[k].b2 = legendre(k, 5);
        }
    }
}


template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::configureTransmission()
{
    // This is a weird function where it is called so many times that we basically have to have two separate paths
    // for when derivates are on/off

    bool compute_deriv = m_input_derivatives.numDerivative() > 0;

    // Sometimes we can get OD's so large that exp(-od) = 0, in which case we can't do log(exp(-od)) to recover od so we have
    // to store OD temporarily
    double od_temp;
    Eigen::VectorXd d_temp;

    if(compute_deriv) {
        for(auto& layer : m_layers) {
            layer->ceiling_beam_transmittance().resize(m_input_derivatives.numDerivative(), false);
            layer->floor_beam_transmittance().resize(m_input_derivatives.numDerivative(), true);

            layer->floor_beam_transmittance().value = 0;

            if (layer->index() == 0) {
                layer->ceiling_beam_transmittance().value = 1;
                od_temp = 0;
                layer->ceiling_beam_transmittance().deriv.setZero();

                d_temp = layer->ceiling_beam_transmittance().deriv;
            } else {
                layer->ceiling_beam_transmittance().value = m_layers[layer->index() - 1]->floor_beam_transmittance().value;
                layer->ceiling_beam_transmittance().deriv = m_layers[layer->index() - 1]->floor_beam_transmittance().deriv;
            }
            for (LayerIndex p = 0; p <= layer->index(); ++p) {
                const auto &dual_thickness = m_layers[p]->dual_thickness();
                layer->floor_beam_transmittance().value += m_chapman_factors(layer->index(), p) * dual_thickness.value;

                const auto seq = Eigen::seq(dual_thickness.layer_start,
                                            dual_thickness.layer_start + dual_thickness.deriv.size() - 1);
                layer->floor_beam_transmittance().deriv(seq) +=
                        m_chapman_factors(layer->index(), p) * dual_thickness.deriv;
            }


            layer->dual_average_secant().value =
                    (layer->floor_beam_transmittance().value - od_temp) /
                    layer->dual_thickness().value;

            layer->dual_average_secant().deriv = (layer->floor_beam_transmittance().deriv + d_temp) / layer->dual_thickness().value;
            const auto seq = Eigen::seq(layer->dual_thickness().layer_start, layer->dual_thickness().layer_start + layer->dual_thickness().deriv.size() - 1);

            if (layer->dual_thickness().deriv.size() > 0) {
                layer->dual_average_secant().deriv(seq) -= layer->dual_thickness().deriv / layer->dual_thickness().value * layer->dual_average_secant().value;
            }

            od_temp = layer->floor_beam_transmittance().value;
            layer->floor_beam_transmittance().value = std::exp(-layer->floor_beam_transmittance().value);

            d_temp = -1 * layer->floor_beam_transmittance().deriv;
            layer->floor_beam_transmittance().deriv = layer->floor_beam_transmittance().value * -1 * layer->floor_beam_transmittance().deriv;
        }
    } else {
        for(auto& layer : m_layers) {
            layer->floor_beam_transmittance().value = 0;

            if (layer->index() == 0) {
                layer->ceiling_beam_transmittance().value = 1;
                od_temp = 0;
            } else {
                layer->ceiling_beam_transmittance().value = m_layers[layer->index() - 1]->floor_beam_transmittance().value;
            }
            for (LayerIndex p = 0; p <= layer->index(); ++p) {
                const auto &dual_thickness = m_layers[p]->dual_thickness();
                layer->floor_beam_transmittance().value += m_chapman_factors(layer->index(), p) * dual_thickness.value;
            }
            layer->dual_average_secant().value =
                    (layer->floor_beam_transmittance().value - od_temp) /
                    layer->dual_thickness().value;

            od_temp = layer->floor_beam_transmittance().value;
            layer->floor_beam_transmittance().value = std::exp(-layer->floor_beam_transmittance().value);
        }
    }

}

template<int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::assignLegendreDerivative(std::vector<LegendreCoefficient<NSTOKES>>& d_legendre,
                                                                                 const Eigen::Matrix<double, Eigen::Dynamic, 6>& species_legendre,
                                                                                 const std::vector<LegendreCoefficient<NSTOKES>>& layer_legendre,
                                                                                 double species_ssa,
                                                                                 double layer_ssa,
                                                                                 double thickness) const {
	for (uint l = 0; l < this->M_NSTR; ++l)
	{
		d_legendre[l].a1 = species_ssa / (thickness * layer_ssa) * (species_legendre(l, 0) - layer_legendre[l].a1);

        if constexpr(NSTOKES > 2) {
            d_legendre[l].a2 = species_ssa / (thickness * layer_ssa) * (species_legendre(l, 1) - layer_legendre[l].a2);
            d_legendre[l].a3 = species_ssa / (thickness * layer_ssa) * (species_legendre(l, 2) - layer_legendre[l].a3);
            d_legendre[l].b1 = species_ssa / (thickness * layer_ssa) * (species_legendre(l, 4) - layer_legendre[l].b1);

        }
        if constexpr(NSTOKES > 3) {
            d_legendre[l].a4 = species_ssa / (thickness * layer_ssa) * (species_legendre(l, 3) - layer_legendre[l].a4);
            d_legendre[l].b2 = species_ssa / (thickness * layer_ssa) * (species_legendre(l, 5) - layer_legendre[l].b1);
        }
	}
}


template<int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::assignLegendreDerivativeTest(std::vector<LegendreCoefficient<NSTOKES>>& d_legendre,
                                                                                     const std::vector<double>& legendre) const {
    for (uint l = 0; l < std::min((size_t)this->M_NSTR, legendre.size()); ++l)
    {
        d_legendre[l].a1 = legendre[l];
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>::configureTest(const PersistentConfiguration<NSTOKES, CNSTR>& config, const std::vector<testing::TestLayer<NSTOKES>>& testlayers)
{
	m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);

    // TODO: if heights are given as part of the layer interface we could calculate real chapman factors
	m_chapman_factors.setConstant(1 / this->M_CSZ);

	m_layers.reserve(this->M_NLYR);
	double od_ceil = 0;
	for(LayerIndex lidx = 0; lidx < this->M_NLYR; ++lidx) {

        // Copy the legendre coefficients to a new vector
        std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>> lephasef(new VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);
        for( int i = 0; i < (int)this->M_NSTR; ++i) {
			(*lephasef)[i] = testlayers[lidx].lephasef[i];
        }

		m_layers.push_back(
			std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(new OpticalLayer<NSTOKES, CNSTR>(config,
										   lidx,
										   testlayers[lidx].ssa,
										   1.0,
										   std::move(lephasef),
										   od_ceil,
										   od_ceil + testlayers[lidx].optical_depth,
				-1, -1,
										   m_input_derivatives)));
		od_ceil += testlayers[lidx].optical_depth;
	}
	// Configure azimuthal dependencies
	for(auto& layer : m_layers) {
		registerAzimuthDependency(*layer);
	}
	registerAzimuthDependency(m_albedo);
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::OpticalLayerArray);

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 1>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 1>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 3>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 3>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 1, 2>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 1, 2>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 3, 2>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 3, 2>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 4, 2>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 4, 2>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 1, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 1, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 3, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 3, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 4, 4>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 4, 4>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 1, 16>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 1, 16>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 3, 16>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 3, 16>;

template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::UP, 4, 16>;
template class sasktran_disco::OpticalLayerArrayIterator<sasktran_disco::Propagating::DOWN, 4, 16>;