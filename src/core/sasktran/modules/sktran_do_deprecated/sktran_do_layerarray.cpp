#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_layerarray.h"

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::computeReflectedIntensities(AEOrder m, const sktran_do_detail::LineOfSight& los) {
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
    int layerStart = input_deriv.layerStartIndex(layer.index());
    int numLayerDeriv = input_deriv.numDerivativeLayer(layer.index());

    Radiance<NSTOKES> diffuse_contrib(input_deriv.numDerivative()), direct_contrib(input_deriv.numDerivative());
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

        Radiance<NSTOKES> stream_contrib(input_deriv.numDerivative());
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
            for(uint k = 0; k < numLayerDeriv; ++k) {
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
            for(uint k = 0; k < numLayerDeriv; ++k) {
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
        for(uint k = 0; k < numLayerDeriv; ++k) {
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
        for(uint k = 0; k < numLayerDeriv; ++k) {
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
sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::OpticalLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
													            double wavelength,
													            const OpticalState<NSTOKES>* opticalstate,
													            const std::vector<LineOfSight>& los,
													            std::unique_ptr<BRDF_Base> brdf,
                                                                RTSDiagnostics* diag,
													            bool include_direct_bounce_in_reflected,
                                                                int wavelidx,
                                                                const GeometryLayerArray<NSTOKES, CNSTR>* geometry_layers)
	: OpticalLayerArrayROP<NSTOKES>(config),
	m_albedo(los, *this->M_MU, this->M_CSZ, std::move(brdf), config.userSpec()->getNumBRDFQuadratureTerms()),
	m_direct_toa(this->M_SOLAR_DIRECT_INTENSITY),
	m_optical_state(opticalstate),
	m_config(config),
	m_include_direct_bounce(include_direct_bounce_in_reflected),
	m_surfaceemission(0.0),
    m_input_derivatives(config.pool().thread_data().input_derivatives())
{
	if(config.thisIsATest()) {
		m_wavel_index = 0;
		auto testcase = dynamic_cast<const SKTRAN_DO_TestSpec<NSTOKES>*>(config.userSpec())->getTestCase();
		configureTest(config, testcase->layers);
		configurePerturbations(wavelength, opticalstate, diag);
		// Configure layer derivatives
		for (auto& layer : m_layers) {
			layer->configureDerivative();
		}
		configureTransmission();

		// resize ground reflection vector
		m_ground_reflection.resize(this->M_NSTR, std::vector<Radiance<NSTOKES>>(
			static_cast<uint>(los.size()), 
			Radiance<NSTOKES>(m_input_derivatives.numDerivative())
		));
        if(m_input_derivatives.numDerivative() > 0) {
            for(auto& l : los) l.wf->resize(m_input_derivatives.numDerivative() * NSTOKES);
        }
		m_reflection_computed.resize(this->M_NSTR, std::vector<bool>(los.size(), false));

		brdf = testcase->getBRDF();
		m_albedo.injectTestingBRDF(std::move(brdf));
		return;
	}
    if(wavelidx == -1) {
        m_wavel_index = opticalstate->wavel_index(wavelength);
    } else {
        m_wavel_index = wavelidx;
    }

	// Setup depth array members
	m_cell_depths = &opticalstate->cumulative_od(m_wavel_index);

	m_layers.reserve(this->M_NLYR);
	m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);
	m_chapman_factors.setZero();

	Eigen::Matrix<double, Eigen::Dynamic, 6> legendre;
	legendre.resize(this->M_NSTR, 6);

    double ceiling_depth = 0;
    double floor_depth = 0;

	for(LayerIndex p = 0; p < this->M_NLYR; ++p) {
		double scat_ext = 0, total_ext = 0;
        double ceil_h = geometry_layers->layer_ceiling()(p);
        double floor_h = geometry_layers->layer_floor()(p);

        m_optical_state->compute_layer_quantities_from_interpolation(total_ext, scat_ext, legendre, geometry_layers->interpolating_matrix(), p, m_wavel_index);
        // Copy the legendre coefficients to a new vector
        std::unique_ptr<VectorDim1<sktran_do_detail::LegendreCoefficient<NSTOKES>>> lephasef(new VectorDim1<sktran_do_detail::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);
        copyLegendre(*lephasef, legendre);

        floor_depth += total_ext;

        scat_ext /= ceil_h - floor_h;
        total_ext /= ceil_h - floor_h;
        scat_ext = std::max(scat_ext, total_ext * this->m_userspec->getSSAEqual1Dither());

		m_layers.push_back(
			std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>(new OpticalLayer<NSTOKES, CNSTR>(config, p, scat_ext, total_ext, std::move(lephasef), ceiling_depth, floor_depth, ceil_h, floor_h, m_input_derivatives))
		);

        ceiling_depth = floor_depth;
	}
    m_chapman_factors = geometry_layers->chapman_factors();

	// Setup perturbations
	configurePerturbations(wavelength, opticalstate, diag);
	if(m_config.perturbation_specs()) {
		for(auto& l : los) l.wf->resize(m_input_derivatives.numDerivative());
	}

	// Post configure the layers, PS beam transmittances and derivative calculations
	for (auto& layer : m_layers) {
		// Start by telling each layer to calculate the derivative of thickness
		layer->configureDerivative();
	}
	configureTransmission();

	// Configure azimuthal dependencies
	for(auto& layer : m_layers) {
		registerAzimuthDependency(*layer);
	}
	registerAzimuthDependency(m_albedo);

	// resize ground reflection vector
	m_ground_reflection.resize(this->M_NSTR, std::vector<Radiance<NSTOKES>>(
		static_cast<uint>(los.size()), 
		Radiance<NSTOKES>(m_input_derivatives.numDerivative())
	));
	m_reflection_computed.resize(this->M_NSTR, std::vector<bool>(los.size(), false));

	// Check if we want to include surface emission
	if (m_config.userSpec()->getSurfaceEmissionValues().size() > 0) {
		// Have emissions
		auto surfaceemission = SurfaceEmission(m_config.userSpec()->getSurfaceEmissionWavelengths(), m_config.userSpec()->getSurfaceEmissionValues());
		m_surfaceemission = surfaceemission.emission(wavelength);
	}
	else {
		m_surfaceemission = 0.0;
	}

	// fill diagnostics
	if(diag == nullptr) return;
	for(uint l = 0; l < this->M_NLYR; ++l) {
		DiscreteAtmosphereDiagnostics& atmo_diag = diag->atmo_diagnostics;
		atmo_diag.layer_optical_depths[l] = layer(l).opticalDepth(Location::INSIDE);
		atmo_diag.layer_ssa[l] = layer(l).ssa();
		atmo_diag.layer_boundary_altitudes[l] = altitudeAt(layer(l).opticalDepth(Location::CEILING));
	}
	diag->atmo_diagnostics.layer_boundary_altitudes[this->M_NLYR] = 0;
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::copyLegendre(VectorDim1<sktran_do_detail::LegendreCoefficient<NSTOKES>>& container, const Eigen::Matrix<double, Eigen::Dynamic, 6>& legendre) {
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
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::configureTransmission()
{
	Dual<double> transmission_od_top(m_input_derivatives.numDerivative());
	Dual<double> transmission_od_bottom(m_input_derivatives.numDerivative());
	// Then calculate the dual of the PS transmission factors
	for (auto& layer : m_layers) {
		// For all layers above/equal to the layer
		transmission_od_bottom.value = 0.0;
		transmission_od_bottom.deriv.setZero();
		for (LayerIndex p = 0; p <= layer->index(); ++p) {
			const auto& dual_thickness = m_layers[p]->dual_thickness();
			const auto seq = Eigen::seq(dual_thickness.layer_start, dual_thickness.layer_start + dual_thickness.deriv.size() - 1);
			transmission_od_bottom.value += m_chapman_factors(layer->index(), p) * dual_thickness.value;
            if (dual_thickness.deriv.size() > 0) {
                transmission_od_bottom.deriv(seq) += m_chapman_factors(layer->index(), p) * dual_thickness.deriv;
            }
		}
		layer->configurePseudoSpherical(transmission_od_top, transmission_od_bottom);
		transmission_od_top = transmission_od_bottom;
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::configurePerturbations(double wavlen, const OpticalState<NSTOKES>* opticalstate, RTSDiagnostics* diag)
{
	// Check that perturbation vector is not null
	if(m_config.perturbation_specs() == nullptr) {
		return;
	}

	uint ptrb_group = 0;
	Eigen::Matrix<double, Eigen::Dynamic, 6> legendre;
	legendre.resize(this->M_NSTR, 6);

    if(!m_input_derivatives.is_geometry_configured()) {
        for (const auto &qty: *m_config.perturbation_specs()) {
            //	// Push back a LayerPerturbation

            // Try to cast it, try Species first
            if (qty.get()->type() == SKTRAN_DO_TestSpec<NSTOKES>::WeightingFunctionSpec::SpeciesWF) {
                auto spec = reinterpret_cast<const typename SKTRAN_DO_TestSpec<NSTOKES>::SpeciesWF *>(qty.get());

                // We map altitude perturbations to layer perturbations
                // Each altitude perturbation is a triangle, that is integrated over each layer to get a layer perturbation
                // The integral is only performed over layers that have a non-zero contribution

                double toa_altitude = m_config.userSpec()->getTopAltitude();
                double ground_altitude = m_config.userSpec()->getBottomAltitude();

                auto const *top_layer_ptr = layerAtAltitude(std::min(spec->altitude + spec->u_width, toa_altitude));
                auto const *bot_layer_ptr = layerAtAltitude(std::max(spec->altitude - spec->l_width, ground_altitude));

                if (top_layer_ptr == nullptr || bot_layer_ptr == nullptr) {
                    throw InternalError("Perturbation slab had a null boundary!");
                }

                const auto &top_layer = *top_layer_ptr;
                const auto &bottom_layer = *bot_layer_ptr;

                double width_above = std::min(toa_altitude, spec->altitude + spec->u_width) - spec->altitude;
                double width_below = spec->altitude - std::max(ground_altitude, spec->altitude - spec->l_width);

                for (LayerIndex p = top_layer.index(); p <= bottom_layer.index(); p++) {
                    const auto &ptrb_layer = layer(p);

                    // roughly we probably expect (layer_height / width) perturbations to be present in each layer
                    // Multiply by 2 for good measure, this just helps the allocs a bit, doesn't have an effect on the calculation
                    size_t size_hint = (ptrb_layer.altitude(Location::CEILING) - ptrb_layer.altitude(Location::FLOOR)) /
                                       (spec->l_width + spec->u_width) * 4;

                    auto &ptrb = m_input_derivatives.addIfNotExists(this->M_NSTR, p, spec->handle, size_hint);

                    ptrb.layer_index = ptrb_layer.index();
                    m_layers[ptrb.layer_index]->takeDerivative(true);

                    // Calculate perturbation quantities
                    // Take optical properties from the center of the optical layer
                    // where the perturbation spec's altitude is.

                    // Find the height of the layer in cm
                    double layer_ceiling = ptrb_layer.altitude(Location::CEILING) * 100;
                    double layer_floor = ptrb_layer.altitude(Location::FLOOR) * 100;

                    double layer_height_cm = (layer_ceiling - layer_floor);

                    // Calculate perturbation quantities
                    double intR_divLayerHeight = triangleFragmentArea(
                            layer_floor, layer_ceiling, spec->altitude * 100,
                            width_above * 100, width_below * 100
                    ) / layer_height_cm;

                    ptrb.d_optical_depth = 1;

                    ptrb.group_and_triangle_fraction.emplace_back(
                            std::pair<uint, double>(ptrb_group, intR_divLayerHeight * layer_height_cm));
                    ptrb.alt_and_widths.emplace_back(
                            std::tuple<double, double, double>(spec->altitude, width_above, width_below));
                    ptrb.extinctions.emplace_back(1);
                }
            } else if (qty.get()->type() == SKTRAN_DO_TestSpec<NSTOKES>::WeightingFunctionSpec::AlbedoWF) {
                auto spec = reinterpret_cast<const typename SKTRAN_DO_TestSpec<NSTOKES>::AlbedoWF *>(qty.get());

                const CLIMATOLOGY_HANDLE albedo_handle = SKCLIMATOLOGY_ALBEDO;

                // albedo can only exist in the lowest layer
                LayerIndex p = m_layers.back()->index();

                auto &ptrb = m_input_derivatives.addIfNotExists(this->M_NSTR, p, albedo_handle);

                ptrb.d_albedo = 1;

                ptrb.group_and_triangle_fraction.emplace_back(std::pair<uint, double>(ptrb_group, 1));
                ptrb.extinctions.emplace_back(1);
            } else if (qty.get()->type() == SKTRAN_DO_TestSpec<NSTOKES>::WeightingFunctionSpec::TestWF) {
                auto spec = reinterpret_cast<const typename SKTRAN_DO_TestSpec<NSTOKES>::TestWF *>(qty.get());

                LayerInputDerivative<NSTOKES> &deriv = m_input_derivatives.addDerivative(this->M_NSTR,
                                                                                         spec->layer_index);

                deriv.d_SSA = spec->ssa;
                deriv.d_optical_depth = spec->optd;

                assignLegendreDerivativeTest(deriv.d_legendre_coeff, spec->legendrecoeff);

                deriv.d_albedo = spec->albedo;

                deriv.group_and_triangle_fraction.emplace_back(std::pair<uint, double>(ptrb_group, 1));
                deriv.extinctions.emplace_back(1);
            } else {
                throw InternalError("An undefined weighting function spec was given");
            }

            ptrb_group = ptrb_group + 1;
        }
        m_input_derivatives.sort(this->M_NLYR);
        m_input_derivatives.set_geometry_configured();
    }
    double cs_scat, cs_ext, num_density;
    for(int k = 0; k < m_input_derivatives.numDerivative(); ++k) {
        auto& deriv = m_input_derivatives.layerDerivatives()[k];
        auto& ptrb_layer = layer(deriv.layer_index);

        if(!(deriv.handle == SKCLIMATOLOGY_UNDEFINED || deriv.handle == SKCLIMATOLOGY_ALBEDO)) {
            for(int i = 0; i < deriv.group_and_triangle_fraction.size(); ++i) {
                double alt = std::get<0>(deriv.alt_and_widths[i]);

                auto type = opticalstate->interpolate_species(deriv.handle, alt, cs_scat, cs_ext, num_density, legendre, m_wavel_index);

                double species_ssa;
                if (cs_ext > 0) {
                    species_ssa = cs_scat / cs_ext;
                }
                else {
                    // Derivitive will be 0 anyways so set this to whatever
                    species_ssa = 1;
                }
                double thickness = ptrb_layer.opticalDepth(Location::FLOOR) - ptrb_layer.opticalDepth(Location::CEILING);
                double layer_ceiling = ptrb_layer.altitude(Location::CEILING) * 100;
                double layer_floor = ptrb_layer.altitude(Location::FLOOR) * 100;

                double layer_height_cm = (layer_ceiling - layer_floor);
                if(i == 0) {
                    deriv.d_SSA = 1 / thickness * (species_ssa - ptrb_layer.ssa());
                    if(type != OpticalState<NSTOKES>::PurelyAbsorbing) {
                        assignLegendreDerivative(deriv.d_legendre_coeff, legendre, ptrb_layer.legendre_coeff(),
                                                 species_ssa, ptrb_layer.ssa(), thickness);
                    } // else Already initialized to 0
                }
                deriv.extinctions[i] = cs_ext;
            }
        }
    }

}

template<int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::assignLegendreDerivative(std::vector<LegendreCoefficient<NSTOKES>>& d_legendre,
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
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::assignLegendreDerivativeTest(std::vector<LegendreCoefficient<NSTOKES>>& d_legendre,
                                                                          const std::vector<double>& legendre) const {
    for (uint l = 0; l < std::min((size_t)this->M_NSTR, legendre.size()); ++l)
    {
        d_legendre[l].a1 = legendre[l];
    }
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayerArray<NSTOKES, CNSTR>::configureTest(const PersistentConfiguration<NSTOKES, CNSTR>& config, const std::vector<testing::TestLayer<NSTOKES>>& testlayers)
{
	m_chapman_factors.resize(this->M_NLYR, this->M_NLYR);

    // TODO: if heights are given as part of the layer interface we could calculate real chapman factors
	m_chapman_factors.setConstant(1 / this->M_CSZ);

	m_layers.reserve(this->M_NLYR);
	double od_ceil = 0;
	for(LayerIndex lidx = 0; lidx < this->M_NLYR; ++lidx) {

        // Copy the legendre coefficients to a new vector
        std::unique_ptr<VectorDim1<sktran_do_detail::LegendreCoefficient<NSTOKES>>> lephasef(new VectorDim1<sktran_do_detail::LegendreCoefficient<NSTOKES>>);

        lephasef->resize(this->M_NSTR);
        for( int i = 0; i < this->M_NSTR; ++i) {
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

INSTANTIATE_TEMPLATE(sktran_do_detail::OpticalLayerArray);

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 1>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 1>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 3>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 3>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 4>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 4>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 1, 2>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 1, 2>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 3, 2>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 3, 2>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 4, 2>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 4, 2>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 1, 4>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 1, 4>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 3, 4>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 3, 4>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 4, 4>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 4, 4>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 1, 16>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 1, 16>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 3, 16>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 3, 16>;

template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::UP, 4, 16>;
template class sktran_do_detail::OpticalLayerArrayIterator<sktran_do_detail::Propagating::DOWN, 4, 16>;