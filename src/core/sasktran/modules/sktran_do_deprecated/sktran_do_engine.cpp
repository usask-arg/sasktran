#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_engine.h"

// The state of the SKTRAN_DO_Engine is managed using a RAII approach. 
// Scope-bound resource management make calculation setup/clean up simple and 
// reliable. Because of this SKTRAN_DO_Engine has a single member, m_config
// which stores all data members who's lifetime is greater than a single 
// radiance calculation.

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::configureModel(const SKTRAN_DO_UserSpec& spec,
									                  const SKTRAN_LineOfSightArray_V21& los,
									                  std::mutex* opticalproperties_mutex,
									                  std::vector<sktran_do_detail::LOSDiagnostics>* los_diag)
{
	using namespace sktran_do_detail;
	m_config.configureUserSpec(spec, los, opticalproperties_mutex, los_diag);
}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::prefillWavelengthTables(const std::vector<double>& wavelengths,
                                               sktran_do_detail::SASKTRANAtmosphereInterface*	atmosphereinterface,
                                               std::mutex* opticalproperties_mutex) {
    m_config.preConfigureWavelengthTables(wavelengths, atmosphereinterface, &m_opticalstate);

    m_geometrylayers = std::make_unique<sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>>(m_config, &m_opticalstate);
}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::configureTest(sktran_do_detail::SKTRAN_DO_TestSpec<NSTOKES>& testspec)
{
	m_config.configureTest(testspec);
}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>
::calculateRadianceSpherical(std::vector<double>& los_radiances,
	double wavelength,
	sktran_do_detail::SASKTRANAtmosphereInterface* opticalstate,
	std::vector<std::vector<double>>* los_wf,
	sktran_do_detail::RTSDiagnostics* diag,
	std::vector<double>* los_transmission,
    std::vector<double>* layer_od,
    std::vector<double>* layer_ssa,
    int wavelidx
    )
{
	using namespace sktran_do_detail;

	double prev_coszen = -999;
	bool recache_lp_coszen;
	double cc_epsilon = m_config.userSpec()->getCCEpsilon();	//< Epsilon in Cauchy convergence criterion


	// Sizing specs
	const uint NSTR = m_config.nstr();
	const uint NLYR = m_config.nlyr();

	uint num_ptrbs;
	if (m_config.perturbation_specs()) {
		num_ptrbs = static_cast<uint>(m_config.perturbation_specs()->size());
	}
	else {
		num_ptrbs = 0;
	}


	// Configure engine for this calculation
	std::vector<LineOfSight> linesofsight;
	std::unique_lock<std::mutex> optical_properties_lock(m_config.opticalPropertiesMutex());
	std::unique_ptr<BRDF_Base> brdf;
	m_config.configureRadianceCalculation(&los_radiances, wavelength, opticalstate, linesofsight, los_wf, brdf, &m_opticalstate); //
    if(m_geometrylayers == nullptr && !m_config.thisIsATest()) {
        m_geometrylayers = std::make_unique<sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>>(m_config, &m_opticalstate);
    }

	configureRTSDiagnostics(diag, linesofsight);

	const std::vector<double>& cos_sza = m_config.spherical_cos_sza();

	std::unique_ptr<SphericalPostProcessing<NSTOKES, CNSTR>> spherical_postprocess = nullptr;
	if (!m_config.ss_only()) {
		spherical_postprocess = std::unique_ptr<SphericalPostProcessing<NSTOKES, CNSTR>>(new SphericalPostProcessing<NSTOKES, CNSTR>(cos_sza, m_config, wavelength, &m_opticalstate, linesofsight, *m_geometrylayers, m_config.userSpec()->getUseUpwellingSpherical()));
	}
	optical_properties_lock.unlock();


	SphericalSolarTransmission<NSTOKES> solar_transmission(*m_config.ray_tracer(), m_opticalstate, wavelength);
	SphericalIntegrator<NSTOKES, CNSTR> integrator(m_opticalstate, solar_transmission, *m_config.coords(), wavelength, spherical_postprocess.get());
	
	auto& all_rays = m_config.traced_rays();

	VectorDim1<Dual<double>> radiance_components;				//< Contribution of fourier expansion components
	VectorDim2<double> lp_coszen;								//< Legendre polynomials evaluated at LOS coszen

	if (los_transmission != nullptr) {
		los_transmission->resize(linesofsight.size(), 1);
	}

	for (LineOfSight& los : linesofsight) {
		auto& ray = all_rays[los.unsorted_index];
		RayOptical rayoptical;
		solar_transmission.add_transmission(ray, rayoptical);

		const double coszen = los.coszenith;
		Dual<double> los_dual_radiance(num_ptrbs);
		// If the LOS zenith is new then we need to start the radiance calculation over
		if (prev_coszen == coszen) {
			recache_lp_coszen = false;
		}
		else {
			recache_lp_coszen = true;
			lp_coszen.clear();
			lp_coszen.reserve(NSTR);
			radiance_components.clear();
			radiance_components.reserve(NSTR);
		}

		// Construct solutions from terms in azimuthal fourier expansion until radiance converges
		AEOrder m = 0;
		std::function<bool(uint)> calcNotConverged;
		uint forced_number_azimuth_terms = m_config.userSpec()->getForcedNumberAzimuthTerms();
		if (forced_number_azimuth_terms) {
			calcNotConverged = [&](uint converged_count) -> bool {
				return m < forced_number_azimuth_terms;
			};
		}
		else {
			calcNotConverged = [&](uint converged_count) -> bool {
				return converged_count < 2 && m < NSTR;
			};
		}
		for (uint converged_count = 0; calcNotConverged(converged_count); ++m) {
			// Solve the DE for order m. Solve will return immediately if the term has already been solved
			if (spherical_postprocess != nullptr) {
				spherical_postprocess->solve(m);
			}

			// Update or continue filling cache
			if (recache_lp_coszen || lp_coszen.size() == m) {
				// Update legendre polynomial cache
				lp_coszen.push_back(VectorDim1<double>(NSTR));
				for (LPOrder l = 0; l < NSTR; ++l) lp_coszen[m][l] = boost::math::legendre_p(l, m, coszen);
				radiance_components.push_back(num_ptrbs);

				radiance_components[m] = integrator.integrate_sources(m, ray, rayoptical, lp_coszen[m]);

			}

			// Accumulate contribution and check for convergence
			//double cos_daz = cos(m * (m_config.saz() - los.azimuth));
			*los.radiance += radiance_components[m].value;
			los_dual_radiance += radiance_components[m];
			if (abs(radiance_components[m].value / (*los.radiance)) < cc_epsilon) ++converged_count;
		}

		// Finish up weighting function calculation
		// TODO implement this log_x mode
		bool return_wf_by_log_x = m_config.userSpec()->getWFFormIsByLogX();
		if (m_config.perturbation_specs()) {
			los.wf->resize(num_ptrbs);
			for (uint k = 0; k < num_ptrbs; ++k)
			{
				(*los.wf)[k] = los_dual_radiance.deriv[k];
			}
		}
		// Note current zenith angle to avoid unnecessary cache updating 
		prev_coszen = los.coszenith;
	}

}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>
::calculateRadiancePlaneParallel(std::vector<double>& los_radiances,
	double wavelength,
	sktran_do_detail::SASKTRANAtmosphereInterface* opticalstate,
	std::vector<std::vector<double>>* los_wf,
	sktran_do_detail::RTSDiagnostics* diag,
	std::vector<double>* los_transmission,
    std::vector<double>* layer_od,
    std::vector<double>* layer_ssa,
    int wavelidx
    )
{
	using namespace sktran_do_detail;

	// Sizing specs
	const uint NSTR = m_config.nstr();
	const uint NLYR = m_config.nlyr();

	// Configure engine for this calculation
	std::vector<LineOfSight> linesofsight;
    if( !m_config.opticalstate_prefilled() ) {
        m_config.opticalPropertiesMutex().lock();
    }
	std::unique_ptr<BRDF_Base> brdf;
	m_config.configureRadianceCalculation(&los_radiances, wavelength, opticalstate, linesofsight, los_wf, brdf, &m_opticalstate); //
	configureRTSDiagnostics(diag, linesofsight);

	// This only needs to be configured once, but has to be done using the geometry information in the optical tables
	// either we do it here for the first time or we do it during prefill wavelength tables
	// Could probably refactor this to move into config.configureRadianceCalculation
    if(m_geometrylayers == nullptr && !m_config.thisIsATest()) {
        m_geometrylayers = std::make_unique<sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>>(m_config, &m_opticalstate);
    }

	// Create stacked homogeneous layers with atmospheric properties to approximate the atmosphere
    // Should pass in wavelidx here but the last wavelength gets messed up if we do? It makes no sense
	OpticalLayerArray<NSTOKES, CNSTR> layers(m_config, wavelength, &m_opticalstate, linesofsight, std::move(brdf), diag, true, wavelidx, m_geometrylayers.get());

	if (!m_config.opticalstate_prefilled()) {
		m_config.opticalPropertiesMutex().unlock();
	}

    if( layer_od != nullptr ) {
        layer_od->resize(NLYR);
        for(int p = 0; p < NLYR; ++p) {
            layer_od->at(p) = layers.layer(p).dual_thickness().value;
        }
    }

    if( layer_ssa != nullptr ) {
        layer_ssa->resize(NLYR);
        for(int p = 0; p < NLYR; ++p) {
            layer_ssa->at(p) = layers.layer(p).ssa();
        }
    }


	uint num_ptrbs;
	if (m_config.perturbation_specs()) {
		num_ptrbs = static_cast<uint>(m_config.perturbation_specs()->size());
	}
	else {
		num_ptrbs = static_cast<uint>(layers.inputDerivatives().numDerivative());
	}

	// Create a DE solver
	RTESolver<NSTOKES, CNSTR> rte(m_config, layers);

	// Local vars used to calculate radiances
	double prev_coszen = 0;										//< last line of sight cosine of zenith angle, used to check if caches need updating
	bool recache_lp_coszen;										//< flag indicating legendre polynomial cache for LOS needs to be updated
	double cc_epsilon = m_config.userSpec()->getCCEpsilon();	//< Epsilon in Cauchy convergence criterion
	VectorDim1<Radiance<NSTOKES>> radiance_components;			//< Contribution of fourier expansion components
	VectorDim2<LegendrePhaseContainer<NSTOKES>> lp_coszen;		//< Legendre polynomials evaluated at LOS coszen

	if (los_transmission != nullptr) {
		los_transmission->resize(linesofsight.size(), 1);
	}

																// Do the radiance calculations for each line of sight
	for (LineOfSight& los : linesofsight) {
		double losod;
		// Use the exact losod for the attenutation calculation
		if (m_config.thisIsATest()) {
			losod = 0.0;
		}
		else {
			losod = layers.opticalDepthAt(los.observeraltitude);
		}

		double layerlosod;
		// And an approximate layerlosod for the multiple scatter integration
		if (losod > 0) {
			const auto& layer = layers.layerAt(losod);

			layerlosod = layer->opticalDepth(Location::FLOOR) - (los.observeraltitude - layer->altitude(Location::FLOOR)) / (layer->altitude(Location::CEILING) - layer->altitude(Location::FLOOR)) * layer->opticalDepth(Location::INSIDE);
		}
		else {
			layerlosod = losod;
		}
		if (diag != nullptr) {
			diag->los_diagnostic[los.unsorted_index].observer_od = losod;
		}

		const double coszen = los.coszenith;
		Radiance<NSTOKES> los_dual_radiance(num_ptrbs);
		// If the LOS zenith is new then we need to start the radiance calculation over
		if (prev_coszen == coszen) {
			recache_lp_coszen = false;
		}
		else {
			recache_lp_coszen = true;
			lp_coszen.clear();
			lp_coszen.reserve(NSTR);
			radiance_components.clear();
			radiance_components.reserve(NSTR);
		}

		// Construct solutions from terms in azimuthal fourier expansion until radiance converges
		AEOrder m = 0;
		std::function<bool(uint)> calcNotConverged;
		uint forced_number_azimuth_terms = m_config.userSpec()->getForcedNumberAzimuthTerms();
		if (forced_number_azimuth_terms) {
			calcNotConverged = [&](uint converged_count) -> bool {
				return m < forced_number_azimuth_terms;
			};
		}
		else {
			calcNotConverged = [&](uint converged_count) -> bool {
				return converged_count < 2 && m < NSTR;
			};
		}
        Dual<double> transmittance(num_ptrbs);
        Radiance<NSTOKES> integral(layers.inputDerivatives().numDerivative());
        Radiance<NSTOKES> integral_converted(num_ptrbs);

        for (uint converged_count = 0; calcNotConverged(converged_count); ++m) {
			// Solve the DE for order m. Solve will return immediately if the term has already been solved
			rte.solve(m);

			// Update or continue filling cache
			if (recache_lp_coszen || lp_coszen.size() == m) {
				// Update legendre polynomial cache
				lp_coszen.push_back(VectorDim1<LegendrePhaseContainer<NSTOKES>>(NSTR));
				for (LPOrder l = 0; l < NSTR; ++l) lp_coszen[m][l].fill(m, l, coszen);
				radiance_components.push_back(layers.inputDerivatives().numDerivative());

				// Get radiance reflected off the ground
				radiance_components[m] = convert_dual_to_wf(layers.reflectedIntensity(m, los),
					layers.inputDerivatives(),
					num_ptrbs);	// Radiance reflected off the ground

				// Calculate radiance recursively upwards through atmosphere
				for (auto layer = layers.template iteratorAcross<Propagating::UP>(); layer.isValid(); ++layer) {
					if (layer.entryOpticalDepth() < losod) {
						// Layer doesn't contribute
						continue;
					}
					layer.attenuationFactor(coszen, losod, los.observeraltitude, transmittance);
					if (m==0 && los_transmission != nullptr) {
						los_transmission->at(los.unsorted_index) *= transmittance.value;
					}
					ParticipatingSourceTerm<NSTOKES, CNSTR>(m, layer.ptr(), layers, coszen, layerlosod, lp_coszen[m], m_config.ss_only()).integrate(integral, m_config.pool().thread_data().postprocessing_cache(layer.layer().index()));

					convert_dual_to_wf(integral, layers.inputDerivatives(), integral_converted);

					radiance_components[m].apply_transmission_factor(transmittance);

					// Add in the source term
					radiance_components[m].value += integral_converted.value;
					radiance_components[m].deriv += integral_converted.deriv;
				}
			}

			// Accumulate contribution and check for convergence
            // This copy is necessary if we want to cache the results for a single zenith and have multiple azimuth
            // lines of sight.  But this scenario is so rare that we should figure out a way to remove this copy
            // if there is only 1 line of sight (the most common scenario)
			Radiance<NSTOKES> component = radiance_components[m];

            component.apply_azimuth_expansion((m_config.saz() - los.azimuth), m);

			accumulateStokes(los.radiance, component);

			los_dual_radiance.value += component.value;
			los_dual_radiance.deriv += component.deriv;

            if(component.converged(*los.radiance, cc_epsilon)) {
                ++converged_count;
            } else {
                converged_count = 0;
            }

        }

		// Finish up weighting function calculation
		// TODO implement this log_x mode
		bool return_wf_by_log_x = m_config.userSpec()->getWFFormIsByLogX();
		if (m_config.perturbation_specs()) {
			uint num_ptrbs = static_cast<uint>(m_config.perturbation_specs()->size());

			assignStokesDeriv(los.wf, los_dual_radiance, num_ptrbs);
		}

		// Note current zenith angle to avoid unnecessary cache updating 
		prev_coszen = los.coszenith;
    }
}

template <>
void SKTRAN_DO_Engine<1>::assignStokesDeriv(std::vector<double>* wf, const sktran_do_detail::Radiance<1>& radiance, int num_ptrb) {
	wf->resize(num_ptrb);
	for (int k = 0; k < num_ptrb; ++k)
	{
		(*wf)[k] = radiance.deriv(k);
	}
}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::assignStokesDeriv(std::vector<double>* wf, const sktran_do_detail::Radiance<NSTOKES>& radiance, int num_ptrb) {
	wf->resize(num_ptrb * NSTOKES);
	for (int k = 0; k < num_ptrb; ++k)
	{
		for (int s = 0; s < NSTOKES; ++s) {
			int linearindex = k * NSTOKES + s;
			(*wf)[linearindex] = radiance.deriv(k, s);
		}

	}
}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::accumulateStokes(double* ptr, const sktran_do_detail::Radiance<NSTOKES>& radiance) {
    if constexpr(NSTOKES == 1) {
        *ptr += radiance.I();
    } else {
        for (int i = 0; i < NSTOKES; ++i) {
            *(ptr+i) += radiance.value(i);
        }
    }

}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>
::calculateRadiance(std::vector<double>& los_radiances,
					double wavelength,
					sktran_do_detail::SASKTRANAtmosphereInterface*	opticalstate,
					std::vector<std::vector<double>>* los_wf, 
					sktran_do_detail::RTSDiagnostics* diag,
					std::vector<double>* los_transmission,
                    std::vector<double>* layer_od,
                    std::vector<double>* layer_ssa,
                    int wavelidx)
{
	if (m_config.use_los_spherical()) {
		calculateRadianceSpherical(los_radiances, wavelength, opticalstate, los_wf, diag, los_transmission, layer_od, layer_ssa, wavelidx);
	} else {
		calculateRadiancePlaneParallel(los_radiances, wavelength, opticalstate, los_wf, diag, los_transmission, layer_od, layer_ssa, wavelidx);
	}
}

// The rest of this file is just development code
#pragma region "Development details"

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::runTest(const sktran_do_detail::testing::TestCase<NSTOKES>& testcase,
													   std::vector<double>& abs_diff)
{
	using namespace sktran_do_detail;
	using namespace testing;

	// Generate test case specifications
	SKTRAN_DO_TestSpec<NSTOKES> testspec;
	testspec.configure(testcase);

    // For test cases we want to disable convergence checking to make sure
    // we get the most accurate results
    testspec.setCauchyCriterion(1e-99);

	// Make an engine and then calculate radiances
	SKTRAN_DO_Engine engine;
	engine.configureTest(testspec);
	std::vector<SKTRAN_StokesScalar> radiances;
	engine.CalculateRadiance(&radiances, 0, 0, nullptr);
	abs_diff.resize(radiances.size());
	for(uint i = 0; i < radiances.size(); ++i) {
		double diff = std::abs(radiances[i] - testcase.correct_radiances->at(i));
		abs_diff[i] = diff;

        // Useful to print output for updating tests
        //if( i % 5 == 0 && i > 0) {
        //    printf("\n");
        //}
        //printf("%.12f, ", radiances[i]);
	}

	return;
}

template <int NSTOKES, int CNSTR>
void SKTRAN_DO_Engine<NSTOKES, CNSTR>::configureRTSDiagnostics(sktran_do_detail::RTSDiagnostics* diag, std::vector<sktran_do_detail::LineOfSight>& linesofsight) const
{
	if(diag == nullptr) return;
	auto nstr = m_config.nstr();
	auto nlyr = m_config.nlyr();
	sktran_do_detail::uint nptrbs = 0;
	if(m_config.perturbation_specs()) {
		nptrbs = static_cast<sktran_do_detail::uint>(m_config.perturbation_specs()->size());
	}
	diag->initialize(nstr, nlyr, static_cast<sktran_do_detail::uint>(linesofsight.size()), nptrbs);
}

INSTANTIATE_TEMPLATE(SKTRAN_DO_Engine);

#if 0
bool SKTRAN_DO_Engine::CompareDISORT(std::vector<SKTRAN_StokesScalar>*	 losradiance,
									 double								 wavelen,
									 size_t								 numordersofscatter,
									 SKTRAN_AtmosphericOpticalState_V21* opticalstate)
{
	// This is all development stuff that probably no one beside me [Liam] cares about
	using namespace sktran_do_detail;
	bool ok = true;
	ok &= CalculateRadiance(losradiance, wavelen, numordersofscatter, opticalstate);
	// Configuration ... this is redundant as it has already been done but it needs to be done to get linesofsight
	std::vector<LineOfSight> los;

	std::unique_ptr<BRDF_Base> brdf;
	std::vector<SKTRAN_StokesScalar> temp;
	m_config.configureRadianceCalculation(&temp, wavelen, opticalstate, los, nullptr, brdf);
	*los[0].radiance = losradiance->at(0);
	NXASSERT((los.size() == 1));
	NXASSERT((m_config.userSpec()->compareWithDISORT()));
	OpticalLayerArray layers(m_config, wavelen, opticalstate, los, std::move(brdf));
	layers.configureAEOrder(0);
	DISORT_FPTR DISORT_F = m_config.userSpec()->getDISORT_FPTR();
	// Configure disort
	int NLYR = m_config.nlyr();
	int NSTR = m_config.nstr();
	std::vector<double> DTAUC(NLYR);
	std::vector<double> SSALB(NLYR);
	std::vector<int>    NMOM(NLYR, NSTR);
	std::vector<double> PMOM((NSTR + 1)*NLYR);
	for(LayerIndex p = 0; p < m_config.nlyr(); ++p) {
		DTAUC[p] = layers[p].opticalDepth(Location::FLOOR) - layers[p].opticalDepth(Location::CEILING);
		SSALB[p] = layers[p].ssa();
		NMOM[p] = NSTR;
		const VectorDim1<double>* lephasef = &layers[p].expansionPhaseFunction(0);
		for(int l = 0; l < NSTR; ++l) {
			PMOM[l + p*(NSTR + 1)] = (*lephasef)[l];
		}
		PMOM[NSTR + p*(NSTR + 1)] = 0;	// temp
	}
	// Get user optical depths (always 0)
	int USRTAU = 1;
	int NTAU = 1;
	std::vector<double> UTAU(NTAU);
	UTAU[0] = 0;
	// Get user directions
	int num_rays = static_cast<int>(los.size());
	int USRANG = true;
	int NUMU = num_rays;
	int NPHI = num_rays;
	std::vector<double> UMU(NUMU);
	std::vector<double> PHI(NPHI);
	for(uint ridx = 0; ridx < los.size(); ++ridx) {
		UMU[ridx] = los[ridx].coszenith;
		if(UMU[ridx] == 1.00) UMU[ridx] -= 0.001;
		PHI[ridx] = los[ridx].azimuth * 180.0 / PI;
	}
	// Additional disort config
	int IBCND = false;
	double FBEAM = m_config.userSpec()->getTopDirectIntensity();
	double FISOT = m_config.userSpec()->getTopDiffuseIntensity();
	int LAMBER = false;															// make false for brdf
	double ALBEDO = 1; // layers.albedo();										SUPER TEMP
	int PLANK = false;
	int ONLYFL = false;
	double ACCUR = 0.01;
	int PRINT[5] = { 1, 0, 1, 0, 1 };
	int MAXCLY = NLYR;
	int MAXULV = NLYR + 1;
	int MAXUMU = 2 * NUMU;
	int MAXPHI = 2 * NPHI;
	int MAXMOM = NSTR;
	std::vector<double> TEMPER(NLYR + 1);
	double WVNMLO;
	double WVNMHI;
	double UMU0 = m_config.csz();
	double PHI0 = m_config.saz() * 180.0 / PI;
	double BTEMP;
	double TTEMP;
	double TEMIS;
	std::vector<double> RFLDIR(MAXULV);
	std::vector<double> RFLDN(MAXULV);
	std::vector<double> FLUP(MAXULV);
	std::vector<double> DFDT(MAXULV);
	std::vector<double> UAVG(MAXULV);
	std::vector<double> UU(MAXUMU * MAXULV * MAXPHI);
	std::vector<double> ALBMED(MAXUMU);
	std::vector<double> TRNMED(MAXUMU);
	char HEADER[127] = "This instance was generated by SKTRAN_DO.";
	// Run disort
	try {

		SuppressConsoleOutputThisScope this_prevents_disort_from_printing_a_bunch_of_stuff;
		DISORT_F(&NLYR, DTAUC.data(), SSALB.data(), NMOM.data(), PMOM.data(),
				 TEMPER.data(), &WVNMLO, &WVNMHI, &USRTAU, &NTAU, UTAU.data(),
				 &NSTR, &USRANG, &NUMU, UMU.data(), &NPHI, PHI.data(), &IBCND,
				 &FBEAM, &UMU0, &PHI0, &FISOT, &LAMBER, &ALBEDO, &BTEMP, &TTEMP,
				 &TEMIS, &PLANK, &ONLYFL, &ACCUR, PRINT, HEADER, &MAXCLY, &MAXULV,
				 &MAXUMU, &MAXPHI, &MAXMOM, RFLDIR.data(), RFLDN.data(), FLUP.data(),
				 DFDT.data(), UAVG.data(), UU.data(), ALBMED.data(), TRNMED.data()
		);

	}
	catch(const std::exception& e) {
		nxLog::Record(NXLOG_ERROR,
					  "SKTRAN_DO_Engine::CompareDISORT, An exception was "
					  "thrown while running disort. Exception message: %s",
					  e.what()
		);
		ok = false;
	}
	double diff;
	double max_diff = 0;
	int LU = 0;
	for(int ridx = 0; ridx < num_rays; ++ridx) {
		//std::cout << "DISORT: \t" << UU[ridx + ridx*(MAXUMU*MAXULV)] << "\n";
		diff = abs(*los[ridx].radiance - UU[ridx + ridx*(MAXUMU*MAXULV)]);
		if(diff > max_diff) max_diff = diff;
	}
	if(max_diff > m_config.userSpec()->getDISORTTolerance()) {
		nxLog::Record(NXLOG_WARNING,
					  "SKTRAN_DO_Engine::CompareDISORT, Maximum difference "
					  "between SKTRAN_DO_Engine was %1.2e", max_diff
		);
	}
	return ok;
}

bool SKTRAN_DO_Engine::generateTestResults(const sktran_do_detail::testing::TestCase& testcase_original)
{
	using namespace sktran_do_detail;
	using namespace testing;
	bool ok = true;
	TestCase testcase = testcase_original;
	testcase.linesofsight.resize(1);
	// loop through linesofsight because disort loader only allows 1 LOS
	uint i = 0;
	for(TestLOS los : testcase_original.linesofsight) {
		testcase.linesofsight[0] = los;

		SKTRAN_DO_TestSpec testspec;
		testspec.configure(testcase);
		testspec.loadDISORT(L"C:/Users/lrb175/ARG-disort/disort-2.0/bin/bin/libdisort.dll", "disort_");

		// Make an engine and then calculate radiances
		SKTRAN_DO_Engine engine;
		engine.configureTest(testspec);
		std::vector<SKTRAN_StokesScalar> radiances;
		ok &= engine.CompareDISORT(&radiances, 0, 0, nullptr);
		std::cout << std::setprecision(std::numeric_limits<double>::digits10) << radiances[0] << ",\n";

		if(!ok) {
			std::cout << "\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EEEEEEEEE RRRRRRRRRRRRRRRRRRRRRR OOOOOOOOOOOOOOOOOOOOOOOOO RRRRRRRRRRRRRRRRRRRRRRRRRRRR~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n\n\n";
		}
		//ok &= (abs(radiances[0] - testcase_original.correct_radiances->at(i++)) < 1e-8);

	}

	return ok;
}
#endif
#pragma endregion