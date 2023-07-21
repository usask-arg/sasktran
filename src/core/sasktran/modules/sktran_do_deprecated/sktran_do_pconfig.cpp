#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_pconfig.h"

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureUserSpec(const SKTRAN_DO_UserSpec& userspec,
													                       const SKTRAN_LineOfSightArray_V21& linesofsight,
                                                                           std::mutex* opticalproperties_mutex,
													                       std::vector<sktran_do_detail::LOSDiagnostics>* los_diag)
{
	if(opticalproperties_mutex) {
		m_opticalpropertiesmutex = opticalproperties_mutex;
	} else {
		m_opticalpropertiesmutex = &m_no_mutex_given_mutex;
	}

	// This function is never called for an internel kernel test
	const_cast<bool&>(m_testing) = false;

	// Setup basic members: M_NSTR, M_NLYR, M_MU, M_WT, M_LP_MU, 
	// m_perturb_quantities, m_solar_position
	configureModelSpecs(&userspec);

	// Setup up m_coordinate and m_rayfactory
	configureRayTracing(linesofsight);
	
	// Configures M_CSZ, M_SAZ, m_lp_csz_storage, M_LP_CSZ, and solar intensities
	configureSolarPosition();

	// Configure solar position and look directions
	m_unsorted_los.clear();
	configureDirections(linesofsight);

	fillLOSDiagnostics(los_diag);

    int num_mem_pool = 1;
    if(use_los_spherical()) {
        num_mem_pool = m_spher_cos_sza.size();
    }

    m_pool.resize(num_mem_pool, sktran_do_detail::MemoryPool<NSTOKES, CNSTR>(this->M_NLYR, this->M_NSTR));
    m_poolindex = 0;
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureTest(const SKTRAN_DO_TestSpec<NSTOKES>& userspec)
{
	m_opticalpropertiesmutex = &m_no_mutex_given_mutex;

	// Manually override all SASKTRAN configuration
	const_cast<bool&>(m_testing) = true;
	auto testcase = userspec.getTestCase();
	const_cast<double&>(this->M_CSZ) = testcase->solar.csz;
	const_cast<double&>(this->M_SAZ) = testcase->solar.saz;
	const_cast<double&>(this->M_SOLAR_DIRECT_INTENSITY) = testcase->solar.intensities.direct;
	m_unsorted_los.clear();
	m_unsorted_los.reserve(testcase->linesofsight.size());
	for(const testing::TestLOS& los : testcase->linesofsight) {
		m_unsorted_los.push_back(LineOfSight());
		m_unsorted_los.back().coszenith = los.coszen;
		m_unsorted_los.back().azimuth = los.az;
	}
	configureModelSpecs(&userspec);
    // Tests will only have 1 CSZ
	m_lp_csz_storage.resize(1);
	m_lp_csz_storage[0] = std::unique_ptr<LegendrePolynomials<NSTOKES>>(new LegendrePolynomials<NSTOKES>(this->M_NSTR, this->M_CSZ));
	this->M_LP_CSZ = m_lp_csz_storage[0].get();
    m_pool.resize(1, sktran_do_detail::MemoryPool<NSTOKES, CNSTR>(this->M_NLYR, this->M_NSTR));
    m_poolindex = 0;
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureModelSpecs(const SKTRAN_DO_UserSpec* userspec)
{
	m_userspec = userspec;
	const_cast<uint&>(this->M_NSTR) = m_userspec->getNumberOfStreams();
	const_cast<uint&>(this->M_NLYR) = m_userspec->getNumberOfLayers();
	this->M_MU = m_userspec->getStreamAbscissae();
	this->M_WT = m_userspec->getStreamWeights();
    configureLP(userspec);
	m_perturb_quantities = m_userspec->perturbations();
	const_cast<bool&>(this->M_USE_PSEUDO_SPHERICAL) = m_userspec->getUsePseudoSpherical();
	const_cast<SKTRAN_DO_UserSpec::LayerConstructionMethod&>(this->M_LAYER_CONSTRUCTION) = m_userspec->getLayerConstructionMethod();

	const_cast<bool&>(this->M_USE_LOS_SPHERICAL) = m_userspec->getUseLOSSpherical();
	const_cast<bool&>(this->M_SS_ONLY) = m_userspec->getSSOnly();
	const_cast<size_t&>(this->M_NUM_SZA) = m_userspec->getNumSZA();
    const_cast<double&>(this->M_SZA_REL_SEP) = m_userspec->getSZARelSep();
    const_cast<bool&>(this->M_USE_GREENS_FUNCTION) = m_userspec->getUseGreensFunction();


	if(m_userspec->getForcedNumberAzimuthTerms() > this->M_NSTR) {
		throw InvalidConfiguration("Forced number of azimuth terms must be less than or equal to the number of streams!");
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureLP(const SKTRAN_DO_UserSpec* userspec)
{
    if constexpr(NSTOKES == 1) {
        this->M_LP_MU = userspec->getAbscissaeLegendreP1();
    } else {
        this->M_LP_MU = (const VectorDim3<LegendrePhaseContainer<NSTOKES>>*) userspec->getAbscissaeLegendreP4();
    }
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureRayTracing(const SKTRAN_LineOfSightArray_V21& linesofsight)
{
	// Setup m_raytracemanager and m_rayfactory

	// loop over all lines of sight and check which ones view the ground
	size_t num_limb = 0;
	size_t num_ground = 0;
	nxVector avg_ground(0, 0, 0);
	for(uint i = 0; i < linesofsight.NumRays(); ++i) {
		nxGeodetic geo;
		const SKTRAN_LineOfSightEntry_V2* los;
		linesofsight.GetRay(i, &los);
		nxVector ground_intersect, exit_pt;
		geo.GetShellHeightLocation(m_userspec->getBottomAltitude(), los->Observer(), los->Look(), &ground_intersect, &exit_pt);
		
		if (ground_intersect.IsValid()) {
			num_ground++;
			avg_ground += ground_intersect;
		}
		else {
			num_limb++;
		}
	}
	SKTRAN_RayTracingRegionManager ray_trace_manager;
	ray_trace_manager.SetSun(*m_userspec->getSolarPosition(linesofsight.MeanMJD()));
	if (num_limb > 0) {
		if (!m_userspec->getUseLOSSpherical()) {
			throw InvalidConfiguration("A non-Nadir line of sight was given when operating in purely plane parallel mode");
		}
		// Use the ray region manager to get the reference point
	}
	else {
		// Use the average ground point to override the reference point
		avg_ground /= (double)num_ground;

		// Configure coordinates from LOS
		nxGeodetic avg_ground_geo;
		avg_ground_geo.FromGeocentricVector(avg_ground);

		ray_trace_manager.SetReferencePoint(avg_ground_geo.GeodeticLatitude(), avg_ground_geo.GeodeticLongitude(), avg_ground_geo.Height(), linesofsight.MeanMJD());
	}

	ray_trace_manager.UpdateUndefinedParametersFromLinesOfSight(linesofsight);
	ray_trace_manager.MakeCoordinateSystem(&m_coordinates, m_userspec->getBottomAltitude(), m_userspec->getTopAltitude());

	double refsza;
	ray_trace_manager.GetSZA(&refsza, &m_minsza, &m_maxsza);

	m_spher_cos_sza.resize(this->M_NUM_SZA);

	if (this->M_NUM_SZA == 1) {
		m_spher_cos_sza[0] = nxmath::cosd(refsza);
	}
	else {
		double cossza_diff = (nxmath::cosd(m_minsza) - nxmath::cosd(m_maxsza)) / (this->M_NUM_SZA - 1);
		for (size_t i = 0; i < this->M_NUM_SZA; ++i) {
			m_spher_cos_sza[i] = nxmath::cosd(m_maxsza) + i * cossza_diff;
		}
	}

	m_raytracer = std::unique_ptr<SphericalRayTracer>(new SphericalRayTracer(m_userspec->getAltitudeGrid(), *m_coordinates));
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureSolarPosition()
{
	// Get cosine of solar zenith
	const_cast<double&>(this->M_CSZ) = m_coordinates->ReferencePoint(m_userspec->getBottomAltitude()).CosSZA();
	const_cast<double&>(this->M_SAZ) = 0; // azimuth is relative to SAZ
	if(this->M_CSZ <= 0) {
        std::cout << this->M_CSZ;
		throw InvalidConfiguration("The sun is below the horizon at the reference point!");
	}
	const_cast<double&>(this->M_SOLAR_DIRECT_INTENSITY) = m_userspec->getTopDirectIntensity();

	if (this->M_USE_LOS_SPHERICAL) {
		m_lp_csz_storage.resize(m_spher_cos_sza.size());
		for (uint i = 0; i < m_spher_cos_sza.size(); ++i) {
			m_lp_csz_storage[i] = std::unique_ptr<LegendrePolynomials<NSTOKES>>(new LegendrePolynomials<NSTOKES>(this->M_NSTR, m_spher_cos_sza[i]));
		}
	}
	else {
		m_lp_csz_storage.resize(1);
		m_lp_csz_storage[0] = std::unique_ptr<LegendrePolynomials<NSTOKES>>(new LegendrePolynomials<NSTOKES>(this->M_NSTR, this->M_CSZ));
	}

	this->M_LP_CSZ = m_lp_csz_storage[0].get();
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::correctLineOfSight(const SKTRAN_LineOfSightEntry_V2* los, nxVector& look, nxVector& obs) {
	nxVector local_reference_point;

	// Check to see if the ray hits the ground
	nxGeodetic geo;
	nxVector ground_intersect, exit_pt;
	geo.GetShellHeightLocation(m_userspec->getBottomAltitude(), los->Observer(), los->Look(), &ground_intersect, &exit_pt);

	if (ground_intersect.IsValid()) {
		local_reference_point = ground_intersect;
	}
	else {
		// Might be limb viewing
		nxVector tangent = geo.FromTangentPointLocation(los->Observer(), los->Look());
		if (tangent.IsValid()) {
			local_reference_point = tangent;
		}
		else {
			// Not limb or ground viewing, odd geometry
			local_reference_point = los->Observer();
		}
	}
	geo.FromGeocentricVector(local_reference_point);
	nxVector up, south, west;
	geo.GetGeodeticWestSouthUp(&west, &south, &up);
	// VZA, altitude, and SZA at the local reference
	double vza = los->Look().AngleTo(-1.0*up);
	double altitude = geo.Height();
	double sza = up.AngleTo(m_coordinates->SunUnit());

	nxVector sun_unit = m_coordinates->SunUnit();
	// Calculate the solar azimuth angle
	nxVector sun_proj = m_coordinates->SunUnit() - (m_coordinates->SunUnit() & up) * up;
	nxVector los_proj = los->Look() - (los->Look() & up) * up;

	double proj = sun_proj.UnitVector() & los_proj.UnitVector();
	if (proj > 1) {
		proj = 1;
	}
	else if (proj < -1) {
		proj = -1;
	}
	double saa = 180 - nxmath::acosd(proj);

	// Calculate the target location in heliodetic coordinates
	double earth_radius = m_coordinates->AltitudeToRadius(m_userspec->getBottomAltitude());
	nxVector target_up = nxVector(0.0, nxmath::sind(sza), nxmath::cosd(sza));
	nxVector target = (altitude + earth_radius) * target_up;

	// Set the geodetic to the observer position so we can get the observer altitude relative to the lowest altitude
	geo.FromGeocentricVector(los->Observer());
	double observer_altitude = geo.Height() - m_userspec->getBottomAltitude();

	// Quadratic equation terms to calculate the new observer distance
	double A = 1;
	double B = -2.0*earth_radius*nxmath::cosd(vza);
	double C = earth_radius * earth_radius - (observer_altitude + earth_radius) * (observer_altitude + earth_radius);

	// Distance along look vector
	double s1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A);
	double s2 = (-B + sqrt(B*B - 4 * A*C)) / (2 * A);

	// Calculate the look vector in heliodetic coordinates
	// We calculate the look vector at SZA=0 and then apply a rotation to the look vector
	nxVector look_vector(nxmath::sind(saa),  nxmath::cosd(saa), 0.0);

	look_vector.RotateAboutXaxis(sza);

	look_vector = -1*nxmath::sind(vza) * look_vector + nxmath::cosd(vza) * target_up;

	obs = target - s1 * look_vector;
	look = -1*look_vector;
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureDirections(const SKTRAN_LineOfSightArray_V21& linesofsight)
{
	// Loop through each line of sight
	const uint N = static_cast<uint>(linesofsight.NumRays());
	m_unsorted_los.resize(N);
	m_traced_rays.resize(N);

	for (uint i = 0; i < N; ++i) {
		// Get ray
		const SKTRAN_LineOfSightEntry_V2* los;
		linesofsight.GetRay(i, &los);

		HELIODETIC_UNITVECTOR look_unit;
		HELIODETIC_VECTOR obs;

		if (m_userspec->getLosAdjustmentMethod() == SKTRAN_DO_UserSpec::LineOfSightAdjustment::translate_observer) {
			look_unit = m_coordinates->GeographicToHelioUnitVector(los->Look());
			obs = m_coordinates->GeographicToHelio(m_coordinates->TranslateGeoidToOsculatingSphere(los->Observer()));
		}
		else if (m_userspec->getLosAdjustmentMethod() == SKTRAN_DO_UserSpec::LineOfSightAdjustment::match_target_angles) {
			nxVector new_observer, new_look;
			correctLineOfSight(los, new_look, new_observer);

			look_unit.SetCoords(new_look.X(), new_look.Y(), new_look.Z());
			obs.SetCoords(new_observer.X(), new_observer.Y(), new_observer.Z());
		}


		ViewingRay ray;
		ray.look_away = look_unit;
		ray.observer.FromVector(obs, m_coordinates.get());
		ray.unsorted_ray_idx = i;
		m_traced_rays[i] = m_raytracer->trace_ray(ray);
		if (use_los_spherical()) {
			m_raytracer->trace_solar_rays(m_traced_rays[i]);
		}

		if(!m_traced_rays[i].ground_is_hit && !m_userspec->getUseLOSSpherical()) {
			std::stringstream ss;
			auto obs = m_traced_rays[i].observer_and_look.observer;
			auto bad_los = m_traced_rays[i].observer_and_look.look_away;
			ss << "observer: [" << obs.Vector().X() << ", " << obs.Vector().Y() << ", " << obs.Vector().Z() << "], ";
			ss << "lineofsight: [" << bad_los.X() << ", " << bad_los.Y() << ", " << bad_los.Z() << "]";
			std::string err_msg = "A non-nadir line of sight was detected when using plane-parallel mode! The following line of sight does not intersect the ground: " + ss.str();
			throw InvalidConfiguration(err_msg);
		}

		// If the ray hits the ground we define the LOS parameters from the ground point,
		// else we use the tangent point.  Although it's unclear how these parameters are
		// used in fully spherical mode
		HELIODETIC_UNITVECTOR local_up;
		if (m_traced_rays[i].ground_is_hit) {
			auto& local = m_traced_rays[i].layers[0].exit;
			local_up = local.UnitVector();
		}
		else {
			// Define the parameters from the reference point, is this correct though?
			local_up = m_coordinates->ReferencePoint(m_coordinates->GroundAltitude()).UnitVector();
		}


		// Calculate cosine of zenith
		HELIODETIC_UNITVECTOR propagating_unit = look_unit;
		propagating_unit.Negate();
		double coszenith = propagating_unit & local_up;
		
		//double coszenith = -los->Observer().UnitVector() & los->Look();
		if(coszenith <= 0 && !m_userspec->getUseLOSSpherical()) {
			std::stringstream ss;
			auto obs = m_traced_rays[i].observer_and_look.observer;
			auto bad_los = m_traced_rays[i].observer_and_look.look_away;
			ss << "observer: [" << obs.Vector().X() << ", " << obs.Vector().Y() << ", " << obs.Vector().Z() << "], ";
			ss << "lineofsight: [" << bad_los.X() << ", " << bad_los.Y() << ", " << bad_los.Z() << "]";
			std::string err_msg = "A bad line of sight was detected! The following line of sight is looking at the ground: " + ss.str();
			throw InvalidConfiguration(err_msg);
		} else if(coszenith > 1.0) {
			coszenith = 1.0;
		}
		
		// Project sun vector onto ground plane
		HELIODETIC_UNITVECTOR sun_unit = m_coordinates->GeographicToHelioUnitVector(m_coordinates->SunUnit());
		HELIODETIC_VECTOR sun_vert_trans(local_up, sun_unit & local_up);
		const HELIODETIC_VECTOR sun(sun_unit, 1);
		HELIODETIC_VECTOR sun_projected = sun - sun_vert_trans;
		
		// Project los vector onto ground plane
		HELIODETIC_VECTOR los_vert_trans(local_up, propagating_unit & local_up);
		const HELIODETIC_VECTOR propagating(propagating_unit, 1);
		HELIODETIC_VECTOR los_projected = propagating - los_vert_trans;

        double diff_azimuth;
        if(los_projected.Magnitude() < 1e-15 && false) {
            // We are looking straight down, so we are symmetric in azimuth for I, however for Q/U there is
            // ambiguity as to what direction we should choose.  For this case we calculate the azimuth of the sun
            // relative to north
            // TODO: is this even correct?
            std::array<HELIODETIC_UNITVECTOR, 3> units;
            m_coordinates->ReferencePoint(0).LocalUnitVectors(&units[0], 3);
            diff_azimuth = acos(units[0] & sun_projected.UnitVector());
        } else {
            // Calculate relative azimuth

            // Resolution to Issue #1
            // In the following dot product it is possible for `proj` to have a magnitude slightly larger than 1.
            // If this is the case then we'll force it to a magnitude of 1 to prevent diff_azimuth from being nan.
            double proj = sun_projected.UnitVector() & los_projected.UnitVector();
            if(proj > 1) {
                proj = 1;
            } else if(proj < -1) {
                proj = -1;
            }
            diff_azimuth = acos(proj);
            // End of resolution to Issue #1
        }

		// Copy results to m_unsorted_los
		m_unsorted_los[i].coszenith = coszenith;
		m_unsorted_los[i].azimuth = PI - diff_azimuth;
		m_unsorted_los[i].observeraltitude = m_traced_rays[i].observer_and_look.observer.Altitude();
		m_unsorted_los[i].cos_scattering_angle = sun_unit & propagating_unit;
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::fillLOSDiagnostics(std::vector<sktran_do_detail::LOSDiagnostics>* los_diag) const
{
	if(los_diag == nullptr) return;
	los_diag->resize(m_unsorted_los.size());
	for(uint i = 0; i < m_unsorted_los.size(); ++i) {
		LOSDiagnostics& diagnostic = los_diag->at(i);
		diagnostic.reference_point = coords()->PointToGeodetic(m_coordinates->ReferencePoint(m_userspec->getBottomAltitude()), coords()->ReferencePointMJD());
		diagnostic.local_look_zentih = acos(m_unsorted_los[i].coszenith);
		diagnostic.local_solar_zenith = acos(this->M_CSZ);
		diagnostic.local_relative_azimuth = m_unsorted_los[i].azimuth;
	}
}


template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::configureRadianceCalculation(std::vector<SKTRAN_StokesScalar>* losradiance,
																  double wavelen,
																  SASKTRANAtmosphereInterface* opticalstate,
																  VectorDim1<LineOfSight>& los,
																  VectorDim2<double>* loswf,
																  std::unique_ptr<sktran_do_detail::BRDF_Base>& brdf,
																  OpticalStateInterface* state)
{
	// Configure sorted line of sight
	los.resize(m_unsorted_los.size());
	losradiance->clear();
	losradiance->resize(m_unsorted_los.size()*NSTOKES);
	if(loswf != nullptr) loswf->resize(m_unsorted_los.size(), std::vector<double>());
	else {
		if(perturbation_specs() != nullptr && perturbation_specs()->size() > 0) {
			throw InternalError("You have configured a weighting function calculation but did not give a place to store the results.");
		}
	}
	std::copy(m_unsorted_los.cbegin(), m_unsorted_los.cend(), los.begin());
	LineOfSight::sort(los, *losradiance, loswf, NSTOKES);
	// If we're testing then we can return
	if(m_testing) {
		return;
	}
    if(!m_opticalstate_prefilled) {
        // Setup optical state
        auto& radiigrid = m_userspec->getAltitudeGrid();
        Eigen::VectorXd altitudes;
        altitudes.resize(radiigrid.size());
        for (uint i = 0; i < altitudes.size(); ++i) {
            altitudes(i) = radiigrid[i];
        }

        state->configure(altitudes, this->M_NSTR);

        state->fill_tables(wavelen, coords()->PointToGeodetic(coords()->ReferencePoint(m_userspec->getBottomAltitude()),
                                                              coords()->ReferencePointMJD()), opticalstate, los);
        if (use_los_spherical()) {
            state->fill_derivative_tables(perturbation_specs());
        }
    }

	// Get albedo
	HELIODETIC_POINT refpt = m_coordinates->ReferencePoint(m_userspec->getBottomAltitude());
	GEODETIC_INSTANT ref_inst = coords()->PointToGeodetic(refpt, coords()->ReferencePointMJD());

	const skBRDF* skbrdf = opticalstate->albedo();

	brdf = std::unique_ptr<Wrapped_skBRDF>(new Wrapped_skBRDF(wavelen, *skbrdf, ref_inst));
	
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>::preConfigureWavelengthTables(const std::vector<double>& wavelengths,
                                                                             SASKTRANAtmosphereInterface* atmosphereinterface,
                                                                             OpticalStateInterface* ostate
                                                                             ) {
    // Have to do some temporary things
    VectorDim1<LineOfSight> los;
    los.resize(m_unsorted_los.size());
    std::copy(m_unsorted_los.cbegin(), m_unsorted_los.cend(), los.begin());
    VectorDim1<SKTRAN_StokesScalar> losradiance;
    losradiance.resize(m_unsorted_los.size()*NSTOKES);
    LineOfSight::sort(los, losradiance, nullptr, NSTOKES);

    // Setup optical state
    auto& radiigrid = m_userspec->getAltitudeGrid();
    Eigen::VectorXd altitudes;
    altitudes.resize(radiigrid.size());
    for (uint i = 0; i < altitudes.size(); ++i) {
        altitudes(i) = radiigrid[i];
    }

    ostate->configure(altitudes, this->M_NSTR);

    ostate->fill_tables(wavelengths, coords()->PointToGeodetic(coords()->ReferencePoint(m_userspec->getBottomAltitude()),
                                                          coords()->ReferencePointMJD()), atmosphereinterface, los);

	if (use_los_spherical()) {
		ostate->fill_derivative_tables(perturbation_specs());
	}

    m_opticalstate_prefilled = true;

    for(auto& lp_csz : m_lp_csz_storage) {
        for(int m = 0; m < this->M_NSTR; ++m) {
            lp_csz->configureAEOrder(m);
        }
    }
}

INSTANTIATE_TEMPLATE(sktran_do_detail::PersistentConfiguration);

