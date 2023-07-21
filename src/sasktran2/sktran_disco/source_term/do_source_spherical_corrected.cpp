#include <sasktran2/do_source.h>
#include <fstream>


namespace sasktran2 {
    template<int NSTOKES, int CNSTR>
    DOSourceSphericallyCorrected<NSTOKES, CNSTR>::DOSourceSphericallyCorrected(const sasktran2::Geometry1D &geometry,
                                                                               const sasktran2::raytracing::RayTracerBase &raytracer) : DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>(geometry, raytracer), m_integrator(false), m_raytracer(raytracer) {

    }

    template<int NSTOKES, int CNSTR>
    void DOSourceSphericallyCorrected<NSTOKES, CNSTR>::initialize_atmosphere(
            const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere) {

        DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::initialize_atmosphere(atmosphere);

        m_integrator.initialize_atmosphere(atmosphere);

    }

    template<int NSTOKES, int CNSTR>
    void DOSourceSphericallyCorrected<NSTOKES, CNSTR>::initialize_geometry(
            const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {
        DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::initialize_geometry(los_rays);

        const auto& do_config = *this->m_thread_storage[0].sza_calculators[0].persistent_config;

        //int nstr = do_config.nstr();
        int nlyr = (int)this->m_diffuse_storage->altitude_grid().grid().size();

        auto* quadrature_angles = do_config.quadrature_cos_angle();

        // angles * numlayers
        int num_rays = (int)quadrature_angles->size() * (int)nlyr;
        m_spherical_correction_rays.resize(num_rays);
        int c = 0;

        sasktran2::viewinggeometry::ViewingRay viewing_ray;
        sasktran2::Location obs_location;

        m_altitude_integration_angles.resize(nlyr);
        m_altitude_integration_weights.resize(nlyr);

        double earth_radius = this->m_geometry.coordinates().earth_radius();

        for(int i = 0; i < quadrature_angles->size(); ++i) {
            for(int j = 0; j < nlyr; ++j) {
                // TODO: Should we always use 0 saa here?
                Eigen::Vector3d location = this->m_geometry.coordinates().solar_coordinate_vector(do_config.csz(), 0, this->m_diffuse_storage->altitude_grid().grid()(j));

                obs_location.position = location;
                viewing_ray.observer = obs_location;

                Eigen::Vector3d look_vector_up_0azi = this->m_geometry.coordinates().look_vector_from_azimuth(location, 0, (*quadrature_angles)[i]);

                viewing_ray.look_away = look_vector_up_0azi;
                m_raytracer.trace_ray(viewing_ray, m_spherical_correction_rays[c]);
                ++c;
            }
        }


        m_integrator.initialize_geometry(m_spherical_correction_rays, this->m_geometry);

        m_do_solutions.resize(this->m_config->num_threads());

        for(auto& thread_soln : m_do_solutions) {
            thread_soln.resize(this->m_config->num_do_streams());
            for(auto& soln : thread_soln) {
                soln.resize(this->m_config->num_do_streams() * (this->m_geometry.altitude_grid().grid().size() - 1), NSTOKES);
            }
        }

        m_spherical_source_interpolator = this->m_diffuse_storage->geometry_interpolator(m_spherical_correction_rays,
                                                                                         false);

        m_ms_source = std::make_unique<DOSourceInterpolatedPostProcessingView<NSTOKES, CNSTR>>(*m_spherical_source_interpolator, *this->m_diffuse_storage);

        m_sources.push_back(m_ms_source.get());

    }

    template<int NSTOKES, int CNSTR>
    void DOSourceSphericallyCorrected<NSTOKES, CNSTR>::calculate(int wavelidx, int threadidx) {
        DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::calculate(wavelidx, threadidx);

        int numrays = (int)m_spherical_correction_rays.size();

        Eigen::Matrix<double, -1, NSTOKES> radiance(numrays, NSTOKES);
        Eigen::Matrix<double, -1, NSTOKES> do_radiance(numrays, NSTOKES);
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> radiance_storage;

        // Calculate spherical correction factors
        for(int i = 0; i < numrays; ++i) {
            radiance_storage.value.setZero();

            m_integrator.integrate(radiance_storage,
                                   m_sources,
                                   wavelidx,
                                   i,
                                   threadidx);

            radiance(i, Eigen::all).array() = radiance_storage.value.array();
        }

        const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

        Eigen::MatrixXd full_rad = m_do_solutions[threadidx][0];

        int numangle = this->m_thread_storage[threadidx].sza_calculators[0].persistent_config->nstr();
        int numheight = (int)this->m_diffuse_storage->altitude_grid().grid().size();


        Eigen::Map<Eigen::MatrixXd> spherical_radiance (radiance(Eigen::all, 0).data(), numheight, numangle);
        Eigen::Map<Eigen::MatrixXd> pp_radiance(full_rad(Eigen::all, 0).data(), numangle, numheight);

        /*
        Eigen::VectorXd weight_vector(numangle);
        for(int i = 0; i < numangle; ++i) {
            weight_vector(i) = (*weights)[i];
        }

        int numanglegrid = (int)this->m_diffuse_storage->cos_angle_grid().grid().size();

        int azioffset = numanglegrid*numheight;

        for(int l = 0; l < numheight; ++l) {
            double pp_weight = pp_radiance(Eigen::all, l).transpose().dot(weight_vector);
            double spher_weight = spherical_radiance(l, Eigen::all).dot(this->m_diffuse_storage->angle_integration_weights(l));
            double correction_factor = spher_weight / pp_weight;

            for(int k = 0; k < numangle; ++k) {
                int start = azioffset*k + l*numanglegrid;
                for(int s = 0; s < NSTOKES; ++s) {
                    this->m_diffuse_storage->source(threadidx)[s].value(
                            Eigen::seq(start, start + numanglegrid - 1)) *= correction_factor;
                }
            }
        }
         */
    }

    template<int NSTOKES, int CNSTR>
    void DOSourceSphericallyCorrected<NSTOKES, CNSTR>::accumulate_solved_azimuth(
            sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR> &optical_layer,
            DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage, int szaidx, sasktran_disco::AEOrder m,
            int threadidx) {
        DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR>::accumulate_solved_azimuth(
                optical_layer, thread_storage, szaidx, m, threadidx
                );

        // Store solutions at midpoint of layers where we are going to trace spherical rays from
        auto& soln = m_do_solutions[threadidx][m];
        soln.setZero();
        for(int l = 0; l < (int)optical_layer.numLayers(); ++l) {
            auto& layer = optical_layer.layer(optical_layer.numLayers() - l - 1);
            auto& layer_solution = layer.solution(m);
            int nstr = this->m_config->num_do_streams();

            const auto& eigval = layer_solution.value.dual_eigval();
            const auto& average_secant = layer.dual_average_secant();
            const sasktran_disco::LayerDual<double>& dual_thickness = layer.dual_thickness();
            double x = layer.dual_thickness().value / 2;


            const auto& transmission = layer.dual_beamTransmittance(sasktran_disco::Location::CEILING, optical_layer.inputDerivatives());

            using MatrixView = Eigen::Map<Eigen::MatrixXd>;
            using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;

            ConstMatrixView homog_plus_matrix(layer_solution.value.dual_homog_plus().value.data(), nstr/2 * NSTOKES, nstr/2 * NSTOKES);
            ConstMatrixView homog_minus_matrix(layer_solution.value.dual_homog_minus().value.data(), nstr/2 * NSTOKES, nstr/2 * NSTOKES);


            auto upwelling_seq = Eigen::seq(l*nstr, l*nstr + nstr/2 - 1);
            auto downwelling_seq = Eigen::seq(l*nstr + nstr/2, l*nstr + nstr - 1);

            // Since we are at layer midpoint, exp(-x) and exp(-(od - x)) are equal and are both exp(-od/2)

            for(int i = 0; i < nstr/2 * NSTOKES; ++i) {
                Eigen::Map<const Eigen::Matrix<double, NSTOKES, -1>> homog_plus(homog_plus_matrix(Eigen::all, i).data(), NSTOKES, nstr/2);
                Eigen::Map<const Eigen::Matrix<double, NSTOKES, -1>> homog_minus(homog_minus_matrix(Eigen::all, i).data(), NSTOKES, nstr/2);

                double atten_factor = exp(-layer.dual_thickness().value / 2 * layer_solution.value.dual_eigval().value(i));

                soln(upwelling_seq, Eigen::all) += atten_factor * homog_plus.transpose() * layer_solution.boundary.L_coeffs.value(i);
                soln(upwelling_seq, Eigen::all) += atten_factor * homog_minus.transpose() * layer_solution.boundary.M_coeffs.value(i);

                soln(downwelling_seq, Eigen::all) += atten_factor * homog_minus.transpose() * layer_solution.boundary.L_coeffs.value(i);
                soln(downwelling_seq, Eigen::all) += atten_factor * homog_plus.transpose() * layer_solution.boundary.M_coeffs.value(i);

                double dp = transmission.value / (eigval.value(i) + average_secant.value) * (exp(-x*average_secant.value) - exp(-dual_thickness.value * average_secant.value) * exp(-(dual_thickness.value - x) * eigval.value(i))); //CMinus
                double dm = transmission.value / (average_secant.value - eigval.value(i)) * (exp(-x*eigval.value(i)) - exp(-x*average_secant.value)); // Cplus

                soln(upwelling_seq, Eigen::all) += dm * homog_plus.transpose() * layer_solution.value.dual_green_A_plus().value(i);
                soln(upwelling_seq, Eigen::all) += dp * homog_minus.transpose() * layer_solution.value.dual_green_A_minus().value(i);

                soln(downwelling_seq, Eigen::all) += dm * homog_minus.transpose() * layer_solution.value.dual_green_A_plus().value(i);
                soln(downwelling_seq, Eigen::all) += dp * homog_plus.transpose() * layer_solution.value.dual_green_A_minus().value(i);
            }
        }
    }

    SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(DOSourceSphericallyCorrected);
}