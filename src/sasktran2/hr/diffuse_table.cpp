#include <sasktran2/hr/diffuse_source.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/solartransmission.h>
#include <sasktran2/do_source.h>

#include <fstream>

namespace sasktran2::hr {

    template<int NSTOKES>
    DiffuseTable<NSTOKES>::DiffuseTable(const sasktran2::raytracing::RayTracerBase& ray_tracer,
                                        const sasktran2::Geometry1D& geometry) :
                                        m_raytracer(ray_tracer),
                                        m_geometry(geometry),
                                        m_integrator(false)
                                        {
    }

    template<int NSTOKES>
    sasktran2::grids::Grid DiffuseTable<NSTOKES>::generate_cos_sza_grid(double min_cos_sza, double max_cos_sza) {
        Eigen::VectorXd cos_sza_grid_values;

        // TODO: should we have separate SZA spacings for DO/HR? Or harmonize the naming?
        if(m_config->num_do_sza() > 1) {
            cos_sza_grid_values.setLinSpaced(m_config->num_do_sza(), min_cos_sza, max_cos_sza);
        } else {
            // TODO: Should this be a reference cos_sza value rather than the middle calculated like this?
            // Maybe the mean tangent point should always be included?
            cos_sza_grid_values.resize(1);
            cos_sza_grid_values.setConstant(m_geometry.coordinates().sun_unit().z());
        }

        return sasktran2::grids::Grid(std::move(cos_sza_grid_values), sasktran2::grids::gridspacing::constant, sasktran2::grids::outofbounds::extend,
                                                sasktran2::grids::interpolation::linear);
    }

    template<int NSTOKES>
    sasktran2::grids::AltitudeGrid DiffuseTable<NSTOKES>::generate_altitude_grid() {
        // TODO: Decouple altitude grid here from the global geometry grid ?

        Eigen::VectorXd alt_values = (m_geometry.altitude_grid().grid()(Eigen::seq(0, Eigen::last - 1)) + m_geometry.altitude_grid().grid()(Eigen::seq(1, Eigen::last))) / 2.0;

        alt_values.setLinSpaced(20, 500, 80000);

        return sasktran2::grids::AltitudeGrid(std::move(alt_values),sasktran2::grids::gridspacing::constant, sasktran2::grids::outofbounds::extend,
                                                        sasktran2::grids::interpolation::linear);
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::construct_diffuse_points() {
        // Create the sphere pairs
        m_unit_sphere_pairs.resize(1);

        // TODO: Get npoints from config, figure out if we want to use more than one type of sphere
        // TODO: How to get number of legendre coefficients?
        m_unit_sphere_pairs[0] = std::make_unique<sasktran2::hr::IncomingOutgoingSpherePair<NSTOKES>>(
                m_config->num_do_streams(),
                std::move(std::make_unique<sasktran2::math::LebedevSphere>(m_config->num_hr_incoming())),
                std::move(std::make_unique<sasktran2::math::LebedevSphere>(m_config->num_hr_outgoing())));


        m_diffuse_points.resize(m_location_interpolator->num_interior_points());

        sasktran2::Location loc;

        for(int i = 0; i < m_diffuse_points.size(); ++i) {
            auto& point = m_diffuse_points[i];

            loc.position = m_location_interpolator->grid_location(m_geometry.coordinates(), i);

            point = std::make_unique<sasktran2::hr::DiffusePoint<NSTOKES>>(*m_unit_sphere_pairs[0], loc);

        }

        // Construct interpolators to the diffuse point locations
        m_diffuse_point_interpolation_weights.resize(m_diffuse_points.size());
        for(int i = 0; i < m_diffuse_points.size(); ++i) {
            auto& point = m_diffuse_points[i];

            m_geometry.assign_interpolation_weights(point->location(), m_diffuse_point_interpolation_weights[i]);
        }

        m_diffuse_incoming_index_map.resize(m_diffuse_points.size());
        m_diffuse_outgoing_index_map.resize(m_diffuse_points.size());

        int start_incoming_idx = 0;
        int start_outgoing_idx = 0;
        for(int i = 0; i < m_diffuse_points.size(); ++i) {
            m_diffuse_incoming_index_map[i] = start_incoming_idx;
            m_diffuse_outgoing_index_map[i] = start_outgoing_idx;
            start_incoming_idx += m_diffuse_points[i]->num_incoming();
            start_outgoing_idx += m_diffuse_points[i]->num_outgoing();
        }

        m_incoming_traced_rays.resize(start_incoming_idx);

        for(auto& storage : m_thread_storage) {
            storage.m_incoming_radiances.resize(start_incoming_idx * NSTOKES, 0, false);
            storage.m_firstorder_radiances.resize(start_incoming_idx * NSTOKES, 0, false);
            storage.m_outgoing_sources.resize(start_outgoing_idx * NSTOKES, 0, false);
            storage.accumulation_matrix.resize(start_incoming_idx * NSTOKES, start_outgoing_idx * NSTOKES);

            storage.point_scattering_matrices.resize(m_diffuse_points.size());
            for(int i = 0; i < m_diffuse_points.size(); ++i) {
                storage.point_scattering_matrices[i].resize(m_diffuse_points[i]->num_outgoing() * NSTOKES, m_diffuse_points[i]->num_incoming() * NSTOKES);
            }
        }
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::trace_incoming_rays() {
        sasktran2::viewinggeometry::ViewingRay viewing_ray;

        for(int i = 0; i < m_diffuse_points.size(); ++i) {
            viewing_ray.observer = m_diffuse_points[i]->location();
            for(int j = 0; j < m_diffuse_points[i]->num_incoming(); ++j) {
                viewing_ray.look_away = m_diffuse_points[i]->incoming_direction(j);

                m_raytracer.trace_ray(viewing_ray,  m_incoming_traced_rays[m_diffuse_incoming_index_map[i] + j]);
            }
        }
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay> &los_rays) {
        // TODO: This is the 1D case, need to add separate logic for 2d/3d

        // find the min/max SZA from the LOS rays and generate the cos_sza_grid
        std::pair<double, double> min_max_cos_sza = sasktran2::raytracing::min_max_cos_sza_of_all_rays(los_rays);

        // create the location interpolator
        m_location_interpolator = std::make_unique<sasktran2::grids::AltitudeSZASourceLocationInterpolator>(generate_altitude_grid(), generate_cos_sza_grid(min_max_cos_sza.first, min_max_cos_sza.second));

        // Construct the actual diffuse points
        construct_diffuse_points();

        // Trace all of the incoming rays
        trace_incoming_rays();

        // Set up the integrator
        m_integrator.initialize_geometry(m_incoming_traced_rays, this->m_geometry);
        // And the initial sources
        for(auto& source : m_initial_owned_sources) {
            source->initialize_geometry(m_incoming_traced_rays);
        }

        int temp;
        generate_source_interpolation_weights(m_incoming_traced_rays, m_diffuse_source_weights, m_total_num_diffuse_weights);

        generate_source_interpolation_weights(los_rays, m_los_source_weights, temp);

        // Have to create a vector of all locations and directions
        std::vector<Eigen::Vector3d> locations, directions;

        for(int i = 0; i < m_diffuse_points.size(); ++i) {
            const auto& point = m_diffuse_points[i];

            for(int j = 0; j < point->num_outgoing(); ++j) {
                locations.push_back(point->location().position);
                directions.push_back(point->sphere_pair().outgoing_sphere().get_quad_position(j));
            }
        }

        if(m_config->initialize_hr_with_do()) {
            m_do_source->storage().create_location_source_interpolator(locations,
                                                                       directions,
                                                                       m_do_to_diffuse_outgoing_interpolator);
        }

      }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES> &atmosphere) {
        m_atmosphere = &atmosphere;

        m_integrator.initialize_atmosphere(atmosphere);

        for(auto& source : m_initial_owned_sources) {
            source->initialize_atmosphere(atmosphere);
        }
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::initialize_config(const sasktran2::Config &config) {
        m_config = &config;

        m_thread_storage.resize(m_config->num_threads());

        m_initial_owned_sources.emplace_back(
                std::make_unique<sasktran2::solartransmission::SingleScatterSource<sasktran2::solartransmission::SolarTransmissionTable, NSTOKES>>(m_geometry, m_raytracer)
        );

        m_initial_sources.push_back(m_initial_owned_sources[0].get());


        if(m_config->initialize_hr_with_do()) {
            m_initial_owned_sources.emplace_back(
                    std::make_unique<sasktran2::DOSourceInterpolatedPostProcessing<NSTOKES, -1>>(m_geometry, m_raytracer)
            );

            m_do_source = static_cast<DOSourceInterpolatedPostProcessing<NSTOKES, -1>*>(m_initial_owned_sources[1].get());
        }


        for(auto& source : m_initial_owned_sources) {
            source->initialize_config(config);
        }
    }

    template<int NSTOKES>
    Eigen::Vector3d DiffuseTable<NSTOKES>::rotate_unit_vector(const Eigen::Vector3d &vector,
                                                   const Eigen::Vector3d &initial_position,
                                                   const Eigen::Vector3d &new_position) const {
        /*
        sasktran2::Location temp;
        temp.position = initial_position;

        double csz_initial, saa_initial, csz_end, saa_end, saa_test, csz_test;
        sasktran2::raytracing::calculate_csz_saz(m_geometry.coordinates().sun_unit(), temp, -1*vector, csz_initial, saa_initial);

        double angle = acos(initial_position.normalized().dot(new_position.normalized()));

        Eigen::Vector3d axis = initial_position.cross(new_position);

        Eigen::AngleAxis<double> rot_transform(angle, axis.normalized());

        temp.position = new_position;

        sasktran2::raytracing::calculate_csz_saz(m_geometry.coordinates().sun_unit(), temp, -1*rot_transform.matrix() * vector, csz_end, saa_end);

        Eigen::AngleAxis<double> azi_transform(-1*(saa_initial - saa_end), new_position.normalized());

        sasktran2::raytracing::calculate_csz_saz(m_geometry.coordinates().sun_unit(), temp, -1*azi_transform.matrix() * rot_transform.matrix() * vector, csz_test, saa_test);

        if(abs(saa_test - saa_initial) > 1e-6) {
            BOOST_LOG_TRIVIAL(warning) << "csz_initial: " << csz_initial << "saa_initial: " << saa_initial << " csz_after:" << csz_test << " saa_after:" << saa_test;
        }

        return azi_transform.matrix() * rot_transform.matrix() * vector;
         */

        sasktran2::Location temp;
        temp.position = initial_position;

        double csz_initial, saa_initial;
        sasktran2::raytracing::calculate_csz_saz(m_geometry.coordinates().sun_unit(), temp, vector, csz_initial, saa_initial);

        return m_geometry.coordinates().look_vector_from_azimuth(new_position, saa_initial, temp.cos_zenith_angle(vector));
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_source_interpolation_weights(
            const std::vector<sasktran2::raytracing::TracedRay> &rays, SInterpolator &interpolator,
            int& total_num_weights
            ) const {
        total_num_weights = 0;

        interpolator.resize(rays.size());

        std::vector<std::pair<int, double>> temp_location_storage;
        std::vector<std::pair<int, double>> temp_direction_storage;

        int num_location, num_direction;

        Eigen::Vector3d rotated_los;

        sasktran2::Location temp_location;

        for(int rayidx = 0; rayidx < rays.size(); ++rayidx) {
            auto& ray_interpolator = interpolator[rayidx];
            const auto& ray = rays[rayidx];

            ray_interpolator.resize(ray.layers.size());

            for(int layeridx = 0; layeridx < ray.layers.size(); ++layeridx) {
                const auto& layer = ray.layers[layeridx];
                auto& layer_interpolator = ray_interpolator[layeridx].second;
                auto& atmosphere_interpolator = ray_interpolator[layeridx].first;

                temp_location.position = (layer.entrance.position + layer.exit.position) / 2.0;

                m_geometry.assign_interpolation_weights(temp_location, atmosphere_interpolator);


                m_location_interpolator->interior_interpolation_weights(m_geometry.coordinates(),
                                                                        temp_location,
                                                                        temp_location_storage,
                                                                        num_location);
                for(int locidx = 0; locidx < num_location; ++locidx) {
                    const auto& contributing_point = m_diffuse_points[temp_location_storage[locidx].first];

                    // Multiply by -1 to get the propagation direction
                    rotated_los = rotate_unit_vector(-1*layer.average_look_away, temp_location.position, contributing_point->location().position);

                    contributing_point->sphere_pair().outgoing_sphere().interpolate(rotated_los,
                                                                                    temp_direction_storage,
                                                                                    num_direction);

                    for(int diridx = 0; diridx < num_direction; ++diridx) {
                        #ifdef SASKTRAN_DEBUG_ASSERTS
                        if(m_diffuse_outgoing_index_map[temp_location_storage[locidx].first] + temp_direction_storage[diridx].first > m_diffuse_outgoing_index_map.back() + m_config->num_hr_outgoing()) {
                            BOOST_LOG_TRIVIAL(error) << "BAD INDEX " << temp_location_storage[locidx].first << " " << temp_direction_storage[diridx].first;
                        }
                        #endif


                        layer_interpolator.emplace_back(std::make_pair(m_diffuse_outgoing_index_map[temp_location_storage[locidx].first] + temp_direction_storage[diridx].first,
                                                                    temp_location_storage[locidx].second * temp_direction_storage[diridx].second
                                                                    ));
                    }

                }
                total_num_weights += num_location * num_direction;
            }

        }
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_scattering_matrices(int wavelidx, int threadidx) {
        auto& storage = m_thread_storage[threadidx];

        for(int i = 0; i < m_diffuse_points.size(); ++i) {
            const auto& point = m_diffuse_points[i];

            point->sphere_pair().calculate_scattering_matrix(m_atmosphere->storage().phase[wavelidx],
                                                             m_diffuse_point_interpolation_weights[i],
                                                             storage.point_scattering_matrices[i].data()
                                                             );
        }
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::generate_accumulation_matrix(int wavelidx, int threadidx) {
        auto& matrix = m_thread_storage[threadidx].accumulation_matrix;



    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::iterate_to_solution(int wavelidx, int threadidx) {
        auto& storage = m_thread_storage[threadidx];

        Eigen::VectorXd old_outgoing_vals;

        for(int spher_iter = 0; spher_iter < m_config->num_hr_spherical_iterations(); ++spher_iter) {
            // Apply the scattering matrices

            old_outgoing_vals = m_thread_storage[threadidx].m_outgoing_sources.value(Eigen::seq(0, Eigen::last, NSTOKES));

            if(spher_iter > 0 || !m_config->initialize_hr_with_do()) {
                for(int i = 0; i < m_diffuse_points.size(); ++i) {
                    const auto& point = m_diffuse_points[i];

                    auto incoming_seq = Eigen::seq(m_diffuse_incoming_index_map[i]*NSTOKES, NSTOKES*(m_diffuse_incoming_index_map[i] + point->num_incoming()) - 1);
                    auto outgoing_seq = Eigen::seq(m_diffuse_outgoing_index_map[i]*NSTOKES, NSTOKES*(m_diffuse_outgoing_index_map[i] + point->num_outgoing()) - 1);

                    storage.m_outgoing_sources.value(outgoing_seq) = storage.point_scattering_matrices[i] * storage.m_incoming_radiances.value(incoming_seq);

                    #ifdef SASKTRAN_DEBUG_ASSERTS
                    if(storage.m_outgoing_sources.value(outgoing_seq).hasNaN()) {
                        BOOST_LOG_TRIVIAL(error) << "NaN in outgoing point: " << i;
                    }
                    #endif
                }
            } else {
                m_thread_storage[threadidx].m_outgoing_sources.value = m_do_to_diffuse_outgoing_interpolator * m_do_source->storage().linear_source(threadidx).value;
            }

            if(spher_iter > 0) {
                double max_change_factor = (m_thread_storage[threadidx].m_outgoing_sources.value(Eigen::seq(0, Eigen::last, NSTOKES)).array() / old_outgoing_vals.array() - 1).abs().maxCoeff();
            }

            // Apply the accumulation matrix
            storage.m_incoming_radiances.value = storage.accumulation_matrix * storage.m_outgoing_sources.value + storage.m_firstorder_radiances.value;
        }


        /*
        std::ofstream file("/Users/dannyz/incoming.txt");
        std::ofstream file2("/Users/dannyz/outgoing.txt");
        if (file.is_open())
        {
            for(int i = 0; i < m_diffuse_points.size(); ++i) {
                for(int j = 0; j < m_config->num_hr_incoming(); ++j) {
                    auto xyz = m_diffuse_points[i]->sphere_pair().incoming_sphere().get_quad_position(j);
                    file << xyz(0) << " " << xyz(1) << " " << xyz(2);
                    for(int s = 0; s < NSTOKES; ++s) {
                        file << " " << storage.m_incoming_radiances.value[(m_diffuse_incoming_index_map[i] + j)*NSTOKES + s];
                    }
                    file << std::endl;
                }

                for(int j = 0; j < m_config->num_hr_outgoing(); ++j) {
                    auto xyz = m_diffuse_points[i]->sphere_pair().outgoing_sphere().get_quad_position(j);
                    file2 << xyz(0) << " " << xyz(1) << " " << xyz(2);
                    for(int s = 0; s < NSTOKES; ++s) {
                        file2 << " " << storage.m_outgoing_sources.value[(m_diffuse_outgoing_index_map[i] + j)*NSTOKES + s];
                    }
                    file2 << std::endl;
                }
            }
        }
         */




    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::calculate(int wavelidx, int threadidx) {
        for(auto& source : m_initial_owned_sources) {
            source->calculate(wavelidx, threadidx);
        }

        // Calculate the first order incoming signal
        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> temp_result;
        std::vector<Eigen::Triplet<double>> triplets;

        triplets.reserve(m_total_num_diffuse_weights);

        for(int rayidx = 0; rayidx < m_incoming_traced_rays.size(); ++rayidx) {
            temp_result.value.setZero();
            // m_integrator.integrate(temp_result, m_initial_sources, wavelidx, rayidx, threadidx);

            m_integrator.integrate_and_emplace_accumulation_triplets(temp_result, m_initial_sources, wavelidx, rayidx, threadidx,
                                                                     m_diffuse_source_weights,
                                                                     triplets
                                                                     );

            m_thread_storage[threadidx].m_firstorder_radiances.value(Eigen::seq(rayidx*NSTOKES, rayidx*NSTOKES + NSTOKES - 1)) = temp_result.value;

            if(temp_result.value.hasNaN()) {
                BOOST_LOG_TRIVIAL(error) << "Incoming Ray: " << rayidx << " has NaN";
            }

        }
        auto& matrix = m_thread_storage[threadidx].accumulation_matrix;

        matrix.setFromTriplets(triplets.begin(), triplets.end());


        // Optional: Intialize total incoming with DO solution

        // Else, start with total first order incoming
        m_thread_storage[threadidx].m_incoming_radiances.value = m_thread_storage[threadidx].m_firstorder_radiances.value;

        // Generate the scattering and accumulation matrices
        generate_scattering_matrices(wavelidx, threadidx);
        generate_accumulation_matrix(wavelidx, threadidx);

        // Iterate to the solution
        iterate_to_solution(wavelidx, threadidx);
    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::integrated_source(int wavelidx, int losidx, int layeridx, int threadidx,
                                                  const sasktran2::raytracing::SphericalLayer &layer,
                                                  const sasktran2::SparseODDualView &shell_od,
                                                  sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {
        auto& storage = m_thread_storage[threadidx];

        auto& interpolator = m_los_source_weights[losidx][layeridx];

        // Start by calculating ssa at the source point
        double omega = 0;
        for(int i = 0; i < interpolator.first.size(); ++i) {
            auto& index_weight = interpolator.first[i];
            omega += m_atmosphere->storage().ssa(index_weight.first, wavelidx) * index_weight.second;
        }

        double source_factor = omega * (1 - exp(-1*shell_od.od));

        for(int i = 0; i < interpolator.second.size(); ++i) {
            auto& index_weight = interpolator.second[i];

            for(int s = 0; s < NSTOKES; ++s) {
                source.value(s) += source_factor * storage.m_outgoing_sources.value(index_weight.first*NSTOKES + s) * index_weight.second;
            }
        }

    }

    template<int NSTOKES>
    void DiffuseTable<NSTOKES>::end_of_ray_source(int wavelidx, int losidx, int threadidx,
                                                  sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &source) const {

    }

    template class DiffuseTable<1>;
    template class DiffuseTable<3>;

}