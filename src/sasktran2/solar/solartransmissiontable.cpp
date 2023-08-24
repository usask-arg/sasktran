#include <sasktran2/solartransmission.h>

namespace sasktran2::solartransmission {

    void
    SolarTransmissionTable::initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay> &integration_rays) {
        // find the min/max SZA from the LOS rays and generate the cos_sza_grid
        std::pair<double, double> min_max_cos_sza = sasktran2::raytracing::min_max_cos_sza_of_all_rays(integration_rays);

        Eigen::VectorXd cos_sza_grid_values;

        // TODO: configure this resolution
        cos_sza_grid_values.setLinSpaced(100, min_max_cos_sza.first, min_max_cos_sza.second);

        Eigen::VectorXd alt_values = m_geometry.altitude_grid().grid();

        // create the location interpolator
        m_location_interpolator = std::make_unique<sasktran2::grids::AltitudeSZASourceLocationInterpolator>(sasktran2::grids::AltitudeGrid(std::move(alt_values),sasktran2::grids::gridspacing::constant, sasktran2::grids::outofbounds::extend,
                                                                                                                                           sasktran2::grids::interpolation::linear),
                                                                                                            sasktran2::grids::Grid(std::move(cos_sza_grid_values), sasktran2::grids::gridspacing::constant, sasktran2::grids::outofbounds::extend,
                                                                                                                                   sasktran2::grids::interpolation::linear));

        // Construct the matrix that calculates OD on the solar transmission table locations
        // i.e. solar_od_on_grid = matrix @ extinction
        m_geometry_matrix.resize(m_location_interpolator->num_interior_points(), m_geometry.size());
        m_geometry_matrix.setZero();

        m_ground_hit_flag.resize(m_location_interpolator->num_interior_points());

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;

        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();

        // common memory
        std::vector<std::pair<int, double>> index_weights;

        raytracing::TracedRay traced_ray;

        for(int i = 0; i < m_location_interpolator->num_interior_points(); ++i) {
            ray_to_sun.observer.position = m_location_interpolator->grid_location(m_geometry.coordinates(), i);

            m_raytracer.trace_ray(ray_to_sun, traced_ray);

            if(!traced_ray.ground_is_hit) {
                assign_dense_matrix_column(i, traced_ray, m_geometry, m_geometry_matrix, index_weights);
                m_ground_hit_flag[i] = false;
            } else {
                m_ground_hit_flag[i] = true;
            }
        }
    }

    void SolarTransmissionTable::generate_interpolation_matrix(
            const std::vector<sasktran2::raytracing::TracedRay> &rays, Eigen::SparseMatrix<double, Eigen::RowMajor> &interpolator,
            std::vector<bool> &ground_hit_flag) const {
        // First calculate the number of points we need to create the matrix for
        // We calculate solar transmission at the boundaries of layers, so it is nlayer+1 for each ray
        int numpoints = 0;
        for(const auto& ray : rays) {
            numpoints += (int)ray.layers.size() + 1;
        }

        // od matrix is such that matrix @ extinction = od
        interpolator.resize(numpoints, m_location_interpolator->num_interior_points());

        // Have to handle rays that hit the ground separately since they have no solar transmission
        ground_hit_flag.resize(numpoints, false);

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        std::vector<std::pair<int, double>> interpolator_storage;
        int num_interp;

        int row = 0;
        for(int i = 0; i < rays.size(); ++i) {
            const auto& ray = rays[i];
            for(int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];

                if( j == 0 ) {
                    // End layer at TOA, need to use layer exit
                    m_location_interpolator->interior_interpolation_weights(m_geometry.coordinates(), layer.exit, interpolator_storage, num_interp);

                    for(int k = 0; k < num_interp; ++k) {
                        tripletList.emplace_back(T(row, interpolator_storage[k].first, interpolator_storage[k].second));
                    }
                    ++row;
                }

                m_location_interpolator->interior_interpolation_weights(m_geometry.coordinates(), layer.entrance, interpolator_storage, num_interp);
                for(int k = 0; k < num_interp; ++k) {
                    tripletList.emplace_back(T(row, interpolator_storage[k].first, interpolator_storage[k].second));
                }
                ++row;
            }
        }
        interpolator.setFromTriplets(tripletList.begin(), tripletList.end());
    }
}