#include <sasktran2/solartransmission.h>

namespace sasktran2::solartransmission {
    void SolarTransmissionExact::generate_geometry_matrix(const std::vector<sasktran2::raytracing::TracedRay> &rays,
                                                          Eigen::MatrixXd &od_matrix,
                                                          std::vector<bool>& ground_hit_flag) const {
        // First calculate the number of points we need to create the matrix for
        // We calculate solar transmission at the boundaries of layers, so it is nlayer+1 for each ray
        int numpoints = 0;
        for(const auto& ray : rays) {
            numpoints += (int)ray.layers.size() + 1;
        }

        // od matrix is such that matrix @ extinction = od
        od_matrix.resize(numpoints, m_geometry.size());
        od_matrix.setZero();

        // Have to handle rays that hit the ground separately since they have no solar transmission
        ground_hit_flag.resize(numpoints, false);

        sasktran2::viewinggeometry::ViewingRay ray_to_sun;

        ray_to_sun.look_away = m_geometry.coordinates().sun_unit();

        // common memory
        std::vector<std::pair<int, double>> index_weights;

        raytracing::TracedRay traced_ray;

        int row = 0;
        for(int i = 0; i < rays.size(); ++i) {
            const auto& ray = rays[i];
            for(int j = 0; j < ray.layers.size(); ++j) {
                const auto& layer = ray.layers[j];

                if( j == 0 ) {
                    // End layer at TOA, need to use layer exit
                    ray_to_sun.observer = layer.exit;
                    m_raytracer.trace_ray(ray_to_sun, traced_ray);

                    if(!traced_ray.ground_is_hit) {
                        assign_dense_matrix_column(row, traced_ray, m_geometry, od_matrix, index_weights);
                    } else {
                        ground_hit_flag[row] = true;
                    }
                    ++row;
                }

                ray_to_sun.observer = layer.entrance;
                m_raytracer.trace_ray(ray_to_sun, traced_ray);

                if(!traced_ray.ground_is_hit) {
                    assign_dense_matrix_column(row, traced_ray, m_geometry, od_matrix, index_weights);
                } else {
                    ground_hit_flag[row] = true;
                }

                ++row;
            }
        }
    }
}