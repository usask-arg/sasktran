#include <sasktran2/raytracing.h>

namespace sasktran2::raytracing {
    void PlaneParallelRayTracer::trace_ray(const sasktran2::viewinggeometry::ViewingRay &ray,
                                           TracedRay &result) const {
        // Set the ray to 0
        result.reset();

        if (ray.observer.radius() - m_earth_radius >= m_alt_grid.grid()(Eigen::last)) {
            // Outside atmosphere, probably
            if(ray.cos_viewing() > 0) {
                // We are looking up, so this is just an empty ray
                result.observer_and_look = ray;
                result.ground_is_hit = false;
                return;
            }

            trace_ray_observer_outside_looking_ground(ray, result);
        }
        else {
            // We are inside the atmosphere
            if (ray.cos_viewing() > 0) {
                // Looking up
                trace_ray_observer_inside_looking_up(ray, result);
            }
            else {
                trace_ray_observer_inside_looking_ground(ray, result);
            }
        }

        // For each layer, we go through and add in the computed solar angles
        for(auto& layer : result.layers) {
            add_solar_parameters(m_geometry.coordinates().sun_unit(), layer);
        }
    }

    void PlaneParallelRayTracer::trace_ray_observer_inside_looking_ground(
            const sasktran2::viewinggeometry::ViewingRay &ray, TracedRay &tracedray) const {
        // We have len(altitudes) - 1 spherical layers, starting from the ground
        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = true;

        result.layers.resize(m_alt_grid.grid().size());

        for (int i = 0; i < m_alt_grid.grid().size() - 1; ++i) {
            complete_layer(result.layers[i], ray, i, ViewingDirection::down);
        }
    }

    void PlaneParallelRayTracer::trace_ray_observer_inside_looking_up(const sasktran2::viewinggeometry::ViewingRay &ray,
                                                                      TracedRay &tracedray) const {
        // Find the index to the first altitude ABOVE the observer
        auto it = std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().end(), ray.observer.radius() - m_earth_radius);
        size_t start_index = std::distance(m_alt_grid.grid().begin(), it);

        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = false;

        result.layers.resize(m_alt_grid.grid().size() - start_index);

        int layer_c = 0;
        for (size_t i = m_alt_grid.grid().size() - 1; i != start_index; --i) {
            auto& layer = result.layers[layer_c];
            complete_layer(layer, ray, i, ViewingDirection::up);
            ++layer_c;
        }
        assert(layer_c == result.layers.size() - 1);
        partial_layer(result.layers[layer_c], ray, start_index, ViewingDirection::up);

        for(int i = 0; i < result.layers.size() - 1; ++i) {
            assert(result.layers[i+1].exit.radius() == result.layers[i].entrance.radius());
        }

        assert(abs(result.layers[0].exit.radius() - m_alt_grid.grid()(Eigen::last) - m_earth_radius) < 1e-8);
        assert(abs(result.layers[result.layers.size() -1].entrance.radius() - ray.observer.radius()) < 1e-8);
    }

    void PlaneParallelRayTracer::trace_ray_observer_outside_looking_ground(
            const sasktran2::viewinggeometry::ViewingRay &ray, TracedRay &tracedray) const {
        auto& result = tracedray;

        result.observer_and_look = ray;
        result.ground_is_hit = true;

        // Find the index to the first altitude ABOVE the observer
        auto it = std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().end(), ray.observer.radius() - m_earth_radius);
        size_t start_index = std::distance(m_alt_grid.grid().begin(), it);

        if(start_index == 0) {
            // Have 0 layers, most likely rounding error
            return;
        }

        result.layers.resize(start_index);

        // Complete layers from the ground
        for (size_t i = 0; i < start_index-1; ++i) {
            auto& layer = result.layers[i];
            complete_layer(layer, ray, i, ViewingDirection::down);
        }
        partial_layer(result.layers[start_index-1], ray, start_index-1, ViewingDirection::down);
    }


    void PlaneParallelRayTracer::complete_layer(SphericalLayer &layer,
                                                const sasktran2::viewinggeometry::ViewingRay &ray, size_t exit_index,
                                                ViewingDirection direction) const {
        layer.type = LayerType::complete;

        double entrance_altitude = m_alt_grid.grid()(exit_index + direction);
        double exit_altitude = m_alt_grid.grid()(exit_index);

        layer.entrance.on_exact_altitude = true;
        layer.entrance.lower_alt_index = int(exit_index + direction);

        layer.exit.on_exact_altitude = true;
        layer.exit.lower_alt_index = int(exit_index);

        double s_entrance = distance_to_altitude(ray, entrance_altitude, direction);
        double s_exit = distance_to_altitude(ray, exit_altitude, direction);

        layer.layer_distance = abs(s_entrance - s_exit);

        layer.entrance.position = ray.observer.position + ray.look_away * s_entrance;
        layer.exit.position = ray.observer.position + ray.look_away * s_exit;

        layer.curvature_factor = 1;

        layer.average_look_away = ray.look_away;

        add_od_quadrature(layer, sasktran2::geometrytype::planeparallel);
    }

    void PlaneParallelRayTracer::partial_layer(SphericalLayer &layer, const sasktran2::viewinggeometry::ViewingRay &ray,
                                               size_t start_index, ViewingDirection direction) const {
        layer.type = LayerType::partial;

        double entrance_altitude = ray.observer.radius() - m_earth_radius;
        double exit_altitude = m_alt_grid.grid()(start_index);

        layer.exit.on_exact_altitude = true;
        layer.exit.lower_alt_index = int(start_index);

        layer.entrance.on_exact_altitude = false;
        layer.entrance.lower_alt_index = direction < 0 ? int(start_index + direction) : int(start_index);

        double s_entrance = distance_to_altitude(ray, entrance_altitude, direction);
        double s_exit = distance_to_altitude(ray, exit_altitude, direction);

        layer.layer_distance = abs(s_entrance - s_exit);

        layer.entrance.position = ray.observer.position + ray.look_away * s_entrance;
        layer.exit.position = ray.observer.position + ray.look_away * s_exit;

        layer.curvature_factor = 1;
        layer.average_look_away = ray.look_away;

        add_od_quadrature(layer, sasktran2::geometrytype::planeparallel);
    }

    double
    PlaneParallelRayTracer::distance_to_altitude(const sasktran2::viewinggeometry::ViewingRay &ray, double altitude,
                                                 ViewingDirection direction) const {
        double mu_viewing = ray.cos_viewing();

        double observer_altitude = ray.observer.radius() - m_earth_radius;

        return direction * (altitude - observer_altitude) / mu_viewing;
    }

}