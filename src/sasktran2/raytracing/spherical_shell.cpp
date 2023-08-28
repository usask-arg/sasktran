#include <sasktran2/raytracing.h>

namespace sasktran2::raytracing {
    void SphericalShellRayTracer::trace_ray(const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& result) const {
        // Set the ray to 0
        result.reset();

        // Calculate the tangent point details in a straight geometry
        double rt = ray.observer.radius() * sqrt(1 - ray.cos_viewing()*ray.cos_viewing());
        double tangent_altitude = rt - m_earth_radius;

        // If the tangent altitude is greater than the TOA altitude then we have an empty ray
        if (tangent_altitude >= m_alt_grid.grid()(Eigen::last)) {
            // Empty ray
            result.observer_and_look = ray;
            result.ground_is_hit = false;
        }

        if (ray.observer.radius() - m_earth_radius >= m_alt_grid.grid()(Eigen::last)) {
            // Outside atmosphere, probably
            if(ray.cos_viewing() > 0) {
                // We are looking up, so this is just an empty ray
                result.observer_and_look = ray;
                result.ground_is_hit = false;
                return;
            }

            // There are two cases, limb viewing and ground viewing
            if (tangent_altitude > m_alt_grid.grid()(0)) {
                // Limb viewing
                trace_ray_observer_outside_limb_viewing(ray, result);
            }
            else {
                // Ground hitting
                trace_ray_observer_outside_ground_viewing(ray, result);
            }
        }
        else {
            // We are inside the atmosphere
            if (ray.cos_viewing() > 0) {
                // Looking up, not limb viewing
                trace_ray_observer_inside_looking_up(ray, result);
            }
            else {
                // Looking downwards, could either be limb viewing or ground hitting
                if (tangent_altitude > m_alt_grid.grid()(0)) {
                    // Limb viewing
                    trace_ray_observer_inside_looking_limb(ray, result);
                }
                else {
                    // Ground hitting
                    trace_ray_observer_inside_looking_ground(ray, result);
                }
            }
        }

        // For each layer, we go through and add in the computed solar angles
        for(auto& layer : result.layers) {
            add_solar_parameters(m_geometry.coordinates().sun_unit(), layer);
        }
    }

    void SphericalShellRayTracer::trace_ray_observer_outside_ground_viewing(const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& tracedray) const {
        // We have len(altitudes) - 1 spherical layers, starting from the ground
        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = true;

        result.layers.resize(m_alt_grid.grid().size());

        for (int i = 0; i < m_alt_grid.grid().size() - 1; ++i) {
            complete_layer(result.layers[i], ray, i, ViewingDirection::down, TangentSide::nearside);
        }
    }

    void SphericalShellRayTracer::trace_ray_observer_outside_limb_viewing(const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& tracedray) const {
        double rt = ray.observer.radius() * sqrt(1 - ray.cos_viewing()*ray.cos_viewing());
        double tangent_altitude = rt - m_earth_radius;

        // Find the index to the first altitude ABOVE the tangent altitude
        auto it = std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().cend(), tangent_altitude);
        size_t above_tangent_idx = std::distance(m_alt_grid.grid().begin(), it);

        auto& result = tracedray;
        result.observer_and_look = ray;
        result.ground_is_hit = false;

        size_t numlayer = 2 * (m_alt_grid.grid().size() - above_tangent_idx);
        result.layers.resize(numlayer);

        if(numlayer == 0) {
            // Empty ray
            return;
        }

        size_t layer_c = 0;
        // We have complete layers from the TOA down to the tangent layer
        for (size_t i = m_alt_grid.grid().size()- 1; i != above_tangent_idx; --i) {
            complete_layer(result.layers[layer_c], ray, i, ViewingDirection::up, TangentSide::farside);
            ++layer_c;
        }

        tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::up, TangentSide::farside);
        ++layer_c;
        tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::down, TangentSide::nearside);
        ++layer_c;

        // And the complete layers from the tangent layer up to the instrument
        for (int i = (int)above_tangent_idx; i < (int)m_alt_grid.grid().size() - 1; ++i) {
            complete_layer(result.layers[layer_c], ray, i, ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        }
        assert(layer_c == result.layers.size());
    }

    void SphericalShellRayTracer::complete_layer(SphericalLayer& layer, const sasktran2::viewinggeometry::ViewingRay& ray, size_t exit_index, ViewingDirection direction, TangentSide side) const
    {
        layer.type = LayerType::complete;

        double entrance_altitude = m_alt_grid.grid()(exit_index + direction);
        double exit_altitude = m_alt_grid.grid()(exit_index);

        layer.entrance.on_exact_altitude = true;
        layer.entrance.lower_alt_index = int(exit_index + direction);

        layer.exit.on_exact_altitude = true;
        layer.exit.lower_alt_index = int(exit_index);

        double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
        double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

        layer.layer_distance = abs(s_entrance - s_exit);

        layer.entrance.position = ray.observer.position + ray.look_away * s_entrance;
        layer.exit.position = ray.observer.position + ray.look_away * s_exit;

        layer.curvature_factor = 1;

        layer.average_look_away = ray.look_away;

        add_od_quadrature(layer);
        add_interpolation_weights(layer, m_geometry);
    }

    void SphericalShellRayTracer::partial_layer(SphericalLayer& layer, const sasktran2::viewinggeometry::ViewingRay& ray, size_t start_index, ViewingDirection direction, TangentSide side) const {
        layer.type = LayerType::partial;

        double entrance_altitude = ray.observer.radius() - m_earth_radius;
        double exit_altitude = m_alt_grid.grid()(start_index);

        layer.exit.on_exact_altitude = true;
        layer.exit.lower_alt_index = int(start_index);

        layer.entrance.on_exact_altitude = false;
        layer.entrance.lower_alt_index = direction < 0 ? int(start_index + direction) : int(start_index);

        double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
        double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

        layer.layer_distance = abs(s_entrance - s_exit);

        layer.entrance.position = ray.observer.position + ray.look_away * s_entrance;
        layer.exit.position = ray.observer.position + ray.look_away * s_exit;

        layer.curvature_factor = 1;
        layer.average_look_away = ray.look_away;

        add_od_quadrature(layer);
        add_interpolation_weights(layer, m_geometry);
    }

    void SphericalShellRayTracer::tangent_layer(SphericalLayer& layer, const sasktran2::viewinggeometry::ViewingRay& ray, size_t upper_index, double tangent_altitude, ViewingDirection direction, TangentSide side) const {
        double entrance_altitude, exit_altitude;

        layer.type = LayerType::tangent;

        if (direction == ViewingDirection::up) {
            entrance_altitude = tangent_altitude;
            exit_altitude = m_alt_grid.grid()(upper_index);

            layer.exit.on_exact_altitude = true;
            layer.exit.lower_alt_index = int(upper_index);

            layer.entrance.on_exact_altitude = false;
            layer.entrance.lower_alt_index = int(upper_index-1);
        }
        else {
            exit_altitude = tangent_altitude;
            entrance_altitude = m_alt_grid.grid()(upper_index);

            layer.entrance.on_exact_altitude = true;
            layer.entrance.lower_alt_index = int(upper_index);

            layer.exit.on_exact_altitude = false;
            layer.exit.lower_alt_index = int(upper_index-1);
        }

        double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
        double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

        layer.layer_distance = abs(s_entrance - s_exit);

        layer.entrance.position = ray.observer.position + ray.look_away * s_entrance;
        layer.exit.position = ray.observer.position + ray.look_away * s_exit;

        layer.curvature_factor = 1;
        layer.average_look_away = ray.look_away;

        add_od_quadrature(layer);
        add_interpolation_weights(layer, m_geometry);
    }


    void SphericalShellRayTracer::partial_tangent_layer(SphericalLayer& layer, const sasktran2::viewinggeometry::ViewingRay& ray, size_t start_index, double tangent_altitude, ViewingDirection direction, TangentSide side) const {
        double entrance_altitude, exit_altitude;

        layer.type = LayerType::tangent;

        if (direction == ViewingDirection::up) {
            BOOST_LOG_TRIVIAL(error) << "Trying to construct a partial tangent layer looking up, this shouldn't be a thing";
            throw std::runtime_error("critical raytracing error");
        }
        else {
            exit_altitude = tangent_altitude;
            entrance_altitude = ray.observer.radius() - m_earth_radius;

            layer.entrance.on_exact_altitude = false;
            layer.entrance.lower_alt_index = int(start_index-1);

            layer.exit.on_exact_altitude = false;
            layer.exit.lower_alt_index = int(start_index-1);
        }

        double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
        double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

        layer.layer_distance = abs(s_entrance - s_exit);

        layer.entrance.position = ray.observer.position + ray.look_away * s_entrance;
        layer.exit.position = ray.observer.position + ray.look_away * s_exit;

        layer.curvature_factor = 1;
        layer.average_look_away = ray.look_away;

        add_od_quadrature(layer);
        add_interpolation_weights(layer, m_geometry);
    }


    void SphericalShellRayTracer::trace_ray_observer_inside_looking_up(const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& tracedray) const {
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
            complete_layer(layer, ray, i, ViewingDirection::up, TangentSide::nearside);
            ++layer_c;
        }
        assert(layer_c == result.layers.size() - 1);
        partial_layer(result.layers[layer_c], ray, start_index, ViewingDirection::up, TangentSide::nearside);

        for(int i = 0; i < result.layers.size() - 1; ++i) {
            assert(result.layers[i+1].exit.radius() == result.layers[i].entrance.radius());
        }

        assert(abs(result.layers[0].exit.radius() - m_alt_grid.grid()(Eigen::last) - m_earth_radius) < 1e-8);
        assert(abs(result.layers[result.layers.size() -1].entrance.radius() - ray.observer.radius()) < 1e-8);
    }

    void SphericalShellRayTracer::trace_ray_observer_inside_looking_ground(const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& tracedray) const {
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
            complete_layer(layer, ray, i, ViewingDirection::down, TangentSide::nearside);
        }
        partial_layer(result.layers[start_index-1], ray, start_index-1, ViewingDirection::down, TangentSide::nearside);
    }
    void SphericalShellRayTracer::trace_ray_observer_inside_looking_limb(const sasktran2::viewinggeometry::ViewingRay& ray, TracedRay& tracedray) const {
        auto& result = tracedray;

        result.observer_and_look = ray;
        result.ground_is_hit = false;

        // Find the index to the first altitude ABOVE the observer
        auto it = std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().end(), ray.observer.radius() - m_earth_radius);
        size_t observer_idx = std::distance(m_alt_grid.grid().begin(), it);

        double rt = ray.observer.radius() * sqrt(1 - ray.cos_viewing()*ray.cos_viewing());
        double tangent_altitude = rt - m_earth_radius;

        // Find the index to the first altitude ABOVE the tangent altitude
        auto it_tangent = std::upper_bound(m_alt_grid.grid().begin(), m_alt_grid.grid().cend(), tangent_altitude);
        size_t above_tangent_idx = std::distance(m_alt_grid.grid().begin(), it_tangent);

        // We have (len(alt) - above_tangent) full layers
        // +2 tangent layers
        // + (observer - above_tangent) layers on the opposite side


        // Above the tangent, far side
        int num_layers = int(m_alt_grid.grid().size() - 1) - int(above_tangent_idx);

        // Tangent layers
        num_layers += 2;

        // Layers on the far side
        if(observer_idx > above_tangent_idx) {
            num_layers += int(observer_idx - above_tangent_idx);
        }

        result.layers.resize(num_layers);

        size_t layer_c = 0;
        // We have complete layers from the TOA down to the tangent layer
        for (size_t i = m_alt_grid.grid().size()- 1; i != above_tangent_idx; --i) {
            complete_layer(result.layers[layer_c], ray, i, ViewingDirection::up, TangentSide::farside);
            ++layer_c;
        }

        // Then always one tangent layer
        tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::up, TangentSide::farside);
        ++layer_c;

        // We have a special case if above_tangent_idx == observer_idx, then we have a partial tangent layer
        if(above_tangent_idx == observer_idx) {
            partial_tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        } else {
            tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        }

        // Complete layers from the tangent to the observer
        for (size_t i = above_tangent_idx; i < observer_idx - 1; ++i) {
            auto& layer = result.layers[layer_c];
            complete_layer(layer, ray, i, ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        }

        // Final partial layer, only if we didn't have a tangent partial layer
        if(observer_idx > above_tangent_idx) {
            partial_layer(result.layers[layer_c], ray, observer_idx-1, ViewingDirection::down, TangentSide::nearside);
            ++layer_c;
        }

        assert(layer_c == result.layers.size());
    }


}