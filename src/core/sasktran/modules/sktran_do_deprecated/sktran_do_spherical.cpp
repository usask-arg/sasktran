#include "modules/sktran_do_deprecated/include/sktran_do_spherical.h"

namespace sktran_do_detail {
	TracedRay SphericalRayTracer::trace_ray(const ViewingRay& ray) const {
		double rt = ray.observer.Radius() * sqrt(1 - ray.cos_viewing()*ray.cos_viewing());
		double tangent_altitude = rt - m_earth_radius;

		if (tangent_altitude > m_altitudes[m_altitudes.size() - 1]) {
			// Empty ray
			TracedRay result;
			result.observer_and_look = ray;
			result.ground_is_hit = false;

			return result;
		}

		if (ray.observer.Altitude() >= m_altitudes[m_altitudes.size() - 1]) {
			// Outside atmosphere

			// There are two cases, limb viewing and ground viewing
			if (tangent_altitude > m_altitudes[0]) {
				// Limb viewing
				return trace_ray_observer_outside_limb_viewing(ray);
			}
			else {
				// Ground hitting
				return trace_ray_observer_outside_ground_viewing(ray);
			}
		}
		else {
			if (ray.cos_viewing() > 0) {
				// Looking up, not limb viewing
				return trace_ray_observer_inside_looking_up(ray);
			}
			else {
				// Looking downwards, could either be limb viewing or ground hitting
				if (tangent_altitude > m_altitudes[0]) {
					// Limb viewing
					return trace_ray_observer_inside_looking_limb(ray);
				}
				else {
					// Ground hitting
					return trace_ray_observer_inside_looking_ground(ray);
				}
			}
		}
	}

	void SphericalRayTracer::trace_solar_rays(TracedRay& ray) const {
		ViewingRay to_sun;
		to_sun.look_away.SetCoords(0, 0, 1); // Sun is always in z direction for heliodetic

		ray.solar_entrance_traced_layers.resize(ray.layers.size());
		ray.solar_exit_traced_layers.resize(ray.layers.size());

		for (int i = 0; i < ray.layers.size(); ++i) {
			const auto& ray_layer = ray.layers[i];
			to_sun.observer = ray_layer.entrance;
			ray.solar_entrance_traced_layers[i] = trace_ray(to_sun).layers;

			to_sun.observer = ray_layer.exit;
			ray.solar_exit_traced_layers[i] = trace_ray(to_sun).layers;
		}
	}

	TracedRay SphericalRayTracer::trace_ray_observer_outside_ground_viewing(const ViewingRay& ray) const {
		// We have len(altitudes) - 1 spherical layers, starting from the ground
		TracedRay result;
		result.observer_and_look = ray;
		result.ground_is_hit = true;
	
		result.layers.resize(m_altitudes.size() - 1);

		for (size_t i = 0; i < m_altitudes.size() - 1; ++i) {
			complete_layer(result.layers[i], ray, i, ViewingDirection::down, TangentSide::nearside);
		}

		return result;
	}

	TracedRay SphericalRayTracer::trace_ray_observer_outside_limb_viewing(const ViewingRay& ray) const {
		double rt = ray.observer.Radius() * sqrt(1 - ray.cos_viewing()*ray.cos_viewing());
		double tangent_altitude = rt - m_earth_radius;

		// Find the index to the first altitude ABOVE the tangent altitude
		auto it = std::upper_bound(m_altitudes.cbegin(), m_altitudes.cend(), tangent_altitude);
		size_t above_tangent_idx = std::distance(m_altitudes.cbegin(), it);

		TracedRay result;
		result.observer_and_look = ray;
		result.ground_is_hit = false;

		size_t numlayer = 2 * (m_altitudes.size() - above_tangent_idx);
		result.layers.resize(numlayer);

		size_t layer_c = 0;
		// We have complete layers from the TOA down to the tangent layer
		for (size_t i = m_altitudes.size() - 1; i != above_tangent_idx; --i) {
			complete_layer(result.layers[layer_c], ray, i, ViewingDirection::up, TangentSide::farside);
			++layer_c;
		}

		tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::up, TangentSide::farside);
		++layer_c;
		tangent_layer(result.layers[layer_c], ray, above_tangent_idx, tangent_altitude, ViewingDirection::down, TangentSide::nearside);
		++layer_c;

		// And the complete layers from the tangent layer up to the instrument
		for (size_t i = above_tangent_idx; i < m_altitudes.size() - 1; ++i) {
			complete_layer(result.layers[layer_c], ray, i, ViewingDirection::down, TangentSide::nearside);
			++layer_c;
		}

		return result;
	}

	void SphericalRayTracer::complete_layer(SphericalLayer& layer, const ViewingRay& ray, size_t exit_index, ViewingDirection direction, TangentSide side) const
	{
		layer.type = LayerType::complete;

		double entrance_altitude = m_altitudes[exit_index + direction];
		double exit_altitude = m_altitudes[exit_index];

		layer.entrance_weights[0].first = exit_index + direction;
		layer.entrance_weights[0].second = 1;
		layer.num_entrance_weights = 1;

		layer.exit_weights[0].first = exit_index;
		layer.exit_weights[0].second = 1;
		layer.num_exit_weights = 1;

		double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
		double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

		layer.layer_distance = abs(s_entrance - s_exit);

		HELIODETIC_VECTOR entrance = ray.observer.Vector();
		entrance.SetCoords(entrance.X() + s_entrance * ray.look_away.X(), entrance.Y() + s_entrance * ray.look_away.Y(), entrance.Z() + s_entrance*ray.look_away.Z());

		HELIODETIC_VECTOR exit = ray.observer.Vector();
		exit.SetCoords(exit.X() + s_exit * ray.look_away.X(), exit.Y() + s_exit * ray.look_away.Y(), exit.Z() + s_exit * ray.look_away.Z());

		layer.entrance.FromVector(entrance, &m_coords);
		layer.exit.FromVector(exit, &m_coords);

		add_od_quadrature(layer, ray);
		add_solar_parameters(layer, ray);
	}

	void SphericalRayTracer::partial_layer(SphericalLayer& layer, const ViewingRay& ray, size_t start_index, ViewingDirection direction, TangentSide side) const {
		layer.type = LayerType::partial;

		double entrance_altitude = ray.observer.Altitude();
		double exit_altitude = m_altitudes[start_index];

		layer.exit_weights[0].first = start_index;
		layer.exit_weights[0].second = 1;
		layer.num_exit_weights = 1;

		double width = abs(m_altitudes[start_index] - m_altitudes[start_index + direction]);
		layer.entrance_weights[0].first = start_index;
		layer.entrance_weights[0].second = 1 - abs(m_altitudes[start_index] - entrance_altitude) / width;

		layer.entrance_weights[1].first = start_index + direction;
		layer.entrance_weights[1].second = 1 - layer.entrance_weights[0].second;

		layer.num_entrance_weights = 2;

		double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
		double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

		layer.layer_distance = abs(s_entrance - s_exit);

		HELIODETIC_VECTOR entrance = ray.observer.Vector();
		entrance.SetCoords(entrance.X() + s_entrance * ray.look_away.X(), entrance.Y() + s_entrance * ray.look_away.Y(), entrance.Z() + s_entrance * ray.look_away.Z());

		HELIODETIC_VECTOR exit = ray.observer.Vector();
		exit.SetCoords(exit.X() + s_exit * ray.look_away.X(), exit.Y() + s_exit * ray.look_away.Y(), exit.Z() + s_exit * ray.look_away.Z());

		layer.entrance.FromVector(entrance, &m_coords);
		layer.exit.FromVector(exit, &m_coords);

		add_od_quadrature(layer, ray);
		add_solar_parameters(layer, ray);
	}

	void SphericalRayTracer::tangent_layer(SphericalLayer& layer, const ViewingRay& ray, size_t upper_index, double tangent_altitude, ViewingDirection direction, TangentSide side) const {
		double entrance_altitude, exit_altitude;

		layer.type = LayerType::tangent;

		if (direction == ViewingDirection::up) {
			entrance_altitude = tangent_altitude;
			exit_altitude = m_altitudes[upper_index];

			layer.exit_weights[0].first = upper_index;
			layer.exit_weights[0].second = 1;
			layer.num_exit_weights = 1;

			double width = abs(m_altitudes[upper_index] - m_altitudes[upper_index - 1]);
			layer.entrance_weights[0].first = upper_index;
			layer.entrance_weights[0].second = abs(m_altitudes[upper_index - 1] - entrance_altitude) / width;

			layer.entrance_weights[1].first = upper_index - 1;
			layer.entrance_weights[1].second = 1 - layer.entrance_weights[0].second;

			layer.num_entrance_weights = 2;
		}
		else {
			exit_altitude = tangent_altitude;
			entrance_altitude = m_altitudes[upper_index];

			layer.entrance_weights[0].first = upper_index;
			layer.entrance_weights[0].second = 1;
			layer.num_entrance_weights = 1;

			double width = abs(m_altitudes[upper_index] - m_altitudes[upper_index - 1]);
			layer.exit_weights[0].first = upper_index;
			layer.exit_weights[0].second = 1 - abs(m_altitudes[upper_index] - exit_altitude) / width;

			layer.exit_weights[1].first = upper_index - 1;
			layer.exit_weights[1].second = 1 - layer.exit_weights[0].second;

			layer.num_exit_weights = 2;
		}

		double s_entrance = distance_to_altitude(ray, entrance_altitude, direction, side);
		double s_exit = distance_to_altitude(ray, exit_altitude, direction, side);

		layer.layer_distance = abs(s_entrance - s_exit);

		HELIODETIC_VECTOR entrance = ray.observer.Vector();
		entrance.SetCoords(entrance.X() + s_entrance * ray.look_away.X(), entrance.Y() + s_entrance * ray.look_away.Y(), entrance.Z() + s_entrance * ray.look_away.Z());

		HELIODETIC_VECTOR exit = ray.observer.Vector();
		exit.SetCoords(exit.X() + s_exit * ray.look_away.X(), exit.Y() + s_exit * ray.look_away.Y(), exit.Z() + s_exit * ray.look_away.Z());

		layer.entrance.FromVector(entrance, &m_coords);
		layer.exit.FromVector(exit, &m_coords);

		add_od_quadrature(layer, ray);
		add_solar_parameters(layer, ray);
	}


	TracedRay SphericalRayTracer::trace_ray_observer_inside_looking_up(const ViewingRay& ray) const {
		// Find the index to the first altitude ABOVE the observer
		auto it = std::upper_bound(m_altitudes.cbegin(), m_altitudes.cend(), ray.observer.Altitude());
		size_t start_index = std::distance(m_altitudes.cbegin(), it);

		TracedRay result;
		result.observer_and_look = ray;
		result.ground_is_hit = false;

		result.layers.resize(m_altitudes.size() - start_index);

		for (size_t i = m_altitudes.size() - 1; i != start_index; --i) {
			auto& layer = result.layers[i - start_index];
			complete_layer(layer, ray, i, ViewingDirection::up, TangentSide::nearside);
		}
		partial_layer(result.layers[0], ray, start_index, ViewingDirection::up, TangentSide::nearside);

		return result;
	}

	TracedRay SphericalRayTracer::trace_ray_observer_inside_looking_ground(const ViewingRay& ray) const {
		TracedRay result;

		result.observer_and_look = ray;
		result.ground_is_hit = true;

		// Find the index to the first altitude ABOVE the observer
		auto it = std::upper_bound(m_altitudes.cbegin(), m_altitudes.cend(), ray.observer.Altitude());
		size_t start_index = std::distance(m_altitudes.cbegin(), it);

		result.layers.resize(start_index + 1);

		// Complete layers from the ground
		for (size_t i = 0; i < start_index; ++i) {
			auto& layer = result.layers[i];
			complete_layer(layer, ray, i, ViewingDirection::down, TangentSide::nearside);
		}
		partial_layer(result.layers[start_index], ray, start_index, ViewingDirection::down, TangentSide::nearside);
		return result;
	}
	TracedRay SphericalRayTracer::trace_ray_observer_inside_looking_limb(const ViewingRay& ray) const {
		TracedRay result;

		nxLog::Record(NXLOG_ERROR, "Trying to trace a ray where the observer is inside looking limb. Which is not yet implemented");
		throw std::runtime_error("Trying to trace a ray where the observer is inside looking limb.");

		return result;
	}


    template <int NSTOKES, int CNSTR>
	void SphericalSolarTransmission<NSTOKES, CNSTR>::add_transmission(TracedRay& ray, RayOptical& rayoptical) {
		ViewingRay to_sun;
		to_sun.look_away.SetCoords(0, 0, 1); // Sun is always in z direction for heliodetic

		rayoptical.solar_transmission_entrance.resize(ray.layers.size(), m_optical_state.num_deriv());
		rayoptical.solar_transmission_exit.resize(ray.layers.size(), m_optical_state.num_deriv());

		for (int i = 0; i < ray.layers.size(); ++i) {
			auto& ray_layer = ray.layers[i];

			for (const auto& layer : ray.solar_entrance_traced_layers[i]) {
				for (size_t j = 0; j < layer.num_entrance_weights; j++) {
					auto index = layer.entrance_weights[j].first;
					auto weight = layer.entrance_weights[j].second;

					rayoptical.solar_transmission_entrance[i].value += layer.od_quad_start * weight * m_optical_state.total_extinction(m_wavel_index)(index);

					const auto& deriv_mapping = m_optical_state.deriv_mapping(m_wavel_index)[index];
					for (const auto& deriv : deriv_mapping) {
						rayoptical.solar_transmission_entrance[i].deriv(deriv.deriv_index) += layer.od_quad_start * weight * deriv.d_totalext;
					}
				}

				for (size_t j = 0; j < layer.num_exit_weights; j++) {
					auto index = layer.exit_weights[j].first;
					auto weight = layer.exit_weights[j].second;

					rayoptical.solar_transmission_entrance[i].value += layer.od_quad_end * weight * m_optical_state.total_extinction(m_wavel_index)(index);

					const auto& deriv_mapping = m_optical_state.deriv_mapping(m_wavel_index)[index];
					for (const auto& deriv : deriv_mapping) {
						rayoptical.solar_transmission_entrance[i].deriv(deriv.deriv_index) += layer.od_quad_end * weight * deriv.d_totalext;
					}
				}
			}
			rayoptical.solar_transmission_entrance[i].value = std::exp(-1.0*rayoptical.solar_transmission_entrance[i].value);
			rayoptical.solar_transmission_entrance[i].deriv = -1.0*rayoptical.solar_transmission_entrance[i].value*rayoptical.solar_transmission_entrance[i].deriv;


			for (const auto& layer : ray.solar_exit_traced_layers[i]) {
				for (size_t j = 0; j < layer.num_entrance_weights; j++) {
					auto index = layer.entrance_weights[j].first;
					auto weight = layer.entrance_weights[j].second;

					rayoptical.solar_transmission_exit[i].value += layer.od_quad_start * weight * m_optical_state.total_extinction(m_wavel_index)(index);

					const auto& deriv_mapping = m_optical_state.deriv_mapping(m_wavel_index)[index];
					for (const auto& deriv : deriv_mapping) {
						rayoptical.solar_transmission_exit[i].deriv(deriv.deriv_index) += layer.od_quad_start * weight * deriv.d_totalext;
					}
				}

				for (size_t j = 0; j < layer.num_exit_weights; j++) {
					auto index = layer.exit_weights[j].first;
					auto weight = layer.exit_weights[j].second;

					rayoptical.solar_transmission_exit[i].value += layer.od_quad_end * weight * m_optical_state.total_extinction(m_wavel_index)(index);

					const auto& deriv_mapping = m_optical_state.deriv_mapping(m_wavel_index)[index];
					for (const auto& deriv : deriv_mapping) {
						rayoptical.solar_transmission_exit[i].deriv(deriv.deriv_index) += layer.od_quad_end * weight * deriv.d_totalext;
					}
				}
			}
			rayoptical.solar_transmission_exit[i].value = std::exp(-1.0*rayoptical.solar_transmission_exit[i].value);
			rayoptical.solar_transmission_exit[i].deriv = -1.0*rayoptical.solar_transmission_exit[i].value*rayoptical.solar_transmission_exit[i].deriv;
		}

	}

    template <int NSTOKES, int CNSTR>
	sktran_do_detail::Dual<double> SphericalIntegrator<NSTOKES, CNSTR>::integrate_sources(AEOrder m, const TracedRay& ray, RayOptical& rayoptical, const VectorDim1<double>& lp_mu) {

		if (m == 0) {
			// Have to resize the optical parameters
			rayoptical.k_entrance.resize(ray.layers.size(), m_optical_state.num_deriv());
			rayoptical.k_exit.resize(ray.layers.size(), m_optical_state.num_deriv());

			rayoptical.k_scat_entrance.resize(ray.layers.size(), m_optical_state.num_deriv());
			rayoptical.k_scat_exit.resize(ray.layers.size(), m_optical_state.num_deriv());

			rayoptical.phase_entrance.resize(ray.layers.size(), m_optical_state.num_deriv());
			rayoptical.phase_exit.resize(ray.layers.size(), m_optical_state.num_deriv());
		}

		size_t phase_idx = ray.observer_and_look.unsorted_ray_idx;

		sktran_do_detail::Dual<double> radiance(m_optical_state.num_deriv());
		sktran_do_detail::Dual<double> source(m_optical_state.num_deriv());
		sktran_do_detail::Dual<double> od(m_optical_state.num_deriv());

		if (ray.ground_is_hit && m == 0) {
			const auto& ground_layer = ray.layers[0];

			double brdf = m_optical_state.exact_brdf(ray.observer_and_look.look_away, ground_layer.exit, m_coords, m_wavel_index);
			radiance = rayoptical.solar_transmission_exit[0] * brdf * ground_layer.exit.CosSZA();

			if (m_ms_source != nullptr) {
				m_ms_source->add_ground_reflectance(m, ray, radiance);
			}

		}
	
		for (size_t i = 0; i < ray.layers.size(); i++) {

			integrate_layer_sources(m, ray.layers[i], ray, rayoptical, i, source, phase_idx, lp_mu);
			od = spherical_layer_od(ray.layers[i], m_optical_state, m_wavel_index);

			// Attenuate the radiances and add sources
			radiance = radiance * sktran_do_detail::dual::exp(-1.0*od) + source;
		}

		return radiance;
	}

    template <int NSTOKES, int CNSTR>
	void SphericalIntegrator<NSTOKES, CNSTR>::integrate_layer_sources(AEOrder m, const SphericalLayer& layer, const TracedRay& ray, RayOptical& rayoptical, size_t layer_idx, sktran_do_detail::Dual<double>& source, size_t phase_idx, const VectorDim1<double>& lp_mu) {
		source.value = 0.0;
		source.deriv.setZero();

		size_t numderiv = m_optical_state.num_deriv();

		const Dual<double>& entrance_transmission = rayoptical.solar_transmission_entrance[layer_idx];
		const Dual<double>& exit_transmission = rayoptical.solar_transmission_exit[layer_idx];

		Dual<double>& phase_entrance = rayoptical.phase_entrance[layer_idx];
		Dual<double>& k_entrance = rayoptical.k_entrance[layer_idx];
		Dual<double>& k_scat_entrance = rayoptical.k_scat_entrance[layer_idx];

		if (m == 0) {
			for (size_t i = 0; i < layer.num_entrance_weights; ++i) {
				auto index = layer.entrance_weights[i].first;
				auto weight = layer.entrance_weights[i].second;

				phase_entrance.value += weight * m_optical_state.phase_function(m_wavel_index)(index, phase_idx);
				k_entrance.value += weight * m_optical_state.total_extinction(m_wavel_index)(index);
				k_scat_entrance.value += weight * m_optical_state.scattering_extinction(m_wavel_index)(index);

				for (const auto& deriv : m_optical_state.deriv_mapping(m_wavel_index)[index]) {
					k_entrance.deriv(deriv.deriv_index) += weight * deriv.d_totalext;
					k_scat_entrance.deriv(deriv.deriv_index) += weight * deriv.d_totalscatext;
					phase_entrance.deriv(deriv.deriv_index) += weight * deriv.d_totalphase(phase_idx);
				}
			}
		}

		Dual<double>& phase_exit = rayoptical.phase_exit[layer_idx];
		Dual<double>& k_exit = rayoptical.k_exit[layer_idx];
		Dual<double>& k_scat_exit = rayoptical.k_scat_exit[layer_idx];

		if (m == 0) {
			for (size_t i = 0; i < layer.num_exit_weights; ++i) {
				auto index = layer.exit_weights[i].first;
				auto weight = layer.exit_weights[i].second;

				phase_exit.value += weight * m_optical_state.phase_function(m_wavel_index)(index, phase_idx);
				k_exit.value += weight * m_optical_state.total_extinction(m_wavel_index)(index);
				k_scat_exit.value += weight * m_optical_state.scattering_extinction(m_wavel_index)(index);

				for (const auto& deriv : m_optical_state.deriv_mapping(m_wavel_index)[index]) {
					k_exit.deriv(deriv.deriv_index) += weight * deriv.d_totalext;
					k_scat_exit.deriv(deriv.deriv_index) += weight * deriv.d_totalscatext;
					phase_exit.deriv(deriv.deriv_index) += weight * deriv.d_totalphase(phase_idx);
				}
			}
		}

		size_t num_splits = layer.type == LayerType::tangent ? m_num_integration_splits : 1;

		double ds = layer.layer_distance;
		double split_distance = ds / float(num_splits);
		HELIODETIC_VECTOR midpoint = layer.entrance.Vector();
		HELIODETIC_VECTOR look_delta;
		look_delta.SetCoords(ray.observer_and_look.look_away.X() * split_distance / 2, ray.observer_and_look.look_away.Y() * split_distance / 2, ray.observer_and_look.look_away.Z() * split_distance / 2);

		source.value = 0.0;
		Dual<double> internal_attenuation_factor(numderiv);
		internal_attenuation_factor.value = 1.0;

		// Containers for the sources
		Dual<double> partic_source(numderiv);
		Dual<double> homog_source(numderiv);
		bool recache = true;


		for (size_t i = 0; i < num_splits; ++i) {
			// Move to the middle of the split
			midpoint += look_delta;

			// Interpolate our quantities to the middle of the split layer
			double altitude = m_coords.RadiusToAltitude(midpoint.Magnitude());
			double w_end = abs((layer.entrance.Altitude() - altitude) / (layer.entrance.Altitude() - layer.exit.Altitude()));
			double w_start = 1 - w_end;

			if (m_ms_source != nullptr) {
				// TODO: Get layer average cos solar zenith of the integration segment?
				double avg_cos_sza = (layer.entrance.CosSZA() + layer.exit.CosSZA()) / 2.0;
				bool is_upwelling = layer.entrance.Altitude() > layer.exit.Altitude();
				double avg_viewing_cos = (layer.entrance.CosZenithAngle(ray.observer_and_look.look_away) + layer.exit.CosZenithAngle(ray.observer_and_look.look_away)) / 2.0;
				m_ms_source->non_integrated_source_unrolled(m, altitude, avg_cos_sza, avg_viewing_cos, is_upwelling, homog_source, partic_source, recache);
				recache = false;
				
				// Apply the azimuth expansion to the DOM sources
				double cos_factor = cos(m*layer.saz);

				partic_source *= cos_factor;
				homog_source *= cos_factor;
			}

			if (((w_start * k_entrance.value + w_end * k_exit.value) * split_distance) < 1e-6) {
				// For really small values we can do a series approximation instead which is more accurate
				m_cached_dist_factor.value = split_distance;
				m_cached_dist_factor.deriv.setZero();
			}
			else {
				m_cached_dist_factor.value = (1 - std::exp(-1.0*(w_start * k_entrance.value + w_end * k_exit.value) * split_distance)) / (w_start * k_entrance.value + w_end * k_exit.value);

				for (uint k = 0; k < numderiv; ++k) {
					m_cached_dist_factor.deriv(k) = (w_start * k_entrance.deriv(k) + w_end * k_exit.deriv(k)) *
						((split_distance * std::exp(-1.0*(w_start * k_entrance.value + w_end * k_exit.value) * split_distance) - m_cached_dist_factor.value) / (w_start * k_entrance.value + w_end * k_exit.value));
				}
			}

			// Add in the SS_Source
			if (m == 0) {
				source.value += internal_attenuation_factor.value * m_cached_dist_factor.value * (w_start*k_scat_entrance.value + w_end*k_scat_exit.value) *
					(w_start*phase_entrance.value + w_end*phase_exit.value) * (w_start * entrance_transmission.value + w_end * exit_transmission.value);


				for (uint k = 0; k < numderiv; ++k) {
					source.deriv(k) += internal_attenuation_factor.deriv(k) * m_cached_dist_factor.value * (w_start*k_scat_entrance.value + w_end * k_scat_exit.value) *
						(w_start*phase_entrance.value + w_end * phase_exit.value) * (w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.deriv(k) * (w_start*k_scat_entrance.value + w_end * k_scat_exit.value) *
						(w_start*phase_entrance.value + w_end * phase_exit.value) * (w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * (w_start*k_scat_entrance.deriv(k) + w_end * k_scat_exit.deriv(k)) *
						(w_start*phase_entrance.value + w_end * phase_exit.value) * (w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * (w_start*k_scat_entrance.value + w_end * k_scat_exit.value) *
						(w_start*phase_entrance.deriv(k) + w_end * phase_exit.deriv(k)) * (w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * (w_start*k_scat_entrance.value + w_end * k_scat_exit.value) *
						(w_start*phase_entrance.value + w_end * phase_exit.value) * (w_start * entrance_transmission.deriv(k) + w_end * exit_transmission.deriv(k));
				}
			}
			// Add in the MS sources
			if (m_ms_source != nullptr) {
				// First the particular source
				source.value += internal_attenuation_factor.value * m_cached_dist_factor.value * partic_source.value * (w_start * k_entrance.value + w_end * k_exit.value) *
					(w_start * entrance_transmission.value + w_end * exit_transmission.value);
				for (uint k = 0; k < numderiv; ++k) {
					source.deriv(k) += internal_attenuation_factor.deriv(k) * m_cached_dist_factor.value * partic_source.value * (w_start * k_entrance.value + w_end * k_exit.value) *
						(w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.deriv(k) * partic_source.value * (w_start * k_entrance.value + w_end * k_exit.value) *
						(w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * partic_source.deriv(k) * (w_start * k_entrance.value + w_end * k_exit.value) *
						(w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * partic_source.value * (w_start * k_entrance.deriv(k) + w_end * k_exit.deriv(k)) *
						(w_start * entrance_transmission.value + w_end * exit_transmission.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * partic_source.value * (w_start * k_entrance.value + w_end * k_exit.value) *
						(w_start * entrance_transmission.deriv(k) + w_end * exit_transmission.deriv(k));
				}

				
				// Then the homogeneous source
				source.value += internal_attenuation_factor.value * m_cached_dist_factor.value * homog_source.value * (w_start * k_entrance.value + w_end * k_exit.value);

				for (uint k = 0; k < numderiv; ++k) {
					source.deriv(k) += internal_attenuation_factor.deriv(k) * m_cached_dist_factor.value * homog_source.value
						* (w_start * k_entrance.value + w_end * k_exit.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.deriv(k) * homog_source.value
						* (w_start * k_entrance.value + w_end * k_exit.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * homog_source.deriv(k)
						* (w_start * k_entrance.value + w_end * k_exit.value);
					source.deriv(k) += internal_attenuation_factor.value * m_cached_dist_factor.value * homog_source.value
						* (w_start * k_entrance.deriv(k) + w_end * k_exit.deriv(k));
				}
			}

			// source += internal_attenuation_factor * m_cached_dist_factor * ((ss_source + partic_source * k_avg) * transmission_avg + homog_source * k_avg);

			internal_attenuation_factor *= std::exp(-1.0*(w_start * k_entrance.value + w_end * k_exit.value) * split_distance);

			for (uint k = 0; k < numderiv; ++k) {
				internal_attenuation_factor.deriv(k) += -1.0*(w_start * k_entrance.deriv(k) + w_end * k_exit.deriv(k)) * split_distance * internal_attenuation_factor.value;
			}
			// Move to the end of the split
			midpoint += look_delta;
		}
	}

    template class SphericalSolarTransmission<1>;
	template class SphericalSolarTransmission<3>;
    template class SphericalSolarTransmission<4>;

    INSTANTIATE_TEMPLATE(SphericalIntegrator);
}
