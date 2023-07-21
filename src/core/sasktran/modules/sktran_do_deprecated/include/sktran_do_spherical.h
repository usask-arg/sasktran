#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"

namespace sktran_do_detail
{
    // A complete layer starts and ends on OpticalState altitude boundaries
    // A partial layer starts or ends (or both) internally in a layer
    // A tangent layer ends of starts at the tnagent point
	enum LayerType {
		complete,
		partial,
		tangent
	};

    // Geometry components of a line of sight
	struct ViewingRay {
		HELIODETIC_POINT observer;
		HELIODETIC_UNITVECTOR look_away;
		size_t unsorted_ray_idx;

		double cos_viewing() const {
			return look_away & observer.UnitVector();
		}

	};

    // Values needed to integrate a layer.  Rather than interpolating anything we directly store weights/indices
    // to the optical table to speed up calculations over multiple wavelengths and allow for easier calculation
    // of derivatives.  Everything here is wavelength independent
	struct SphericalLayer {
		HELIODETIC_POINT entrance;
		HELIODETIC_POINT exit;
		double layer_distance;     // Total distance of the ray within the layer

        // Quadrature parameters
		double od_quad_start;
		double od_quad_end;

        // Solar azimuth
		double saz;

		// Properties relating to the optical table
        // Faster to store 2 weights and the total number of weights instead of allocating
		std::array<std::pair<size_t, double>, 2> entrance_weights; // Indicies and interpolation weights to the entrance altitude in the optical table
		std::array<std::pair<size_t, double>, 2> exit_weights;     // Indicies and interpolation weights to the exit altitude in the optical table
		size_t num_entrance_weights;							   // Number of interpolation points for the entrance
		size_t num_exit_weights;                                   // Number of interpolation points for the exit

		LayerType type;
	};

    // Calculates the OD quadrature terms (analytic formula assuming linear variation in altitude of extinction)
    // and updates them inside a spherical layer.  Note these are wavelength independent and purely geometry quantities.
	inline void add_od_quadrature(SphericalLayer& layer, const ViewingRay& ray) {
		double r0 = layer.entrance.Radius();
		double r1 = layer.exit.Radius();
		double dr = r1 - r0;

		if (abs(dr) < 0.001) {
			// Tiny layer, we can just use the average of the extinctions
			layer.od_quad_start = layer.layer_distance / 2;
			layer.od_quad_end = layer.layer_distance / 2;
		}
		else {
			double costheta0 = layer.entrance.CosZenithAngle(ray.look_away);
			double costheta1 = layer.exit.CosZenithAngle(ray.look_away);

			double t0 = r0 * costheta0;
			double t1 = r1 * costheta1;
			double rt = r0 * sqrt(1.0 - costheta0 * costheta0);

			double dt1;
			double dt2;

			if (t1 >= t0) {
				dt1 = t1 - t0;
				if (abs(rt) < 10) {
					dt2 = 0.5*((r1*t1 - r0 * t0));
				}
				else {
					dt2 = 0.5*((r1*t1 - r0 * t0) + rt * rt*log((r1 + t1) / (r0 + t0)));
				}
			}
			else {
				dt1 = t0 - t1;
				if (abs(rt) < 10) {
					dt2 = 0.5*((r0 + t0 - r1 * t1));
				}
				else {
					dt2 = 0.5*((r0*t0 - r1 * t1) + rt * rt*log((r0 + t0) / (r1 + t1)));
				}
			}
			layer.od_quad_start = (r1*dt1 - dt2) / dr;
			layer.od_quad_end = -1 * (r0*dt1 - dt2) / dr;
		}
	}

    // Result after tracing a ray, still geometry independent.
	struct TracedRay {
		ViewingRay observer_and_look;

		bool ground_is_hit;					// True if the ground is hit

		std::vector<SphericalLayer> layers;
		std::vector<std::vector<SphericalLayer>> solar_entrance_traced_layers;
		std::vector<std::vector<SphericalLayer>> solar_exit_traced_layers;
	};

    // Optical components of a ray, optical parameters along the ray.
	struct RayOptical {
        // Solar transmission at entrance/exit points
		std::vector<Dual<double>> solar_transmission_entrance;
		std::vector<Dual<double>> solar_transmission_exit;

		// TODO: NSTOKES phase matrix
		std::vector<Dual<double>> phase_entrance;
		std::vector<Dual<double>> phase_exit;

        // Extinctions at entrance/exit
		std::vector<Dual<double>> k_entrance;
		std::vector<Dual<double>> k_exit;

        // Scattering extinctions at entrance/exit
		std::vector<Dual<double>> k_scat_entrance;
		std::vector<Dual<double>> k_scat_exit;
	};

    // Uses the geometry quantities in layer and the OpticalState to calculation the optical depth along the
    // layer
    template <int NSTOKES, int CNSTR=-1>
	inline sktran_do_detail::Dual<double> spherical_layer_od(const SphericalLayer& layer, const sktran_do_detail::OpticalState<NSTOKES>& optical_state,
                                                             size_t wavel_index) {

		Dual<double> result(optical_state.num_deriv());

		for (size_t i = 0; i < layer.num_entrance_weights; i++) {
			auto index = layer.entrance_weights[i].first;
			auto weight = layer.entrance_weights[i].second;

			result.value += weight * layer.od_quad_start * optical_state.total_extinction(wavel_index)(index);

			for (const auto& deriv : optical_state.deriv_mapping(wavel_index)[index]) {
				result.deriv(deriv.deriv_index) += weight * layer.od_quad_start * deriv.d_totalext;
			}

		}
		for (size_t i = 0; i < layer.num_exit_weights; i++) {
			auto index = layer.exit_weights[i].first;
			auto weight = layer.exit_weights[i].second;

			result.value += weight * layer.od_quad_end * optical_state.total_extinction(wavel_index)(index);

			for (const auto& deriv : optical_state.deriv_mapping(wavel_index)[index]) {
				result.deriv(deriv.deriv_index) += weight * layer.od_quad_end * deriv.d_totalext;
			}
		}		

		return result;
	};

    // Traces rays in a spherical atmosphere
	class SphericalRayTracer {
	public:
		SphericalRayTracer(const std::vector<double>& altitudes, const SKTRAN_CoordinateTransform_V2& coords) :
			m_altitudes(altitudes),
			m_coords(coords)
		{
			m_earth_radius = coords.AltitudeToRadius(0.0);
		}

		TracedRay trace_ray(const ViewingRay& ray) const;
		void trace_solar_rays(TracedRay& ray) const;

	private:
		const std::vector<double>& m_altitudes;
		double m_earth_radius;
		const SKTRAN_CoordinateTransform_V2& m_coords;

		enum ViewingDirection { up = -1, down = 1 };
		enum TangentSide {farside = -1, nearside = 1};

		TracedRay trace_ray_observer_outside_ground_viewing(const ViewingRay& ray) const;
		TracedRay trace_ray_observer_inside_looking_up(const ViewingRay& ray) const;
		TracedRay trace_ray_observer_inside_looking_ground(const ViewingRay& ray) const;
		TracedRay trace_ray_observer_inside_looking_limb(const ViewingRay& ray) const;
		TracedRay trace_ray_observer_outside_limb_viewing(const ViewingRay& ray) const;

		void complete_layer(SphericalLayer& layer, const ViewingRay& ray, size_t exit_index, ViewingDirection direction, TangentSide side) const;
		void partial_layer(SphericalLayer& layer, const ViewingRay& ray, size_t start_index, ViewingDirection direction, TangentSide side) const;
		void tangent_layer(SphericalLayer& layer, const ViewingRay& ray, size_t upper_index, double tangent_altitude, ViewingDirection direction, TangentSide side) const;

		double distance_to_altitude(const ViewingRay& ray, double altitude, ViewingDirection direction, TangentSide side) const {
			double cos_zenith = abs(ray.cos_viewing());
			double ro = ray.observer.Radius();
			double re = m_earth_radius + altitude;


			double rtsq = ro * ro*(1 - cos_zenith * cos_zenith);

			// Distance to the tangent point
			double tangent_distance = side * direction * ro * cos_zenith;

			// Distance from the tangent point to the requested altitude
			double dist_from_tangent;
			if (rtsq > re*re) {
				// Either a rounding error or a problem
				if (abs(rtsq - re * re) < 0.1) {
					dist_from_tangent = 0.0;
				}
				else {
					throw("Error, requesting distance to a shell that does not exist");
				}
			}
			else {
				dist_from_tangent = side * direction * sqrt(re*re - rtsq);
			}


			if (side == TangentSide::nearside) {
				return tangent_distance - dist_from_tangent;
			}
			else {
				return tangent_distance + dist_from_tangent;
			}
		}

		void add_solar_parameters(SphericalLayer& layer, const ViewingRay& ray) const {
            // Need to calculate SAZ in the layer in order to do the azimuth expansion?
            // Not 100% sure on this, should we use SAZ of the layer or SAZ at the CSZ the multiple scatter source
            // is calculated at?

			HELIODETIC_VECTOR midpoint = (layer.exit.Vector() + layer.entrance.Vector());
			midpoint *= 0.5;

			HELIODETIC_UNITVECTOR sun_unit = m_coords.GeographicToHelioUnitVector(m_coords.SunUnit());
			HELIODETIC_UNITVECTOR local_up = midpoint.UnitVector();
			HELIODETIC_VECTOR sun_vert_trans(local_up, sun_unit & local_up);
			HELIODETIC_VECTOR los_vert_trans(local_up, ray.look_away & local_up);

			const HELIODETIC_VECTOR sun(sun_unit, 1);
			HELIODETIC_VECTOR sun_projected = sun - sun_vert_trans;

			const HELIODETIC_VECTOR propagating(ray.look_away, 1);
			HELIODETIC_VECTOR los_projected = propagating - los_vert_trans;

			// TODO: We should be handling the special cases here where either the sun or los direction is
			// directly perpinduclar to the tangent plane

			double proj = sun_projected.UnitVector() & los_projected.UnitVector();
			if (proj > 1) {
				proj = 1;
			}
			else if (proj < -1) {
				proj = -1;
			}

			layer.saz = acos(proj);
		}
	};

    // TOOD: Doesn't actually depend on NSTOKES, should define an interface for OpticalState that implements the spherical
    // solar transmission functionality and remove this template
    template <int NSTOKES, int CNSTR=-1>
	class SphericalSolarTransmission {
	public:
		SphericalSolarTransmission(const SphericalRayTracer& ray_tracer, const sktran_do_detail::OpticalState<NSTOKES>& optical_state,
                                   double wavelength) :
			m_ray_tracer(ray_tracer),
			m_optical_state(optical_state)
		{
            m_wavel_index = optical_state.wavel_index(wavelength);
		}

		void add_transmission(TracedRay& ray, RayOptical& rayoptical);
	private:
		const SphericalRayTracer& m_ray_tracer;
		const sktran_do_detail::OpticalState<NSTOKES>& m_optical_state;

        size_t m_wavel_index;
	};

    template <int NSTOKES, int CNSTR=-1>
	class SphericalPostProcessing;

    template <int NSTOKES, int CNSTR=-1>
	class SphericalIntegrator {
	public:
		SphericalIntegrator(const sktran_do_detail::OpticalState<NSTOKES>& optical_state, const SphericalSolarTransmission<NSTOKES>& solar_transmission,
			const SKTRAN_CoordinateTransform_V2& coords, double wavelength, SphericalPostProcessing<NSTOKES, CNSTR>* ms_source = nullptr) :
			m_optical_state(optical_state),
			m_solar_transmission(solar_transmission),
			m_coords(coords),
			m_num_integration_splits(5),
			m_ms_source(ms_source)
		{
			m_cached_dist_factor.resize(m_optical_state.num_deriv());
            m_wavel_index = optical_state.wavel_index(wavelength);
		}

		sktran_do_detail::Dual<double> integrate_sources(AEOrder m, const TracedRay& ray, RayOptical& rayoptical, const VectorDim1<double>& lp_mu);

	private:
		const sktran_do_detail::OpticalState<NSTOKES>& m_optical_state;
		const SphericalSolarTransmission<NSTOKES>& m_solar_transmission;
		const SKTRAN_CoordinateTransform_V2& m_coords;
		SphericalPostProcessing<NSTOKES, CNSTR>* m_ms_source;
		const size_t m_num_integration_splits;
        size_t m_wavel_index;

		Dual<double> m_cached_dist_factor;

		void integrate_layer_sources(AEOrder m, const SphericalLayer& layer, const TracedRay& ray, RayOptical& rayoptical, size_t layer_idx, sktran_do_detail::Dual<double>& source, size_t phase_idx, const VectorDim1<double>& lp_mu);
	};
}