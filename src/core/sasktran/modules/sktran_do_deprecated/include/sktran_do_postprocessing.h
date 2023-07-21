#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"

namespace sktran_do_detail
{
	// Performs post-processing on a RTE solution. Essentially integrates an 
	// interpolation of the RTE solutions in a given layer to arbitrary zenith 
	// angles.
    template <int NSTOKES, int CNSTR=-1>
	using STIntegralProperties = ReadOnlyProperties<BasicProperties<NSTOKES>, SolarProperties<NSTOKES>>;
    template <int NSTOKES, int CNSTR=-1>
	class ParticipatingSourceTerm: public STIntegralProperties<NSTOKES>
	{
	using HomogType = typename std::conditional<NSTOKES != 5, double, std::complex<double>>::type;

	public: // Interface
		ParticipatingSourceTerm(AEOrder m, const OpticalLayer<NSTOKES, CNSTR>* layer, const OpticalLayerArray<NSTOKES, CNSTR>& layers, double mu, double obsod, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu, bool ssonly=false,
                                bool includesinglescatter=true);

		// Integrates the participating source term across the entire layer
		inline void integrate(Radiance<NSTOKES>& result, PostProcessingCache<NSTOKES>& cache, const sktran_do_detail::Radiance<NSTOKES>* manual_ss_source=nullptr) const
		{
			if (M_OBSOD < M_EXIT_OPTICAL_DEPTH) {
				integrate(M_EXIT_OPTICAL_DEPTH, result, cache, manual_ss_source);
			}
			else {
				integrate(M_OBSOD, result, cache, manual_ss_source);
			}
		}
		// Integrate the participating source term from the entry point to the 
		// layer to the terminal point at optical_depth. An exception is thrown
		// if optical_depth is not valid.
		void integrate(double optical_depth, Radiance<NSTOKES>& result, PostProcessingCache<NSTOKES>& cache, const sktran_do_detail::Radiance<NSTOKES>* manual_ss_source) const;
	private: // Functions
		enum class SolutionType { PARTICULAR, HOMOGENEOUS };

		// See Eq (39)
		void h_plus(double x,
				    double thickness,
				    SolutionIndex j, 
					LayerDual<HomogType>& xform) const;

		// See Eq (40)
		void h_minus(double x,
			         double thickness,
			         SolutionIndex j,
			         LayerDual<HomogType>& xform) const;

		// See Eq (41)
		void E(double x,
			   double thickness,
			   const Dual<double>& transmission,
			   Dual<double>& xform) const;

		inline bool verifyIntegrationUpperBound(double optical_depth) const {
			if(M_COSZEN > 0) { // upwards
				return (optical_depth < M_EXIT_OPTICAL_DEPTH + M_OPTICAL_THICKNESS && optical_depth >= M_EXIT_OPTICAL_DEPTH);
			} else {
				return (optical_depth <= M_EXIT_OPTICAL_DEPTH && optical_depth > M_EXIT_OPTICAL_DEPTH - M_OPTICAL_THICKNESS);
			}
		}

	private: // Members
		const OpticalLayer<NSTOKES, CNSTR>& m_layer;
		const OpticalLayerArray<NSTOKES, CNSTR>& m_layers;
		const LayerSolution<NSTOKES, CNSTR>& m_solution;
		const AEOrder M_FEORDER;
		const double M_EXIT_OPTICAL_DEPTH;
		const double M_OPTICAL_THICKNESS;
		const double M_COSZEN;
		const double M_OBSOD;
		const std::vector<LegendrePhaseContainer<NSTOKES>>& M_LP_COSZEN;
		bool M_SSONLY;
        bool M_INCLUDESS;
	};

    template <int NSTOKES, int CNSTR>
	class SphericalPostProcessing: BasicProperties<NSTOKES> {
	public:
		SphericalPostProcessing(const std::vector<double>& cos_sza, PersistentConfiguration<NSTOKES, CNSTR>& config,
			double wavelength,
			const OpticalState<NSTOKES>* opticalstate,
			const std::vector<LineOfSight>& los,
            const GeometryLayerArray<NSTOKES, CNSTR>& geometrylayers,
            bool use_upwelling_spher = true
			) :
			BasicProperties<NSTOKES>(config),
			m_cos_sza(cos_sza),
			m_config(config),
			m_opticalstate(opticalstate),
			m_los(los),
			m_use_upwelling_spher(use_upwelling_spher),
            m_geometrylayers(geometrylayers)
		{
			const auto* coords = m_config.coords();
			m_opticallayers.reserve(m_cos_sza.size());
			m_rte.reserve(m_cos_sza.size());
			for (size_t i = 0; i < m_cos_sza.size(); ++i) {
				// Get albedo
				HELIODETIC_POINT refpt = coords->ReferencePoint(m_config.userSpec()->getBottomAltitude());
				GEODETIC_INSTANT ref_inst = coords->PointToGeodetic(refpt, coords->ReferencePointMJD());

				const skBRDF* skbrdf = opticalstate->brdf_object();

				m_brdfs.emplace_back(std::unique_ptr<Wrapped_skBRDF>(new Wrapped_skBRDF(wavelength, *skbrdf, ref_inst)));

				config.override_csz(i);
				m_opticallayers.emplace_back(config, wavelength, opticalstate, los, std::move(m_brdfs.back()), nullptr, false, -1, &m_geometrylayers);

				m_rte.emplace_back(config, m_opticallayers.back());
			}

			if (m_config.perturbation_specs()) {
				m_numwf = m_config.perturbation_specs()->size();
			}
			else {
				m_numwf = 0;
			}

			// Generate caches
			m_cached_layer_homog_source.resize(m_opticallayers[0].inputDerivatives().numDerivative());
			m_cached_layer_part_source.resize(m_opticallayers[0].inputDerivatives().numDerivative());

			m_cached_lp_mu.resize(this->M_NSTR);
		}

		void solve(AEOrder m) {
			for (auto& rte : m_rte) {
				rte.solve(m);
			}
		}

		void non_integrated_source_unrolled(AEOrder m, double altitude, double csz, double cos_viewing, bool isupwelling, Dual<double>& homog_source, Dual<double>& partic_source, bool recalculate_caches = true);

		void add_ground_reflectance(AEOrder m, const TracedRay& ray, Dual<double>& radiance);

	private:
		const std::vector<double>& m_cos_sza;
		PersistentConfiguration<NSTOKES, CNSTR>& m_config;
        const GeometryLayerArray<NSTOKES, CNSTR>& m_geometrylayers;
		const OpticalState<NSTOKES>* m_opticalstate;
		const std::vector<LineOfSight>& m_los;

		std::vector<OpticalLayerArray<NSTOKES, CNSTR>> m_opticallayers;
		std::vector<RTESolver<NSTOKES, CNSTR>> m_rte;
		std::vector<std::unique_ptr<BRDF_Base>> m_brdfs;

		Dual<double> m_cached_layer_part_source;
		Dual<double> m_cached_layer_homog_source;

		VectorLayerDual<double> m_cached_lpsum_plus;
		VectorLayerDual<double> m_cached_lpsum_minus;

		VectorDim1<LegendrePhaseContainer<NSTOKES>> m_cached_lp_mu;

		size_t m_numwf;

		bool m_use_upwelling_spher;

		void interpolation_index(double csz, std::array<size_t, 2>& indicies, std::array<double, 2>& weights) const {
			if (m_cos_sza.size() == 1) {
				// Special case,
				indicies[0] = 0;
				indicies[1] = 0;
				weights[0] = 1.0;
				weights[1] = 0.0;

				return;
			}

			if (csz >= m_cos_sza[m_cos_sza.size() - 1]) {
				indicies[0] = m_cos_sza.size() - 1;
				indicies[1] = 0;
				weights[0] = 1.0;
				weights[1] = 0.0;

				return;
			}
			if (csz <= m_cos_sza[0]) {
				indicies[0] = 0;
				indicies[1] = 0;
				weights[0] = 1.0;
				weights[1] = 0.0;

				return;
			}

			nxLinearInterpolate::FindBoundingIndicesAscending(m_cos_sza.cbegin(), m_cos_sza.cend(), csz, &indicies[0],
				&indicies[1], &weights[0], &weights[1]);

			weights[1] = (csz - weights[0]) / (weights[1] - weights[0]);
			weights[0] = (1 - weights[1]);
		}
	};
}