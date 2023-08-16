#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco
{
	// Provides a concise optical description of a atmospheric layer, caches 
	// some calculations, stores DO solutions to the RTE for this layer.
    template <int NSTOKES, int CNSTR=-1>
	using OpticalLayerROP = ReadOnlyProperties<BasicProperties<NSTOKES>, SolarProperties<NSTOKES>, UserSpecProperties>;
    template <int NSTOKES, int CNSTR=-1>
	class OpticalLayer:
		public OpticalLayerROP<NSTOKES>,
		public AzimuthDependencyCascade
	{
		using ScalarDerivativeResult = sasktran_disco::LayerFundamentalDerivative<NSTOKES>;
		using HomogType = typename std::conditional<NSTOKES != 5, double, std::complex<double>>::type;
	public:
        OpticalLayer(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                     LayerIndex index,
                     double altitude_ceiling,
                     double altitude_floor,
                     const InputDerivatives<NSTOKES>& input_derivs,
                     const ThreadData<NSTOKES, CNSTR>& thread_data
                     );

        void set_optical(double scat_ext, double tot_ext,
                         const VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>& phasef_expansion,
                         double optical_depth_ceiling,
                         double optical_depth_floor
                         );

		OpticalLayer(const PersistentConfiguration<NSTOKES, CNSTR>& config,
					 LayerIndex index,
					 double scat_ext,
					 double tot_ext,
					 std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>> phasef_expansion,
					 double optical_depth_ceiling,
					 double optical_depth_floor,
					 double altitude_ceiling,
					 double altitude_floor,
					 const InputDerivatives<NSTOKES>& input_derivs);

		void takeDerivative(bool take_deriv);
		void configureDerivative();
		void configurePseudoSpherical(const Dual<double>& od_top, const Dual<double>& od_bottom);

		// Returns the value of the scattering phase function originating from
		// stream with index originating_stream. lp_out is a vector of 
		// associated Legendre polynomials evaluated at the outgoing zenith 
		// angle.
		void vectordual_scatPhaseF(AEOrder m,
			const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_out,
			const InputDerivatives<NSTOKES>& in_deriv,
			VectorLayerDual<double>& result_positive_stream,
			VectorLayerDual<double>& result_negative_stream) const;

		// Returns the value of the scattering phase function originating from 
		// origin to outgoing. This is a very efficient call.
		inline void inplace_scatPhaseFAndDerivative(AEOrder m, StreamIndex outgoing, StreamIndex origin, TripleProductDerivativeHolder<NSTOKES>& dual) const {
			m_legendre_sum[m].inplace_dual(outgoing, origin, dual);
		}


		// Return the value of the single scatter source term. The result 
		// of the singleScatST evaluated at the negative of the coszenith
		// is stored in result_negative_coszenith (it is free to calculate).
		void singleScatST(AEOrder m,
						  const VectorDim1<LegendrePhaseContainer<NSTOKES>>& outgoing_coszenith,
						  InhomogeneousSourceHolder<NSTOKES>& result,
						  InhomogeneousSourceHolder<NSTOKES>& result_negative_coszenith) const;

		// Returns the single scatter albedo in this layer
		inline double ssa() const
		{
			return M_SSA;
		}

		inline double opticalDepth(Location loc) const {
			switch (loc) {
				case Location::FLOOR: return M_OPTICALDEPTH_FLOOR; break;
				case Location::CEILING: return M_OPTICALDEPTH_CEILING; break;
				case Location::INSIDE: return M_OPTICAL_THICKNESS; break;
				default:
					abort();
			}
		}

		inline double altitude(Location loc) const {
			switch (loc) {
			case Location::FLOOR: return M_ALT_FLOOR; break;
			case Location::CEILING: return M_ALT_CEILING; break;
			default:
				abort();
			}
		}
		
		inline double beamTransmittance(Location loc, double x = -1) const {
			switch (loc) {
				case Location::FLOOR: return m_dual_bt_floor.value; break;
				case Location::CEILING: return m_dual_bt_ceiling.value; break;
				case Location::INSIDE:
					if (x < 0)
						abort();
					return m_dual_bt_ceiling.value * std::exp(-x * m_average_secant.value);
				default:
					abort();
			}
		}

		inline double d_beamTransmittance(Location loc, const LayerInputDerivative<NSTOKES>& deriv, uint derivindex, double x = -1) const {
			LayerIndex d_idx = deriv.layer_index;
			LayerIndex cur_layer = M_INDEX;

			switch (loc) {
			case Location::FLOOR:
				return m_dual_bt_floor.deriv(derivindex);
				break;
			case Location::CEILING:
				return m_dual_bt_ceiling.deriv(derivindex);
				break;
			case Location::INSIDE:
				if (x < 0)
					abort();
				if (d_idx < cur_layer) {
					// transmittance is bt_ceiling * exp(-x * secant)
					// dx = 0.0
					double dx = 0.0;
					return std::exp(-x * m_average_secant.value) *(m_dual_bt_ceiling.deriv(derivindex) - m_average_secant.value * m_dual_bt_ceiling.value * dx - x * m_average_secant.deriv(derivindex) * m_dual_bt_ceiling.value);
				}
				else if (d_idx == cur_layer) {
					// transmittance is = bt_ceiling * exp(-x * secant)
					// dx = layer_fraction * deriv.d_optical_depth

					double layerfraction = x / (M_OPTICAL_THICKNESS);
					double dx = layerfraction * deriv.d_optical_depth;

					return std::exp(-x * m_average_secant.value) *(m_dual_bt_ceiling.deriv(derivindex) - m_average_secant.value * m_dual_bt_ceiling.value * dx - x*m_average_secant.deriv(derivindex) * m_dual_bt_ceiling.value);
				}
				else {
					return 0.0;
				}
				break;
			default:
				abort();
			}
		}

		inline Dual<double> dual_beamTransmittance(Location loc, const InputDerivatives<NSTOKES>& deriv, double x = -1) const {
			Dual<double> result(deriv.numDerivative());

			result.value = beamTransmittance(loc, x);

			for (uint i = 0; i < deriv.numDerivative(); ++i) {
				result.deriv[i] = d_beamTransmittance(loc, deriv.layerDerivatives()[i], i, x);
			}

			return result;
		}

		// Returns the transmittance factor for the given separation constant.
		inline HomogType streamTransmittance(Location loc, AEOrder m, SolutionIndex j) const {
			if (loc != Location::INSIDE) {
				abort();
			}
			return exp(-std::abs(m_solutions[m].value.eigval(j)) * opticalDepth(Location::INSIDE));
		}

		inline HomogType d_streamTransmittance(Location loc, AEOrder m, SolutionIndex j, uint derivIndex, const LayerInputDerivative<NSTOKES>& deriv) const {
			if (loc != Location::INSIDE) {
				abort();
			}
			double d_od = deriv.d_optical_depth;

			return -(m_solutions[m].value.dual_eigval().deriv(derivIndex, j) * M_OPTICAL_THICKNESS + m_solutions[m].value.eigval(j) * d_od) * streamTransmittance(loc, m, j);
		}

		inline Dual<HomogType> dual_streamTransmittance(Location loc, AEOrder m, SolutionIndex j, const InputDerivatives<NSTOKES>& deriv) const {
			// TODO: Layerdual instead?

			size_t derivStart = deriv.layerStartIndex(M_INDEX);
			if (loc != Location::INSIDE) {
				abort();
			}
			Dual<HomogType> result(deriv.numDerivative());
			result.value = streamTransmittance(loc, m, j);
			// Have no cross derivatives
			for (uint i = 0; i < deriv.numDerivativeLayer(M_INDEX); ++i) {
				result.deriv[i + derivStart] = d_streamTransmittance(loc, m, j, i, deriv.layerDerivatives()[derivStart + i]);
			}

			return result;
		}

		// Returns the index of this layer
		inline LayerIndex index() const
		{
			return M_INDEX;
		}

		// Returns the expansion of the phase function for this layer
		inline const VectorDim1<LegendreCoefficient<NSTOKES>>& expansionPhaseFunction(AEOrder m) const {
			return *m_lephasef;
		}

		inline bool computeDerivative() const {
			return m_compute_deriv;
		}
		
		inline const LayerSolution<NSTOKES, CNSTR>& solution(AEOrder m) const {
			return m_solutions[m];
		}

		inline LayerSolution<NSTOKES, CNSTR>& solution(AEOrder m) {
			return m_solutions[m];
		}

		inline double scatExt() const {
			return M_SCAT_EXT;
		}
		
		inline double totalExt() const {
			return M_TOT_EXT;
		}

		inline const LayerDual<double>& dual_thickness() const {
			return m_dual_thickness;
		}

        inline const LayerDual<double>& dual_ssa() const {
            return m_dual_ssa;
        }

		inline const Dual<double>& dual_average_secant() const {
			return m_average_secant;
		}

        inline const Dual<double>& ceiling_beam_transmittanc() const {
            return m_dual_bt_ceiling;
        }

        inline Dual<double>& ceiling_beam_transmittance() const {
            return m_dual_bt_ceiling;
        }

        inline Dual<double>& floor_beam_transmittance() const {
            return m_dual_bt_floor;
        }

        inline Dual<double>& dual_average_secant() {
            return m_average_secant;
        }

		inline const std::vector<LegendreCoefficient<NSTOKES>>& legendre_coeff() const { return *m_lephasef; }

        void integrate_source(AEOrder m,
                              double mu,
                              double obsod,
                              const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu,
                              Radiance<NSTOKES>& result,
                              const sasktran_disco::Radiance<NSTOKES>* manual_ss_source=nullptr,
                              bool include_ss=true
                              ) const;

        // See Eq (39)
        void h_plus(AEOrder m,
                    double mu,
                    double x,
                    double thickness,
                    SolutionIndex j,
                    LayerDual<HomogType>& xform) const;

        // See Eq (40)
        void h_minus(AEOrder m,
                     double mu,
                     double x,
                     double thickness,
                     SolutionIndex j,
                     LayerDual<HomogType>& xform) const;

        // See Eq (41)
        void E(double mu,
                double x,
                double thickness,
                const Dual<double>& transmission,
                Dual<double>& xform) const;

	protected:
		// Optical properties
		double								        M_SSA;
		double								        M_SCAT_EXT;
		double								        M_TOT_EXT;
		double								        M_OPTICALDEPTH_FLOOR;
		double								        M_OPTICALDEPTH_CEILING;
		double								        M_OPTICAL_THICKNESS;
		const double								M_ALT_CEILING;
		const double								M_ALT_FLOOR;
        std::unique_ptr<VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>>	m_lephasef;
		const LayerIndex							M_INDEX;
		LegendreSumMatrix<NSTOKES>					m_legendre_sum;

        LayerCache<NSTOKES>&                        m_layercache;

        PostProcessingCache<NSTOKES>&               m_postprocessing_cache;

		bool										m_compute_deriv;
		std::vector<LayerSolution<NSTOKES, CNSTR>>&	m_solutions;
		const InputDerivatives<NSTOKES>&			m_input_derivs;
		LayerDual<double>&						    m_dual_thickness;
        LayerDual<double>&                          m_dual_ssa;
		Dual<double>&								m_average_secant;
		Dual<double>&								m_dual_bt_floor;
		Dual<double>&								m_dual_bt_ceiling;

        // These are used to avoid reallocs, should probably move into ThreadData class along with
        // the above Duals?
        TripleProductDerivativeHolder<NSTOKES>& m_triple_product_holder_0;
        TripleProductDerivativeHolder<NSTOKES>& m_triple_product_holder_1;
        LPTripleProduct<NSTOKES>& m_triple_product;
	};

	
}