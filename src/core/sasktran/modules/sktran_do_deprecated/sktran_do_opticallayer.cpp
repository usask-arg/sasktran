#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_opticallayer.h"

template <int NSTOKES, int CNSTR>
sktran_do_detail::OpticalLayer<NSTOKES, CNSTR>
::OpticalLayer(const PersistentConfiguration<NSTOKES, CNSTR>& config,
			   LayerIndex index,
			   double scat_ext,
			   double tot_ext,
			   std::unique_ptr<VectorDim1<sktran_do_detail::LegendreCoefficient<NSTOKES>>> phasef_expansion,
			   double optical_depth_ceiling,
			   double optical_depth_floor,
			   double altitude_ceiling,
	           double altitude_floor,
	           const InputDerivatives<NSTOKES>& input_derivatives):
	OpticalLayerROP<NSTOKES>(config),
	M_SSA(scat_ext/tot_ext),
	M_SCAT_EXT(scat_ext),
	M_TOT_EXT(tot_ext),
	M_OPTICALDEPTH_CEILING(optical_depth_ceiling),
	M_OPTICALDEPTH_FLOOR(optical_depth_floor),
	M_OPTICAL_THICKNESS(optical_depth_floor - optical_depth_ceiling),
	M_ALT_CEILING(altitude_ceiling),
	M_ALT_FLOOR(altitude_floor),
	m_lephasef(std::move(phasef_expansion)),
	M_INDEX(index),
	m_legendre_sum(config.nstr(), scat_ext / tot_ext, *this->M_LP_MU, *m_lephasef, config.pool().thread_data().legendre_sum_storage(index)),
	m_solutions(config.pool().thread_data().rte_solution(index)),
	m_compute_deriv(false),
	m_input_derivs(input_derivatives),
    m_layercache(config.pool().thread_data().layer_cache(index)),
    m_triple_product_holder(m_layercache.triple_product_holder),
    m_triple_product(m_layercache.triple_product),
    m_dual_thickness(m_layercache.dual_thickness),
    m_average_secant(m_layercache.average_secant),
    m_dual_bt_ceiling(m_layercache.dual_bt_ceiling),
    m_dual_bt_floor(m_layercache.dual_bt_floor)
{
	// Check that SSA is not approximately zero. If so then emit warning once and then dither SSA.
	double ssa_dither = this->m_userspec->getSSAEqual1Dither();
	if(1 - M_SSA < ssa_dither) {
		const_cast<double&>(M_SSA) = 1 - ssa_dither;
		m_legendre_sum.adjustSSA(M_SSA);
	}
	// Configure azimuthal expansion 
	registerAzimuthDependency(m_legendre_sum);
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayer<NSTOKES, CNSTR>::configureDerivative() {
	m_dual_thickness.resize(m_input_derivs.numDerivativeLayer(M_INDEX));
	m_dual_thickness.layer_index = M_INDEX;
	m_dual_thickness.layer_start = (uint)m_input_derivs.layerStartIndex(M_INDEX);

	m_dual_thickness.value = M_OPTICAL_THICKNESS;
	for (uint i = 0; i < m_input_derivs.numDerivativeLayer(M_INDEX); ++i) {
		m_dual_thickness.deriv(i) = m_input_derivs.layerDerivatives()[m_dual_thickness.layer_start + i].d_optical_depth;
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayer<NSTOKES, CNSTR>::configurePseudoSpherical(const Dual<double>& od_top, const Dual<double>& od_bottom) {
	m_dual_bt_ceiling.resize(od_top.deriv.size(), false);
	m_dual_bt_floor.resize(od_bottom.deriv.size(), false);
	m_average_secant.resize(od_top.deriv.size(), false);

	m_dual_bt_ceiling.value = std::exp(-1*od_top.value);
	m_dual_bt_ceiling.deriv = m_dual_bt_ceiling.value * -1 * od_top.deriv;

	m_dual_bt_floor.value = std::exp(-1 * od_bottom.value);
	m_dual_bt_floor.deriv = m_dual_bt_floor.value * -1 * od_bottom.deriv;

	m_average_secant.value = (od_bottom.value - od_top.value) / m_dual_thickness.value;
	m_average_secant.deriv = (od_bottom.deriv - od_top.deriv) / m_dual_thickness.value;

	const auto seq = Eigen::seq(m_dual_thickness.layer_start, m_dual_thickness.layer_start + m_dual_thickness.deriv.size() - 1);

    if (m_dual_thickness.deriv.size() > 0) {
        m_average_secant.deriv(seq) += m_dual_thickness.deriv * od_top.value / (m_dual_thickness.value * m_dual_thickness.value) -
            m_dual_thickness.deriv * od_bottom.value / (m_dual_thickness.value * m_dual_thickness.value);
    }

}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayer<NSTOKES, CNSTR>::takeDerivative(bool take_deriv)
{
	m_compute_deriv = take_deriv;
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayer<NSTOKES, CNSTR>::vectordual_scatPhaseF(AEOrder m,
                                                                    const std::vector<sktran_do_detail::LegendrePhaseContainer<NSTOKES>>& lp_out,
                                                                    const InputDerivatives<NSTOKES>& in_deriv,
                                                                    VectorLayerDual<double>& result_positive_stream,
                                                                    VectorLayerDual<double>& result_negative_stream) const
{
    uint derivStart = (uint)in_deriv.layerStartIndex(M_INDEX);
    uint numDeriv = (uint)in_deriv.numDerivativeLayer(M_INDEX);

    // NSTOKES * NSTOKES matrix for each N
    result_positive_stream.resize(this->M_NSTR / 2 * NSTOKES * NSTOKES, numDeriv, M_INDEX, derivStart);
    result_negative_stream.resize(this->M_NSTR / 2 * NSTOKES * NSTOKES, numDeriv, M_INDEX, derivStart);

    typename std::conditional<NSTOKES==1, double, Eigen::Matrix<double, NSTOKES, NSTOKES>>::type derivtemp;

    for (StreamIndex j = 0; j < this->M_NSTR / 2; ++j)
    {
        m_triple_product.calculate(m, *m_lephasef, lp_out, (*this->M_LP_MU)[m][j]);
        m_triple_product.negations_derivative_emplace(0, m_triple_product_holder);

        // Add scattering parameters
        m_triple_product_holder.ssa = M_SSA;
        m_triple_product_holder.value *= 0.5 * M_SSA;

        if constexpr(NSTOKES == 1) {
            m_triple_product_holder.d_by_legendre_coeff *= 0.5 * M_SSA;
        }
        if constexpr(NSTOKES > 2) {
            m_triple_product_holder.a1deriv *= 0.5 * M_SSA;
            m_triple_product_holder.a2deriv *= 0.5 * M_SSA;
            m_triple_product_holder.a3deriv *= 0.5 * M_SSA;
            m_triple_product_holder.b1deriv *= 0.5 * M_SSA;
        }
        if constexpr(NSTOKES > 3) {
            m_triple_product_holder.b2deriv *= 0.5 * M_SSA;
            m_triple_product_holder.a4deriv *= 0.5 * M_SSA;
        }
        if constexpr(NSTOKES ==1) {
            result_positive_stream.value(j) = m_triple_product_holder.value;
        } else {
            for(int s1 = 0; s1 < NSTOKES; ++s1) {
                for(int s2 = 0; s2 < NSTOKES; ++s2) {
                    int linearindex = j * NSTOKES * NSTOKES + s1 + NSTOKES * s2;
                    result_positive_stream.value(linearindex) = m_triple_product_holder.value(s1, s2);
                }
            }
        }

        for (uint i = 0; i < numDeriv; ++i) {
            m_triple_product_holder.reduce(in_deriv[i + derivStart], derivtemp);

            if constexpr(NSTOKES == 1) {
                result_positive_stream.deriv(i, j) = derivtemp;
            } else {
                for(int s1 = 0; s1 < NSTOKES; ++s1) {
                    for(int s2 = 0; s2 < NSTOKES; ++s2) {
                        int linearindex = j * NSTOKES * NSTOKES + s1 + NSTOKES * s2;
                        result_positive_stream.deriv(i,linearindex) = derivtemp(s1, s2);
                    }
                }
            }
        }

        m_triple_product.negations_derivative_emplace(1, m_triple_product_holder);

        // Add scattering parameters
        m_triple_product_holder.ssa = M_SSA;
        m_triple_product_holder.value *= 0.5 * M_SSA;

        if constexpr(NSTOKES == 1) {
            m_triple_product_holder.d_by_legendre_coeff *= 0.5 * M_SSA;
        }

        if constexpr(NSTOKES > 2) {
            m_triple_product_holder.a1deriv *= 0.5 * M_SSA;
            m_triple_product_holder.a2deriv *= 0.5 * M_SSA;
            m_triple_product_holder.a3deriv *= 0.5 * M_SSA;
            m_triple_product_holder.b1deriv *= 0.5 * M_SSA;
        }
        if constexpr(NSTOKES > 3) {
            m_triple_product_holder.a4deriv *= 0.5 * M_SSA;
            m_triple_product_holder.b2deriv *= 0.5 * M_SSA;
        }

        if constexpr(NSTOKES == 1) {
            result_negative_stream.value(j) = m_triple_product_holder.value;

        } else {
            for(int s1 = 0; s1 < NSTOKES; ++s1) {
                for(int s2 = 0; s2 < NSTOKES; ++s2) {
                    int linearindex = j * NSTOKES * NSTOKES + s1 + NSTOKES * s2;
                    result_negative_stream.value(linearindex) = m_triple_product_holder.value(s1, s2);
                }
            }
        }


        for (uint i = 0; i < numDeriv; ++i) {
            m_triple_product_holder.reduce(in_deriv[i + derivStart], derivtemp);

            if constexpr(NSTOKES == 1) {
                result_negative_stream.deriv(i, j) = derivtemp;

            } else {
                for(int s1 = 0; s1 < NSTOKES; ++s1) {
                    for(int s2 = 0; s2 < NSTOKES; ++s2) {
                        int linearindex = j * NSTOKES * NSTOKES + s1 + NSTOKES * s2;
                        result_negative_stream.deriv(i,linearindex) = derivtemp(s1, s2);
                    }
                }
            }


        }
    }

}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::OpticalLayer<NSTOKES, CNSTR>
::singleScatST(AEOrder m,
               const sktran_do_detail::VectorDim1<LegendrePhaseContainer<NSTOKES>>& lp_out,
               sktran_do_detail::InhomogeneousSourceHolder<NSTOKES>& result,
               sktran_do_detail::InhomogeneousSourceHolder<NSTOKES>& result_negative_coszenith) const
{
    m_triple_product.calculate(m, *m_lephasef, lp_out, (*this->M_LP_CSZ)[m]);
    m_triple_product.negations_derivative_emplace(1, m_triple_product_holder);

    if constexpr(NSTOKES == 1) {
        result.value = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * M_SSA * (2 - kronDelta(m, 0)) * m_triple_product_holder.value;
        result.d_by_legendre_coeff = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * M_SSA * (2 - kronDelta(m, 0)) * m_triple_product_holder.d_by_legendre_coeff;

        m_triple_product.negations_derivative_emplace(2, m_triple_product_holder);

        result_negative_coszenith.value = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * M_SSA * (2 - kronDelta(m, 0)) * m_triple_product_holder.value;
        result_negative_coszenith.d_by_legendre_coeff = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * M_SSA * (2 - kronDelta(m, 0)) * m_triple_product_holder.d_by_legendre_coeff;

        // And the derivative with respect to SSA is equal to the value divided by SSA
        result.d_by_ssa = result.value / M_SSA;
        result_negative_coszenith.d_by_ssa= result_negative_coszenith.value / M_SSA;
    } else {
        double factor = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * M_SSA * (2 - kronDelta(m, 0));

        result.value = factor * m_triple_product_holder.value(Eigen::all, 0);
        result.d_by_a1 = factor * m_triple_product_holder.a1deriv;
        result.d_by_b1_first = factor * m_triple_product_holder.b1deriv(Eigen::all, 0);
        result.d_by_b1_second = factor * m_triple_product_holder.b1deriv(Eigen::all, 1);
        result.d_by_ssa = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * (2 - kronDelta(m, 0)) * m_triple_product_holder.value(Eigen::all, 0);

        m_triple_product.negations_derivative_emplace(2, m_triple_product_holder);

        result_negative_coszenith.value = factor * m_triple_product_holder.value(Eigen::all, 0);
        result_negative_coszenith.d_by_a1 = factor * m_triple_product_holder.a1deriv;
        result_negative_coszenith.d_by_b1_first = factor * m_triple_product_holder.b1deriv(Eigen::all, 0);
        result_negative_coszenith.d_by_b1_second = factor * m_triple_product_holder.b1deriv(Eigen::all, 1);
        result_negative_coszenith.d_by_ssa = (1.0 / (4 * PI)) * this->M_SOLAR_DIRECT_INTENSITY * (2 - kronDelta(m, 0)) * m_triple_product_holder.value(Eigen::all, 0);
    }

}

INSTANTIATE_TEMPLATE(sktran_do_detail::OpticalLayer);
