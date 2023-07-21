#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_postprocessing.h"

template <int NSTOKES, int CNSTR>
sktran_do_detail::ParticipatingSourceTerm<NSTOKES, CNSTR>::ParticipatingSourceTerm(AEOrder m, const OpticalLayer<NSTOKES, CNSTR>* layer, const OpticalLayerArray<NSTOKES, CNSTR>& layers, double mu, double obsod, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp_mu, bool ssonly, bool includess)
	: STIntegralProperties<NSTOKES>(*layer),
	m_layer(*layer),
	m_layers(layers),
	m_solution(m_layer.solution(m)),
	M_FEORDER(m),
	M_COSZEN(mu),
	M_LP_COSZEN(lp_mu),
	M_OBSOD(obsod),
	M_EXIT_OPTICAL_DEPTH(mu > 0 ? m_layer.opticalDepth(Location::CEILING) : m_layer.opticalDepth(Location::FLOOR)),
	M_OPTICAL_THICKNESS(m_layer.opticalDepth(Location::INSIDE)),
	M_SSONLY(ssonly),
    M_INCLUDESS(includess)
{
	if(M_COSZEN < 0) {
		throw InternalError("ParticipatingSourceTerm: downward's calculations have never been tested so don't trust it without testing!");
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::ParticipatingSourceTerm<NSTOKES, CNSTR>::integrate(double upper_bound, Radiance<NSTOKES>& result, PostProcessingCache<NSTOKES>& cache, const sktran_do_detail::Radiance<NSTOKES>* manual_ss_source) const
{
    // Note that for the source function integration we switch convention for upwelling/downwelling
    // from the solution convention
    bool use_green_function = m_solution.value.use_green_function();

    LayerIndex p = m_layer.index();
    size_t numderiv = m_layers.inputDerivatives().numDerivativeLayer(p);
    size_t numtotalderiv = m_layers.inputDerivatives().numDerivative();
    uint layerStart = (uint)m_layers.inputDerivatives().layerStartIndex(p);

    cache.resize(this->M_NSTR, p, numderiv, layerStart, numtotalderiv);

    // Check upper bound is valid if we are not integrating across the entire layer
    if (M_EXIT_OPTICAL_DEPTH != upper_bound && !verifyIntegrationUpperBound(upper_bound)) {
        throw InternalError("Participating source term integral upperbound is not valid");
    }

    double x = upper_bound - M_EXIT_OPTICAL_DEPTH;
    const auto& transmission = x == 0 ? m_layer.dual_beamTransmittance(Location::CEILING, m_layers.inputDerivatives()) : m_layer.dual_beamTransmittance(Location::INSIDE, m_layers.inputDerivatives(), x);

    const auto& average_secant = m_layer.dual_average_secant();
    const LayerDual<double>& dual_thickness = m_layer.dual_thickness();

    // Calculate legendre sums multiplied by stream weights
    VectorLayerDual<double>& dual_lpsum_plus = cache.dual_lpsum_plus;
    VectorLayerDual<double>& dual_lpsum_minus = cache.dual_lpsum_minus;

    m_layer.vectordual_scatPhaseF(M_FEORDER, M_LP_COSZEN, m_layers.inputDerivatives(), dual_lpsum_minus, dual_lpsum_plus);

    // dual_lpsum_plus and minus are now linearly storage of N * NSTOKES * NSTOKES, have to multiply the blocks by
    // the weights
    for(int i = 0; i < this->M_NSTR/2; ++i) {
        for(int j = 0; j < NSTOKES*NSTOKES; ++j) {
            int linearindex = i*NSTOKES*NSTOKES + j;
            dual_lpsum_minus.value(linearindex) *= (*this->M_WT)[i];
            dual_lpsum_plus.value(linearindex) *= (*this->M_WT)[i];

            dual_lpsum_minus.deriv(Eigen::all, linearindex) *= (*this->M_WT)[i];
            dual_lpsum_plus.deriv(Eigen::all, linearindex) *= (*this->M_WT)[i];
        }
    }

    const auto& dual_particular_plus = m_solution.value.dual_particular_plus();
    const auto& dual_particular_minus = m_solution.value.dual_particular_minus();

    const auto& dual_Aplus = m_solution.value.dual_green_A_plus();
    const auto& dual_Aminus = m_solution.value.dual_green_A_minus();

    // Eq (44)
    auto& V = cache.V;

    // Eq (43)
    auto& Qtemp = cache.Qtemp;
    auto& temp = cache.temp;
    auto& Q = cache.Q;

    if (M_INCLUDESS) {
        m_layer.singleScatST(M_FEORDER, M_LP_COSZEN, Qtemp, temp);

        Q.value = Qtemp.value;
        Q.deriv.setZero();

        typedef typename std::conditional<NSTOKES==1, double, Eigen::Vector<double, NSTOKES>>::type DerivType;

        DerivType derivtemp;
        for(uint k = 0; k < numderiv; ++k) {
            Qtemp.reduce(m_layers.inputDerivatives().layerDerivatives()[layerStart + k], derivtemp);
            if constexpr(NSTOKES==1) {
                Q.deriv(k + layerStart) = derivtemp;
            } else {
                Q.deriv(k + layerStart, Eigen::all) = derivtemp;
            }
        }
    } else {
        if( manual_ss_source == nullptr) {
            Q.setzero();
        }
        else {
            // Can just reassign Q to the manual source
            Q.value = (*manual_ss_source).value;
            Q.deriv = (*manual_ss_source).deriv;
        }
    }

    auto& hp = cache.hp;
    auto& hm = cache.hm;
    auto& J = cache.J;

    auto& Y_plus = cache.Y_plus;
    auto& Y_minus = cache.Y_minus;

    using MatrixView = Eigen::Map<Eigen::MatrixXd>;
    using ConstMatrixView = Eigen::Map<const Eigen::MatrixXd>;
    MatrixView Y_plus_matrix(Y_plus.value.data(), NSTOKES, this->M_NSTR/2 *NSTOKES);
    MatrixView Y_minus_matrix(Y_minus.value.data(), NSTOKES, this->M_NSTR/2 *NSTOKES);

    ConstMatrixView lpsum_plus_matrix(dual_lpsum_plus.value.data(), NSTOKES, this->M_NSTR/2 * NSTOKES);
    ConstMatrixView lpsum_minus_matrix(dual_lpsum_minus.value.data(), NSTOKES, this->M_NSTR/2 * NSTOKES);

    ConstMatrixView homog_plus_matrix(m_solution.value.dual_homog_plus().value.data(), this->M_NSTR/2 * NSTOKES, this->M_NSTR/2 *NSTOKES);
    ConstMatrixView homog_minus_matrix(m_solution.value.dual_homog_minus().value.data(), this->M_NSTR/2 * NSTOKES, this->M_NSTR/2 *NSTOKES);

    Y_plus_matrix = lpsum_plus_matrix * homog_plus_matrix + lpsum_minus_matrix * homog_minus_matrix;
    Y_minus_matrix = lpsum_plus_matrix * homog_minus_matrix + lpsum_minus_matrix * homog_plus_matrix;

    // Assign derivatives of Y_plus and Y_minus
    // Only have layer derivatives
    std::vector<Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>> Y_plus_deriv;
    std::vector<Eigen::Map<Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>>> Y_minus_deriv;

    Y_plus_deriv.reserve(numderiv);
    Y_minus_deriv.reserve(numderiv);

    for(uint k = 0; k < numderiv; ++k) {
        Y_plus_deriv.emplace_back(&Y_plus.deriv(k, 0), NSTOKES, this->M_NSTR/2 * NSTOKES, Eigen::InnerStride<>(numderiv));
        Y_minus_deriv.emplace_back(&Y_minus.deriv(k, 0), NSTOKES, this->M_NSTR/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> lpsum_plus_deriv(
                &dual_lpsum_plus.deriv(k, 0), NSTOKES, this->M_NSTR/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> lpsum_minus_deriv(
                &dual_lpsum_minus.deriv(k, 0), NSTOKES, this->M_NSTR/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> homog_minus_deriv(
                &m_solution.value.dual_homog_minus().deriv(k, 0), this->M_NSTR/2 * NSTOKES, this->M_NSTR/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

        Eigen::Map<const Eigen::MatrixXd, 0, Eigen::InnerStride<Eigen::Dynamic>> homog_plus_deriv(
                &m_solution.value.dual_homog_plus().deriv(k, 0), this->M_NSTR/2 * NSTOKES, this->M_NSTR/2 * NSTOKES, Eigen::InnerStride<>(numderiv));

        Y_plus_deriv[k].noalias() = lpsum_plus_deriv * homog_plus_matrix + lpsum_plus_matrix * homog_plus_deriv +
                lpsum_minus_deriv * homog_minus_matrix + lpsum_minus_matrix * homog_minus_deriv;

        Y_minus_deriv[k].noalias() = lpsum_plus_deriv * homog_minus_matrix + lpsum_plus_matrix * homog_minus_deriv +
                lpsum_minus_deriv * homog_plus_matrix + lpsum_minus_matrix * homog_plus_deriv;
    }

    auto& Dm = cache.Dm;
    auto& Dp = cache.Dp;
    auto& Eform = cache.Eform;

    E(x, M_OPTICAL_THICKNESS, transmission, Eform);

    const auto& dual_L = m_solution.boundary.L_coeffs;
    const auto& dual_M = m_solution.boundary.M_coeffs;

    J.setzero();
    V.setzero();

    for (SolutionIndex i = 0; i < this->M_NSTR / 2 * NSTOKES; ++i) {
        h_plus(x, M_OPTICAL_THICKNESS, i, hp);
        h_minus(x, M_OPTICAL_THICKNESS, i, hm);

        if constexpr(NSTOKES==1) {
            J.value += Y_plus_matrix(0, i) * hp.value * dual_L.value(i);
            J.value += Y_minus_matrix(0, i) * hm.value * dual_M.value(i);
        } else {
            J.value += Y_plus_matrix(Eigen::all, i) * hp.value * dual_L.value(i);
            J.value += Y_minus_matrix(Eigen::all, i) * hm.value * dual_M.value(i);
        }

        // Y and Hp/Hm only have layer derivatives
        for(uint k = 0; k < numderiv; ++k) {
            J.deriv(layerStart + k, Eigen::all).noalias() += Y_plus_deriv[k](Eigen::all, i) * hp.value * dual_L.value(i);
            J.deriv(layerStart + k, Eigen::all).noalias() += Y_minus_deriv[k](Eigen::all, i) * hm.value * dual_M.value(i);

            J.deriv(layerStart + k, Eigen::all).noalias() += Y_plus_matrix(Eigen::all, i) * hp.deriv(k) * dual_L.value(i);
            J.deriv(layerStart + k, Eigen::all).noalias() += Y_minus_matrix(Eigen::all, i) * hm.deriv(k) * dual_M.value(i);
        }
        // But L/M have full derivatives
        for(uint k = 0; k < numtotalderiv; ++k) {
            J.deriv(k, Eigen::all).noalias() += Y_plus_matrix(Eigen::all, i) * hp.value * dual_L.deriv(k, i);
            J.deriv(k, Eigen::all).noalias() += Y_minus_matrix(Eigen::all, i) * hm.value * dual_M.deriv(k, i);
        }

        const auto& eigval = m_solution.value.dual_eigval();

        double expfactor = exp(-1*(dual_thickness.value - x) * average_secant.value);

        Dp.value = (-1*transmission.value * expfactor * hm.value + Eform.value) / (average_secant.value + m_solution.value.dual_eigval().value(i));
        Dm.value = (transmission.value * hp.value - Eform.value) / (average_secant.value - m_solution.value.dual_eigval().value(i));

        Dp.deriv.setZero();
        Dm.deriv.setZero();

        Dp.deriv.noalias() += average_secant.deriv * (transmission.value * hm.value * dual_thickness.value * expfactor - Dp.value) / (average_secant.value + eigval.value(i));
        Dp.deriv.noalias() += Eform.deriv / (average_secant.value + eigval.value(i));
        Dp.deriv.noalias() += transmission.deriv * (-hm.value * expfactor) / (average_secant.value + eigval.value(i));

        Dm.deriv.noalias() += average_secant.deriv * (-Dm.value) / (average_secant.value - eigval.value(i));
        Dm.deriv.noalias() += transmission.deriv * (hp.value) / (average_secant.value - eigval.value(i));
        Dm.deriv.noalias() += Eform.deriv * -1.0 / (average_secant.value - eigval.value(i));

        for(uint k = 0; k < numderiv; ++k) {
            Dp.deriv(layerStart + k) += eigval.deriv(k, i) * (-Dp.value) / (average_secant.value + eigval.value(i));
            Dp.deriv(layerStart + k) += hm.deriv(k) * (-transmission.value * expfactor) / (average_secant.value + eigval.value(i));
            Dp.deriv(layerStart + k) += (dual_thickness.value - x) / dual_thickness.value * dual_thickness.deriv(k) * (transmission.value * average_secant.value * hm.value * expfactor) / (average_secant.value + eigval.value(i));

            Dm.deriv(layerStart + k) += hp.deriv(k) * transmission.value / (average_secant.value - eigval.value(i));
            Dm.deriv(layerStart + k) += eigval.deriv(k, i) * (Dm.value) / (average_secant.value - eigval.value(i));
        }

        if constexpr(NSTOKES==1) {
            V.value += dual_Aplus.value(i) * Y_plus_matrix(0, i) * Dm.value + dual_Aminus.value(i) * Y_minus_matrix(0, i) * Dp.value;
        } else {
            V.value += dual_Aplus.value(i) * Y_plus_matrix(Eigen::all, i) * Dm.value + dual_Aminus.value(i) * Y_minus_matrix(Eigen::all, i) * Dp.value;
        }

        // Y and A only have layer derivatives
        for( uint k = 0; k < numderiv; ++k) {
            V.deriv(layerStart + k, Eigen::all).noalias() += dual_Aplus.value(i) * Y_plus_deriv[k](Eigen::all, i) * Dm.value;
            V.deriv(layerStart + k, Eigen::all).noalias() += dual_Aminus.value(i) * Y_minus_deriv[k](Eigen::all, i) * Dp.value;

			V.deriv(layerStart + k, Eigen::all).noalias() += dual_Aplus.deriv(k, i) * Y_plus_matrix(Eigen::all, i) * Dm.value;
			V.deriv(layerStart + k, Eigen::all).noalias() += dual_Aminus.deriv(k, i) * Y_minus_matrix(Eigen::all, i) * Dp.value;
        }
        // But D has full derivatives
        for( uint k = 0; k < numtotalderiv; ++k) {
            V.deriv(k, Eigen::all).noalias() += dual_Aplus.value(i) * Y_plus_matrix(Eigen::all, i) * Dm.deriv(k);
            V.deriv(k, Eigen::all).noalias() += dual_Aminus.value(i) * Y_minus_matrix(Eigen::all, i) * Dp.deriv(k);
        }
    }

    if (M_SSONLY) {
        result.value = Q.value * Eform.value;

        result.deriv.noalias() = Q.deriv * Eform.value;
        for (uint k = 0; k < numtotalderiv; ++k) {
            if constexpr(NSTOKES==1) {
                result.deriv(k) += Q.value * Eform.deriv(k);
            } else {
                result.deriv(k, Eigen::all).noalias() += Q.value * Eform.deriv(k);
            }
        }
    }
    else {
        result.value = J.value + V.value + Q.value * Eform.value;

        result.deriv.noalias() = J.deriv + V.deriv + Q.deriv * Eform.value;
		for (uint k = 0; k < numtotalderiv; ++k) {
            if constexpr(NSTOKES==1) {
                result.deriv(k) += Q.value * Eform.deriv(k);
            } else {
                result.deriv(k, Eigen::all).noalias() += Q.value * Eform.deriv(k);
            }
		}
    }

}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::ParticipatingSourceTerm<NSTOKES, CNSTR>::h_plus(double x,
													            double thickness,
													            SolutionIndex j,
													            LayerDual<HomogType>& xform) const
{
	double mu = std::abs(M_COSZEN);
	uint p = m_layer.index();

	auto v_eigval = m_solution.value.dual_eigval().value(j);
	const auto& d_eigval = m_solution.value.dual_eigval().deriv(Eigen::all, j);

	const LayerDual<double>& dual_thickness = m_layer.dual_thickness();

	// Exit depth is x, 0 is exit at the top of the layer
	// h_layer = h_low + (1 - x/thickness) * (h_upper - h_lower)
	// 0 = -dx * lh / thickness - x/(thickness*thickness) * dthickness * lh
	// dx = x/thickness * dthickness

	double layerfraction = 1 - x / thickness;
	HomogType den = (1.0 + mu * v_eigval);

	if(std::abs(den) > 0.0001) {
		HomogType e1 = std::exp(-x * v_eigval);
		HomogType e2 = std::exp(-1.0*dual_thickness.value * v_eigval) * std::exp(-(dual_thickness.value - x) / mu);

		xform.value = (e1 - e2) / (den);
	
		if (xform.deriv.size() > 0) {
			// Start with the numerator derivative
			xform.deriv.noalias() = (-1.0*e1*(d_eigval*x + dual_thickness.deriv * (1 - layerfraction) * v_eigval) + e2 * ((v_eigval + layerfraction / mu) * dual_thickness.deriv + dual_thickness.value * d_eigval)) / den;

			// Then add in the denominator derivative
			xform.deriv -= 1.0 / (den)* xform.value * mu * d_eigval;
		}
	} else { // limit as eigval -> -1/mu, take two terms in the taylor series expansion
		HomogType e1 = std::exp(-x * v_eigval);
		xform.value = (dual_thickness.value - x) / mu * e1 * (HomogType(1) - (dual_thickness.value - x) * (v_eigval + HomogType(1) / mu));
		if (xform.deriv.size() > 0) {
			xform.deriv.noalias() = -1.0*x*d_eigval * xform.value;
			xform.deriv += 1 / mu * e1 * dual_thickness.deriv * (HomogType(1) - (dual_thickness.value - x) * (v_eigval + HomogType(1) / mu));
			xform.deriv += 1 / mu * e1 * (dual_thickness.value - x) * -1.0 * (v_eigval + 1 / mu) * dual_thickness.deriv;
		}

		// xform = dual_thickness / mu;
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::ParticipatingSourceTerm<NSTOKES, CNSTR>::h_minus(double x,
														         double thickness,
														         SolutionIndex j,
														         LayerDual<HomogType>& xform) const
{
	double mu = std::abs(M_COSZEN);
	uint p = m_layer.index();

	auto v_eigval = m_solution.value.dual_eigval().value(j);
	const auto& d_eigval = m_solution.value.dual_eigval().deriv(Eigen::all, j);

	const LayerDual<double>& dual_thickness = m_layer.dual_thickness();
	double layerfraction = 1 - x / thickness;

	HomogType den = (1.0 - mu * v_eigval);

	if (std::abs(den) > 0.0001) {
		HomogType e1 = std::exp(-1.0*v_eigval*(dual_thickness.value - x));
		HomogType e2 = std::exp(-1/mu * (dual_thickness.value - x));

		xform.value = (e1 - e2) / den;

		if (xform.deriv.size() > 0)
		{
			xform.deriv.noalias() = (-1.0*e1*(d_eigval*(dual_thickness.value - x) + dual_thickness.deriv * layerfraction * v_eigval) + e2 * dual_thickness.deriv * layerfraction / mu) / den;
			xform.deriv += 1.0 / den * xform.value * mu * d_eigval;
		}

	}
	else { // limit as eigval -> -1/mu, take two terms in the taylor series expansion
		HomogType e1 = std::exp(-1.0*v_eigval*(dual_thickness.value - x));

		xform.value = e1 * (dual_thickness.value - x) / mu * (HomogType(1) - (dual_thickness.value - x) * (v_eigval - HomogType(1) / mu));
		if (xform.deriv.size() > 0)
		{
			xform.deriv.noalias() = -1.0*(v_eigval * dual_thickness.deriv + dual_thickness.value * d_eigval) * xform.value;
			xform.deriv.noalias() += e1 * dual_thickness.deriv / mu * (HomogType(1) - (dual_thickness.value - x) * (v_eigval - HomogType(1) / mu));
			xform.deriv.noalias() += e1 * (dual_thickness.value - x) / mu * -1.0 * dual_thickness.deriv * (v_eigval - 1 / mu);
		}
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::ParticipatingSourceTerm<NSTOKES, CNSTR>::E(double x,
												           double thickness,
												           const Dual<double>& transmission,
												           Dual<double>& xform) const
{
	double mu = std::abs(M_COSZEN);
	uint p = m_layer.index();
	uint layerStart = (uint)m_layers.inputDerivatives().layerStartIndex(p);
	uint numDeriv = (uint)m_layers.inputDerivatives().numDerivativeLayer(p);

	const LayerDual<double>& dual_thickness = m_layer.dual_thickness();
	const Dual<double>& avg_secant = m_layer.dual_average_secant();
	double layerfraction = 1 - x / thickness;

	// Start by assigning cross derivative terms of xform
	double e1 = exp(-x * avg_secant.value);
	double e2 = exp(-1.0*dual_thickness.value * avg_secant.value) * exp(-1.0*(dual_thickness.value - x) / mu);
	double den = (1.0 + mu * avg_secant.value);

	xform.value = transmission.value / den * (e1 - e2);
	if (xform.deriv.size() > 0) {
		xform.deriv.noalias() = transmission.deriv / den * (e1 - e2);
		xform.deriv += transmission.value / den * e1 * (-1.0*x * avg_secant.deriv);
		xform.deriv += transmission.value / den * e2 * dual_thickness.value * avg_secant.deriv;
		xform.deriv -= xform.value * avg_secant.deriv / den * mu;
	}

	// Then add in the terms that contain cross derivatives
	if (numDeriv > 0) {
		xform.deriv(Eigen::seq(layerStart, layerStart + numDeriv - 1), Eigen::all) += transmission.value / den * (e2 * dual_thickness.deriv * (avg_secant.value + layerfraction / mu) - e1 * dual_thickness.deriv * (1 - layerfraction) * avg_secant.value);
	}
}

template <int NSTOKES, int CNSTR>
void sktran_do_detail::SphericalPostProcessing<NSTOKES, CNSTR>::add_ground_reflectance(AEOrder m, const TracedRay& ray, Dual<double>& radiance) {
	// Start by determining which CSZ solutions contribute
	std::array<size_t, 2> csz_indicies;
	std::array<double, 2> csz_weights;

	const auto& ground = ray.layers[0].exit;

	interpolation_index(ground.CosSZA(), csz_indicies, csz_weights);

	sktran_do_detail::Radiance<1> temp(m_opticallayers[csz_indicies[0]].inputDerivatives().numDerivative());

	for (size_t i = 0; i < 2; ++i) {
		double w = csz_weights[i];
		if (w == 0) {
			continue;
		}
        if constexpr(NSTOKES == 1) {
            temp = convert_dual_to_wf(m_opticallayers[csz_indicies[i]].reflectedIntensity(m, m_los[0]),
                                      m_opticallayers[csz_indicies[i]].inputDerivatives(), m_numwf);
        } else {
            // TODO: ???
        }

		radiance.value += w*temp.value;
		radiance.deriv += w*temp.deriv;
	}
	radiance *= cos(m*ray.layers[0].saz);
}

template <>
void sktran_do_detail::SphericalPostProcessing<3>::add_ground_reflectance(AEOrder m, const TracedRay& ray, Dual<double>& radiance) {
	// TODO: implement?
}

template <>
void sktran_do_detail::SphericalPostProcessing<4>::add_ground_reflectance(AEOrder m, const TracedRay& ray, Dual<double>& radiance) {
	// TODO: implement?
}

template <>
void sktran_do_detail::SphericalPostProcessing<3>::non_integrated_source_unrolled(AEOrder m, double altitude, double csz, double cos_viewing, bool is_upwelling, Dual<double>& homog_source, Dual<double>& partic_source, bool recalculate_caches) {
	// TODO: Have to implement 
}

template <>
void sktran_do_detail::SphericalPostProcessing<4>::non_integrated_source_unrolled(AEOrder m, double altitude, double csz, double cos_viewing, bool is_upwelling, Dual<double>& homog_source, Dual<double>& partic_source, bool recalculate_caches) {
	// TODO: Have to implement complex solutions
}

// This just needs bo be completely redone and rethought out
template <int NSTOKES, int CNSTR>
void sktran_do_detail::SphericalPostProcessing<NSTOKES, CNSTR>::non_integrated_source_unrolled(AEOrder m, double altitude, double csz, double cos_viewing, bool is_upwelling, Dual<double>& homog_source, Dual<double>& partic_source, bool recalculate_caches) {
	// Start by determining which CSZ solutions contribute
	std::array<size_t, 2> csz_indicies;
	std::array<double, 2> csz_weights;

	interpolation_index(csz, csz_indicies, csz_weights);

	partic_source.value = 0.0;
	partic_source.deriv.setZero();
	homog_source.value = 0.0;
	homog_source.deriv.setZero();

	if (recalculate_caches)
	{
		for (LPOrder l = 0; l < this->M_NSTR; ++l) m_cached_lp_mu[l].fill(m, l, cos_viewing);
	}

	double downwellingscale = 1.0;
	double upwellingscale = 1.0;

	if (m_use_upwelling_spher)
	{
		double earthrad = m_config.coords()->AltitudeToRadius(m_config.coords()->GroundAltitude());

		double sincrit = earthrad / (earthrad + altitude);
		double theta = asin(sincrit);

		downwellingscale = theta / (nxmath::Pi / 2);
		upwellingscale = (nxmath::Pi - theta) / (nxmath::Pi / 2);
	}

	for (size_t i = 0; i < 2; ++i) {
		if (csz_weights[i] == 0.0) {
			continue;
		}

		double w = csz_weights[i];
		const auto& layers = m_opticallayers[csz_indicies[i]];
		const auto& rte = m_rte[csz_indicies[i]];

		Dual<double>& layer_partic_source = m_cached_layer_part_source;
		Dual<double>& layer_homog_source = m_cached_layer_homog_source;

		layer_partic_source.value = 0.0;
		layer_partic_source.deriv.setZero();
		layer_homog_source.deriv.setZero();
		layer_homog_source.value = 0.0;

		// Find the layer corresponding to this altitude
		double optical_depth = layers.opticalDepthAt(altitude);
		const auto layer = layers.layerAt(optical_depth);

        const auto& average_secant = layer->dual_average_secant();
        const auto& thickness = layer->dual_thickness();
        const auto& transmission = layer->dual_beamTransmittance(Location::CEILING, layers.inputDerivatives());

		// We need the partial optical depth within the layer, ideally we would use the optical depth
		// we just calculated but that is linearized with respect to input extinctions and
		// our layer terms are linearized with respect to layer quantities
		// So we will just approximate the partial layer optical depth as
		// od = od_ceiling + layer_thickness / layer_width * (h - h_floor);
		// x = layer_thickness / layer_width * (h - h_floor)

		// x is the partial layer optical depth relative to the top of the layer
		// LayerDual x = dual_thickness * layerfraction
		double layer_fraction = (layer->altitude(Location::CEILING) - altitude) / (layer->altitude(Location::CEILING) - layer->altitude(Location::FLOOR));

		const auto& solution = layer->solution(m);

		// Configure layer quantities
		LayerIndex p = layer->index();
		size_t numderiv = layers.inputDerivatives().numDerivativeLayer(p);
		size_t numtotalderiv = layers.inputDerivatives().numDerivative();
		uint layerStart = (uint)layers.inputDerivatives().layerStartIndex(p);

		// Calculate the source terms
		// Calculate legendre sums multiplied by stream weights
		VectorLayerDual<double>& dual_lpsum_plus = m_cached_lpsum_plus;
		VectorLayerDual<double>& dual_lpsum_minus = m_cached_lpsum_minus;
		if (recalculate_caches) {
			layer->vectordual_scatPhaseF(m, m_cached_lp_mu, layers.inputDerivatives(), dual_lpsum_minus, dual_lpsum_plus);
			dual_lpsum_plus *= *this->M_WT;
			dual_lpsum_minus *= *this->M_WT;
		}

        if(!solution.value.use_green_function()) {
            // Downwelling layers swap these terms
            const auto& dual_particular_plus = is_upwelling ? solution.value.dual_particular_plus() : solution.value.dual_particular_minus();
            const auto& dual_particular_minus = is_upwelling ? solution.value.dual_particular_minus() : solution.value.dual_particular_plus();

            double x = thickness.value * layer_fraction;

            // Eq (44)
            for (int j = 0; j < this->M_NSTR / 2; ++j) {
                layer_homog_source.value += upwellingscale * dual_lpsum_plus.value(j) * dual_particular_plus.value(j) * transmission.value * exp(-x*average_secant.value);
                layer_homog_source.value += downwellingscale * dual_lpsum_minus.value(j) * dual_particular_minus.value(j) * transmission.value * exp(-x*average_secant.value);

                for (int k = 0; k < numtotalderiv; ++k) {
                    layer_homog_source.deriv(k) += (dual_particular_plus.deriv(k, j) * dual_lpsum_plus.value(j) * upwellingscale * transmission.value * exp(-x*average_secant.value) +
                                                     dual_particular_minus.deriv(k, j) * dual_lpsum_minus.value(j) * downwellingscale * transmission.value * exp(-x*average_secant.value));
                    layer_homog_source.deriv(k) += upwellingscale * dual_lpsum_plus.value(j) * dual_particular_plus.value(j) * transmission.deriv(k) * exp(-x*average_secant.value);
                    layer_homog_source.deriv(k) += downwellingscale * dual_lpsum_minus.value(j) * dual_particular_minus.value(j) * transmission.deriv(k) * exp(-x*average_secant.value);
                    layer_homog_source.deriv(k) += upwellingscale * dual_lpsum_plus.value(j) * dual_particular_plus.value(j) * transmission.value * exp(-x*average_secant.value) * -x * average_secant.deriv(k);
                    layer_homog_source.deriv(k) += downwellingscale * dual_lpsum_minus.value(j) * dual_particular_minus.value(j) * transmission.value * exp(-x*average_secant.value) * -x * average_secant.deriv(k);

                }
                for (int k = 0; k < numderiv; ++k) {
                    layer_homog_source.deriv(k + layerStart) += (dual_lpsum_minus.deriv(k, j) * dual_particular_minus.value(j) * downwellingscale * transmission.value * exp(-x*average_secant.value) +
                                                                  dual_lpsum_plus.deriv(k, j) * dual_particular_plus.value(j) * upwellingscale * transmission.value * exp(-x*average_secant.value));

                    layer_homog_source.deriv(k) += upwellingscale * dual_lpsum_plus.value(j) * dual_particular_plus.value(j) * transmission.value * exp(-x*average_secant.value) * -average_secant.value * thickness.deriv(k) * layer_fraction;
                    layer_homog_source.deriv(k) += downwellingscale * dual_lpsum_minus.value(j) * dual_particular_minus.value(j) * transmission.value * exp(-x*average_secant.value) * -average_secant.value * thickness.deriv(k) * layer_fraction;

                }
            }
        }

        Dual<double> Dm(layers.inputDerivatives().numDerivative());
        Dual<double> Dp(layers.inputDerivatives().numDerivative());

		// Now calculate the homogeneous source
		LayerDual<double> Y_plus(numderiv, p, layerStart), Y_minus(numderiv, p, layerStart);
		LayerDual<double> hp(numderiv, p, layerStart), hm(numderiv, p, layerStart);

		for (SolutionIndex j = 0; j < this->M_NSTR / 2; ++j) {
			const auto& dual_homog_plus = solution.value.dual_homog_plus();
			const auto& dual_homog_minus = solution.value.dual_homog_minus();
			size_t start_idx = j * this->M_NSTR / 2;
			size_t end_idx = (j + 1)*(this->M_NSTR / 2);
			auto v_eigval = solution.value.dual_eigval().value(j);
			const auto& d_eigval = solution.value.dual_eigval().deriv(Eigen::all, j);

			Y_plus.value = 0.0;
			Y_plus.deriv.setZero();
			Y_minus.value = 0.0;
			Y_minus.deriv.setZero();

			for (size_t k = start_idx; k < end_idx; ++k) {
				// source = (Yp dot L * hp) + (Ym dot M * hp) for this solution
				// Ypm = lpsum_p dot homog_pm + lpsum_m dot homog_mp
				Y_plus.value += upwellingscale * dual_lpsum_plus.value(k - start_idx) * dual_homog_plus.value(k);
				Y_plus.value += downwellingscale * dual_lpsum_minus.value(k - start_idx) * dual_homog_minus.value(k);

				Y_minus.value += upwellingscale * dual_lpsum_plus.value(k - start_idx) * dual_homog_minus.value(k);
				Y_minus.value += downwellingscale * dual_lpsum_minus.value(k - start_idx) * dual_homog_plus.value(k);

				for (size_t n = 0; n < numderiv; ++n) {
					Y_plus.deriv(n) += upwellingscale * (dual_lpsum_plus.deriv(n, k - start_idx) * dual_homog_plus.value(k) +
						dual_homog_plus.deriv(n, k) * dual_lpsum_plus.value(k - start_idx));
					Y_plus.deriv(n) += downwellingscale * (dual_lpsum_minus.deriv(n, k - start_idx) * dual_homog_minus.value(k) +
						dual_homog_minus.deriv(n, k) * dual_lpsum_minus.value(k - start_idx));

					Y_minus.deriv(n) += upwellingscale * (dual_lpsum_plus.deriv(n, k - start_idx) * dual_homog_minus.value(k) +
						dual_homog_minus.deriv(n, k) * dual_lpsum_plus.value(k - start_idx));
					Y_minus.deriv(n) += downwellingscale * (dual_lpsum_minus.deriv(n, k - start_idx) * dual_homog_plus.value(k) +
						dual_homog_plus.deriv(n, k) * dual_lpsum_minus.value(k - start_idx));
				}
			}

			hp.value = std::exp(-1.0*v_eigval * layer->dual_thickness().value * layer_fraction);
			hm.value = std::exp(-v_eigval * (layer->dual_thickness().value - layer->dual_thickness().value * layer_fraction));
			for (size_t n = 0; n < numderiv; ++n) {
				hp.deriv(n) = -1.0*hp.value*layer_fraction*(v_eigval * layer->dual_thickness().deriv(n) + d_eigval(n) * layer->dual_thickness().value);
				hm.deriv(n) = -1.0*hm.value*(1 - layer_fraction)*(v_eigval * layer->dual_thickness().deriv(n) + d_eigval(n) * layer->dual_thickness().value);
			}

			const auto& hpp = is_upwelling ? hp : hm;
			const auto& hmm = is_upwelling ? hm : hp;

			layer_homog_source.value += Y_plus.value * solution.boundary.L_coeffs.value(j) * hpp.value;
			layer_homog_source.value += Y_minus.value * solution.boundary.M_coeffs.value(j) * hmm.value;

			for (size_t n = 0; n < numderiv; ++n) {
				layer_homog_source.deriv(n + layerStart) +=  (Y_plus.deriv(n) * solution.boundary.L_coeffs.value(j) * hpp.value +
					Y_plus.value * solution.boundary.L_coeffs.value(j) * hpp.deriv(n));

				layer_homog_source.deriv(n + layerStart) +=  (Y_minus.deriv(n) * solution.boundary.M_coeffs.value(j) * hmm.value +
					Y_minus.value * solution.boundary.M_coeffs.value(j) * hmm.deriv(n));
			}

			for (size_t n = 0; n < numtotalderiv; ++n) {
				layer_homog_source.deriv(n) +=  (Y_plus.value * solution.boundary.L_coeffs.deriv(n, j) * hpp.value);
				layer_homog_source.deriv(n) += (Y_minus.value * solution.boundary.M_coeffs.deriv(n, j) * hmm.value);
			}

            if(solution.value.use_green_function()) {
                // We do not include solar transmission because a spherical solar transmission is added on later
                // by the other integrator
                // Derives from greens function C +/- that are multiplied by exp(x lambda) to remove the solar
                // dependence
                // Dp = 1 / (lambda + k) * (1 - exp(-(Delta-x) (lambda+k))
                // Dm = 1 / (lambda - k) * (exp(-x (k-lambda)) - 1)

                // Or
                // Dp = 1 / (lambda + k) * (1 - exp(-Delta lambda) Hp)
                // Dm = 1 / (lambda - k) * (exp(-x lambda) Hm - 1)

                // Remember x = Delta * layer_fraction
                // Then (Delta-x) = Delta (1 - layerfraction)
                double x = thickness.value * layer_fraction;

                Dp.value = transmission.value / (v_eigval + average_secant.value) * (exp(-x*average_secant.value) - exp(-thickness.value * average_secant.value) * exp(-(thickness.value - x) * v_eigval));
                Dm.value = transmission.value / (average_secant.value - v_eigval) * (exp(-x*v_eigval) - exp(-x*average_secant.value));

                Dp.deriv.setZero();
                Dm.deriv.setZero();

                // Calculate the derivatives
                // eigval, thickness, and x have layer derivatives while average_secant has full derivatives
                for(uint n = 0; n < numderiv; ++n) {
                    Dp.deriv(n + layerStart) += thickness.deriv(n) * transmission.value * exp(-thickness.value * average_secant.value) * exp(-(thickness.value - x) * v_eigval);
                    Dp.deriv(n + layerStart) += thickness.deriv(n) * layer_fraction * transmission.value / (v_eigval + average_secant.value) * (-average_secant.value * exp(-x*average_secant.value) - v_eigval * exp(-thickness.value * average_secant.value) * exp(-(thickness.value - x) * v_eigval));
                    Dp.deriv(n + layerStart) -= d_eigval(n) * Dp.value / (v_eigval + average_secant.value);
                    Dp.deriv(n + layerStart) += d_eigval(n) * transmission.value / (v_eigval + average_secant.value) * exp(-thickness.value * average_secant.value) * exp(-(thickness.value - x) * v_eigval) * (thickness.value - x);

                    Dm.deriv(n + layerStart) += layer_fraction * thickness.deriv(n) * transmission.value / (average_secant.value - v_eigval) * (-1 * v_eigval * exp(-x*v_eigval) + average_secant.value * exp(-x*average_secant.value));
                    Dm.deriv(n + layerStart) += d_eigval(n) * Dm.value / (average_secant.value - v_eigval);
                    Dm.deriv(n + layerStart) += d_eigval(n) * transmission.value / (average_secant.value - v_eigval) * (-x) * exp(-x*v_eigval);
                }

                for(uint n = 0; n < numtotalderiv; ++n) {
                    Dp.deriv(n) += average_secant.deriv(n) / (v_eigval + average_secant.value) * transmission.value * (thickness.value * exp(-thickness.value * average_secant.value) * exp(-(thickness.value - x) * v_eigval) -x * exp(-x*average_secant.value));
                    Dp.deriv(n) -= average_secant.deriv(n) * Dp.value / (v_eigval + average_secant.value);
                    Dp.deriv(n) += transmission.deriv(n) * Dp.value / transmission.value;

                    Dm.deriv(n) -= average_secant.deriv(n) * Dm.value / (average_secant.value - v_eigval);
                    Dm.deriv(n) += average_secant.deriv(n) * transmission.value / (average_secant.value - v_eigval) * x * exp(-x*average_secant.value);
                    Dm.deriv(n) += transmission.deriv(n) * Dm.value / transmission.value;
                }

                // Add the terms to the particular source
                layer_homog_source.value += solution.value.dual_green_A_plus().value(j) * Y_plus.value * Dm.value;
                layer_homog_source.value += solution.value.dual_green_A_minus().value(j) * Y_minus.value * Dp.value;


                for(uint n = 0; n < numderiv; ++n) {
                    layer_homog_source.deriv(n + layerStart) += solution.value.dual_green_A_plus().value(j) * Y_plus.deriv(n) * Dm.value;
                    layer_homog_source.deriv(n + layerStart) += solution.value.dual_green_A_minus().value(j) * Y_minus.deriv(n) * Dp.value;

                    layer_homog_source.deriv(n + layerStart) +=
                            solution.value.dual_green_A_plus().deriv(n, j) * Y_plus.value * Dm.value;

                    layer_homog_source.deriv(n + layerStart) +=
                            solution.value.dual_green_A_minus().deriv(n, j) * Y_minus.value * Dp.value;
                }

                for(uint n = 0; n < numtotalderiv; ++n) {
                    layer_homog_source.deriv(n) +=
                            solution.value.dual_green_A_plus().value(j) * Y_plus.value * Dm.deriv(n);

                    layer_homog_source.deriv(n) +=
                            solution.value.dual_green_A_minus().value(j) * Y_minus.value * Dp.deriv(n);
                }
            }

		}
		homog_source.value += w*layer_homog_source.value;
		partic_source.value += w*layer_partic_source.value;
		for (uint l = 0; l < layers.inputDerivatives().numDerivative(); ++l)
		{
			const auto& qty = layers.inputDerivatives().layerDerivatives()[l];
			for (uint k = 0; k < qty.group_and_triangle_fraction.size(); ++k) {
				homog_source.deriv[qty.group_and_triangle_fraction[k].first] += w*layer_homog_source.deriv[l] * qty.group_and_triangle_fraction[k].second * qty.extinctions[k];
				partic_source.deriv[qty.group_and_triangle_fraction[k].first] += w*layer_partic_source.deriv[l] * qty.group_and_triangle_fraction[k].second * qty.extinctions[k];

			}
		}
	}
}

INSTANTIATE_TEMPLATE(sktran_do_detail::ParticipatingSourceTerm);
INSTANTIATE_TEMPLATE(sktran_do_detail::SphericalPostProcessing);
