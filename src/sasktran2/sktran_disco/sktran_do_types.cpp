#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_types.h"

template struct sasktran_disco::Dual<double>;
template struct sasktran_disco::Dual<std::complex<double>>;

template struct sasktran_disco::LayerDual<double>;
template struct sasktran_disco::LayerDual<std::complex<double>>;

template struct sasktran_disco::VectorDual<double>;
template struct sasktran_disco::VectorDual<std::complex<double>>;

template struct sasktran_disco::VectorLayerDual<double>;
template struct sasktran_disco::VectorLayerDual<std::complex<double>>;


template <int NSTOKES, int CNSTR>
void sasktran_disco::LPTripleProduct<NSTOKES, CNSTR>::calculate(const std::vector<LegendreCoefficient<NSTOKES>>& lephase, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2) {
    assert(lephase.size() == lp1.size() && lp1.size() == lp2.size());

    // Specialized version for NSTOKES == 1
    if constexpr (NSTOKES == 1) {
        m_aux.first.value = 0.0;
        m_aux.first.d_by_legendre_coeff.setZero();

        m_aux.second.value = 0.0;
        m_aux.second.d_by_legendre_coeff.setZero();

        for (int l = m_association_order; l < (int)m_nstr; ++l) {
            const auto& coeff = lephase[l];
            const auto& lp1l = lp1[l];
            const auto& lp2l = lp2[l];

            int negation_factor = 1;
            if ((l - m_association_order) % 2 != 0) {
                negation_factor *= -1;
            }

            double lp_factor = lp1l.P() * lp2l.P();
            double full_factor = coeff.a1 * lp_factor;

            m_aux.first.value += full_factor;
            m_aux.first.d_by_legendre_coeff[l] += lp_factor;

            m_aux.second.value += full_factor * negation_factor;
            m_aux.second.d_by_legendre_coeff[l] += lp_factor * negation_factor;
        }
    }
    else {
        // Fall back to direct calculation
        m_aux.first.calculate(lephase, lp1, lp2, false, m_association_order);
        m_aux.second.calculate(lephase, lp1, lp2, true, m_association_order);
    }
}


template<int NSTOKES, int CNSTR>
void sasktran_disco::LPTripleProduct<NSTOKES, CNSTR>::negations_derivative_emplace(uint num, sasktran_disco::TripleProductDerivativeHolder<NSTOKES>& holder) {
    // Implementation for NSTOKES == 1

    if constexpr (NSTOKES == 1) {
        uint nstr = (uint)m_aux.first.d_by_legendre_coeff.size();

        if (num % 2 == 0) {
            holder.value = m_aux.first.value;
            holder.d_by_legendre_coeff = m_aux.first.d_by_legendre_coeff;
        }
        else {
            holder.value = m_aux.second.value;
            holder.d_by_legendre_coeff = m_aux.second.d_by_legendre_coeff;
        }
    }
    else if constexpr (NSTOKES == 3) {
        if (num % 2 == 0) {
            holder.value = m_aux.first.value;
            holder.a1deriv = m_aux.first.a1deriv;
            holder.a2deriv = m_aux.first.a2deriv;
            holder.a3deriv = m_aux.first.a3deriv;
            holder.b1deriv = m_aux.first.b1deriv;
        }
        else {
            holder.value = m_aux.second.value;
            holder.a1deriv = m_aux.second.a1deriv;
            holder.a2deriv = m_aux.second.a2deriv;
            holder.a3deriv = m_aux.second.a3deriv;
            holder.b1deriv = m_aux.second.b1deriv;
        }
    }
    else {
        BOOST_LOG_TRIVIAL(error) << "LPTripleProduct<NSTOKES, CNSTR>::negations_derivative_emplace has no implementation for NSTOKES == " << NSTOKES;
    }
    
}


template <int NSTOKES, int CNSTR>
void sasktran_disco::LPTripleProduct<NSTOKES, CNSTR>::calculate_and_emplace(
        sasktran_disco::AEOrder m,
        const std::vector<LegendreCoefficient<NSTOKES>>& lephase,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1,
        const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2,
        sasktran_disco::TripleProductDerivativeHolder<NSTOKES>& holder_0_negation,
        sasktran_disco::TripleProductDerivativeHolder<NSTOKES>& holder_1_negation,
        double ssa) {
    holder_1_negation.ssa = ssa;
    holder_0_negation.ssa = ssa;

    // Implementation for NSTOKES == 1, direct calculation and emplacement
    if constexpr (NSTOKES == 1) {
        holder_0_negation.value = 0.0;
        holder_1_negation.value = 0.0;
        holder_0_negation.d_by_legendre_coeff.setZero();
        holder_1_negation.d_by_legendre_coeff.setZero();

        // Negation = 1
        for (int l = m; l < (int)m_nstr; l += 2) {
            double lp_factor = 0.5 * ssa * lp1[l].P() * lp2[l].P();
            double full_factor = lephase[l].a1 * lp_factor;

            holder_0_negation.value += full_factor;
            holder_0_negation.d_by_legendre_coeff[l] = lp_factor;

            holder_1_negation.value += full_factor;
            holder_1_negation.d_by_legendre_coeff[l] = lp_factor;
        }

        // Negation = -1
        for (int l = m + 1; l < (int)m_nstr; l += 2) {
            double lp_factor = 0.5 * ssa * lp1[l].P() * lp2[l].P();
            double full_factor = lephase[l].a1 * lp_factor;

            holder_0_negation.value += full_factor;
            holder_0_negation.d_by_legendre_coeff[l] = lp_factor;

            holder_1_negation.value -= full_factor;
            holder_1_negation.d_by_legendre_coeff[l] -= lp_factor;
        }
    } else if constexpr (NSTOKES == 3) {
        // Implementation for NSTOKES == 3, standard calculation 

        this->calculate(m, lephase, lp1, lp2);
        this->negations_derivative_emplace(0, holder_0_negation);
        this->negations_derivative_emplace(1, holder_1_negation);

        holder_0_negation.value *= 0.5 * ssa;
        holder_1_negation.value *= 0.5 * ssa;

        holder_0_negation.a1deriv *= 0.5 * ssa;
        holder_0_negation.a2deriv *= 0.5 * ssa;
        holder_0_negation.a3deriv *= 0.5 * ssa;
        holder_0_negation.b1deriv *= 0.5 * ssa;

        holder_1_negation.a1deriv *= 0.5 * ssa;
        holder_1_negation.a2deriv *= 0.5 * ssa;
        holder_1_negation.a3deriv *= 0.5 * ssa;
        holder_1_negation.b1deriv *= 0.5 * ssa;
    } else {
        BOOST_LOG_TRIVIAL(error) << "LPTripleProduct<NSTOKES, CNSTR>::calculate_and_emplace has no implementation for NSTOKES == " << NSTOKES;
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::LPTripleProduct)

template<int NSTOKES, int CNSTR>
void sasktran_disco::LayerInputDerivative<NSTOKES, CNSTR>::setZero() {
    d_SSA = 0.0;            // Change in layer SSA
    d_optical_depth = 0.0;  // Change in layer optical depth
    d_albedo = 0.0;         // Change in layer albedo ( only for surface)
    for (uint i = 0; i < d_legendre_coeff.size(); ++i) {
        d_legendre_coeff[i].a1 = 0.0;  // Change in layer legendre coeff

        if constexpr (NSTOKES == 3) {
            d_legendre_coeff[i].a2 = 0.0;
            d_legendre_coeff[i].a3 = 0.0;
            d_legendre_coeff[i].b1 = 0.0;
        }
    }
}

SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::LayerInputDerivative);
