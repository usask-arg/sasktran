#include "modules/sktran_do_deprecated/include/sktran_do.h"
#include "modules/sktran_do_deprecated/include/sktran_do_types.h"

template class sktran_do_detail::Dual<double>;
template class sktran_do_detail::Dual<std::complex<double>>;

template class sktran_do_detail::LayerDual<double>;
template class sktran_do_detail::LayerDual<std::complex<double>>;

template class sktran_do_detail::VectorDual<double>;
template class sktran_do_detail::VectorDual<std::complex<double>>;

template class sktran_do_detail::VectorLayerDual<double>;
template class sktran_do_detail::VectorLayerDual<std::complex<double>>;

template <int NSTOKES, int CNSTR>
void sktran_do_detail::LPTripleProduct<NSTOKES, CNSTR>::calculate(const std::vector<LegendreCoefficient<NSTOKES>>& lephase, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp1, const std::vector<LegendrePhaseContainer<NSTOKES>>& lp2) {
    assert(lephase.size() == lp1.size() && lp1.size() == lp2.size());

    m_aux.first.calculate(lephase, lp1, lp2, false, m_association_order);
    m_aux.second.calculate(lephase, lp1, lp2, true, m_association_order);
}

template<int NSTOKES, int CNSTR>
void sktran_do_detaill::LPTripleProduct<NSTOKES, CNSTR>::negations_derivative_emplace(uint num, sktran_do_detail::TripleProductDerivativeHolder<NSTOKES, CNSTR>& holder) {
    if constexpr(NSTOKES == 1) {
        uint nstr = (uint)m_aux.first.d_by_legendre_coeff.size();

        if (num % 2 == 0) {
            holder.value = m_aux.first.value;
            holder.d_by_legendre_coeff = m_aux.first.d_by_legendre_coeff;
        }
        else {
            holder.value = m_aux.second.value;
            holder.d_by_legendre_coeff =  m_aux.second.d_by_legendre_coeff;
        }
    }

    if constexpr(NSTOKES == 3) {
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

    if constexpr(NSTOKES == 4) {
        if (num % 2 == 0) {
            holder.value = m_aux.first.value;
            holder.a1deriv = m_aux.first.a1deriv;
            holder.a2deriv = m_aux.first.a2deriv;
            holder.a3deriv = m_aux.first.a3deriv;
            holder.a4deriv = m_aux.first.a4deriv;
            holder.b1deriv = m_aux.first.b1deriv;
            holder.b2deriv = m_aux.first.b2deriv;
        }
        else {
            holder.value = m_aux.second.value;
            holder.a1deriv = m_aux.second.a1deriv;
            holder.a2deriv = m_aux.second.a2deriv;
            holder.a3deriv = m_aux.second.a3deriv;
            holder.a4deriv = m_aux.second.a4deriv;
            holder.b1deriv = m_aux.second.b1deriv;
            holder.b2deriv = m_aux.second.b2deriv;
        }
    }
}


template class sktran_do_detail::LPTripleProduct<1>;
template class sktran_do_detail::LPTripleProduct<3>;
template class sktran_do_detail::LPTripleProduct<4>;


template <>
void sktran_do_detail::LayerInputDerivative<1>::setZero() {
    d_SSA = 0.0;            // Change in layer SSA
    d_optical_depth = 0.0;  // Change in layer optical depth
    d_albedo = 0.0;         // Change in layer albedo ( only for surface)
    for (uint i = 0; i < d_legendre_coeff.size(); ++i) {
        d_legendre_coeff[i].a1 = 0.0;  // Change in layer legendre coeff
    }
    handle = SKCLIMATOLOGY_UNDEFINED;
}

template <>
void sktran_do_detail::LayerInputDerivative<3>::setZero() {
    d_SSA = 0.0;            // Change in layer SSA
    d_optical_depth = 0.0;  // Change in layer optical depth
    d_albedo = 0.0;         // Change in layer albedo ( only for surface)
    for (uint i = 0; i < d_legendre_coeff.size(); ++i) {
        d_legendre_coeff[i].a1 = 0.0;
        d_legendre_coeff[i].a2 = 0.0;
        d_legendre_coeff[i].a3 = 0.0;
        d_legendre_coeff[i].b1 = 0.0;
    }
    handle = SKCLIMATOLOGY_UNDEFINED;

}

template <>
void sktran_do_detail::LayerInputDerivative<4>::setZero() {
    d_SSA = 0.0;            // Change in layer SSA
    d_optical_depth = 0.0;  // Change in layer optical depth
    d_albedo = 0.0;         // Change in layer albedo ( only for surface)
    for (uint i = 0; i < d_legendre_coeff.size(); ++i) {
        d_legendre_coeff[i].a1 = 0.0;
        d_legendre_coeff[i].a2 = 0.0;
        d_legendre_coeff[i].a3 = 0.0;
        d_legendre_coeff[i].a4 = 0.0;
        d_legendre_coeff[i].b1 = 0.0;
        d_legendre_coeff[i].b2 = 0.0;
    }
    handle = SKCLIMATOLOGY_UNDEFINED;

}

template class sktran_do_detail::LayerInputDerivative<1>;
template class sktran_do_detail::LayerInputDerivative<3>;
template class sktran_do_detail::LayerInputDerivative<4>;