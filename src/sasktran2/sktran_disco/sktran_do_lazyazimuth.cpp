#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_lazyazimuth.h"

template <>
void sasktran_disco::LegendrePolynomials<1>::calculateAEOrder(AEOrder m, std::vector<sasktran_disco::LegendrePhaseContainer<1>>& lepolys)
{
    auto calculator = sasktran_disco::WignerDCalculator(m, 0);

    double theta = acos(m_value);

	// Evaluate the associated Legendre polynomial ( order m, degree l=[0:M_NSTR) ) at m_value.
	lepolys.resize(this->M_NSTR);
	for(LPOrder l = 0; l < this->M_NSTR; ++l) {
        lepolys[l].value = calculator.d(theta, l);
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::LegendrePolynomials<NSTOKES, CNSTR>::calculateAEOrder(AEOrder m, std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>>& lepolys)
{
    // TODO: breaks for NSTOKES=2?
    auto calculatorP = sasktran_disco::WignerDCalculator(m, 0);
    auto calculatorneg = sasktran_disco::WignerDCalculator(m, -2);
    auto calculatorpos = sasktran_disco::WignerDCalculator(m, 2);

    double theta = acos(m_value);

    // Evaluate the associated Legendre polynomial ( order m, degree l=[0:M_NSTR) ) at m_value.
    lepolys.resize(this->M_NSTR);
    for(LPOrder l = 0; l < this->M_NSTR; ++l) {
        lepolys[l].P() = calculatorP.d(theta, l);
        lepolys[l].R() = -0.5 * (calculatorpos.d(theta, l) + calculatorneg.d(theta, l));
        lepolys[l].T() = -0.5 * (calculatorpos.d(theta, l) - calculatorneg.d(theta, l));
    }
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::LegendreSumMatrix<NSTOKES, CNSTR>::calculateAEOrder(AEOrder m, LegendreSumMatrixStorage<NSTOKES>& sum_matrix)
{

	// Compute scattering matrix from streams to streams.
	sum_matrix.M_SSA = M_SSA;
	sum_matrix.resize(this->M_NSTR);

	const auto& le_phasef = *M_LPE_PHASEF;

	LPTripleProduct<NSTOKES> triple_product(m, (uint)le_phasef.size());
    sasktran_disco::TripleProductDerivativeHolder<NSTOKES> holder(this->M_NSTR);

	const uint N = this->M_NSTR / 2;
	for(StreamIndex i = 0; i < N; ++i) {
		const auto& lp_out = M_LP_MU[m][i];
		for(StreamIndex j = 0; j <= i; ++j) {
			const uint linear_index_0 = sum_matrix.linear_index(i, j);
			const uint linear_index_1 = sum_matrix.linear_index(i, j + N);

			const auto& lp_in = M_LP_MU[m][j];
			triple_product.calculate(le_phasef, lp_out, lp_in);

            triple_product.negations_derivative_emplace(0, holder);
			assign(linear_index_0, holder, sum_matrix);

            triple_product.negations_derivative_emplace(1, holder);
            assign(linear_index_1, holder, sum_matrix);
		}
	}
}

template <>
void sasktran_disco::LegendreSumMatrix<1>::calculateAEOrder(AEOrder m, LegendreSumMatrixStorage<1>& sum_matrix)
{
    // Compute scattering matrix from streams to streams.
    sum_matrix.M_SSA = M_SSA;
    sum_matrix.resize(this->M_NSTR);

    const auto& le_phasef = *M_LPE_PHASEF;

    const uint N = this->M_NSTR / 2;
    for(StreamIndex i = 0; i < N; ++i) {
        const auto& lp_out = M_LP_MU[m][i];
        for(StreamIndex j = 0; j <= i; ++j) {
            const uint linear_index_0 = sum_matrix.linear_index(i, j);
            const uint linear_index_1 = sum_matrix.linear_index(i, j + N);

            const auto& lp_in = M_LP_MU[m][j];

            auto& first = sum_matrix.storage[linear_index_0];
            auto& second = sum_matrix.storage[linear_index_1];

            sum_matrix.triple_product->calculate_and_emplace(m,
                                                             le_phasef,
                                                             lp_out,
                                                             lp_in,
                                                             first,
                                                             second,
                                                             M_SSA);
        }
    }
}


template <>
void sasktran_disco::LegendreSumMatrix<3>::calculateAEOrder(AEOrder m, LegendreSumMatrixStorage<3>& sum_matrix)
{

    // Compute scattering matrix from streams to streams.
    sum_matrix.M_SSA = M_SSA;
    sum_matrix.resize(this->M_NSTR);

    const auto& le_phasef = *M_LPE_PHASEF;

    const uint N = this->M_NSTR / 2;
    for(StreamIndex i = 0; i < N; ++i) {
        const auto& lp_out = M_LP_MU[m][i];
        for(StreamIndex j = 0; j <= i; ++j) {
            const uint linear_index_0 = sum_matrix.linear_index(i, j);
            const uint linear_index_1 = sum_matrix.linear_index(i, j + N);

            const auto& lp_in = M_LP_MU[m][j];
            sum_matrix.triple_product->calculate(m, le_phasef, lp_out, lp_in);

            sum_matrix.triple_product->negations_derivative_emplace(0, sum_matrix.storage[linear_index_0]);

            sum_matrix.storage[linear_index_0].value *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_0].a1deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_0].a2deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_0].a3deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_0].b1deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_0].ssa = M_SSA;

            sum_matrix.triple_product->negations_derivative_emplace(1, sum_matrix.storage[linear_index_1]);

            sum_matrix.storage[linear_index_1].value *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_1].a1deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_1].a2deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_1].a3deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_1].b1deriv *= 0.5 * M_SSA;
            sum_matrix.storage[linear_index_1].ssa = M_SSA;
        }
    }
}

template <>
void sasktran_disco::LegendreSumMatrix<1>::assign(int linear_index, const sasktran_disco::TripleProductDerivativeHolder<1>& val, LegendreSumMatrixStorage<1>& sum_matrix) {
    sum_matrix.storage[linear_index].value = val.value * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].d_by_legendre_coeff = val.d_by_legendre_coeff * 0.5 * M_SSA;

    sum_matrix.storage[linear_index].ssa = M_SSA;
}


template <>
void sasktran_disco::LegendreSumMatrix<3>::assign(int linear_index, const sasktran_disco::TripleProductDerivativeHolder<3>& val, LegendreSumMatrixStorage<3>& sum_matrix) {
    sum_matrix.storage[linear_index].value = val.value * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a1deriv = val.a1deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a2deriv = val.a2deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a3deriv = val.a3deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].b1deriv = val.b1deriv * 0.5 * M_SSA;

    sum_matrix.storage[linear_index].ssa = M_SSA;

}

template <>
void sasktran_disco::LegendreSumMatrix<4>::assign(int linear_index, const sasktran_disco::TripleProductDerivativeHolder<4>& val, LegendreSumMatrixStorage<4>& sum_matrix) {
    sum_matrix.storage[linear_index].value = val.value * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a1deriv = val.a1deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a2deriv = val.a2deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a3deriv = val.a3deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].a4deriv = val.a4deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].b1deriv = val.b1deriv * 0.5 * M_SSA;
    sum_matrix.storage[linear_index].b2deriv = val.b2deriv * 0.5 * M_SSA;

    sum_matrix.storage[linear_index].ssa = M_SSA;

}

void sasktran_disco::AlbedoExpansion::calculateAEOrder(AEOrder m, sasktran_disco::Albedo& albedo)
{
	// Configure the BRDF object for the current order of the azimuth expansion.
	albedo.configure(m, M_LOS, M_MU, M_CSZ, m_brdf.get(), m_nterms);
}

template class sasktran_disco::LegendrePolynomials<1>;
template class sasktran_disco::LegendrePolynomials<3>;

template class sasktran_disco::LegendreSumMatrix<1>;
template class sasktran_disco::LegendreSumMatrix<3>;
