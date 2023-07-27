#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_types.h"

template <int NSTOKES, int CNSTR>
sasktran_disco::LayerFundamentalDerivative<NSTOKES, CNSTR>::LayerFundamentalDerivative(uint nstr)
{
	d_by_legendre_coeff.resize(nstr);
	d_by_legendre_coeff.data().setZero();
	d_by_SSA = 0.0;
	d_by_opticalDepth = 0.0;
	value = 0.0;
}

template <>
void sasktran_disco::LayerFundamentalDerivative<1>::reduce(const LayerInputDerivative<1>& layer_deriv, double& output) const
{
	output = 0.0;

	output += layer_deriv.d_SSA * d_by_SSA;
	output += layer_deriv.d_optical_depth * d_by_opticalDepth;

	for (uint i = 0; i < d_by_legendre_coeff.data().size(); ++i)
	{
		output += layer_deriv.d_legendre_coeff[i].a1 * d_by_legendre_coeff.data()(i);
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::LayerFundamentalDerivative<NSTOKES, CNSTR>::reduce(const LayerInputDerivative<NSTOKES>& layer_deriv, double& output) const
{
    // TODO: To be implemented
    // probably will just delete this class anyway
}


template class sasktran_disco::LayerFundamentalDerivative<1>;
template class sasktran_disco::LayerFundamentalDerivative<3>;
template class sasktran_disco::LayerFundamentalDerivative<4>;