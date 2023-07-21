#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_polarization_types.h"

void sasktran_disco::TripleProductDerivativeHolder<1>::reduce(const LayerInputDerivative<1>& layer_deriv, double& deriv) const {
    // TripleProduct derivative with respect to thickness is always 0

    // TripleProduct derivative with respect to SSA is always value / ssa
    // Have a special case if ssa = 0 then value is also 0 and the derivative is 0

    // TODO: Is that actually correct? The derivative probably isn't 0 if the SSA is 0
    deriv = ssa != 0.0 ? layer_deriv.d_SSA * value / ssa : 0.0;
    for( int l = 0; l < nstr; ++l) {
        deriv += layer_deriv.d_legendre_coeff[l].a1 * d_by_legendre_coeff[l];
    }
}

void sasktran_disco::TripleProductDerivativeHolder<3>::reduce(const LayerInputDerivative<3>& layer_deriv, Eigen::Matrix<double, 3, 3>& deriv) const {
    // TripleProduct derivative with respect to thickness is always 0

    // TripleProduct derivative with respect to SSA is always  value / ssa
    if (ssa != 0.0) {
        deriv = value / ssa * layer_deriv.d_SSA;
    }
    else {
        deriv.setZero();
    }
    for (int l = 0; l < nstr; ++l) {
        const auto& dl = layer_deriv.d_legendre_coeff[l];
        // Unravel our internal sparse derivative storage
        deriv(0, 0) += dl.a1 * a1deriv(l);
        deriv(0, 1) += dl.b1 * b1deriv(l, 0);
        deriv(0, 2) += dl.b1 * b1deriv(l, 1);

        deriv(1, 0) += dl.b1 * b1deriv(l, 2);
        deriv(1, 1) += dl.a2 * a2deriv(l, 0) * dl.a3 * a3deriv(l, 0);
        deriv(1, 2) += dl.a2 * a2deriv(l, 1) * dl.a3 * a3deriv(l, 1);

        deriv(2, 0) += dl.b1 * b1deriv(l, 3);
        deriv(2, 1) += dl.a2 * a2deriv(l, 2) + dl.a3 * a3deriv(l, 2);
        deriv(2, 2) += dl.a2 * a2deriv(l, 3) + dl.a3 * a3deriv(l, 3);
    }
}

void sasktran_disco::TripleProductDerivativeHolder<4>::reduce(const LayerInputDerivative<4>& layer_deriv, Eigen::Matrix<double, 4, 4>& deriv) const {
    // TripleProduct derivative with respect to thickness is always 0

    // TripleProduct derivative with respect to SSA is always  value / ssa
    if( ssa != 0.0 ) {
        deriv = value / ssa * layer_deriv.d_SSA;
    } else {
        deriv.setZero();
    }
    for( int l = 0; l < nstr; ++l) {
        const auto& dl = layer_deriv.d_legendre_coeff[l];
        // Unravel our internal sparse derivative storage
        deriv(0, 0) += dl.a1 * a1deriv(l);
        deriv(0, 1) += dl.b1 * b1deriv(l, 0);
        deriv(0, 2) += dl.b1 * b1deriv(l, 1);

        deriv(1, 0) += dl.b1 * b1deriv(l, 2);
        deriv(1, 1) += dl.a2 * a2deriv(l, 0) * dl.a3 * a3deriv(l, 0);
        deriv(1, 2) += dl.a2 * a2deriv(l, 1) * dl.a3 * a3deriv(l, 1);
        deriv(1, 3) += dl.b2 * b2deriv(l, 0);

        deriv(2, 0) += dl.b1 * b1deriv(l, 3);
        deriv(2, 1) += dl.a2 * a2deriv(l, 2) + dl.a3 * a3deriv(l, 2);
        deriv(2, 2) += dl.a2 * a2deriv(l, 3) + dl.a3 * a3deriv(l, 3);
        deriv(2, 3) += dl.b2 * b2deriv(l, 1);

        deriv(3, 1) += dl.b2 * b2deriv(l, 2);
        deriv(3, 2) += dl.b2 * b2deriv(l, 3);
        deriv(3, 3) += dl.a4 * a4deriv(l);
    }
}

template class sasktran_disco::TripleProductDerivativeHolder<1>;
template class sasktran_disco::TripleProductDerivativeHolder<3>;
template class sasktran_disco::TripleProductDerivativeHolder<4>;

void sasktran_disco::InhomogeneousSourceHolder<1>::reduce(const LayerInputDerivative<1>& layer_deriv, double& deriv) const {
    // TripleProduct derivative with respect to thickness is always 0

    // TripleProduct derivative with respect to SSA is always value / ssa
    // Have a special case if ssa = 0 then value is also 0 and the derivative is 0
    deriv = d_by_ssa * layer_deriv.d_SSA;
    for (int l = 0; l < nstr; ++l) {
        deriv += layer_deriv.d_legendre_coeff[l].a1 * d_by_legendre_coeff[l];
    }
}


template <int NSTOKES, int CNSTR>
void sasktran_disco::InhomogeneousSourceHolder<NSTOKES, CNSTR>::reduce(const LayerInputDerivative<NSTOKES>& layer_deriv, Eigen::Vector<double, NSTOKES>& deriv) const {
    deriv = d_by_ssa * layer_deriv.d_SSA;
    for (int l = 0; l < nstr; ++l) {
        deriv(0) += layer_deriv.d_legendre_coeff[l].a1 * d_by_a1(l);
        deriv(1) += layer_deriv.d_legendre_coeff[l].b1 * d_by_b1_first(l);
        // TODO: this breaks for NSTOKES=2 right?
        deriv(2) += layer_deriv.d_legendre_coeff[l].b1 * d_by_b1_second(l);
    }
}

template class sasktran_disco::InhomogeneousSourceHolder<1>;
template class sasktran_disco::InhomogeneousSourceHolder<3>;
template class sasktran_disco::InhomogeneousSourceHolder<4>;


void sasktran_disco::Radiance<1>::apply_transmission_factor(const Dual<double>& transmission) {
    deriv = deriv * transmission.value + value * transmission.deriv;
    value = value * transmission.value;
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::Radiance<NSTOKES, CNSTR>::apply_transmission_factor(const Dual<double>& transmission) {
    deriv *= transmission.value;
    
    if (transmission.deriv.size() > 0) {
        for (int s = 0; s < NSTOKES; ++s) {
            for (int k = 0; k < transmission.deriv.size(); ++k) {
                deriv(k, s) += value(s) * transmission.deriv(k);
            }
        }
    }
    value = value * transmission.value;
}

bool sasktran_disco::Radiance<1>::converged(double I, double epsilon) {
    return abs(this->I() / I) < epsilon;
}

template <int NSTOKES, int CNSTR>
bool sasktran_disco::Radiance<NSTOKES, CNSTR>::converged(double I, double epsilon) {
    bool converged = true;

    for(uint i = 0; i < NSTOKES; ++i) {
        converged = converged && abs(value(i) / I) < epsilon;
    }

    return converged;
}


void sasktran_disco::Radiance<1>::apply_azimuth_expansion(double angle, int m) {
    double cos_daz = cos(m*angle);

    value *= cos_daz;
    deriv *= cos_daz;
}

template <>
void sasktran_disco::Radiance<3>::apply_azimuth_expansion(double angle, int m) {
    double cos_daz = cos(m*angle);
    double sin_daz = sin(m*angle);

    value(0) *= cos_daz;
    value(1) *= cos_daz;
    value(2) *= sin_daz;

    deriv(Eigen::all, 0) *= cos_daz;
    deriv(Eigen::all, 1) *= cos_daz;
    deriv(Eigen::all, 2) *= sin_daz;
}

template <>
void sasktran_disco::Radiance<4>::apply_azimuth_expansion(double angle, int m) {
    double cos_daz = cos(m*angle);
    double sin_daz = sin(m*angle);

    value(0) *= cos_daz;
    value(1) *= cos_daz;
    value(2) *= sin_daz;
    value(3) *= sin_daz;

    deriv(Eigen::all, 0) *= cos_daz;
    deriv(Eigen::all, 1) *= cos_daz;
    deriv(Eigen::all, 2) *= sin_daz;
    deriv(Eigen::all, 3) *= sin_daz;
}

template struct sasktran_disco::Radiance<1>;
template struct sasktran_disco::Radiance<3>;
template struct sasktran_disco::Radiance<4>;
