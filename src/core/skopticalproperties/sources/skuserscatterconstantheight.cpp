#include <skopticalproperties21.h>
#include <boost/math/special_functions/legendre.hpp>

bool skOpticalProperties_UserDefinedScatterConstantHeight::CalculateCrossSections(double wavenumber, double* absxs, double* extxs, double* scattxs) {
	std::array<size_t, 2> indicies;
	std::array<double, 2> weights;

	InterpolationWeights(1e7 / wavenumber, indicies, weights);

	*absxs = 0.0;
	*scattxs = 0.0;

	for (int j = 0; j < 2; ++j) {
		*absxs += weights[j] * m_xs_abs.at(indicies[j]);
		*scattxs += weights[j] * m_xs_scat.at(indicies[j]);
	}
	*extxs = *absxs + *scattxs;

	return true;
}
bool skOpticalProperties_UserDefinedScatterConstantHeight::CalculatePhaseMatrix(double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix) {
	phasematrix->SetTo(0.0);

    std::pair<double, size_t> cosscatterandindex;
    cosscatterandindex.first = cosscatterangle;

    double temp;
    CalculateP11(wavenumber, cosscatterandindex, temp);

    phasematrix->At(1, 1) = temp;

	return true;
}
bool skOpticalProperties_UserDefinedScatterConstantHeight::CalculateP11(double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11) {
	
	p11 = 0.0;

	std::array<size_t, 2> indicies;
	std::array<double, 2> weights;

	InterpolationWeights(1e7 / wavenumber, indicies, weights);

	for (int idx = 0; idx < (int)m_a1.YSize(); idx++) {
		double legendre = boost::math::legendre_p(idx, cosscatterandindex.first);

		p11 += m_a1.at(indicies[0], idx) * weights[0] * legendre;
		p11 += m_a1.at(indicies[1], idx) * weights[1] * legendre;
	}

	return true;
}
bool skOpticalProperties_UserDefinedScatterConstantHeight::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) {
	std::array<size_t, 2> indicies;
	std::array<double, 2> weights;

	InterpolationWeights(1e7 / wavenumber, indicies, weights);


	for (int i = 0; i < usermaxcoeff; ++i) {
		coeff[i] = 0.0;
		if (i < m_a1.YSize()) {
			for (int j = 0; j < 2; ++j) {
				coeff[i] += weights[j] * m_a1.at(indicies[j], i);
			}
		}
	}

	opticalmaxcoeff = std::min(usermaxcoeff, (int)m_a1.YSize());

	return true;
}

bool skOpticalProperties_UserDefinedScatterConstantHeight::LegendreCoefficientsPolarized(double wavenumber, double* a1, double* a2, double* a3, double* a4, double* b1, double* b2, int usermaxcoeff, int& opticalmaxcoeff) {
	std::array<size_t, 2> indicies;
	std::array<double, 2> weights;

	InterpolationWeights(1e7 / wavenumber, indicies, weights);

	for (int i = 0; i < usermaxcoeff; ++i) {
		a1[i] = 0.0;
		a2[i] = 0.0;
		a3[i] = 0.0;
		a4[i] = 0.0;
		b1[i] = 0.0;
		b2[i] = 0.0;
		if (i < m_a1.YSize()) {
			for (int j = 0; j < 2; ++j) {
				a1[i] += weights[j] * m_a1.at(indicies[j], i);
				a2[i] += weights[j] * m_a2.at(indicies[j], i);
				a3[i] += weights[j] * m_a3.at(indicies[j], i);
				a4[i] += weights[j] * m_a4.at(indicies[j], i);
				b1[i] += weights[j] * m_b1.at(indicies[j], i);
				b2[i] += weights[j] * m_b2.at(indicies[j], i);
			}
		}
	}

	opticalmaxcoeff = std::min(usermaxcoeff, (int)m_a1.YSize());

	return true;
}

void skOpticalProperties_UserDefinedScatterConstantHeight::InterpolationWeights(double wavelength, std::array<size_t, 2>& indicies, std::array<double, 2>& weights) {
	nxLinearInterpolate::FindBoundingIndicesAscending(m_wavelengths.begin(),
		m_wavelengths.end(),
		wavelength,
		&indicies[0],
		&indicies[1],
		&weights[0],
		&weights[1]);

	if (wavelength > weights[1])
	{
		weights[0] = 0.0;
		weights[1] = 1.0;
	}
	else if (wavelength < weights[0])
	{
		weights[0] = 1.0;
		weights[1] = 0.0;
	}
	else
	{
		weights[0] = (weights[1] - wavelength) / (weights[1] - weights[0]);
		weights[1] = 1.0 - weights[0];
	}
}
