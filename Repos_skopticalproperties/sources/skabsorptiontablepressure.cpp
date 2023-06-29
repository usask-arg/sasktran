#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <skopticalproperties21.h>

std::array<double, 2> skAbsorptionPressureTemperatureEntry::InterpolationWeights(const std::array<size_t, 2>& indicies,
																				 const std::array<double, 2>& values,
																				 double value)
{
	// Calculates the interpolation weights given indicies and values at those indicies
	// Values outside are truncated 

	std::array<double, 2> weights;

	if (value > values[1])
	{
		weights[0] = 0.0;
		weights[1] = 1.0;
	}
	else if (value < values[0])
	{
		weights[0] = 1.0;
		weights[1] = 0.0;
	}
	else
	{
		weights[0] = (values[1] - value) / (values[1] - values[0]);
		weights[1] = 1.0 - weights[0];
	}

	return weights;
}

double skAbsorptionPressureTemperatureEntry::Interpolate(const std::array<size_t, 2>& lowTemperatureIndex,
                                                         const std::array<double, 2>& lowTemperatureWeight,
                                                         const std::array<size_t, 2>& highTemperatureIndex,
                                                         const std::array<double, 2>& highTemperatureWeight,
                                                         const std::array<size_t, 2>& pressureIndex,
                                                         const std::array<double, 2>& pressureWeight,
                                                         double wavelength) {
    // Last the indicies and weights for the wavelengths
    std::array<size_t, 2> wavelIndex;
    std::array<double, 2> wavelWeight;
    nxLinearInterpolate::FindBoundingIndicesAscending(m_wavelength.begin(),
                                                      m_wavelength.end(),
                                                      wavelength,
                                                      &wavelIndex[0],
                                                      &wavelIndex[1],
                                                      &wavelWeight[0],
                                                      &wavelWeight[1]);

    wavelWeight = this->InterpolationWeights(wavelIndex, wavelWeight, wavelength);

    double xs = 0;
    for (size_t pidx = 0; pidx < 2; pidx++)
    {
        double pweight = pressureWeight[pidx];
        size_t pindex = pressureIndex[pidx];
        for (size_t tidx = 0; tidx < 2; tidx++)
        {
            double tweight;
            size_t tindex;
            if (pidx == 0)
            {
                tweight = lowTemperatureWeight[tidx];
                tindex = lowTemperatureIndex[tidx];
            }
            else
            {
                tweight = highTemperatureWeight[tidx];
                tindex = highTemperatureIndex[tidx];
            }
            for (size_t widx = 0; widx < 2; widx++)
            {
                double wweight = wavelWeight[widx];
                size_t windex = wavelIndex[widx];

                xs += m_xs.at(windex, tindex, pindex) * wweight * tweight * pweight;
            }
        }
    }
    return xs;
}

void skAbsorptionPressureTemperatureEntry::calcInterpolationWeights(double temperature,
                                                                    double pressure,
                                                                    std::array<size_t, 2>& lowTemperatureIndex,
                                                                    std::array<double, 2>& lowTemperatureWeight,
                                                                    std::array<size_t, 2>& highTemperatureIndex,
                                                                    std::array<double, 2>& highTemperatureWeight,
                                                                    std::array<size_t, 2>& pressureIndex,
                                                                    std::array<double, 2>& pressureWeight) {
    // Calculate the indicies and weights for pressure
    nxLinearInterpolate::FindBoundingIndicesAscending(m_pressure.begin(),
                                                      m_pressure.end(),
                                                      pressure,
                                                      &pressureIndex[0],
                                                      &pressureIndex[1],
                                                      &pressureWeight[0],
                                                      &pressureWeight[1]);

    pressureWeight = this->InterpolationWeights(pressureIndex, pressureWeight, pressure);

    const double* start = &m_temperature.at(0, pressureIndex[0]);
    nxLinearInterpolate::FindBoundingIndicesAscending(start,
                                                      start + m_temperature.XSize(),
                                                      temperature, &lowTemperatureIndex[0],
                                                      &lowTemperatureIndex[1],
                                                      &lowTemperatureWeight[0],
                                                      &lowTemperatureWeight[1]);

    lowTemperatureWeight = this->InterpolationWeights(lowTemperatureIndex, lowTemperatureWeight, temperature);

    start = &m_temperature.at(0, pressureIndex[1]);
    nxLinearInterpolate::FindBoundingIndicesAscending(start,
                                                      start + m_temperature.XSize(),
                                                      temperature, &highTemperatureIndex[0],
                                                      &highTemperatureIndex[1],
                                                      &highTemperatureWeight[0],
                                                      &highTemperatureWeight[1]);

    highTemperatureWeight = this->InterpolationWeights(highTemperatureIndex, highTemperatureWeight, temperature);

}


double skAbsorptionPressureTemperatureEntry::Interpolate(double temperature, double pressure, double wavelength)
{
	// Interpolates the entry linearly in temperature, pressure, and wavelength
	std::array<size_t, 2> pressureIndex;
	std::array<double, 2> pressureWeight;

	std::array<size_t, 2> lowTemperatureIndex;
	std::array<double, 2> lowTemperatureWeight;

	std::array<size_t, 2> highTemperatureIndex;
	std::array<double, 2> highTemperatureWeight;

    calcInterpolationWeights(temperature,
                             pressure,
                             lowTemperatureIndex,
                             lowTemperatureWeight,
                             highTemperatureIndex,
                             highTemperatureWeight,
                             pressureIndex,
                             pressureWeight);

    return Interpolate(lowTemperatureIndex, lowTemperatureWeight,
                       highTemperatureIndex, highTemperatureWeight,
                       pressureIndex, pressureWeight, wavelength);
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::AddEntry(const nx3dArray<double>& xs,
																 const nx2dArray<double>& temperature,
																 const nx1dArray<double>& pressure,
																 const nx1dArray<double>& wavelength,
																 double broadenervmr)
{
	bool ok = true;

	m_xsentry.emplace_back(skAbsorptionPressureTemperatureEntry(xs, temperature, pressure, wavelength));
	m_broadenervmr.push_back(broadenervmr);

	return ok;
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::CalculateCrossSections(double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool ok = true;
	double xs;

	if (m_xsentry.size() == 1)
	{
		// Don't need to worry about interpolating in broadener vmr
		xs = m_xsentry[0].Interpolate(m_temperature, m_pressure, 1e7 / wavenumber);
	}
	else
	{
		// Calculate the broadener interpolation weights
		std::array<size_t, 2> vmrIndex;
		std::array<double, 2> vmrWeight;
		ok = ok && nxLinearInterpolate::FindBoundingIndicesAscending(m_broadenervmr.begin(),
																	 m_broadenervmr.end(),
																	 m_currentbroadenervmr,
																	 &vmrIndex[0],
																	 &vmrIndex[1],
																	 &vmrWeight[0],
																	 &vmrWeight[1]);
		vmrWeight = skAbsorptionPressureTemperatureEntry::InterpolationWeights(vmrIndex, vmrWeight, m_currentbroadenervmr);

		xs = 0;
		for (size_t vmridx = 0; vmridx < 2; vmridx++)
		{
			xs += m_xsentry[vmrIndex[vmridx]].Interpolate(m_temperature, m_pressure, 1e7 / wavenumber) * vmrWeight[vmridx];
		}
	}


	*absxs = xs;
	*extxs = xs;
	*scattxs = 0.0;

	return ok;
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::CalculateCrossSectionsArray(const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs) {
    bool ok = true;

    std::array<size_t, 2> pressureIndex;
    std::array<double, 2> pressureWeight;

    std::array<size_t, 2> lowTemperatureIndex;
    std::array<double, 2> lowTemperatureWeight;

    std::array<size_t, 2> highTemperatureIndex;
    std::array<double, 2> highTemperatureWeight;

    m_xsentry[0].calcInterpolationWeights(m_temperature,
                                          m_pressure,
                                          lowTemperatureIndex,
                                          lowTemperatureWeight,
                                          highTemperatureIndex,
                                          highTemperatureWeight,
                                          pressureIndex,
                                          pressureWeight);

    if (m_xsentry.size() == 1)
    {
        // Don't need to worry about interpolating in broadener vmr
        #pragma omp parallel for schedule(guided, 1)
        for(int wavidx = 0; wavidx < numwavenumber; ++wavidx) {
            extxs[wavidx] = m_xsentry[0].Interpolate(lowTemperatureIndex,
                                                     lowTemperatureWeight,
                                                     highTemperatureIndex,
                                                     highTemperatureWeight,
                                                     pressureIndex,
                                                     pressureWeight,
                                                     1e7 / wavenumber[wavidx]);
            absxs[wavidx] = extxs[wavidx];
            scattxs[wavidx] = 0.0;
        }
    }
    else
    {
        // Calculate the broadener interpolation weights
        std::array<size_t, 2> vmrIndex;
        std::array<double, 2> vmrWeight;
        ok = ok && nxLinearInterpolate::FindBoundingIndicesAscending(m_broadenervmr.begin(),
                                                                     m_broadenervmr.end(),
                                                                     m_currentbroadenervmr,
                                                                     &vmrIndex[0],
                                                                     &vmrIndex[1],
                                                                     &vmrWeight[0],
                                                                     &vmrWeight[1]);
        vmrWeight = skAbsorptionPressureTemperatureEntry::InterpolationWeights(vmrIndex, vmrWeight, m_currentbroadenervmr);


        #pragma omp parallel for schedule(guided, 1)
        for(int wavidx = 0; wavidx < numwavenumber; ++wavidx) {
            extxs[wavidx] = 0.0;
            scattxs[wavidx] = 0.0;
            for (size_t vmridx = 0; vmridx < 2; vmridx++)
            {
                extxs[wavidx] += m_xsentry[vmrIndex[vmridx]].Interpolate(lowTemperatureIndex,
                                                                         lowTemperatureWeight,
                                                                         highTemperatureIndex,
                                                                         highTemperatureWeight,
                                                                         pressureIndex,
                                                                         pressureWeight,
                                                              1e7 / wavenumber[wavidx]) * vmrWeight[vmridx];
            }
            absxs[wavidx] = extxs[wavidx];
        }

    }
    return ok;
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::SetLocation(const GEODETIC_INSTANT& pt, bool* crosssectionschanged)
{
	bool ok = true;

	ok = ok && m_backgroundatmosphere->GetParameter(SKCLIMATOLOGY_TEMPERATURE_K, pt, &m_temperature, true);
	ok = ok && m_backgroundatmosphere->GetParameter(SKCLIMATOLOGY_PRESSURE_PA, pt, &m_pressure, false);

	if (m_broadener != nullptr)
	{
		// Get the broadener number density
		ok = ok && m_broadener->GetParameter(m_broadenerhandle, pt, &m_currentbroadenervmr, true);
		// now convert to VMR with the ideal gas law
		double airdens = m_pressure / nxcgs::KBOLTZMAN  / m_temperature * 10;
		m_currentbroadenervmr /= airdens;
	}
	else
	{
		m_currentbroadenervmr = 0;
	}

	return ok;
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::SetAtmosphericState(skClimatology* neutralatmosphere)
{
	bool ok = true;

	m_backgroundatmosphere = neutralatmosphere;
	m_backgroundatmosphere->AddRef();

	return ok;
}

skOpticalProperties_UserDefinedAbsorptionPressure::skOpticalProperties_UserDefinedAbsorptionPressure()
{
	m_quietwavelengthtruncation = true;
	m_backgroundatmosphere = nullptr;
	m_temperature = -999;
	m_pressure = -999;
	m_currentbroadenervmr = -999;
	m_broadener = nullptr;
}

skOpticalProperties_UserDefinedAbsorptionPressure::~skOpticalProperties_UserDefinedAbsorptionPressure()
{
	if (m_backgroundatmosphere != nullptr) 
		m_backgroundatmosphere->Release();
	if (m_broadener != nullptr)
		m_broadener->Release();
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::SetBroadener(skClimatology* broadener)
{
	m_broadener = broadener;
	if (m_broadener != nullptr)
	{
		m_broadener->AddRef();
	}
	return true;
}

bool skOpticalProperties_UserDefinedAbsorptionPressure::SetBroadenerHandle(const CLIMATOLOGY_HANDLE& handle)
{
	m_broadenerhandle = handle;
	return true;
}
#pragma clang diagnostic pop