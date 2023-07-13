

class skOpticalProperties_UserDefinedScatterConstantHeight : public skOpticalProperties
{
private:
	nx2dArray<double> m_a1; // Legendre phase moments
	nx2dArray<double> m_a2;
	nx2dArray<double> m_a3;
	nx2dArray<double> m_a4;
	nx2dArray<double> m_b1;
	nx2dArray<double> m_b2;

	nx1dArray<double> m_wavelengths; // Wavelengths in nm
	nx1dArray<double> m_xs_scat; // Scattering cross section in cm2
	nx1dArray<double> m_xs_abs; // Absorption cross section in cm2

	void InterpolationWeights(double wavelength, std::array<size_t, 2>& indicies, std::array<double, 2>& weights);

public:
	skOpticalProperties_UserDefinedScatterConstantHeight() {};
	virtual ~skOpticalProperties_UserDefinedScatterConstantHeight() {};

	virtual bool						SetAtmosphericState(skClimatology* neutralatmosphere) { return true; };
	virtual bool						SetLocation(const GEODETIC_INSTANT& pt, bool* crosssectionschanged) { *crosssectionschanged = false; return true; };
	virtual bool						InternalClimatology_UpdateCache(const GEODETIC_INSTANT& pt) { return true; };
	virtual bool						CalculateCrossSections(double wavenumber, double* absxs, double* extxs, double* scattxs);
	virtual bool						CalculatePhaseMatrix(double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix);
	virtual bool						CalculateP11(double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11);
	virtual bool						LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff);
	virtual bool                        LegendreCoefficientsPolarized(double wavenumber, double* a1, double* a2, double* a3, double* a4, double* b1, double* b2, int usermaxcoeff, int& opticalmaxcoeff);
	virtual bool						IsScatterer() const { return true; };
	virtual bool						IsAbsorber() const { return true; };
	virtual bool						IsHeightDependent() const { return false; };

	void SetWavelengths(const nx1dArray<double>& wavelengths) { m_wavelengths = wavelengths; };
	void SetLegendreMoments(const nx2dArray<double>& lm) { m_a1 = lm; };
	void SetLegendrea1(const nx2dArray<double>& lm) { m_a1 = lm; };
	void SetLegendrea2(const nx2dArray<double>& lm) { m_a2 = lm; };
	void SetLegendrea3(const nx2dArray<double>& lm) { m_a3 = lm; };
	void SetLegendrea4(const nx2dArray<double>& lm) { m_a4 = lm; };
	void SetLegendreb1(const nx2dArray<double>& lm) { m_b1 = lm; };
	void SetLegendreb2(const nx2dArray<double>& lm) { m_b2 = lm; };

	void SetScatteringCrossSection(const nx1dArray<double>& xs_scat) { m_xs_scat = xs_scat; };
	void SetAbsorptionCrossSection(const nx1dArray<double>& xs_abs) { m_xs_abs = xs_abs; };

};
