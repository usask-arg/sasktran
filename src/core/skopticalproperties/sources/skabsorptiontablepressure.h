#pragma once

// Implements a user defined absorption 

class skAbsorptionPressureTemperatureEntry
{
private:
	nx3dArray<double> m_xs;
	nx2dArray<double> m_temperature;
	nx1dArray<double> m_pressure;
	nx1dArray<double> m_wavelength;

public:
	static std::array<double, 2> InterpolationWeights(const std::array<size_t, 2>& indicies,
											          const std::array<double, 2>& values,
											          double value);

public:
	skAbsorptionPressureTemperatureEntry(const nx3dArray<double>& xs,
		                                 const nx2dArray<double>& temperature,
		                                 const nx1dArray<double>& pressure,
		                                 const nx1dArray<double>& wavelength) : m_xs{ xs }, m_temperature{ temperature }, m_pressure{ pressure }, m_wavelength{ wavelength } {}

	double Interpolate(double temperature, double pressure, double wavelength);

    void calcInterpolationWeights(double temperature,
                                  double pressure,
                                  std::array<size_t, 2>& lowTemperatureIndex,
                                  std::array<double, 2>& lowTemperatureWeight,
                                  std::array<size_t, 2>& highTemperatureIndex,
                                  std::array<double, 2>& highTemperatureWeight,
                                  std::array<size_t, 2>& pressureIndex,
                                  std::array<double, 2>& pressureWeight);

    double Interpolate(const std::array<size_t, 2>& lowTemperatureIndex,
                       const std::array<double, 2>& lowTemperatureWeight,
                       const std::array<size_t, 2>& highTemperatureIndex,
                       const std::array<double, 2>& highTemperatureWeight,
                       const std::array<size_t, 2>& pressureIndex,
                       const std::array<double, 2>& pressureWeight,
                       double wavelength);

};


class skOpticalProperties_UserDefinedAbsorptionPressure : public skOpticalProperties
{
private:
	skClimatology*													m_backgroundatmosphere;
	bool															m_quietwavelengthtruncation; // If true then dont fail a cross-section request outside the range of wavelengths
	RefractiveIndexDryAirSTP										m_refractiveindexair;		 // used to convert wavelength in air to wavelength in vacuum.
	bool															m_wavelengthinvacuum;		 // True if the wavelengthtables are in vacuum

	// Cross section is a function of pressure, temperature, interfering species, and wavelength
	// Temperature is also a function of pressure
	std::vector<skAbsorptionPressureTemperatureEntry>               m_xsentry;
	std::vector<double>												m_broadenervmr;
	double															m_temperature;
	double															m_pressure;
	CLIMATOLOGY_HANDLE												m_broadenerhandle;
	skClimatology*													m_broadener;
	bool															m_broadenerincluded;
	double															m_currentbroadenervmr;

public:
	skOpticalProperties_UserDefinedAbsorptionPressure();
	virtual									   ~skOpticalProperties_UserDefinedAbsorptionPressure() override;
	void										SetVacuumWavelengths(bool isvacuum) { m_wavelengthinvacuum = isvacuum; }		// True if the wavelengthtables are in vacuum
	void										SetQuietWavelengthTruncation(bool bequiet) { m_quietwavelengthtruncation = bequiet; }

public:
	virtual bool								SetAtmosphericState(skClimatology* neutralatmosphere)  override;
	virtual bool								SetLocation(const GEODETIC_INSTANT& pt, bool* crosssectionschanged)  override;
	virtual bool								CalculateCrossSections(double wavenumber, double* absxs, double* extxs, double* scattxs) override;
    virtual bool						        CalculateCrossSectionsArray(const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs) override;
    virtual bool								IsScatterer() const  override { return false; }
	virtual bool								IsAbsorber() const  override { return true; }
	virtual bool								SetBroadener(skClimatology* broadener);
	virtual bool								SetBroadenerHandle(const CLIMATOLOGY_HANDLE& handle);
	virtual bool								InternalClimatology_UpdateCache(const GEODETIC_INSTANT& /*pt*/)  override { return true; }
	bool										AddEntry(const nx3dArray<double>& xs,
														 const nx2dArray<double>& temperature,
														 const nx1dArray<double>& pressure,
														 const nx1dArray<double>& wavelength,
														 double interferingvmr);
};

