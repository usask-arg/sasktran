/**
 * SASKTRAN TIR Optical Properties Table
 */

/**
 * SKTRAN_TIR_TableOpticalProperties
 * 2018-10-09
 *
 * Optical property table for the TIR engine which can be 1D, 2D, or 3D
 * depending on the type of unit sphere used. Basically, this is identical
 * to SKTRAN_TableOpticalProperties_3D_UnitSphere but anything related to
 * scattering has been removed to save memory and computation time.
 */
class SKTRAN_TIR_TableOpticalProperties : public nxUnknown
{
protected:
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coordinates;

	const SKTRAN_UnitSphere_V2*								m_unitsphere;
	const SKTRAN_GridDefOpticalPropertiesRadii_V21*			m_alts;
	std::vector<std::vector<std::vector<double>>>			m_absextinction;	// arrays of [wavelength][horizontal location][altitude]
	std::vector<std::vector<std::vector<double>>>			m_emission;
	std::vector<std::vector<std::vector<double>>>			m_groundemission;
	std::vector<std::vector<std::vector<double>>>			m_dk_dT;			// derivative of absorption extinction with respect to temperature
	std::vector<std::vector<std::vector<double>>>			m_dB_dT;			// derivative of Planck function with respect to temperature
	bool													m_calctemperaturederivatives;
	std::vector<std::vector<double>>						m_airnumberdensity;		// array of [horizontal location][altitude]
	double													m_groundemissivity;
	double													m_mjd;
	mutable std::vector<size_t>								m_speedhelper;

	std::vector<double>										m_wavelengtharray;
	bool													m_firsttime;
	bool													m_forcecacheupdates;

	std::map<CLIMATOLOGY_HANDLE, std::vector<std::vector<std::vector<double>>>> m_speciesxs;
	std::map <CLIMATOLOGY_HANDLE, std::vector<std::vector<double>>> m_speciesnumberdensity_previousrun;

private:
	virtual double InterpTable(const std::vector<std::vector<double>>& table,
							   const HELIODETIC_POINT& loc) const;
	virtual bool FillTablesAtIndexMultiWavel(size_t altidx,
											 size_t locidx,
											 SKTRAN_TIR_AtmosphericOpticalState& opticalstate);
	virtual bool FillGroundEmissionTableAtIndexMultiWavel(size_t locidx,
														  SKTRAN_TIR_AtmosphericOpticalState& opticalstate);

	double PlanckBlackbody(double wavelen_nm,
						   double temperature_K);

	SKTRAN_GridIndex TableSubToInd(size_t wavidx,
								   SKTRAN_GridIndex locidx,
								   SKTRAN_GridIndex altidx) const;
	void ReleaseResources();

protected:
	virtual bool CalcSphereIndices(const HELIODETIC_POINT& loc,
								   double* weights,
								   size_t* indices,
								   size_t& numindex) const;
	virtual bool CalcAltIndices(const HELIODETIC_POINT& loc,
								double* weights,
								size_t* indices,
								size_t& numindex) const;

public:
	SKTRAN_TIR_TableOpticalProperties();
	~SKTRAN_TIR_TableOpticalProperties();
	void SetCoords(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coordinates) { m_coordinates = coordinates; }
	const SKTRAN_CoordinateTransform_V2* CoordinatesPtr() const { return m_coordinates.get(); }
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2> Coordinates() const { return m_coordinates; }

public:
	bool ConfigureOptical(SKTRAN_TIR_AtmosphericOpticalState& opticalstate);
	bool ConfigureOpticalFromCache(SKTRAN_TIR_AtmosphericOpticalState& opticalstate);
	bool ConfigureGeometry(const SKTRAN_SpecsUser_Base* specs);
	bool SetAltitudes(const SKTRAN_GridDefOpticalPropertiesRadii_V21& alts);
	bool SetUnitSphere(const SKTRAN_UnitSphere_V2& unitsphere);
	bool SetWavelengths(std::vector<double>& wavelengths);
	void SetForceCacheUpdates(bool force) { m_forcecacheupdates = force; }
	bool GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayStorage_Base* r,
												size_t startPtIndex,
												double* sigma0,
												double* sigma1,
												size_t wavelidx) const;
	bool GetEffectiveExtinctionPerCMWithHeight1(const HELIODETIC_POINT& startpoint,
												const HELIODETIC_POINT& endpoint,
												double* sigma0,
												double* sigma1,
												size_t wavelidx) const;
	bool TotalExtinctionPerCM(const HELIODETIC_POINT& point,
							  std::vector<double>& extinction) const;
	double TotalExtinctionPerCMAtWavel(const HELIODETIC_POINT& point,
									   size_t wavelidx) const;
	bool SpeciesCrossSectionCM2(const HELIODETIC_POINT& point,
								const CLIMATOLOGY_HANDLE& species,
								std::vector<double>& xs) const;
	double SpeciesCrossSectionCM2AtWavel(const HELIODETIC_POINT& point,
										 const CLIMATOLOGY_HANDLE& species,
										 size_t wavelidx) const;

	double AirNumberDensityAtPoint(const HELIODETIC_POINT& point) const;

	// from emission table
	bool SourceTermAtPoint(const HELIODETIC_POINT& point,
						   std::vector<double>& source) const;
	bool GroundSourceAtPoint(const HELIODETIC_POINT& point,
							 std::vector<double>& source) const;
	double SourceTermAtPointAndWavel(const HELIODETIC_POINT& point,
									 size_t wavelidx) const;
	double GroundSourceAtPointAndWavel(const HELIODETIC_POINT& point,
									   size_t wavelidx) const;

	// for temperature derivatives
	bool PlanckFunctionTemperatureDerivativeAtPoint(const HELIODETIC_POINT& point,
													std::vector<double>& derivative) const;
	bool AbsorptionTemperatureDerivativeAtPoint(const HELIODETIC_POINT& point,
												std::vector<double>& derivative) const;
	double PlanckFunctionTemperatureDerivativeAtPointAndWavel(const HELIODETIC_POINT& point,
															  size_t wavelidx) const;
	double AbsorptionTemperatureDerivativeAtPointAndWavel(const HELIODETIC_POINT& point,
														  size_t wavelidx) const;
};
