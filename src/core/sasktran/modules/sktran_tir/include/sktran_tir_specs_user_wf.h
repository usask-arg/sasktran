/**
 * SASKTRAN TIR User Weighting Function Specifications
 */

/**
 * SKTRAN_TIR_Specs_User_wf
 *
 * Allows the user to configure properties of the TIR weighting functions.
 * Available settings are:
 *
 *   WFSpecies - Species to calculate weighting functions for. Supports
 *     calculating the number density jacobians of absorbing species.
 *   WeightingFunctionHeights - If specified, this vector sets manual heights
 *     to calculate weighting functions. If not specified, evenly spaced
 *     heights are used with a spacing given by WFRes.
 *   WeightingFunctionWidths - If specified, this vector sets manual widths
 *     for each perturbation. This is the vertical distance from the height
 *     of the perturbation to the point where the perturbation is zero. Must
 *     be the same length as WeightingFunctionHeights
 *   WFRes - Sets a constant spacing in meters between perturbation altitudes.
 *   MaxWFHeight - Sets the maximum perturbation altitude, in meters. If manual
 *     heights are not specified, the perturbations will be spaced evenly from
 *     0.5 * WFRes to MaxWFHeight with a spacing of WFRes
 *   WFInterpWidth - Sets a constant width for all perturbations, in meters.
 */
class SKTRAN_TIR_Specs_User_wf
{
private:
	std::vector<CLIMATOLOGY_HANDLE>		m_wfspecies;
	bool								m_dotemperaturewf;
	std::vector<double>					m_manualwfheights;
	std::vector<double>					m_manualwfwidths;
	double								m_manualwfresolution;
	double								m_maxwfheight;
	double								m_wfinterpwidth;
	SpeciesWFUnitTIR					m_wfunit;

private:
	void ConfigureDefaults();

public:
	SKTRAN_TIR_Specs_User_wf();

	void AddWeightingFunctionSpecies(const CLIMATOLOGY_HANDLE& species);
	void SetWeightingFunctionHeight(const std::vector<double>& heights) { m_manualwfheights = heights; }
	void SetWeightingFunctionWidth(const std::vector<double>& widths) { m_manualwfwidths = widths; }
	void SetWFRes(double val) { m_manualwfresolution = val; }
	void SetMaxWFHeight(double val) { m_maxwfheight = val; }
	void SetWFInterpWidth(double val) { m_wfinterpwidth = val; }
	void SetSpeciesWFUnit(SpeciesWFUnitTIR unit) { m_wfunit = unit; }
	void ClearWFSpecies() { m_wfspecies.clear(); }
	void SetDoTemperatureWF(bool dotemperaturewf) { m_dotemperaturewf = dotemperaturewf; }

	size_t NumWFSpecies() const { return m_wfspecies.size(); }
	const std::vector<CLIMATOLOGY_HANDLE>& GetWFSpecies() const { return m_wfspecies; }
	const std::vector<double>& GetWeightingFunctionHeight() const { return m_manualwfheights; }
	const std::vector<double>& GetWeightingFunctionWidth() const { return m_manualwfwidths; }
	double GetWFRes() const { return m_manualwfresolution; }
	double GetMaxWFHeight() const { return m_maxwfheight; }
	double GetWFInterpWidth() const { return m_wfinterpwidth; }
	SpeciesWFUnitTIR GetSpeciesWFUnit() const { return m_wfunit; }
	const bool GetDoTemperatureWF() const { return m_dotemperaturewf; }

	friend class SKTRAN_TIR_Specs_Internal_wf;
	friend class SKTRAN_TIR_Engine;
};
