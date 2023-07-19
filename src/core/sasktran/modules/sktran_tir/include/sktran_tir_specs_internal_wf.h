/**
 * SASKTRAN TIR Internal Weighting Function Specs
 */

/**
 * SKTRAN_TIR_Specs_Internal_wf
 *
 * Factory class which reads in the user's weighting function specifications
 * and creates a list of perturbations for each species and altitude that
 * weighting functions are to be calculated at. Supports perturbations which
 * vary uniformly with altitude.
 */
class SKTRAN_TIR_Specs_Internal_wf
{
private:
	std::vector<CLIMATOLOGY_HANDLE>							m_wfspecies;
	bool													m_dotemperaturewf;
	std::vector<double>										m_manualwfheights;
	std::vector<double>										m_manualwfwidths;
	double													m_manualwfresolution;
	double													m_maxwfheight;
	double													m_wfinterpwidth;
	std::vector<nxVector>									m_normals;
	std::vector<double>										m_alts;
	std::vector<SKTRAN_TIR_Perturbation_Absorption_Linear>	m_wfperturbs;
	
private:
	bool MakeOneDimUniformPerturbations();

public:
	bool Configure(const SKTRAN_TIR_Specs_User& specs);

	bool MakePerturbationList();

	bool DoWfCalculation() const { return m_wfspecies.size() > 0; }

	bool DoTemperatureWFCalculation() const { return m_dotemperaturewf; }

	bool AddWfGeometryToRayTracer(SKTRAN_RayTracer_Straight_Generic& raytracer,
								  std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords) const;

	const std::vector<CLIMATOLOGY_HANDLE> WFSpecies() const { return m_wfspecies; }

	std::vector<SKTRAN_TIR_Perturbation_Absorption_Linear>& PertList() { return m_wfperturbs; }
	const std::vector<SKTRAN_TIR_Perturbation_Absorption_Linear>& PertList() const { return m_wfperturbs; }
};
