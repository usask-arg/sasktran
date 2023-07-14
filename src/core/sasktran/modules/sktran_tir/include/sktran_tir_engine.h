/**
 * SaskTran Thermal Infra-Red Engine
 */

/**
 * SKTRAN_TIR_Engine
 * 2018-09-13
 */
class SKTRAN_TIR_Engine
{
private:
	SKTRAN_TIR_Specs_Internal_Core								m_internalspecs;
	std::vector<double>											m_wavelengths;

	IntegratorPtrTIR											m_integrator;
	SKTRAN_TIR_Thread_Manager									m_threadmanager;

	OpticalTablePtrTIR											m_opticalpropertiestable;

	RayFactoryPtrTIR											m_linesofsightrayfactory;
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>		m_coordinates;
	SKTRAN_TIR_RayTracingRegionManager							m_raymanager;
	SKTRAN_TIR_LinesOfSightTable								m_linesofsighttable;
	bool														m_docurvedrays;

	bool														m_calcwf;
	bool														m_dotemperaturewf;
	std::vector<std::vector<std::vector<std::vector<double>>>>	m_wf;	// shape of weighting functions is [wavelength][line_of_sight][species][perturbation_altitude]
	std::vector<CLIMATOLOGY_HANDLE>								m_wfspecies;
	std::vector<SKTRAN_TIR_Perturbation_Absorption_Linear>		m_wfperturbs;


private:
	void ReleaseResources();

	// Weighting function calculations
	bool CalculateWeightingFunctionsForRay(const size_t rayidx,
										   SKTRAN_RayOptical_Base* rayptr,
										   const size_t wavelidx);

	bool CalculateRadianceMultiThreaded(std::vector<std::vector<SKTRAN_StokesScalar>>* losradiance,
										SKTRAN_TIR_AtmosphericOpticalState* opticalstate);

public:
	SKTRAN_TIR_Engine();
	~SKTRAN_TIR_Engine();

	/*!
	 * Configures the engine with the user's specifications, sets the lines of
	 * sight, and sets the number of threads to perform the computation with.
	 * @param[in] modelspecifications User specifications.
	 * @param[in] linesofsight The lines of sight to calculate radiance at.
	 * @param[in] numthreads The number of threads. If no value is given or if
	 * a value of 0 is specified, the engine will set the number of threads to
	 * the maximum allowed by the user's system.
	 *
	 * @pre None.
	 * @post Radiance calculations can be made.
	 @ @return true if there no issues.
	 */
	bool ConfigureModel(SKTRAN_SpecsUser_Base& modelspecifications,
						const std::vector<double>& wavelen,
						const SKTRAN_LineOfSightArray_V21& linesofsight,
						size_t numthreads = 0);

	/*!
	 * Calculates the radiance at the specified wavelengths, along the lines of
	 * sight specified by calling ConfigureModel.
	 * @param[out] losradiance The scalar radiance result, which has units of
	 * photon / (s cm^2 sr nm)
	 * @param[in] wavelen The wavelengths to calculate radiance at; must be in
	 * increasing order with units of nm.
	 * @param[in] opticalstate The current optical state.
	 * @param[in] usecachedcrosssections If true, tells the engine to use cross
	 * sections from the previous radiance calculation. This option can only be
	 * used if a radiance calculation has already been performed and is
	 * intended for use in iterative retrieval computations where the model is
	 * ran successively with the only change between runs being the number
	 * density of atmospheric species whose weighting functions are to be
	 * calculated. Therefore, is this option is set to true, the user must be
	 * careful to only change the number density of species which they also
	 * requested weighting functions to be computed for.
	 *
	 * @pre Engine must be configured by calling ConfigureModel.
	 * @post losradiance and losvector will be resized to
	 * [# wavelengths][# lines of sight]. Weighting functions from this
	 * calculation can be retrieved by calling GetWF, and weighting functions
	 * from any previous calculation can no longer be retrieved.
	 * @return true if there were no issues.
	 */
	bool CalculateRadiance(std::vector<std::vector<SKTRAN_StokesScalar>>* losradiance,
						   SKTRAN_TIR_AtmosphericOpticalState* opticalstate,
						   bool usecachedcrosssections = false);

	/*!
	 * Retrieve the weighting functions from the most recent radiance calculation.
	 * The shape of the returned array is [wavelength][line of sight][species][perturbation altitude]
	 *
	 * @pre Performed the calculation by calling CalculateRadiance without error.
	 * @post None.
	 * @return Array of shape [wavelength][line of sight][species][perturbation altitude]
	 */
	const std::vector<std::vector<std::vector<std::vector<double>>>>& GetWF() { return m_wf; }

	virtual GEODETIC_INSTANT ReferencePoint() const;
	const SKTRAN_CoordinateTransform_V2* Coordinates() { return m_coordinates.get(); }
	
	double WFHeightAt(size_t idx) { return m_wfperturbs[idx].PerturbationLocation(*m_coordinates).Altitude(); }
	size_t NumWF() const { return m_wfperturbs.size(); }
	std::vector<double> WFHeights() const;
	
	const SKTRAN_TIR_Specs_Internal_Core& InternalSpecs() const { return m_internalspecs; };
	SKTRAN_TIR_Specs_Internal_Core* InternalSpecsVar() { return &m_internalspecs; }
	SKTRAN_TIR_LinesOfSightTable& LinesOfSight() { return m_linesofsighttable; }
	SKTRAN_TIR_RayTracingRegionManager* RayTracingRegionManager() { return &m_raymanager; }
};
