/**
 * SASKTRAN TIR Internal Integrator Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_Internal_Integrator::SKTRAN_TIR_Specs_Internal_Integrator
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_Integrator::SKTRAN_TIR_Specs_Internal_Integrator()
{
	m_opttype = OpticalPropertiesIntegratorTypeTIR::adaptive;
	m_srctype = SourceTermOrderTIR::order0;
	m_maxopticaldepth = 0.0;
	m_minextinctionratio = 0.0;
	m_docurvedrays = false;
	m_extinctiontype = LayerExtinctionTypeTIR::linearwithheight;
	m_wfunit = SpeciesWFUnitTIR::numberdensity;
}

/**
 * SKTRAN_TIR_Specs_Internal_Integrator::~SKTRAN_TIR_Specs_Internal_Integrator
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_Integrator::~SKTRAN_TIR_Specs_Internal_Integrator()
{

}

/**
 * SKTRAN_TIR_Specs_Internal_Integrator::ConfigureDefaults
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_Integrator::ConfigureDefaults()
{
	bool ok = true;

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_Integrator::Configure
 * 2019-04-22
 *
 * Copies user settings to the internal specifications object.
 *
 * @param[in] specs The SKTRAN_TIR_Specs_User object containing all user settings
 */
bool SKTRAN_TIR_Specs_Internal_Integrator::Configure(
	const SKTRAN_TIR_Specs_User& specs)
{
	bool ok = true;
	
	m_opttype = specs.IntegratorSpecsConst().GetOpticalPropertiesType();
	m_srctype = specs.IntegratorSpecsConst().GetSourceTermOrder();
	m_maxopticaldepth = specs.IntegratorSpecsConst().GetMaxOpticalDepth();
	m_minextinctionratio = specs.IntegratorSpecsConst().GetMinExtinctionRatio();
	m_docurvedrays = specs.RayTracingSpecsConst().GetUseCurve();
	m_extinctiontype = specs.IntegratorSpecsConst().GetLayerExtinctionType();
	m_wfunit = specs.WeightingFunctionSpecsConst().GetSpeciesWFUnit();

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_Integrator::CreateIntegrator
 * 2019-04-22
 *
 * Creates the line of sight integrator for TIR calculations.
 *
 * @param[in] optproptable The optical properties table, which the created integrator keeps a reference to
 * @param[out] integrator The line of sight integrator is allocated in memory and ready for use
 *
 * @pre This object's Configure method has been called with the desired user settings
 */
bool SKTRAN_TIR_Specs_Internal_Integrator::CreateIntegrator(
	const SKTRAN_TIR_TableOpticalProperties& optproptable,
	IntegratorPtrTIR& integrator)
{
	bool ok = true;

	std::unique_ptr<SKTRAN_TIR_Integrator> integrator_local(new SKTRAN_TIR_Integrator);

	integrator_local->AddRef();
	ok = ok && integrator_local->SetOpticalProps(&optproptable);
	integrator_local->SetMaxOpticalDepthOfCell(m_maxopticaldepth);
	integrator_local->SetMinExtinctionRatioOfCell(m_minextinctionratio);
	integrator_local->SetLayerExtinctionType(m_extinctiontype);
	integrator_local->SetOpticalPropertiesIntegrationType(m_opttype);
	integrator_local->SetSourceTermOrder(m_srctype);
	integrator_local->SetSpeciesWFUnit(m_wfunit);
	integrator = std::move(integrator_local);

	return ok;
}
