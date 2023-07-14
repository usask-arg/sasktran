/**
 * SASKTRAN TIR User Specifications
 * 2019-04-16
 */

/**
 * SKTRAN_TIR_Specs_User
 *
 * Interface for allowing the user to configure the model, passed to
 * SKTRAN_TIR_Engine::ConfigureModel.
 */
class SKTRAN_TIR_Specs_User : public SKTRAN_Specifications_Base
{
private:
	SKTRAN_TIR_Specs_User_RayTracer m_raytracerspecs;
	SKTRAN_TIR_Specs_User_Integrator m_integratorspecs;
	SKTRAN_TIR_Specs_User_OpticalPropertiesTable m_optpropspecs;
	SKTRAN_TIR_Specs_User_wf m_wfspecs;
	bool m_calcwf;

public:
	SKTRAN_TIR_Specs_User() { }

	SKTRAN_TIR_Specs_User_RayTracer& RayTracingSpecs() { return m_raytracerspecs; }
	const SKTRAN_TIR_Specs_User_RayTracer& RayTracingSpecsConst() const { return m_raytracerspecs; }
	SKTRAN_TIR_Specs_User_Integrator& IntegratorSpecs() { return m_integratorspecs; }
	const SKTRAN_TIR_Specs_User_Integrator& IntegratorSpecsConst() const { return m_integratorspecs; }
	SKTRAN_TIR_Specs_User_OpticalPropertiesTable& OpticalPropertiesSpecs() { return m_optpropspecs; }
	const SKTRAN_TIR_Specs_User_OpticalPropertiesTable& OpticalPropertiesSpecsConst() const { return m_optpropspecs; }
	SKTRAN_TIR_Specs_User_wf& WeightingFunctionSpecs() { return m_wfspecs; }
	const SKTRAN_TIR_Specs_User_wf& WeightingFunctionSpecsConst() const { return m_wfspecs; }
};
