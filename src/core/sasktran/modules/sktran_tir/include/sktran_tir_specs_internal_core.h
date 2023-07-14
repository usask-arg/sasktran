/**
 * SASKTRAN TIR Internal Core Specifications
 * 2019-04-16
 */

/**
* SKTRAN_TIR_Specs_Internal_Core
*
* Factory class is responsible for creating and configuring objects required 
* by the engine.  A SKTRAN_TIR_Specs_User object is passed in first to
* initialize any user settings, then base class pointers are bassed in and
* set to created derived objects.
*/
class SKTRAN_TIR_Specs_Internal_Core
{
private:
	SKTRAN_TIR_RayTracingRegionManager						m_raymanager;
	const SKTRAN_LineOfSightArray_V21*						m_linesofsight;
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coords;
	bool													m_in3dmode;
	bool													m_calcwf;
	SKTRAN_TIR_Specs_Internal_RayTracer						m_raytracerspecs;
	SKTRAN_TIR_Specs_Internal_Integrator					m_integratorspecs;
	SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable		m_opttablespecs;
	SKTRAN_TIR_Specs_Internal_wf							m_wfspecs;

private:
	bool ReleaseResources();
	bool RotateStartAndEnd(const HELIODETIC_UNITVECTOR& look,
						   const HELIODETIC_UNITVECTOR& in,
						   HELIODETIC_UNITVECTOR& out,
						   double theta);
	bool CalcReferencePoints(nxVector& in,
							 nxVector& out);
	const SKTRAN_CoordinateTransform_V2* CoordinatePtr() const { return m_coords.get(); }

public:
	SKTRAN_TIR_Specs_Internal_Core();
	virtual ~SKTRAN_TIR_Specs_Internal_Core();

	SKTRAN_TIR_Specs_Internal_RayTracer& RayTracerSpecs() { return m_raytracerspecs; }
	SKTRAN_TIR_Specs_Internal_Integrator& IntegratorSpecs() { return m_integratorspecs; }
	SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable& OpticalPropertiesSpecs() { return m_opttablespecs; }
	const SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable& OpticalPropertiesSpecs() const { return m_opttablespecs; }
	SKTRAN_TIR_Specs_Internal_wf& WeightingFunctionSpecs() { return m_wfspecs; }
	const SKTRAN_TIR_Specs_Internal_wf& WeightingFunctionSpecs() const { return m_wfspecs; }
	
	bool Is3dMode() { return m_in3dmode; }
	bool CalcWf() const { return m_wfspecs.DoWfCalculation(); }
	SKTRAN_TIR_RayTracingRegionManager& RayManager() { return m_raymanager; }

public:
	bool Configure(const SKTRAN_Specifications_Base& specs,
				   const SKTRAN_LineOfSightArray_V21& linesofsight);
	bool CreateCoordinates(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* coords,
						   const SKTRAN_LineOfSightArray_V21& linesofsight,
						   double toaHeight,
						   nxGeodetic::GEOID_MODEL geoidmodel,
						   bool userdefinedgeoidmodel);
};
