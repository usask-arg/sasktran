//#include "sktran_hr_internals.h"


/** @defgroup HR_SPECS High Resolution User Specs
 *	@ingroup sktranhr
 *  The user specifications for the high resolution engine.  Controls
 *  specs of individual objects in the engine, including,
 *    - Ray Tracing
 *    - Diffuse Profile Management
 *    - Integration Techniques
 *    - Optical properties
 *
 *  Upon creation of this object, defaults are set for each of the above listed
 *  areas.  These defaults are outlined in detail in their corresponding classes.
 **/

/*-----------------------------------------------------------------------------
 *					SKTRAN_Specs_User		2013-08-23*/
/**  @ingroup HR_SPECS
 *   Container for more specific user specifications, contains specs objects
 *   for the
 *      - Ray Tracing
 *      - Diffuse Profile Management
 *      - Integration Techniques
 *		- Optical Properties ( 1D, 3D )
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_User : public SKTRAN_Specifications_Base
{
	private:
		SKTRAN_HR_Specs_User_RayTracer					m_raytracerspecs;				//!< The raytracing specifications
		SKTRAN_HR_Specs_User_Diffuse					m_diffusespecs;					//!< Specifications for diffuse profile management
		SKTRAN_HR_Specs_User_Integrator					m_integratorspecs;				//!< Specifications for different integration techniques
		SKTRAN_HR_Specs_User_OpticalPropertiesTable		m_optpropspecs;					//!< Specifications for different optical properties tables
		SKTRAN_HR_Specs_User_wf							m_wfspecs;						//!< Specifications for weighting function calculations
		bool											m_calcwf;
		int 											m_scatterorder;

	public:
		SKTRAN_HR_Specs_User() { m_scatterorder = 50; }

		SKTRAN_HR_Specs_User_RayTracer&					RayTracingSpecs() { return m_raytracerspecs; }
		const SKTRAN_HR_Specs_User_RayTracer&			RayTracingSpecsConst() const { return m_raytracerspecs; }
		SKTRAN_HR_Specs_User_Diffuse&					DiffuseSpecs() { return m_diffusespecs; }
		const SKTRAN_HR_Specs_User_Diffuse&				DiffuseSpecsConst() const { return m_diffusespecs; }
		SKTRAN_HR_Specs_User_Integrator&				IntegratorSpecs() { return m_integratorspecs; }
		const SKTRAN_HR_Specs_User_Integrator&			IntegratorSpecsConst() const { return m_integratorspecs; }
		SKTRAN_HR_Specs_User_OpticalPropertiesTable&	OpticalPropertiesSpecs() { return m_optpropspecs; }
		const SKTRAN_HR_Specs_User_OpticalPropertiesTable& OpticalPropertiesSpecsConst() const { return m_optpropspecs; }
		SKTRAN_HR_Specs_User_wf&						WeightingFunctionSpecs() { return m_wfspecs; }
		const SKTRAN_HR_Specs_User_wf&					WeightingFunctionSpecsConst() const { return m_wfspecs; }
		
		void 											SetScatterOrder(int scatterorder) { m_scatterorder = scatterorder; }
		int 											ScatterOrder() const { return m_scatterorder; }
};
