//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Integrator		2014-10-30*/
/**  Factory class for creating integrators used within the HR model
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_Internal_Integrator
{
	private:
		SKTRAN_HR_Integrator_Type		m_integratortype;
		double							m_maxadaptiveopticaldepth;
		double							m_maxadaptiverayopticaldepth;
		double							m_maxextinctiongradient;
        bool                            m_usesolartransmission;
        bool                            m_useemissions;

	private:
		virtual bool					ConfigureDefaults();
	public:
										SKTRAN_HR_Specs_Internal_Integrator();
		virtual							~SKTRAN_HR_Specs_Internal_Integrator();

		virtual bool					Configure								( const SKTRAN_HR_Specs_User& specs );
		virtual bool					CreateIntegrator						( const SKTRAN_TableOpticalProperties_Base& optproptable, OptIntegratorPtr& optint, SrcIntegratorPtr& srcint );
        bool                            GetUseSolarTransmission                 ( ) const;
        bool                            GetUseEmissions                         ( ) const;
};
