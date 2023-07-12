//#include "sktran_hr_internals.h"

class SKTRAN_HR_Specs_User_Integrator
{
	private:
		SKTRAN_HR_Integrator_Type	m_integratortype;
		double						m_adaptivemaxopticaldepth;
		double						m_adaptivemaxrayopticaldepth;
		double						m_maxextinctiongradient;
        bool                        m_useemissions;
        bool                        m_usesolartransmission;

	private:
		void						ConfigureDefaults();
	public:
									SKTRAN_HR_Specs_User_Integrator() { ConfigureDefaults(); }
		void						UseLegacySasktran21Technique( bool legacy );

		void						SetMaxOpticalDepth       ( double opticaldepth ) { m_adaptivemaxopticaldepth = opticaldepth; }
		void						SetMaxRayOpticalDepth(double opticaldepth) { m_adaptivemaxrayopticaldepth = opticaldepth; }
		void						SetMaxExtinctionGradient ( double gradient     ) { m_maxextinctiongradient = gradient; }
        void                        SetIntegratorType        ( SKTRAN_HR_Integrator_Type t ) {m_integratortype = t;}
        void                        SetUseEmissions          ( bool   use          ) { m_useemissions = use;}
        void                        SetUseSolarTransmission  ( bool   use          ) { m_usesolartransmission = use;}

		double						GetMaxExtinctionGradient ( ) const { return m_maxextinctiongradient; }
		double						GetMaxOpticalDepth       ( ) const { return m_adaptivemaxopticaldepth; }
		double						GetMaxRayOpticalDepth    ( ) const { return m_adaptivemaxrayopticaldepth; }
		SKTRAN_HR_Integrator_Type	GetIntegratorType        ( ) const { return m_integratortype; }
        bool                        GetUseEmissions          ( ) const { return m_useemissions;}
        bool                        GetUseSolarTransmission  ( ) const { return m_usesolartransmission;}
};
