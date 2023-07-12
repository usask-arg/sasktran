#include "include/sktran_hr_internals.h"

void SKTRAN_HR_Specs_User_Integrator::ConfigureDefaults()
{
	m_integratortype          = SKTRAN_HR_IntegratorType_Adaptive;
	m_adaptivemaxopticaldepth = 1e6;
	m_maxextinctiongradient   = 0.9;
    m_usesolartransmission    = true;
    m_useemissions            = false;
}

void SKTRAN_HR_Specs_User_Integrator::UseLegacySasktran21Technique( bool legacy )
{
	if( legacy )
	{
		m_integratortype = SKTRAN_HR_IntegratorType_Straight;
	}
	else
	{
		m_integratortype = SKTRAN_HR_IntegratorType_Adaptive;
	}
}
