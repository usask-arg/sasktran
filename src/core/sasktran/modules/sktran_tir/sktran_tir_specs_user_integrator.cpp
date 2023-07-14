/**
 * SASKTRAN TIR User Integrator Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_User_Integrator::ConfigureDefaults
 * 2019-04-22
 */
void SKTRAN_TIR_Specs_User_Integrator::ConfigureDefaults()
{
	m_opttype = OpticalPropertiesIntegratorTypeTIR::adaptive;
	m_srctype = SourceTermOrderTIR::order0;
	m_maxopticaldepth = 0.1;
	m_minextinctionratio = 0.9;
	m_extinctiontype = LayerExtinctionTypeTIR::linearwithheight;
}
