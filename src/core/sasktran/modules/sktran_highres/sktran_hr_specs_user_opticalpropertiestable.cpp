#include "include/sktran_hr_internals.h"

bool SKTRAN_HR_Specs_User_OpticalPropertiesTable::ConfigureDefaults()
{
	m_optproptype = SKTRAN_HR_OpticalPropertiesTableType_1d;
	m_maxPolarizationOrder = 0;
	m_polHigherOrderBehaviour = SKTRAN_HR_PolHOType::unpolarized;
	m_atmosphereHasDelta = SKTRAN_HR_AtmosphereHasDelta::no;
	m_numprofiles = 30;
	m_heightres = 500;
	m_scatres = 0.5;

	m_numcones = 10;
	m_coneanglesep = 1;
	m_profilepercone = 10;

	// Invalid values mean use normal and reference from the line of sight plane
	m_normal.SetInvalid();
	m_reference.SetInvalid();

	m_forcecacheupdates = false;

	return true;
}
