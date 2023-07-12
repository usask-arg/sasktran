#include "include/sktran_hr_internals.h"

SKTRAN_HR_Specs_User_RayTracer::SKTRAN_HR_Specs_User_RayTracer()
{
	ConfigureDefaults();
}

SKTRAN_HR_Specs_User_RayTracer::~SKTRAN_HR_Specs_User_RayTracer()
{

}

bool SKTRAN_HR_Specs_User_RayTracer::ConfigureDefaults()
{
	//m_linesofsighttype	= SKTRAN_HR_RayTracer_Shells;
	m_linesofsighttype		= SKTRAN_HR_RayTracer_Straight_Generic;
	m_diffusetype			= SKTRAN_HR_RayTracer_Straight_Generic;
	m_solartype				= SKTRAN_HR_RayTracer_Straight_Generic;
	m_shellspacing			= 1000.0;
	m_setshellspacing		= false;

	size_t numshells = 101;
	double spacing = 1000.0;
	m_manualshells.resize(numshells);
	for (size_t idx = 0; idx < numshells; idx++)
	{
		m_manualshells[idx] = idx*spacing;
	}

	m_setmanualshells		= false;
	m_solarshellspacing		= 1000.0;
	m_setsolarshellspacing	= false;
	m_manualsolarshells		= m_manualshells;
	m_setmanualsolarshells	= false;
	m_curvedseparation		= 1000.0;
	m_usecurve				= true;
	m_groundshiftalt		= 0.0;
    m_toaHeight				= 100000.0;
	m_nadir_referencepoint_onground = false;

	return true;
}

bool SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters() const
{
	if (m_setmanualshells && m_setshellspacing)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters, shell spacing and manual shells have both been specified, the manual shells will be ignored.");
	}
	if ((m_setmanualshells) && (!m_setshellspacing))
	{
		double minshell = *std::min_element(m_manualshells.begin(), m_manualshells.end());
		double maxshell = *std::max_element(m_manualshells.begin(), m_manualshells.end());
		if ((minshell > 0.0) || (maxshell < m_toaHeight - m_groundshiftalt))
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters, custom shells do not span from the surface to top of atmosphere.");
		}
		if (minshell < m_groundshiftalt)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters, custom shells go below 0. The lowest shell height should be exactly 0.");
		}
	}

	if (m_setmanualsolarshells && m_setsolarshellspacing)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters, solar shell spacing and manual solar shells have both been specified, the manual shells will be ignored.");
	}
	if ((m_setmanualsolarshells) && (!m_setsolarshellspacing))
	{
		double minshell = *std::min_element(m_manualsolarshells.begin(), m_manualsolarshells.end());
		double maxshell = *std::max_element(m_manualsolarshells.begin(), m_manualsolarshells.end());		
		if ((minshell > 0.0) || (maxshell < m_toaHeight - m_groundshiftalt))
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters, custom solar shells do not span from the surface to top of atmosphere.");
		}
		if (minshell < m_groundshiftalt)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::CheckShellParameters, custom solar shells go below 0. The lowest shell height should be exactly 0.");
		}
	}
	return true;
}

bool SKTRAN_HR_Specs_User_RayTracer::SetManualShells(std::vector<double> customshells)
{
	m_manualshells = customshells;
	m_setmanualshells = true;
	if (!std::is_sorted(m_manualshells.begin(), m_manualshells.end()))
	{
		std::sort(m_manualshells.begin(), m_manualshells.end());
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::SetCustomShells, custom shells were out of order, they have been sorted.");
	}
	// round the highest shell up to the nearest integer: if it is <= top-of-atmosphere and rounds up to the nearest 0.001 it will crash
	if (m_manualshells.back() != std::ceil(m_manualshells.back()))
	{
		m_manualshells.back() = std::ceil(m_manualshells.back());
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::SetCustomSolarShells, the top shell has been rounded up to the nearest integer.");
	}
	return true;
}

bool SKTRAN_HR_Specs_User_RayTracer::SetManualSolarShells(std::vector<double> customshells)
{
	m_manualsolarshells = customshells;
	m_setmanualsolarshells = true;
	if (!std::is_sorted(m_manualsolarshells.begin(), m_manualsolarshells.end()))
	{
		std::sort(m_manualsolarshells.begin(), m_manualsolarshells.end());
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::SetCustomSolarShells, custom solar shells were out of order, they have been sorted.");
	}
	// round the highest shell up to the nearest integer: if it is <= top-of-atmosphere and rounds up to the nearest 0.001 it will crash
	if (m_manualsolarshells.back() != std::ceil(m_manualsolarshells.back()))
	{
		m_manualsolarshells.back() = std::ceil(m_manualsolarshells.back());
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_User_RayTracer::SetCustomSolarShells, the top shell has been rounded up to the nearest integer.");
	}
	return true;
}