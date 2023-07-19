/**
 * SASKTRAN TIR User Ray Tracer Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_User_RayTracer::SKTRAN_TIR_Specs_User_RayTracer
 * 2019-04-22
 */
SKTRAN_TIR_Specs_User_RayTracer::SKTRAN_TIR_Specs_User_RayTracer()
{
	ConfigureDefaults();
}

/**
 * SKTRAN_TIR_Specs_User_RayTracer::~SKTRAN_TIR_Specs_User_RayTracer
 * 2019-04-22
 */
SKTRAN_TIR_Specs_User_RayTracer::~SKTRAN_TIR_Specs_User_RayTracer()
{

}

/**
 * SKTRAN_TIR_Specs_User_RayTracer::ConfigureDefaults
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_User_RayTracer::ConfigureDefaults()
{
	m_linesofsighttype = RayTracerTypeTIR::straight;
	m_shellspacing = 1000.0;
	m_setshellspacing = false;

	size_t numshells = 101;
	double spacing = 1000.0;
	m_manualshells.resize(numshells);
	for (size_t idx = 0; idx < numshells; idx++)
	{
		m_manualshells[idx] = idx * spacing;
	}

	m_setmanualshells = false;
	m_usecurve = true;
	m_groundshiftalt = 0.0;
	m_toaHeight = 100000.0;

	m_setmanualgeoidmodel = false;
	m_geoidmodel = nxGeodetic::WGS84;

	return true;
}

/**
 * SKTRAN_TIR_Specs_User_RayTracer::CheckShellParameters
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_User_RayTracer::CheckShellParameters() const
{
	if (m_setmanualshells && m_setshellspacing)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Specs_User_RayTracer::CheckShellParameters, shell spacing and manual shells have both been specified, the manual shells will be ignored.");
	}
	if ((m_setmanualshells) && (!m_setshellspacing))
	{
		double minshell = *std::min_element(m_manualshells.begin(), m_manualshells.end());
		double maxshell = *std::max_element(m_manualshells.begin(), m_manualshells.end());
		if ((minshell > m_groundshiftalt) || (maxshell < m_toaHeight))
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Specs_User_RayTracer::CheckShellParameters, custom shells do not span from the surface to top of atmosphere.");
		}
	}
	return true;
}

/**
 * SKTRAN_TIR_Specs_User_RayTracer::SetManualShells
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_User_RayTracer::SetManualShells(
	std::vector<double> customshells)
{
	m_manualshells = customshells;
	m_setmanualshells = true;
	if (!std::is_sorted(m_manualshells.begin(), m_manualshells.end()))
	{
		std::sort(m_manualshells.begin(), m_manualshells.end());
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Specs_User_RayTracer::SetCustomShells, custom shells were out of order, they have been sorted.");
	}
	// round the highest shell up to the nearest integer: if it is <= top-of-atmosphere and rounds up to the nearest 0.001 it will crash
	if (m_manualshells.back() != std::ceil(m_manualshells.back()))
	{
		m_manualshells.back() = std::ceil(m_manualshells.back());
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Specs_User_RayTracer::SetManualShells, the top shell has been rounded up to the nearest interer.");
	}
	return true;
}
