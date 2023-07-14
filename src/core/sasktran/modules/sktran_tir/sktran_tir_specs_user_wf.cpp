/**
 * SASKTRAN TIR User Weighting Function Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_User_wf::SKTRAN_TIR_Specs_User_wf
 * 2019-04-22
 */
SKTRAN_TIR_Specs_User_wf::SKTRAN_TIR_Specs_User_wf()
{
	ConfigureDefaults();
}

/**
 * SKTRAN_TIR_Specs_User_wf::ConfigureDefaults
 * 2019-04-22
 */
void SKTRAN_TIR_Specs_User_wf::ConfigureDefaults()
{
	m_wfspecies.clear();
	m_dotemperaturewf = false;
	m_wfinterpwidth = 1000.0;
	m_maxwfheight = 99500.0;
	m_manualwfresolution = 1000.0;
	m_wfunit = SpeciesWFUnitTIR::numberdensity;
}

/**
 * SKTRAN_TIR_Specs_User_wf::AddWeightingFunctionSpecies
 * 2019-04-22
 */
void SKTRAN_TIR_Specs_User_wf::AddWeightingFunctionSpecies(
	const CLIMATOLOGY_HANDLE& species)
{
	// Make sure that species is not duplicated in the list
	for (size_t i = 0; i < m_wfspecies.size(); i++)
	{
		if (species == m_wfspecies[i])
		{
			return;
		}
	}
	m_wfspecies.push_back(species);
}
