/**
 * SASKTRAN TIR User Optical Properties Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_User_OpticalPropertiesTable::ConfigureDefaults
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_User_OpticalPropertiesTable::ConfigureDefaults()
{
	m_tabledim = AtmosphereDimensionTIR::dim1;
	m_numprofiles = 30;
	m_heightres = 500;

	// Invalid values mean use normal and reference from the line of sight plane
	m_normal.SetInvalid();
	m_reference.SetInvalid();

	m_forcecacheupdates = false;

	m_groundemissivity = 1.0;

	return true;
}
