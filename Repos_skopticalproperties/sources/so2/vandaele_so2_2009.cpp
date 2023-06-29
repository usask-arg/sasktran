#include <skopticalproperties21.h>
#include "so2_vandaele_2009_xsection.h"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_NO2_Vandaele1998::skOpticalProperties_NO2_Vandaele1998		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_SO2_Vandaele2009::skOpticalProperties_SO2_Vandaele2009()
{
	SetQuietWavelengthTruncation( true);			// quitely set cross-sections outside the wavelength range to zero.
	SetVacuumWavelengths(false);					// we are using air measurements
	SetInstrumentPSF_FWHMWavenumber(2.0);			// Vandaele instrument has constant resolution of 2 cm-1
	SetInstrumentPointSpacingWavenumber( 0.5);		// Vandaele instrument has constant spacing of 0.5 wavenumber in air
	ConfigureEntries();								// Add the cross-section data
	Set_Temperature(298.0);							// set the default temperature
}


/*-----------------------------------------------------------------------------
 *					AddEntry		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_SO2_Vandaele2009::ConfigureEntries()
{

	AddEntry( 298.0, &g_vandaele_so2_2009_298K[0], (int)sizeof(double)*2, &g_vandaele_so2_2009_298K[1], (int)sizeof(double)*2, (int)(N_ELEMENTS(g_vandaele_so2_2009_298K)/2));
	AddEntry( 318.0, &g_vandaele_so2_2009_318K[0], (int)sizeof(double)*2, &g_vandaele_so2_2009_318K[1], (int)sizeof(double)*2, (int)(N_ELEMENTS(g_vandaele_so2_2009_318K)/2));
	AddEntry( 338.0, &g_vandaele_so2_2009_338K[0], (int)sizeof(double)*2, &g_vandaele_so2_2009_338K[1], (int)sizeof(double)*2, (int)(N_ELEMENTS(g_vandaele_so2_2009_338K)/2));
	AddEntry( 358.0, &g_vandaele_so2_2009_358K[0], (int)sizeof(double)*2, &g_vandaele_so2_2009_358K[1], (int)sizeof(double)*2, (int)(N_ELEMENTS(g_vandaele_so2_2009_358K)/2));
	return true;
}



#include "so2_freeman_1984_xsection.h"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SO2_Freeman1984::skOpticalProperties_SO2_Freeman1984		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_SO2_Freeman1984::skOpticalProperties_SO2_Freeman1984()
{
	SetQuietWavelengthTruncation( true);			// quitely set cross-sections outside the wavelength range to zero.
	SetVacuumWavelengths(false);					// we are using air measurements
	SetInstrumentPSF_FWHM(0.002);
	SetInstrumentPointSpacing(0.0004);
	ConfigureEntries();								// Add the cross-section data
	Set_Temperature(298.0);							// set the default temperature
}


/*-----------------------------------------------------------------------------
 *					AddEntry		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_SO2_Freeman1984::ConfigureEntries()
{

	AddEntry( 213.0, &g_freeman_so2_1984_213K[0], (int)sizeof(double)*2, &g_freeman_so2_1984_213K[1], (int)sizeof(double)*2, (int)(N_ELEMENTS(g_freeman_so2_1984_213K)/2));
	return true;
}

#include "so2_rufus_2003_xsection.h"

/*---------------------------------------------------------------------------
 * skOpticalProperties_SO2_Rufus2003::skOpticalProperties_SO2_Rufus20032019-05-14 */
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_SO2_Rufus2003::skOpticalProperties_SO2_Rufus2003()
{
	SetQuietWavelengthTruncation(true);			// quitely set cross-sections outside the wavelength range to zero.
	SetVacuumWavelengths(false);					// we are using air measurements
	SetInstrumentPSF_FWHM(0.0005);
	SetInstrumentPointSpacing(0.0006);
	ConfigureEntries();								// Add the cross-section data
	Set_Temperature(295.0);							// set the default temperature
}


/*-----------------------------------------------------------------------------
 *					AddEntry		2009-11-6*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool skOpticalProperties_SO2_Rufus2003::ConfigureEntries()
{

	AddEntry(295.0, &g_rufus_so2_2003_295K[0], (int)sizeof(double) * 2, &g_rufus_so2_2003_295K[1], (int)sizeof(double) * 2, (int)(N_ELEMENTS(g_rufus_so2_2003_295K) / 2));
	return true;
}


#include "so2_bogumil_2003_xsection.h"

/*---------------------------------------------------------------------------
 * skOpticalProperties_SO2_Bogumil2003::skOpticalProperties_SO2_Bogumil2003 -05-14 */
 /** **/
 /*---------------------------------------------------------------------------*/

skOpticalProperties_SO2_Bogumil2003::skOpticalProperties_SO2_Bogumil2003()
{
	SetQuietWavelengthTruncation(true);			// quitely set cross-sections outside the wavelength range to zero.
	SetVacuumWavelengths(false);				// we are using air measurements
	SetInstrumentPSF_FWHM(0.24);
	SetInstrumentPointSpacing(0.12);
	ConfigureEntries();								// Add the cross-section data
	Set_Temperature(273.0);							// set the default temperature
}


/*-----------------------------------------------------------------------------
 *					AddEntry		2009-11-6*/
 /** **/
 /*---------------------------------------------------------------------------*/

bool skOpticalProperties_SO2_Bogumil2003::ConfigureEntries()
{

	AddEntry(203.0, &g_bogumil_so2_2003_203K[0], (int)sizeof(double) * 2, &g_bogumil_so2_2003_203K[1], (int)sizeof(double) * 2, (int)(N_ELEMENTS(g_bogumil_so2_2003_203K) / 2));
	AddEntry(223.0, &g_bogumil_so2_2003_223K[0], (int)sizeof(double) * 2, &g_bogumil_so2_2003_223K[1], (int)sizeof(double) * 2, (int)(N_ELEMENTS(g_bogumil_so2_2003_223K) / 2));
	AddEntry(243.0, &g_bogumil_so2_2003_243K[0], (int)sizeof(double) * 2, &g_bogumil_so2_2003_243K[1], (int)sizeof(double) * 2, (int)(N_ELEMENTS(g_bogumil_so2_2003_243K) / 2));
	AddEntry(273.0, &g_bogumil_so2_2003_273K[0], (int)sizeof(double) * 2, &g_bogumil_so2_2003_273K[1], (int)sizeof(double) * 2, (int)(N_ELEMENTS(g_bogumil_so2_2003_273K) / 2));
	AddEntry(293.0, &g_bogumil_so2_2003_293K[0], (int)sizeof(double) * 2, &g_bogumil_so2_2003_293K[1], (int)sizeof(double) * 2, (int)(N_ELEMENTS(g_bogumil_so2_2003_293K) / 2));
	return true;
}





