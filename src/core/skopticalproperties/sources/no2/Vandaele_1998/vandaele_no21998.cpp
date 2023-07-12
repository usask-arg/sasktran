#include <skopticalproperties21.h>
#include "vandaelele_xsection.h"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_NO2_Vandaele1998::skOpticalProperties_NO2_Vandaele1998		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_NO2_Vandaele1998::skOpticalProperties_NO2_Vandaele1998()
{
	SetQuietWavelengthTruncation( true);			// quitely set cross-sections outside the wavelength range to zero.
	SetVacuumWavelengths(false);					// we are using air measurements
	SetInstrumentPSF_FWHMWavenumber(2.0);			// Vandaele instrument has constant resolution of 2 cm-1
	SetInstrumentPointSpacingWavenumber( 0.96457);	// Vandaele instrument has constant spacing of 0.96457 wavenumber in air
	ConfigureEntries();								// Add the cross-section data
	Set_Temperature(220.0);							// set the default temperature
}


/*-----------------------------------------------------------------------------
 *					AddEntry		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_NO2_Vandaele1998::ConfigureEntries()
{
	size_t	numpoints = N_ELEMENTS(g_vandaele97_no2)/3;

	AddEntry( 220.0, &g_vandaele97_no2[0], (int)sizeof(double)*3, &g_vandaele97_no2[1], (int)sizeof(double)*3, (int)numpoints );
	AddEntry( 294.0, &g_vandaele97_no2[0], (int)sizeof(double)*3, &g_vandaele97_no2[2], (int)sizeof(double)*3, (int)numpoints );
	return true;
}



