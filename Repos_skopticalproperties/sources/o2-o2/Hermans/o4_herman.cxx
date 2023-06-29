#include <skopticalproperties21.h>


#include "o4_herman.hpp"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O4_Fally2000::skOpticalProperties_O4_Fally2000		2009-11-6*/
 /** **/
 /*---------------------------------------------------------------------------*/

skOpticalProperties_O4_Fally2000::skOpticalProperties_O4_Fally2000()
{

	SetVacuumWavelengths(false);				// The wavelength tables are for air not Vacuum 
	SetQuietWavelengthTruncation(true);			// quietly set cross-sections outside the wavelength range to zero.

	ClearEntries();

	size_t N = N_ELEMENTS(o4_raw_xsc) / 2;
	m_wavenm.SetSize(N);
	m_o4xsc.SetSize(N);

	size_t index = 2*N - 1;
	for (size_t i = 0; i < N; i++)
	{
		m_o4xsc [i] = o4_raw_xsc[index--];
		m_wavenm[i] = 1.0E7/o4_raw_xsc[index--];
	}
	AddUserEntry( 293.0, m_wavenm, m_o4xsc);
	Set_Temperature(293.0);
	SetInstrumentPSF_FWHMWavenumber(2.0);			// FWHM resolution is 2 cm-1 between 
	SetInstrumentPointSpacingWavenumber(1.0); // Nominal spacing is 1.0 nm

}

