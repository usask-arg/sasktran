#include <skopticalproperties21.h>


#include "o4_thalman2013_203K_335.749-600.802nm.hpp"
#include "o4_thalman2013_233K_335.646-600.033nm.hpp"
#include "o4_thalman2013_253K_335.147-595.957nm.hpp"
#include "o4_thalman2013_273K_335.696-599.976nm.hpp"
#include "o4_thalman2013_293K_335.749-600.802nm.hpp"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O4_Hermans2000::skOpticalProperties_O4_Hermans2000		2009-11-6*/
 /** **/
 /*---------------------------------------------------------------------------*/

skOpticalProperties_O4_Thalman2013::skOpticalProperties_O4_Thalman2013()
{

	SetVacuumWavelengths(false);				// The wavelength tables are for air not Vacuum 
	SetQuietWavelengthTruncation(true);			// quietly set cross-sections outside the wavelength range to zero.

	ClearEntries();
	AddEntry( 203.0, &thalman2013_o4_203K[0], 2*sizeof(thalman2013_o4_203K[0]), &thalman2013_o4_203K[1], 2*sizeof(thalman2013_o4_203K[0]),  N_ELEMENTS(thalman2013_o4_203K)/2 );
	AddEntry( 233.0, &thalman2013_o4_233K[0], 2*sizeof(thalman2013_o4_233K[0]), &thalman2013_o4_233K[1], 2*sizeof(thalman2013_o4_233K[0]),  N_ELEMENTS(thalman2013_o4_233K)/2 );
	AddEntry( 253.0, &thalman2013_o4_253K[0], 2*sizeof(thalman2013_o4_253K[0]), &thalman2013_o4_253K[1], 2*sizeof(thalman2013_o4_253K[0]),  N_ELEMENTS(thalman2013_o4_253K)/2 );
	AddEntry( 273.0, &thalman2013_o4_273K[0], 2*sizeof(thalman2013_o4_273K[0]), &thalman2013_o4_273K[1], 2*sizeof(thalman2013_o4_273K[0]),  N_ELEMENTS(thalman2013_o4_273K)/2 );
	AddEntry( 293.0, &thalman2013_o4_293K[0], 2*sizeof(thalman2013_o4_293K[0]), &thalman2013_o4_293K[1], 2*sizeof(thalman2013_o4_293K[0]),  N_ELEMENTS(thalman2013_o4_293K)/2 );
	Set_Temperature				(293.0);
	SetInstrumentPSF_FWHM		( 0.4);
	SetInstrumentPointSpacing	( 0.06);

}

