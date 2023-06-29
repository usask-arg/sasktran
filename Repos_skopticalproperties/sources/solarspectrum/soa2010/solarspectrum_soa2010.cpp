#include <skopticalproperties21.h>

#include "soa2010.h"


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::skSolarSpectrum_SAO2010		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_SAO2010::skSolarSpectrum_SAO2010()
{
	AttachToTable	( &g_soa2010spectrum[0] , N_ELEMENTS(g_soa2010spectrum));
	SetInstrumentPSF_FWHM			( 0.04);			// The SAO2010 solar spectrum is specified as a resolution of 0.04 nm FWHM
	SetInstrumentPointSpacing		( 0.01);			// The SAO2010 solar spectrum has a spacing of 0.01 nm between points
}


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::~skSolarSpectrum_SAO2010		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_SAO2010::~skSolarSpectrum_SAO2010()
{
}


		
/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::NanometerResolutionFWHM		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum_SAO2010::NanometerResolutionFWHM	(double wavelen_nm_vacuum) const
{
	return GetInstrumentPSF_FWHM(wavelen_nm_vacuum);
}

double skSolarSpectrum_SAO2010::SampleSpacing(double wavelen_nm_vacuum) const
{
	return GetInstrumentPointSpacing(wavelen_nm_vacuum);
}

