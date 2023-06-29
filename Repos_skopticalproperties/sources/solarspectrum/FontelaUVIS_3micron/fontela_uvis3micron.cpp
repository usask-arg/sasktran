#include <skopticalproperties21.h>

#include "fontela_UVIS3micron.h"


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_FontelaUVIS3Micron::skSolarSpectrum_FontelaUVIS3Micron		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_FontelaUVIS3Micron::skSolarSpectrum_FontelaUVIS3Micron()
{
	AttachToTable	( &g_fontelaspectrum[0] , N_ELEMENTS(g_fontelaspectrum));
	SetInstrumentPSF_FWHM			( 0.1);				// The Fontela spectrum is convolved to 0.1 nm FWHM
	SetInstrumentPointSpacing		( 0.02);			// The Fontela spectrum has a spacing of 0.01 nm between points
}

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::~skSolarSpectrum_SAO2010		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_FontelaUVIS3Micron::~skSolarSpectrum_FontelaUVIS3Micron()
{
}
		
/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_SAO2010::NanometerResolutionFWHM		 2015- 1- 6*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum_FontelaUVIS3Micron::NanometerResolutionFWHM	(double wavelen_nm_vacuum) const
{
	return GetInstrumentPSF_FWHM(wavelen_nm_vacuum);
}

double skSolarSpectrum_FontelaUVIS3Micron::SampleSpacing(double wavelen_nm_vacuum) const
{
	return GetInstrumentPointSpacing(wavelen_nm_vacuum);
}
