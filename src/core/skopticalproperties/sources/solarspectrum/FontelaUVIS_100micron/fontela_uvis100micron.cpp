#include <skopticalproperties21.h>

#include "fontela_UVIS100micron.h"


/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_FontelaUVIS100Micron::skSolarSpectrum_FontelaUVIS100Micron		 2017- 8- 11*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_FontelaUVIS100Micron::skSolarSpectrum_FontelaUVIS100Micron()
{
	AttachToTable	( &g_fontelaspectrum100[0] , N_ELEMENTS(g_fontelaspectrum100));
	SetInstrumentPSF_FWHM			( 1.0);				// The Fontela spectrum is convolved to 1.0 nm FWHM
	SetInstrumentPointSpacing		( 0.2);				// The Fontela spectrum has a spacing of 0.2 nm between points
}

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_FontelaUVIS100Micron::~skSolarSpectrum_FontelaUVIS100Micron		 2017- 8- 11*/
/** **/
/*---------------------------------------------------------------------------*/

skSolarSpectrum_FontelaUVIS100Micron::~skSolarSpectrum_FontelaUVIS100Micron()
{
}

/*-----------------------------------------------------------------------------
 *					skSolarSpectrum_FontelaUVIS100Micron::NanometerResolutionFWHM		 2017- 8- 11*/
/** **/
/*---------------------------------------------------------------------------*/

double skSolarSpectrum_FontelaUVIS100Micron::NanometerResolutionFWHM(double wavelen_nm_vacuum) const
{
	return GetInstrumentPSF_FWHM(wavelen_nm_vacuum);
}


/*---------------------------------------------------------------------------
 *                   skSolarSpectrum_FontelaUVIS100Micron::SampleSpacing    2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/
double skSolarSpectrum_FontelaUVIS100Micron::SampleSpacing(double wavelen_nm_vacuum) const
{
	return GetInstrumentPointSpacing(wavelen_nm_vacuum);
}

