#include <skopticalproperties21.h>
#include <boost/algorithm/string.hpp>
#include "sources/sasktranif_opticalimpl/skoptprop_stubs.h"
#include "sources/sasktranif_opticalimpl/skemission_stubs.h"


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::ISKSolarSpectrum_Stub_Base		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKSolarSpectrum_Stub_Base::ISKSolarSpectrum_Stub_Base( skSolarSpectrum* solar)
{
	m_solar = solar;
	if (m_solar != nullptr) m_solar->AddRef();
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::~ISKSolarSpectrum_Stub_Base		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKSolarSpectrum_Stub_Base::~ISKSolarSpectrum_Stub_Base( )
{
	if (m_solar != nullptr) m_solar->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::RawObjectPointer		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

//nxUnknown*	ISKSolarSpectrum_Stub_Base::RawObjectPointer()
//{ 
//	return m_solar;
//}

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::Irradiance		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::Irradiance( double wavelen_nm_vacuum, double* irradiance )
{
	*irradiance =  m_solar->Irradiance( wavelen_nm_vacuum);
	return true;
}

bool ISKSolarSpectrum_Stub_Base::IrradianceArray( const double* wavelen_nm_vacuum, double* irradiance, int numpoints)
{
	for (int i = 0; i < numpoints; i++)
	{
		irradiance[i] = m_solar->Irradiance( wavelen_nm_vacuum[i]);
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::IrradianceAt1AU		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::IrradianceAt1AU( double wavelen_nm_vacuum, double* irradiance )
{
	*irradiance = m_solar->IrradianceAt1AU( wavelen_nm_vacuum);
	return true;
}

bool ISKSolarSpectrum_Stub_Base::IrradianceAt1AUArray( const double* wavelen_nm_vacuum, double* irradiance, int numpoints )
{
	for (int i = 0; i < numpoints; i++)
	{
		irradiance[i] = m_solar->IrradianceAt1AU( wavelen_nm_vacuum[i]);
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::SetSolarDistanceFromMjd		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::SetSolarDistanceFromMjd( double mjd )
{
	return m_solar->SetSolarDistanceFromMjd( mjd );

}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::SetSolarDistanceFromAU		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::SetSolarDistanceFromAU( double au )
{
	return m_solar->SetSolarDistanceFromAU(au);
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::MinValidWavelength		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::MinValidWavelength( double* minwavelength_nm) 
{
	*minwavelength_nm = m_solar->MinValidWavelength();
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::MaxValidWavelength		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::MaxValidWavelength( double* maxwavelength_nm) 
{
	*maxwavelength_nm = m_solar->MaxValidWavelength();
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::NanometerResolutionFWHM		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::NanometerResolutionFWHM( double wavelen_nm_vacuum, double* resolution_nm_fwhm)
{
	*resolution_nm_fwhm = m_solar->NanometerResolutionFWHM( wavelen_nm_vacuum );
	return true;
}


/*---------------------------------------------------------------------------
 *    ISKSolarSpectrum_Stub_Base::NanometerResolutionFWHMArray    2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::NanometerResolutionFWHMArray		( const double* wavelen_nm_vacuum, double* resolution_nm_fwhm, int numpoints)
{
	for (int i = 0; i < numpoints; i++)
	{
		resolution_nm_fwhm[i] = m_solar->NanometerResolutionFWHM( wavelen_nm_vacuum[i]);
	}
	return true;

}

/*-----------------------------------------------------------------------------
 *					ISKSolarSpectrum_Stub_Base::NanometerResolutionFWHM		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKSolarSpectrum_Stub_Base::SampleSpacing( double wavelen_nm_vacuum, double* sample_spacing)
{
	*sample_spacing = m_solar->SampleSpacing( wavelen_nm_vacuum);
	return true;
}


/*---------------------------------------------------------------------------
 *         ISKSolarSpectrum_Stub_Base::SampleSpacingArray         2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/
bool ISKSolarSpectrum_Stub_Base::SampleSpacingArray( const double* wavelen_nm_vacuum, double* sample_spacing, int numpoints)
{
	for (int i = 0; i < numpoints; i++)
	{
		sample_spacing[i] = m_solar->SampleSpacing( wavelen_nm_vacuum[i]);
	}
	return true;


}



