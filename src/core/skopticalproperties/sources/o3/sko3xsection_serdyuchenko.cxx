#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					_skxsectarray_serdy		2013-6-13*/
/** \internal
**/
/*---------------------------------------------------------------------------*/

 struct _skxsectarray_serdy
 {
	 double nm;
	 double sigma293;
	 double sigma283;
	 double sigma273;
	 double sigma263;
	 double sigma253;
	 double sigma243;
	 double sigma233;
	 double sigma223;
	 double sigma213;
	 double sigma203;
	 double sigma193;
 };

#include "serdyuchenko/serdyuchenkogorshelev5digits.hpp"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_SerdyuchenkoV1::skOpticalProperties_O3_SerdyuchenkoV1		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O3_SerdyuchenkoV1::skOpticalProperties_O3_SerdyuchenkoV1()
{

	SetVacuumWavelengths( true );			// The wavelength tables are for Vacuum , not air
	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.

	ClearEntries();
	AddEntry( 193.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma193, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 203.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma203, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 213.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma213, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 223.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma223, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 233.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma233, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 243.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma243, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 253.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma253, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 263.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma263, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 273.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma273, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 283.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma283, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	AddEntry( 293.0, &o3skserdy[0].nm,  sizeof(o3skserdy[0]),  &o3skserdy[0].sigma293, sizeof(o3skserdy[0]), N_ELEMENTS(o3skserdy) );
	Set_Temperature(243.0);
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_SerdyuchenkoV1::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_SerdyuchenkoV1*	clone;
	bool						ok;

	clone = new skOpticalProperties_O3_SerdyuchenkoV1;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_SerdyuchenkoV1::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		clone->skOpticalProperties_UserDefinedAbsorption::DeepCopy(*this);
		clone->skWavelengthToPSF_SerdyuchenkoV1::DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_SerdyuchenkoV1::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_SerdyuchenkoV1::skWavelengthToPSF_SerdyuchenkoV1		2012-7-27*/
/**  - Spectral Resolution(HWHM):
 *			- 0.01 nm below 290 nm},
 *			- 1 cm-1 between 290 nm and 350 nm},
 *			- 0.01 nm between 350 nm and 450 nm},
 *			- 1 cm-1 between 450 nm and 1100 nm,
 *			.
**/
/*---------------------------------------------------------------------------*/

skWavelengthToPSF_SerdyuchenkoV1::skWavelengthToPSF_SerdyuchenkoV1()
{
	m_psfwavenum.SetInstrumentPSF_FWHMWavenumber(2.0);			// FWHM resolution is 2 cm-1 between 
}

/*-----------------------------------------------------------------------------
 *					skWavelengthToPSF_SerdyuchenkoV1::GetInstrumentPSF_FWHM		2012-7-27*/
/** **/
/*---------------------------------------------------------------------------*/

double skWavelengthToPSF_SerdyuchenkoV1::GetInstrumentPSF_FWHM( double nm ) const
{
	double fwhm = 0.02;	// if (nm <= 290.0) fwhm = 0.02

	if       ( ( nm >  290.0) && (nm < 350.0) ) fwhm = m_psfwavenum.GetInstrumentPSF_FWHM(nm);
	else if  ( ( nm >= 350.0) && (nm < 450.0) ) fwhm = 0.02;
	else if  ( ( nm >= 450.0) )                 fwhm = m_psfwavenum.GetInstrumentPSF_FWHM(nm);
//  else fwhm = 0.02;
	return fwhm;
}
