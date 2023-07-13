#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					_skxsectarray_sciav4		2013-6-13*/
/** \internal
**/
/*---------------------------------------------------------------------------*/

 struct _skxsectarray_sciav4
 {
	 double nm;
	 double sigma203;
	 double sigma223;
	 double sigma243;
	 double sigma273;
	 double sigma293;
 };

#include "sciabogumil/scia_o3_temp_crosssection_v4.hpp"

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_O3BassPaur_OSIRISRes::skOpticalProperties_O3GomeBurrows_OSIRISRes		2003-11-28
 *-------------------------------------------------------------------------*/

skOpticalProperties_O3_SciaBogumilV4::skOpticalProperties_O3_SciaBogumilV4()
{
	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.

	SetInstrumentPointSpacing( 0.1 );		// nominal spacing of points in nm, its about right
	SetVacuumWavelengths( true );

	//! Spectral Resolution:	0.32 nm below 311.7 nm,
	//							0.21 nm between 311.0 nm and 397.4 nm,
	//							0.52 nm between 397.4 nm and 598.2 nm,
	//							0.47 nm between 598.2 nm and 780.1 nm,
	//!                         0.62 nm between 780.1 nm and 1053.5 nm,
	//							1.45 nm above 1053.5 nm (FWHM) 

	AddPSFEntry(  220.0, 0.32 );
	AddPSFEntry(  311.0, 0.21 );
	AddPSFEntry(  397.4, 0.52 );
	AddPSFEntry(  598.2, 0.47 );
	AddPSFEntry(  780.1, 0.62 );
	AddPSFEntry( 1052.5, 1.45 );


	Set_Temperature(243.0);
	AddEntry( 203.0, &sciabogumilv4[0].nm, sizeof(sciabogumilv4[0]), &sciabogumilv4[0].sigma203, sizeof(sciabogumilv4[0]), N_ELEMENTS(sciabogumilv4) );
	AddEntry( 223.0, &sciabogumilv4[0].nm, sizeof(sciabogumilv4[0]), &sciabogumilv4[0].sigma223, sizeof(sciabogumilv4[0]), N_ELEMENTS(sciabogumilv4) );
	AddEntry( 243.0, &sciabogumilv4[0].nm, sizeof(sciabogumilv4[0]), &sciabogumilv4[0].sigma243, sizeof(sciabogumilv4[0]), N_ELEMENTS(sciabogumilv4) );
	AddEntry( 273.0, &sciabogumilv4[0].nm, sizeof(sciabogumilv4[0]), &sciabogumilv4[0].sigma273, sizeof(sciabogumilv4[0]), N_ELEMENTS(sciabogumilv4) );
	AddEntry( 293.0, &sciabogumilv4[0].nm, sizeof(sciabogumilv4[0]), &sciabogumilv4[0].sigma293, sizeof(sciabogumilv4[0]), N_ELEMENTS(sciabogumilv4) );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_SciaBogumilV4::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_SciaBogumilV4::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_SciaBogumilV4*	clone;
	bool				ok;

	clone = new skOpticalProperties_O3_SciaBogumilV4;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_OSIRISRes_SciaBogumil::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok =       clone->skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables(*this);
		ok = ok && clone->skWavelengthToPSF_TableArray::DeepCopy(*this);

		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_OSIRISRes_SciaBogumil::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/

