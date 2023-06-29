#include <skopticalproperties21.h>


#include "daumont_brion_malicet_o3/o3dbmat218_new.hpp"
#include "daumont_brion_malicet_o3/o3dbmat228_new.hpp"
#include "daumont_brion_malicet_o3/o3dbmat243_new.hpp"
#include "daumont_brion_malicet_o3/o3dbmat273_new.hpp"
#include "daumont_brion_malicet_o3/o3dbmat295_new.hpp"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_DaumontBrionMalicet::skOpticalProperties_O3_DaumontBrionMalicet		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O3_DaumontBrionMalicet::skOpticalProperties_O3_DaumontBrionMalicet()
{
	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.
	SetInstrumentPSF_FWHM    ( 0.02 );
	SetInstrumentPointSpacing( 0.01  );		// Spacing of data points in nm

	ClearEntries();
	AddEntry( 218.0, &o3dbmat218[0].nm, sizeof(o3dbmat218[0]), &o3dbmat218[0].xsect, sizeof(o3dbmat218[0]), N_ELEMENTS(o3dbmat218) );
	AddEntry( 228.0, &o3dbmat228[0].nm, sizeof(o3dbmat228[0]), &o3dbmat228[0].xsect, sizeof(o3dbmat228[0]), N_ELEMENTS(o3dbmat228) );
	AddEntry( 243.0, &o3dbmat243[0].nm, sizeof(o3dbmat243[0]), &o3dbmat243[0].xsect, sizeof(o3dbmat243[0]), N_ELEMENTS(o3dbmat243) );
	AddEntry( 273.0, &o3dbmat273[0].nm, sizeof(o3dbmat273[0]), &o3dbmat273[0].xsect, sizeof(o3dbmat273[0]), N_ELEMENTS(o3dbmat273) );
	AddEntry( 295.0, &o3dbmat295[0].nm, sizeof(o3dbmat295[0]), &o3dbmat295[0].xsect, sizeof(o3dbmat295[0]), N_ELEMENTS(o3dbmat295) );
//	Set_Temperature(241.0);
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_DaumontBrionMalicet::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_DaumontBrionMalicet*	clone;
	bool						ok;

	clone = new skOpticalProperties_O3_DaumontBrionMalicet;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		clone->skOpticalProperties_UserDefinedAbsorption::DeepCopy(*this);
		clone->skWavelengthToPSF_TableConstant::DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}

*/
