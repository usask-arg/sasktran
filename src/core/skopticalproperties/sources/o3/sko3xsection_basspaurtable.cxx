#include <skopticalproperties21.h>

#include "bass_paur/bp_203clc.hpp"
#include "bass_paur/bp_223clc.hpp"
#include "bass_paur/bp_246clc.hpp"
#include "bass_paur/bp_273clc.hpp"
#include "bass_paur/bp_276clc.hpp"
#include "bass_paur/bp_280clc.hpp"

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaur::skOpticalProperties_O3_BassPaur		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O3_BassPaur::skOpticalProperties_O3_BassPaur()
{

	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.
	SetInstrumentPSF_FWHM    ( 0.025 );		// FWHM point spread function
	SetInstrumentPointSpacing( 0.05  );		// Spacing of data points in nm

	ClearEntries();
	AddEntry( 203.0, &o3bpat203[0].nm, sizeof(o3bpat203[0]), &o3bpat203[0].xsect, sizeof(o3bpat203[0]), N_ELEMENTS(o3bpat203) );
	AddEntry( 223.0, &o3bpat223[0].nm, sizeof(o3bpat223[0]), &o3bpat223[0].xsect, sizeof(o3bpat223[0]), N_ELEMENTS(o3bpat223) );
	AddEntry( 246.0, &o3bpat246[0].nm, sizeof(o3bpat246[0]), &o3bpat246[0].xsect, sizeof(o3bpat246[0]), N_ELEMENTS(o3bpat246) );
	AddEntry( 273.0, &o3bpat273[0].nm, sizeof(o3bpat273[0]), &o3bpat273[0].xsect, sizeof(o3bpat273[0]), N_ELEMENTS(o3bpat273) );
	AddEntry( 276.0, &o3bpat276[0].nm, sizeof(o3bpat276[0]), &o3bpat276[0].xsect, sizeof(o3bpat276[0]), N_ELEMENTS(o3bpat276) );
	AddEntry( 280.0, &o3bpat280[0].nm, sizeof(o3bpat280[0]), &o3bpat280[0].xsect, sizeof(o3bpat280[0]), N_ELEMENTS(o3bpat280) );
	Set_Temperature(241.0);
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_BassPaur::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_BassPaur::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_BassPaur*	clone;
	bool						ok;

	clone = new skOpticalProperties_O3_BassPaur;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_BassPaur::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		clone->skOpticalProperties_UserDefinedAbsorption::DeepCopy(*this);
		clone->skWavelengthToPSF_TableConstant::DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_BassPaur::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}

*/
