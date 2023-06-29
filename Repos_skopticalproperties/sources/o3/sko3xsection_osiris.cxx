#include <skopticalproperties21.h>

#include "osiris_v507/o3_203_osiris.hpp"
#include "osiris_v507/o3_223_osiris.hpp"
#include "osiris_v507/o3_243_osiris.hpp"
#include "osiris_v507/o3_273_osiris.hpp"
#include "osiris_v507/o3_293_osiris.hpp"


/*---------------------------------------------------------------------------
 *'					skOpticalProperties_O3_OSIRISRes::skOpticalProperties_O3_OSIRISRes		2003-11-28
 *-------------------------------------------------------------------------*/

skOpticalProperties_O3_OSIRISRes::skOpticalProperties_O3_OSIRISRes()
{
	SetQuietWavelengthTruncation( true);
	Set_Temperature(241.0);
	AddEntry( 203.0, &o3at203[0].nm, sizeof(o3at203[0]), &o3at203[0].xsect, sizeof(o3at203[0]), N_ELEMENTS(o3at203) );
	AddEntry( 223.0, &o3at223[0].nm, sizeof(o3at223[0]), &o3at223[0].xsect, sizeof(o3at223[0]), N_ELEMENTS(o3at223) );
	AddEntry( 243.0, &o3at243[0].nm, sizeof(o3at243[0]), &o3at243[0].xsect, sizeof(o3at243[0]), N_ELEMENTS(o3at243) );
	AddEntry( 273.0, &o3at273[0].nm, sizeof(o3at273[0]), &o3at273[0].xsect, sizeof(o3at273[0]), N_ELEMENTS(o3at273) );
	AddEntry( 293.0, &o3at293[0].nm, sizeof(o3at293[0]), &o3at293[0].xsect, sizeof(o3at293[0]), N_ELEMENTS(o3at293) );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_OSIRISRes::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_OSIRISRes*	clone;
	bool				ok;

	clone = new skOpticalProperties_O3_OSIRISRes;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}

*/
