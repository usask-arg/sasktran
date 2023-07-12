#include <skopticalproperties21.h>

/* Cross-sections of NO2 */


/*-----------------------------------------------------------------------------
 *					sk_dummy_no2xsectarray		2013-6-13*/
/** \internal
 **/
/*---------------------------------------------------------------------------*/

struct sk_dummy_no2xsectarray
{
	double	nm;
	double	xsect;
};

#include "../o3/osiris_v507/no2_220_osiris.hpp"
#include "../o3/osiris_v507/no2_294_osiris.hpp"

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_NO2_OSIRISRes::skOpticalProperties_NO2_OSIRISRes		2003-11-28
 *-------------------------------------------------------------------------*/

skOpticalProperties_NO2_OSIRISRes::skOpticalProperties_NO2_OSIRISRes()
{
	SetQuietWavelengthTruncation( true);
	NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** skOpticalProperties_NO2_OSIRISRes, The NO2 cross-sections should be updated from 2 temperatures to 4\n"));
	Set_Temperature( 221.0);
	AddEntry( 221.0, &no2at221[0].nm, sizeof(no2at221[0]), &no2at221[0].xsect, sizeof(no2at221[0]), N_ELEMENTS(no2at221) );
	AddEntry( 293.0, &no2at293[0].nm, sizeof(no2at293[0]), &no2at293[0].xsect, sizeof(no2at293[0]), N_ELEMENTS(no2at293) );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_NO2_OSIRISRes::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_NO2_OSIRISRes*	clone;
	bool				ok;

	clone = new skOpticalProperties_NO2_OSIRISRes;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_NO2_OSIRISRes::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_NO2_OSIRISRes::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/
