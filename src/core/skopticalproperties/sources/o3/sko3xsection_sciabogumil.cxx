#include <skopticalproperties21.h>

#include "sciabogumil/sciabogumil_203k.hpp"
#include "sciabogumil/sciabogumil_223k.hpp"
#include "sciabogumil/sciabogumil_243k.hpp"
#include "sciabogumil/sciabogumil_273k.hpp"
#include "sciabogumil/sciabogumil_293k.hpp"

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_O3BassPaur_OSIRISRes::skOpticalProperties_O3GomeBurrows_OSIRISRes		2003-11-28
 *-------------------------------------------------------------------------*/

skOpticalProperties_O3_SciaBogumilV3::skOpticalProperties_O3_SciaBogumilV3()
{
	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.

	SetInstrumentPointSpacing( 0.1 );		// nominal spacing of points in nm, its about right
	SetVacuumWavelengths( true );

	AddPSFEntry(  220.0, 0.32 );
	AddPSFEntry(  310.8, 0.21 );
	AddPSFEntry(  402.0, 0.52 );
	AddPSFEntry(  597.8, 0.47 );
	AddPSFEntry(  781.6, 0.62 );
	AddPSFEntry( 1052.5, 1.45 );


	Set_Temperature(243.0);
	AddEntry( 203.0, &sciabogumil_203[0].nm, sizeof(sciabogumil_203[0]), &sciabogumil_203[0].xsect, sizeof(sciabogumil_203[0]), N_ELEMENTS(sciabogumil_203) );
	AddEntry( 223.0, &sciabogumil_223[0].nm, sizeof(sciabogumil_223[0]), &sciabogumil_223[0].xsect, sizeof(sciabogumil_223[0]), N_ELEMENTS(sciabogumil_223) );
	AddEntry( 243.0, &sciabogumil_243[0].nm, sizeof(sciabogumil_243[0]), &sciabogumil_243[0].xsect, sizeof(sciabogumil_243[0]), N_ELEMENTS(sciabogumil_243) );
	AddEntry( 273.0, &sciabogumil_273[0].nm, sizeof(sciabogumil_273[0]), &sciabogumil_273[0].xsect, sizeof(sciabogumil_273[0]), N_ELEMENTS(sciabogumil_273) );
	AddEntry( 293.0, &sciabogumil_293[0].nm, sizeof(sciabogumil_293[0]), &sciabogumil_293[0].xsect, sizeof(sciabogumil_293[0]), N_ELEMENTS(sciabogumil_293) );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_SciaBogumilV3::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_SciaBogumilV3::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_SciaBogumilV3*	clone;
	bool				ok;

	clone = new skOpticalProperties_O3_SciaBogumilV3;
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

