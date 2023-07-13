#include <skopticalproperties21.h>

/*---------------------------------------------------------------------------
 *'					O3 Cross-section data		2003-11-28
 *	These are the tabulated versions of O3 cross-sections as discussed by Burrows et al. 1998.
 *
 * References
 * ------------
 * J. P. Burrows, A. Dehn, B. Deters, S. Himmelmann, A. Richter, S. Voigt, and J. Orphal: 
 * "Atmospheric Remote-Sensing Reference Data from GOME: 2. Temperature-Dependent Absorption Cross Sections of O3 in the 231-794 nm Range", 
 * Journal of Quantitative Spectroscopy and Radiative Transfer 61, 509-517, 1999. 
 *-------------------------------------------------------------------------*/


#include "gomeburrows/gomeburrows_202k.hpp"
#include "gomeburrows/gomeburrows_221k.hpp"
#include "gomeburrows/gomeburrows_241k.hpp"
#include "gomeburrows/gomeburrows_273k.hpp"
#include "gomeburrows/gomeburrows_293k.hpp"

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_O3_HiRes::skOpticalProperties_O3_HiRes		2003-11-28
 *-------------------------------------------------------------------------*/

skOpticalProperties_O3_GomeBurrows::skOpticalProperties_O3_GomeBurrows()
{
	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.
	SetInstrumentPointSpacing		( 0.15 );		// nominal spacing of points in nm, its about right

	AddPSFEntry( 231.0, 0.20 );		//1A.  Map the variation of PSF with wavelength. Note that we cant get 
	AddPSFEntry( 307.0, 0.20 );		//1B.  an exact mapping as there are some overlaps in the spectral
	AddPSFEntry( 311.0, 0.17 );		//2.   of Gome. But it is pretty good.
	AddPSFEntry( 405.0, 0.29 );		//3.
	AddPSFEntry( 600.0, 0.33 );		//4.

	Set_Temperature(241.0);
	AddEntry( 202.0, &o3at202[0].nm, sizeof(o3at202[0]), &o3at202[0].xsect, sizeof(o3at202[0]), N_ELEMENTS(o3at202) );
	AddEntry( 221.0, &o3at221[0].nm, sizeof(o3at221[0]), &o3at221[0].xsect, sizeof(o3at221[0]), N_ELEMENTS(o3at221) );
	AddEntry( 241.0, &o3at241[0].nm, sizeof(o3at241[0]), &o3at241[0].xsect, sizeof(o3at241[0]), N_ELEMENTS(o3at241) );
	AddEntry( 273.0, &o3at273[0].nm, sizeof(o3at273[0]), &o3at273[0].xsect, sizeof(o3at273[0]), N_ELEMENTS(o3at273) );
	AddEntry( 293.0, &o3at293[0].nm, sizeof(o3at293[0]), &o3at293[0].xsect, sizeof(o3at293[0]), N_ELEMENTS(o3at293) );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_HiRes::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_GomeBurrows::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_GomeBurrows*	clone;
	bool				ok;

	clone = new skOpticalProperties_O3_GomeBurrows;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_HiRes::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables(*this);
		ok = ok && clone->skWavelengthToPSF_TableArray::DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3_HiRes::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}

*/

