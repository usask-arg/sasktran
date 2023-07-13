#include <skopticalproperties21.h>



#include "ftsvoigt/o3_203l.hpp"
#include "ftsvoigt/o3_203h.hpp"
#include "ftsvoigt/o3_223l.hpp"
#include "ftsvoigt/o3_223h.hpp"
#include "ftsvoigt/o3_246l.hpp"
#include "ftsvoigt/o3_246h.hpp"
#include "ftsvoigt/o3_280l.hpp"
#include "ftsvoigt/o3_280h.hpp"
#include "ftsvoigt/o3_293l.hpp"
#include "ftsvoigt/o3_293h.hpp"




/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_FTSVoigt::skOpticalProperties_O3_FTSVoigt		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_O3_FTSVoigt::skOpticalProperties_O3_FTSVoigt()
{
	SetQuietWavelengthTruncation( true);	// quietly set cross-sections outside the wavelength range to zero.
	SetVacuumWavelengths(true);
	SetInstrumentPSF_FWHMWavenumber(5.0);		// Voigt instrument has constant resolution of 5 cm-1
	SetInstrumentPointSpacingWavenumber( 1.929);
	Set_Temperature(241.0);
	m_usinglowpressure = false;					// Set to false to make UseLowPressureEntries work
	UseLow100mBPressureEntries();				// Default to using low pressure entries at the 100 mb level
}


/*-----------------------------------------------------------------------------
 *					AddEntry		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_FTSVoigt::UseLow100mBPressureEntries()
{
	if (!m_usinglowpressure)
	{
		ClearEntries();
		AddEntry( 203.0, &o3_203l[0].nm, sizeof(o3_203l[0]), &o3_203l[0].xsect, sizeof(o3_203l[0]), N_ELEMENTS(o3_203l) );
		AddEntry( 223.0, &o3_223l[0].nm, sizeof(o3_223l[0]), &o3_223l[0].xsect, sizeof(o3_223l[0]), N_ELEMENTS(o3_223l) );
		AddEntry( 246.0, &o3_246l[0].nm, sizeof(o3_246l[0]), &o3_246l[0].xsect, sizeof(o3_246l[0]), N_ELEMENTS(o3_246l) );
		AddEntry( 280.0, &o3_280l[0].nm, sizeof(o3_280l[0]), &o3_280l[0].xsect, sizeof(o3_280l[0]), N_ELEMENTS(o3_280l) );
		AddEntry( 293.0, &o3_293l[0].nm, sizeof(o3_293l[0]), &o3_293l[0].xsect, sizeof(o3_293l[0]), N_ELEMENTS(o3_293l) );
		m_usinglowpressure = true;
	}
	return true;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_FTSVoigt::UseHigh1000mBPressureEntries		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_FTSVoigt::UseHigh1000mBPressureEntries()
{
	if (m_usinglowpressure)
	{
		ClearEntries();
		AddEntry( 203.0, &o3_203h[0].nm, sizeof(o3_203h[0]), &o3_203h[0].xsect, sizeof(o3_203h[0]), N_ELEMENTS(o3_203h) );
		AddEntry( 223.0, &o3_223h[0].nm, sizeof(o3_223h[0]), &o3_223h[0].xsect, sizeof(o3_223h[0]), N_ELEMENTS(o3_223h) );
		AddEntry( 246.0, &o3_246h[0].nm, sizeof(o3_246h[0]), &o3_246h[0].xsect, sizeof(o3_246h[0]), N_ELEMENTS(o3_246h) );
		AddEntry( 280.0, &o3_280h[0].nm, sizeof(o3_280h[0]), &o3_280h[0].xsect, sizeof(o3_280h[0]), N_ELEMENTS(o3_280h) );
		AddEntry( 293.0, &o3_293h[0].nm, sizeof(o3_293h[0]), &o3_293h[0].xsect, sizeof(o3_293h[0]), N_ELEMENTS(o3_293h) );
		m_usinglowpressure = false;
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3_FTSVoigt::DeepCopy		2009-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_O3_FTSVoigt::DeepCopy( const skOpticalProperties_O3_FTSVoigt& other )
{
	bool	ok = true;

	if (other.m_usinglowpressure) ok  = UseLow100mBPressureEntries();
	else                          ok = UseHigh1000mBPressureEntries();

	ok = ok && skOpticalProperties_UserDefinedAbsorption::DeepCopyWithoutTables(other);
	ok = ok && skWavelengthToPSF_TableConstantWavenumber::DeepCopy(other);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_O3::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_O3_FTSVoigt::CreateClone( skOpticalProperties** userclone) const
{
	skOpticalProperties_O3_FTSVoigt*	clone;
	bool						ok;

	clone = new skOpticalProperties_O3_FTSVoigt;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		clone->DeepCopy(*this);
		if(!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_O3::CreateClone, Error copying this object over to clone");
		}
	}
	*userclone = clone;
	return ok;
}
*/

