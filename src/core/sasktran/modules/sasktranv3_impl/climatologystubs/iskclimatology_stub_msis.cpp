
#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"



/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_MSIS		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_MSIS::ISKClimatology_Stub_MSIS( skClimatology_MSIS90* msis)
	                                :ISKClimatology_Stub_Base(msis)
{
	m_msis = msis;
	MakeSetFunctions();
}


/*---------------------------------------------------------------------------
 *              ~ISKClimatology_Stub_MSIS              2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_MSIS:: ~ISKClimatology_Stub_MSIS()
{
}


/*---------------------------------------------------------------------------
 *     ISKClimatology_Stub_MSIS::MakeSetFunctions     2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_MSIS::MakeSetFunctions()
{
	AddSetScalarFunction("MaxHeightKMS",
		[&,this](double h)
		{
			m_msis->CachedMSIS().SetMaxHeightKMS(h);
			return true;
		}
	);
	AddSetScalarFunction("HeightSpacingKMS",
		[&,this](double h)
		{
			m_msis->CachedMSIS().SetHeightSpacingKMS(h);
			return true;
		}
	);

	AddSetStringFunction( "AddSpecies",
		[&,this](const char* handlename)
		{
			bool ok;
			CLIMATOLOGY_HANDLE* handleptr;
			handleptr = FindGlobalClimatologyHandle ( handlename);
			ok = (handleptr != nullptr);
			if (ok)
			{
				ok = m_msis->CachedMSIS().AddSpecies( *handleptr );
			}
			return ok;
		}
	);
	AddSetScalarFunction("F10.7",
		[&,this](double f107)
		{
			m_msis->CachedMSIS().SetF10p7(f107);
			return true;
		}
	);
	AddSetScalarFunction("F10.7Avg",
		[&,this](double f107avg)
		{
			m_msis->CachedMSIS().SetF10p7avg(f107avg);
			return true;
		}
	);
	AddSetVectorFunction( "Ap", 
		[&,this](const double* ap, int n)
		{
			return m_msis->CachedMSIS().SetAp( ap, n);
		}
	);

}

