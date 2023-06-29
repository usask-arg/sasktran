
#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_LinearCombination::ISKClimatology_Stub_LinearCombination		 2015- 6- 30*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_LinearCombination::ISKClimatology_Stub_LinearCombination( skClimatologyLinearCombination* ptclimate)
	                         :ISKClimatology_Stub_Base(ptclimate)
{
	m_ptclimate = ptclimate;
	MakeSetFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_LinearCombination::~ISKClimatology_Stub_LinearCombination		 2015- 6- 30*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_LinearCombination::~ISKClimatology_Stub_LinearCombination()
{
}
/*---------------------------------------------------------------------------
 *        ISKClimatology_Stub_LinearCombination::MakeSetFunctions         2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_LinearCombination::MakeSetFunctions()
{
	this->AddSetObjectFunction( "SetFirstClimatology", 
		[&, this](nxUnknown* obj)
		{
			bool ok;
			skClimatology*	climate = dynamic_cast<skClimatology*>(obj);
			ok = climate != nullptr;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKClimatology_Stub_LinearCombination::SetFirstClimatology, The object passed in is not an instance of skClimatology");
			}
			else
			{
				ok = m_ptclimate->SetFirstClimatology(climate);
			}
			return ok;
		}
	);
	this->AddSetObjectFunction( "SetSecondClimatology", 
		[&, this](nxUnknown* obj)
		{
			bool ok;
			skClimatology*	climate = dynamic_cast<skClimatology*>(obj);
			ok = climate != nullptr;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKClimatology_Stub_LinearCombination::SetSecondClimatology, The object passed in is not an instance of skClimatology");
			}
			else
			{
				ok = m_ptclimate->SetSecondClimatology(climate);
			}
			return ok;
		}
	);
	this->AddSetVectorFunction( "SetHeightProfileCoeffsOfFirstClimatology",
		[&, this](const double* constdata, int npts)
		{
			bool ok;
			double*				data = (double*)(intptr_t)constdata;
			nx1dArray<double>	h;
			nx1dArray<double>	coeffs;
			size_t				dims    = npts/2;
			size_t				strides = 2*sizeof(double);
			ok =       h.nxArrayLinear::Attach     ( 1, &dims, data,   NULL, &strides	);
			ok = ok && coeffs.nxArrayLinear::Attach( 1, &dims, data+1, NULL, &strides	);
			ok = ok && m_ptclimate->SetHeightProfileCoeffsOfFirstProperty( h, coeffs);
			return ok;
		}
	);
}
