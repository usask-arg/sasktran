#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_LinearCombination::ISKBrdf_Stub_LinearCombination		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_LinearCombination::ISKBrdf_Stub_LinearCombination(SKTRAN_BRDF_LinearCombination* brdf)
	: ISKBrdf_Stub_Base(brdf)
{
	m_brdf = brdf;
}


/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_LinearCombination::~ISKBrdf_Stub_LinearCombination		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_LinearCombination::~ISKBrdf_Stub_LinearCombination()
{}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_LinearCombination::SetPropertyArray		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_LinearCombination::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	nxString					name(propertyname);
	bool						ok = m_brdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_LinearCombination::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "KernelWeights")
		{
			int numkernels = m_brdf->NumKernels();
			ok = ok && (numpoints == numkernels);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_LinearCombination::SetProperty(BRDFParameters), BRDF has %d kernels, you passed in %d parameters", (int)numkernels, (int)numpoints);
			}
			else
			{
				ok = m_brdf->SetKernelWeights(value, numpoints);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_LinearCombination::SetProperty(KernelWeights), failed to set kernel weights");
				}
			}
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyArray(propertyname, value, numpoints);
		}
	}
	return ok;

}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_LinearCombination::SetPropertyObject		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/
bool ISKBrdf_Stub_LinearCombination::SetPropertyObject(const char* propertyname, nxUnknown* object)
{
	nxString					name(propertyname);
	bool						ok = m_brdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_LinearCombination::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "AddKernel")
		{
			skBRDF* brdf_stub_ptr = dynamic_cast<skBRDF*>(object);
			ok = ok && brdf_stub_ptr;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_LinearCombiation::SetPropertyObject, the added object must be a ISKBrdf_Stub_Base object");
			}
			else
			{
				ok = ok && m_brdf->AddKernel(brdf_stub_ptr);
			}
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyObject(propertyname, object);
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
*					ISKBrdf_Stub_LinearCombination::SetPropertyScalar		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/
bool ISKBrdf_Stub_LinearCombination::SetPropertyScalar(const char* propertyname, double value)
{
	nxString					name(propertyname);
	bool						ok = m_brdf != nullptr;

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "ISKBrdf_Stub_LinearCombination::SetProperty, the internal C++ object is undefined. Thats not good");
	}
	else
	{
		if (name == "RemoveKernel")
		{
			m_brdf->RemoveKernel((int)value);
		}
		else
		{
			ok = ISKBrdf_Stub_Base::SetPropertyScalar(propertyname, value);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombination::SKTRAN_BRDF_LinearCombination		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_BRDF_LinearCombination::SKTRAN_BRDF_LinearCombination()
{}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombination::NumKernels		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

int SKTRAN_BRDF_LinearCombination::NumKernels()
{
	return (int)m_kernels.size();
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombination::AddKernel		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_LinearCombination::AddKernel(skBRDF* kernel)
{
	m_kernels.emplace_back(kernel);
	return true;
}


/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombination::SetKernelWeights		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_LinearCombination::SetKernelWeights(const double* f, int length)
{
	m_f.assign(f, f + length);
	return true;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombination::RemoveKernel		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_LinearCombination::RemoveKernel(int index)
{
	int nkernels = NumKernels();
	bool ok = (nkernels > 0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LinearCombination::RemoveKernel, there are no kernels to remove");
	}
	else
	{
		ok = ok && (index >= 0) && (index < nkernels);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LinearCombination::RemoveKernel, there are currently %d kernels and the given index is %d. The index must be between 0 and %d", nkernels, index, nkernels - 1);
		}
		else
		{
			m_kernels.erase(m_kernels.begin() + index);
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_BRDF_LinearCombinationBase::BRDF		 2017-07-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_BRDF_LinearCombinationBase::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf) const
{
	bool ok = (m_f.size() == m_kernels.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LinearCombniation::BRDF, the number of kernels (%d) is not equal to the number of weights (%d). Cannot calculate BRDF until you define values for all kernels", (int)m_kernels.size(), (int)m_f.size());
		*brdf = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		for (int i = 0; i < m_f.size(); i++)
			ok = ok && NXFINITE(m_f[i]);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LinearCombination::BRDF, one or more of the linear kernel weights is NaN. Cannot calculate BRDF until you define values for all kernels");
			*brdf = std::numeric_limits<double>::quiet_NaN();
		}
		else
		{
			// force appropriate cosine values
			CheckCosines(MU_in, MU_out, COSDPHI, "SKTRAN_BRDF_LinearCombinationBase::BRDF");

			// loop through kernels
			double kernel_value;
			double total_brdf = 0.0;
			for (int i = 0; i < m_kernels.size(); i++)
			{
				ok = ok && m_kernels[i]->BRDF(wavelennm, pt, MU_in, MU_out, COSDPHI, &kernel_value);
				total_brdf += m_f[i] * kernel_value;
			}

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_BRDF_LinearCombinaion::BRDF, one or more of the kernels failed");
				*brdf = std::numeric_limits<double>::quiet_NaN();
			}
			else
			{
				*brdf = total_brdf;
			}
		}
	}
	return ok;
}



