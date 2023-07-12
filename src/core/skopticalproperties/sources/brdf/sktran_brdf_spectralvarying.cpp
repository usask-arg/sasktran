#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>


SKTRAN_BRDF_SpectralVarying::SKTRAN_BRDF_SpectralVarying(const std::vector<double>& wavelengths,  const nx1dArray<skBRDF*>& brdfs)
{
	Assign(wavelengths, brdfs);
}

SKTRAN_BRDF_SpectralVarying::~SKTRAN_BRDF_SpectralVarying()
{
	for (size_t idx = 0; idx < m_brdfs.size(); idx++)
	{
		m_brdfs.at(idx)->Release();
	}
}

bool SKTRAN_BRDF_SpectralVarying::Assign(const std::vector<double>& wavelengths, const nx1dArray<skBRDF*>& brdfs)
{
	bool ok = true;

	m_wavelengths = wavelengths;

	m_brdfs = brdfs;

	for (size_t idx = 0; idx < m_brdfs.size(); idx++)
	{
		m_brdfs.at(idx)->AddRef();
	}

	return ok;
}

void SKTRAN_BRDF_SpectralVarying::Interpolate(double wavelength, std::array<size_t, 2>& indicies, std::array<double, 2>& weights, int& numweights) const
{
	double lowerx, upperx;

	nxLinearInterpolate::FindBoundingIndicesAscending(m_wavelengths, wavelength, &indicies[0], &indicies[1], &lowerx, &upperx);
	numweights = 2;

	if (wavelength > upperx)
	{
		weights[0] = 0;
		weights[1] = 1;
	}
	else if (wavelength < lowerx)
	{
		weights[0] = 1;
		weights[1] = 0;
	}
	else
	{
		weights[0] = (upperx - wavelength) / (upperx - lowerx);
		weights[1] = 1 - weights[0];
	}
}

bool SKTRAN_BRDF_SpectralVarying::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const
{
	bool ok = true;

	double tempbrdf;

	std::array<double, 2> wavelweights;
	std::array<size_t, 2> wavelindicies;
	int numwavel;

	Interpolate(wavelennm, wavelindicies, wavelweights, numwavel);

	*brdf = 0.0;
	for (int widx = 0; widx < numwavel; widx++)
	{
		ok = ok && m_brdfs.At(wavelindicies[widx])->BRDF(wavelennm, pt, MU_in, MU_out, DPHI, &tempbrdf);
		*brdf += tempbrdf * wavelweights[widx];
	}
	return ok;
}

ISKBrdf_Stub_SpectralVarying::ISKBrdf_Stub_SpectralVarying(SKTRAN_BRDF_SpectralVarying* brdf) : ISKBrdf_Stub_Base(brdf)
{
	m_brdf = brdf;

	MakeVectorSetFunctions();
	MakeObjectSetFunctions();
	MakeScalarSetFunctions();
}

ISKBrdf_Stub_SpectralVarying::~ISKBrdf_Stub_SpectralVarying()
{
	for (size_t idx = 0; idx < m_brdfs.size(); idx++)
	{
		m_brdfs.at(idx)->Release();
	}
}

bool ISKBrdf_Stub_SpectralVarying::SetPropertyScalar(const char* propertyname, double value)
{
	return ISKModuleBase_Stub::SetPropertyScalar(propertyname, value);
}

bool ISKBrdf_Stub_SpectralVarying::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	return ISKModuleBase_Stub::SetPropertyArray(propertyname, value, numpoints);
}

bool ISKBrdf_Stub_SpectralVarying::SetPropertyObject(const char* propertyname, nxUnknown* object)
{
	return ISKModuleBase_Stub::SetPropertyObject(propertyname, object);
}

void ISKBrdf_Stub_SpectralVarying::MakeVectorSetFunctions()
{
	AddSetVectorFunction("wavelengths",
		[&, this](const double* values, int n)
	{
		m_wavelengths.assign(values, values + n);

		m_brdfs.SetSize(m_wavelengths.size());

		return true;
	}
	);
}

void ISKBrdf_Stub_SpectralVarying::MakeObjectSetFunctions()
{
	AddSetObjectFunction("brdf",
		[&, this](nxUnknown* obj)
	{
		skBRDF* brdf = dynamic_cast<skBRDF*>(obj);

		m_brdfs.at(m_brdfindex) = brdf;
		brdf->AddRef();

		return true;
	}
	);
}

void ISKBrdf_Stub_SpectralVarying::MakeScalarSetFunctions()
{
	AddSetScalarFunction("brdfindex",
		[&, this](double d)
	{
		bool ok = true;
		int specifier = (int)ceil(d - 0.5);

		m_brdfindex = specifier;

		return ok;
	}
	);

	AddSetScalarFunction("initialize",
		[&, this](double d)
	{
		return m_brdf->Assign(m_wavelengths, m_brdfs);
	}
	);
}