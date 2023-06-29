#include <skopticalproperties21.h>
#include <sources/sasktranif_opticalimpl/skbrdf_stubs.h>

SKTRAN_BRDF_UserDefinedLatLon::SKTRAN_BRDF_UserDefinedLatLon()
{
}

SKTRAN_BRDF_UserDefinedLatLon::SKTRAN_BRDF_UserDefinedLatLon(const std::vector<double>& latitudes, const std::vector<double>& longitudes, const nx2dArray<skBRDF*>& brdfs)
{
	Assign(latitudes, longitudes, brdfs);
}

SKTRAN_BRDF_UserDefinedLatLon::~SKTRAN_BRDF_UserDefinedLatLon()
{
	for (size_t idx = 0; idx < m_brdfs.XSize(); idx++)
	{
		for (size_t idy = 0; idy < m_brdfs.YSize(); idy++)
		{
			if (m_brdfs.at(idx, idy) != nullptr)
			{
				m_brdfs.at(idx, idy)->Release();
			}
		}
	}
}

bool SKTRAN_BRDF_UserDefinedLatLon::Assign(const std::vector<double>& latitudes, const std::vector<double>& longitudes, const nx2dArray<skBRDF*>& brdfs)
{
	bool ok = true;

	m_latitudes = latitudes;
	m_longitudes = longitudes;

	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR, "Error copying grids SKTRAN_BRDF_UserDefinedLatLon");
	}

	m_brdfs = brdfs;

	for (size_t idx = 0; idx < m_brdfs.XSize(); idx++)
	{
		for (size_t idy = 0; idy < m_brdfs.YSize(); idy++)
		{
			m_brdfs.at(idx, idy)->AddRef();
		}
	}

	return ok;
}

void SKTRAN_BRDF_UserDefinedLatLon::LatitudeInterpolate(double latitude, std::array<size_t, 2>& indicies, std::array<double, 2>& weights, int& numweights) const
{
	double lowerx, upperx;

	nxLinearInterpolate::FindBoundingIndicesAscending(m_latitudes, latitude, &indicies[0], &indicies[1], &lowerx, &upperx);
	numweights = 2;

	if (latitude > upperx)
	{
		weights[0] = 0;
		weights[1] = 1;
	}
	else if (latitude < lowerx)
	{
		weights[0] = 1;
		weights[1] = 0;
	}
	else
	{
		weights[0] = (upperx - latitude) / (upperx - lowerx);
		weights[1] = 1 - weights[0];
	}
}

void SKTRAN_BRDF_UserDefinedLatLon::LongitudeInterpolate(double longitude, std::array<size_t, 2>& indicies, std::array<double, 2>& weights, int& numweights) const
{
	// Interpolate longitude cyclically

	double lowerx, upperx;
	double modLongitude = fmod(longitude, 360.0);

	nxLinearInterpolate::FindBoundingIndicesAscendingCyclic(m_longitudes, modLongitude, 360.0, &indicies[0], &indicies[1], &lowerx, &upperx);

	lowerx = m_longitudes[indicies[0]];
	upperx = m_longitudes[indicies[1]];

	if (lowerx > upperx)
	{
		std::swap(lowerx, upperx);
		std::swap(indicies[0], indicies[1]);
	}

	while ((upperx < modLongitude) && (lowerx < modLongitude))
	{
		lowerx += 360;
	}
	while ((lowerx > modLongitude) && (upperx > modLongitude))
	{
		upperx -= 360;
	}


	numweights = 2;

	weights[0] = (upperx - modLongitude) / (upperx - lowerx);
	weights[1] = 1 - weights[0];
}

bool SKTRAN_BRDF_UserDefinedLatLon::BRDF(double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double DPHI, double* brdf) const
{
	bool ok = true;

	double tempbrdf;

	double latitude = pt.latitude;
	double longitude = pt.longitude;

	std::array<double, 2> latweights;
	std::array<size_t, 2> latindicies;
	int numlat;

	std::array<double, 2> lonweights;
	std::array<size_t, 2> lonindicies;
	int numlon;

	LatitudeInterpolate(latitude, latindicies, latweights, numlat);
	LongitudeInterpolate(longitude, lonindicies, lonweights, numlon);
	
	*brdf = 0.0;
	for (int latidx = 0; latidx < numlat; latidx++)
	{
		for (int lonidx = 0; lonidx < numlon; lonidx++)
		{
			ok = ok && m_brdfs.At(lonindicies[lonidx], latindicies[latidx])->BRDF(wavelennm, pt, MU_in, MU_out, DPHI, &tempbrdf);
			*brdf += tempbrdf * lonweights[lonidx] * latweights[latidx];
		}
	}
	return ok;
}

ISKBrdf_Stub_UserDefinedLatLon::ISKBrdf_Stub_UserDefinedLatLon(SKTRAN_BRDF_UserDefinedLatLon* brdf) : ISKBrdf_Stub_Base(brdf)
{
	m_brdf = brdf;

	MakeVectorSetFunctions();
	MakeObjectSetFunctions();
	MakeScalarSetFunctions();
}

ISKBrdf_Stub_UserDefinedLatLon::~ISKBrdf_Stub_UserDefinedLatLon()
{
	for (size_t idx = 0; idx < m_brdfs.XSize(); idx++)
	{
		for (size_t idy = 0; idy < m_brdfs.YSize(); idy++)
		{
			if (m_brdfs.at(idx, idy) != nullptr)
			{
				m_brdfs.at(idx, idy)->Release();
			}
		}
	}
}

bool ISKBrdf_Stub_UserDefinedLatLon::SetPropertyScalar(const char* propertyname, double value)
{
	return ISKModuleBase_Stub::SetPropertyScalar(propertyname, value);
}

bool ISKBrdf_Stub_UserDefinedLatLon::SetPropertyArray(const char* propertyname, const double* value, int numpoints)
{
	return ISKModuleBase_Stub::SetPropertyArray(propertyname, value, numpoints);
}

bool ISKBrdf_Stub_UserDefinedLatLon::SetPropertyObject(const char* propertyname, nxUnknown* object)
{
	return ISKModuleBase_Stub::SetPropertyObject(propertyname, object);
}

void ISKBrdf_Stub_UserDefinedLatLon::MakeVectorSetFunctions()
{
	AddSetVectorFunction("latitudes",
		[&, this](const double* values, int n)
		{
			m_latitudes.assign(values, values + n);

			m_brdfs.SetSize(m_longitudes.size(), m_latitudes.size());

			return true;
		}
	);

	AddSetVectorFunction("longitudes",
		[&, this](const double* values, int n)
		{
			m_longitudes.assign(values, values + n);

			m_brdfs.SetSize(m_longitudes.size(), m_latitudes.size());

			return true;
		}
	);
}

void ISKBrdf_Stub_UserDefinedLatLon::MakeObjectSetFunctions()
{
	AddSetObjectFunction("brdf",
		[&, this](nxUnknown* obj)
		{
			skBRDF* brdf = dynamic_cast<skBRDF*>(obj);

			size_t xlen = m_brdfs.XSize();
			
			size_t nx = m_brdfindex % xlen;
			size_t ny = m_brdfindex / xlen;

			m_brdfs.at(nx, ny) = brdf;
			brdf->AddRef();

			return true;
		}
	);
}

void ISKBrdf_Stub_UserDefinedLatLon::MakeScalarSetFunctions()
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
			return m_brdf->Assign(m_latitudes, m_longitudes, m_brdfs);
		}
	);
}