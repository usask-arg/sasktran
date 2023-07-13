#include "include/sktran_montecarlo_internals.h"

SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere()
{
	m_clim = nullptr;
	m_optProp = nullptr;
	m_firsttime = true;
	omp_init_lock(&m_optPropLock);
}

SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::~SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere()
{
	ReleaseResources();
}

void SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::ReleaseResources()
{
	//if (m_baseTable != nullptr) m_baseTable->Release();
	m_baseTable = nullptr;
	m_baseTableMC = nullptr;
	if (m_clim != nullptr) m_clim->Release();
	m_clim = nullptr;
	if (m_optProp != nullptr) m_optProp->Release();
	m_optProp = nullptr;
	omp_destroy_lock(&m_optPropLock);
}


bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::ConfigureGeometry(const SKTRAN_TableOpticalProperties_Base* baseTable)
{
	bool ok = true;

	m_baseTable = baseTable;
	ok = ok && m_baseTable != nullptr;
	//if (ok) m_baseTable->AddRef(); // this causes problems because of cyclic dependencies (if the base table is deleted, this table should be deleted as well anyways)

	m_baseTableMC = dynamic_cast<const SKTRAN_TableOpticalProperties_3D_UnitSphere_MC*>(baseTable);
	ok = ok && m_baseTableMC != nullptr;

	if (ok)
	{
		size_t		numloc = m_baseTableMC->m_unitsphere->NumUnitVectors();
		size_t		numalt = m_baseTableMC->m_alts->NumAltitudes();
		size_t		numscat = m_baseTableMC->m_scatteranglegrid->NumAngles();
		size_t		numwavel = m_baseTableMC->m_wavelengthgrid == NULL ? 1 : m_baseTableMC->m_wavelengthgrid->NumWavelengths();

		if (0 == numwavel)
		{
			numwavel = 1;
		}

		m_inelasticextinction.resize(numwavel);
		for (size_t wavelidx = 0; wavelidx < numwavel; wavelidx++)
		{
			m_inelasticextinction[wavelidx].resize(numloc);
			for (size_t locidx = 0; locidx < numloc; locidx++)
			{
				m_inelasticextinction[wavelidx][locidx].resize(numalt);
			}
		}

		// the inelastic scatter tables are one dimensional since it doesn't depend on wavelength or position
		m_inelasticscatprops->Allocate(numscat);
		m_scattercdf.resize(numscat);
	}

	return ok;
}


bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::ConfigureOptical(double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate)
{
	bool ok = true;

	// hold onto a copy of the inelastic climatology and optical property
	CLIMATOLOGY_HANDLE handle = SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3;
	
	if (m_clim != nullptr) m_clim->Release();
	ok = ok && opticalstate.GetSpeciesClimatology(handle, &m_clim);
	ok = ok && m_clim != nullptr;
	if (ok) m_clim->AddRef();

	skOpticalProperties* optProp;
	ok = ok && opticalstate.GetSpeciesOpticalProperties(handle, &optProp);
	ok = ok && optProp->IsInelasticScatterer();
	if (ok)
	{
		if (m_optProp != nullptr) m_optProp->Release();
		m_optProp = dynamic_cast<skOpticalProperties_Inelastic*>(optProp);
		ok = ok && m_optProp != nullptr;
		if (ok) m_optProp->AddRef();
	}

	// allocate wavelength cdf lookup space
	if (ok)
	{
		size_t numlines = m_optProp->NumInelasticLines();
		m_wnlookup.resize(numlines);
		m_cdflookup.resize(numlines);
	}

	if (ok)
	{
		// fill inelastic extinction table
		m_wavelindex = std::find_if(std::begin(m_baseTableMC->m_wavenumber), std::end(m_baseTableMC->m_wavenumber), [&](double a) {return abs(a - 1e7 / wavelen) < 1e-10;}) - std::begin(m_baseTableMC->m_wavenumber);

		GEODETIC_INSTANT		geopoint;
		HELIODETIC_VECTOR		heliovector;
		nxVector				geovector;

		double cosangle;
		skRTPhaseMatrix pmatrix;

		m_mjd = opticalstate.GetTimeAndLocation().mjd;
		m_wavelen = wavelen;
		opticalstate.SetWavelength(m_wavelen);

		geopoint.mjd = m_mjd;
		for (size_t locidx = 0; locidx < m_baseTableMC->m_unitsphere->NumUnitVectors(); locidx++)
		{	// unit sphere is in heliodetic coordinates, so convert to geographic
			heliovector.SetCoords(m_baseTableMC->m_unitsphere->UnitVectorAt(locidx).X(),
				m_baseTableMC->m_unitsphere->UnitVectorAt(locidx).Y(),
				m_baseTableMC->m_unitsphere->UnitVectorAt(locidx).Z());

			geovector = m_baseTableMC->CoordinatesPtr()->HelioVectorToGeographic(heliovector);
			geopoint.latitude = geovector.Latitude();
			geopoint.longitude = geovector.Longitude();
			geopoint.heightm = 0.0;
			ok = ok && opticalstate.SetTimeAndLocation(geopoint, m_baseTableMC->m_forcecacheupdates);
			for (size_t altidx = 0; altidx < m_baseTableMC->m_alts->NumAltitudes(); altidx++)
			{
				geopoint.heightm = m_baseTableMC->m_alts->At(altidx);
				ok = ok && opticalstate.SetTimeAndLocation(geopoint, false);
				if (0 == m_baseTableMC->m_wavelengtharray.size())
					ok = ok && FillTablesAtIndex(altidx, locidx, opticalstate);
				else
					ok = ok && FillTablesAtIndexMultiWavel(altidx, locidx, opticalstate);
			}
		}
		m_firsttime = false;
	
		// fill the scatter table (1D table since it doesn't depend on position or wavelength)
		for (SKTRAN_GridIndex scatidx = 0; scatidx < m_baseTableMC->m_scatteranglegrid->NumAngles(); scatidx++)
		{
			cosangle = m_baseTableMC->m_scatteranglegrid->At(scatidx);
			ok = ok && GetInelasticScatteringMatrix(opticalstate, m_wavelen, cosangle, &pmatrix);
			m_inelasticscatprops->StorePolarizationPropsCM2(scatidx, pmatrix, opticalstate); // storing the normalized matrix here, no extinction information
		}
	}

	ok = ok && makeScatterCdf(opticalstate);

	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::makeScatterCdf(SKTRAN_AtmosphericOpticalState_V21& opticalstate)
{
	bool ok = true;
	std::vector<double>::iterator it = m_scattercdf.begin();
	double oldVal, newVal;
	skRTPhaseMatrix pmatrix;
	*(it++) = 0.0;

	size_t scatidx = 0;
	newVal = m_inelasticscatprops->PhaseMatrixAccess(scatidx);

	for (scatidx = 1; scatidx < m_baseTableMC->m_scatteranglegrid->NumAngles(); scatidx++)
	{
		oldVal = newVal;
		newVal = m_inelasticscatprops->PhaseMatrixAccess(scatidx);
		*(it++) = *(it - 1) + 0.5 * (oldVal + newVal);
	}

	return ok;
}

double SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::InelasticExtinctionPerCM(const double& wavelength, const HELIODETIC_POINT& point) const
{
	return InterpTable(m_inelasticextinction, wavelength, point);
}


double SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::InterpTable(const std::vector< std::vector< std::vector<double> > >& table, double wavelength, const HELIODETIC_POINT& loc) const
{
	double result = 0;

	double					wlweights[2];			// assume two wavelength indices
	double					altweights[2];			// assume only two altitude indices
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					wlindices[2];
	size_t					altindices[2];
	size_t					locindices[3];
	size_t						numwl;
	size_t						numloc;
	size_t						numalt;

	m_baseTableMC->CalcAltIndices(loc, altweights, altindices, numalt);
	m_baseTableMC->CalcSphereIndices(loc, locweights, locindices, numloc);
	m_baseTableMC->CalcWavelengthIndices(wavelength, wlweights, wlindices, numwl);

	for (size_t wlidx = 0; wlidx < numwl; wlidx++)
	{
		for (size_t locidx = 0; locidx < numloc; locidx++)
		{
			for (size_t altidx = 0; altidx < numalt; altidx++)
			{
				result += wlweights[wlidx] * locweights[locidx] * altweights[altidx] * table[wlindices[wlidx]][locindices[locidx]][altindices[altidx]];
			}
		}
	}
	return result;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetUniquePointWeights(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[24], double gridweights[24], size_t& numNonZero) const
{
	bool ok = true;

	double					wlweights[2];			// assume two wavelength indices
	double					altweights[2];			// assume only two altitude indices
	double					scatweights[2];
	double					locweights[3];			// assumone three location indicies, for triangular interp
	size_t					wlindex[2];
	size_t					altindex[2];
	size_t					locindex[3];
	SKTRAN_GridIndex		scatindex[2];
	size_t	                numwl, numalt, numloc, numscat;
	const double            minval = 0.0;

	ok = ok && m_baseTableMC->CalcWavelengthIndices(wavelength, wlweights, wlindex, numwl);
	ok = ok && m_baseTableMC->CalcAltIndices(point, altweights, altindex, numalt);
	ok = ok && m_baseTableMC->CalcScatterIndices(cosAngle, scatweights, scatindex, numscat);
	ok = ok && m_baseTableMC->CalcSphereIndices(point, locweights, locindex, numloc);
	//if(NULL!=pmatrix) *pmatrix = (lowerAltZeroAngle[loIndex]*(loAltWeight*loWeight) + lowerAltZeroAngle[hiIndex]*(loAltWeight*hiWeight)) + (upperAltZeroAngle[loIndex]*(hiAltWeight*loWeight) + upperAltZeroAngle[hiIndex]*(hiAltWeight*hiWeight));

	numNonZero = 0;
	for (int wlidx = 0; wlidx < numwl; ++wlidx)
	{
		for (int locidx = 0; locidx < numloc; ++locidx)
		{
			for (int altidx = 0; altidx < numalt; ++altidx)
			{
				for (int scatidx = 0; scatidx < numscat; ++scatidx)
				{
					if (minval < locweights[locidx] * altweights[altidx] * scatweights[scatidx])
					{
						gridindices[numNonZero] = TableSubToInd(wlindex[wlidx], locindex[locidx], altindex[altidx], scatindex[scatidx]);
						gridweights[numNonZero] = wlweights[wlidx] * locweights[locidx] * altweights[altidx] * scatweights[scatidx];
						++numNonZero;
					}
				}
			}
		}
	}

	return ok;
}

inline SKTRAN_GridIndex SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::TableSubToInd(size_t wavidx, SKTRAN_GridIndex locidx, SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx) const
{
	return ((wavidx
		*m_baseTableMC->m_unitsphere->NumUnitVectors() + locidx)
		*m_baseTableMC->m_alts->NumAltitudes() + altidx)
		*m_baseTableMC->m_scatteranglegrid->NumAngles() + angidx;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::FillTablesAtIndex(size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate)
{
	bool ok = true;
	
	//double					cosangle;
	skRTPhaseMatrix			pmatrix;

	ok = ok && GetInelasticExtinction(opticalstate, m_wavelen, &m_inelasticextinction[m_wavelindex][locidx][altidx]);

	//for (SKTRAN_GridIndex scatidx = 0; scatidx < m_baseTableMC->m_scatteranglegrid->NumAngles(); scatidx++)
	//{
	//	cosangle = m_baseTableMC->m_scatteranglegrid->At(scatidx);
	//	ok = ok && GetInelasticScatteringMatrix(opticalstate, m_wavelen, cosangle, &pmatrix);
	//	m_inelasticscatprops->StorePolarizationPropsCM2(TableSubToInd(m_wavelindex, locidx, altidx, scatidx), pmatrix, opticalstate);
	//}


	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::FillTablesAtIndex, Error configuring the Optical State at altidx[%i] locidx[%i]", (int)altidx, (int)locidx);
	}


	return ok;
}


bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::FillTablesAtIndexMultiWavel(size_t altidx, size_t locidx, SKTRAN_AtmosphericOpticalState_V21& opticalstate)
{
	bool             ok = true;

	double					wavelength;	

	if (m_firsttime)
	{
		for (size_t wavelidx = 0; wavelidx < m_baseTableMC->m_wavenumber.size(); wavelidx++)
		{
			wavelength = m_baseTableMC->m_wavelengtharray[wavelidx];
			ok = ok && GetInelasticExtinction(opticalstate, wavelength, &m_inelasticextinction[wavelidx][locidx][altidx]);
		}

	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::FillTablesAtIndex, Error configuring the Optical State at altidx[%i] locidx[%i]", (int)altidx, (int)locidx);
	}

	return ok;
}


//inline SKTRAN_GridIndex SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::TableSubToInd(size_t wavidx, SKTRAN_GridIndex locidx, SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx) const
//{
//	return ((wavidx
//		*m_baseTableMC->m_unitsphere->NumUnitVectors() + locidx)
//		*m_baseTableMC->m_alts->NumAltitudes() + altidx)
//		*m_baseTableMC->m_scatteranglegrid->NumAngles() + angidx;
//}


bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetInelasticExtinction(SKTRAN_AtmosphericOpticalState_V21& opticalstate, double wavelength, double* kinel) const
{
	bool ok = true;
	double inelxs;
	double numdensity;

	GEODETIC_INSTANT geopoint = opticalstate.GetTimeAndLocation();
	CLIMATOLOGY_HANDLE handle = SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3;
	
	const size_t numlines = m_optProp->NumInelasticLines();
	ok = ok && m_optProp->SetLocation(geopoint, nullptr);
	ok = ok && m_optProp->CalculateInelasticCrossSections_Reverse(1e7 / wavelength, &inelxs);
	ok = ok && m_clim->GetParameter(handle, geopoint, &numdensity, false);

	*kinel = inelxs * numdensity;

	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetInelasticScatteringMatrix(SKTRAN_AtmosphericOpticalState_V21& opticalstate, double wavelength, double cosscatterangle, skRTPhaseMatrix* pmatrix) const
{
	bool ok = true;

	// location does not have be set, the normalized phase matrix is independent of wavelength and temperature
	ok = ok && m_optProp->CalculateInelasticPhaseMatrix(1e7 / wavelength, cosscatterangle, pmatrix);

	// factor of kscatt / 4pi is applied in SKTRAN_AtmosphericOpticalState_V21::VectorPhaseMatrix()
	// just doing the 4pi for now, the kscatt will be applied later
	*pmatrix *= 0.25 / nxmath::Pi;

	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::SetPolarizationProperties(std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polprops)
{
	m_inelasticscatprops = std::move(polprops);
	return m_inelasticscatprops != nullptr;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetScatteringMatrixCM2(double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix) const
{
	bool ok = true;

	size_t numindex;
	size_t indices[2];
	double weights[2];

	ok = ok && m_baseTableMC->CalcScatterIndices(cosAngle, weights, indices, numindex);
	ok = ok && m_inelasticscatprops->GetScatteringMatrixCM2(indices, weights, numindex, pmatrix);
	pmatrix *= InelasticExtinctionPerCM(wavelength, point); // multiply by the scattering extinction here so we don't have to store scatter matrices at every wavelength and position
	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetScatteringCoefficientCM2(double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const
{
	bool ok = true;

	size_t numindex;
	size_t indices[2];
	double weights[2];

	ok = ok && m_baseTableMC->CalcScatterIndices(cosangle, weights, indices, numindex);
	ok = ok && m_inelasticscatprops->GetPhaseFunctionCM2(indices, weights, numindex, *scatcoeff);
	*scatcoeff *= InelasticExtinctionPerCM(wavelength, point); // multiply by the scattering extinction here so we don't have to store scatter matrices at every wavelength and position

	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetIncomingWavelength(const double & outgoingwavelength, const HELIODETIC_POINT & point, const double & randNum, double & incomingwavelength) const
{
	bool ok = true;

	size_t numlines = m_optProp->NumInelasticLines();
	double wavenumout = 1e7 / outgoingwavelength;
	double wavenumin;

	GEODETIC_INSTANT geopoint = m_baseTableMC->Coordinates()->PointToGeodetic(point, m_mjd);

	omp_set_lock(&m_optPropLock); // lock to set the location for the skOpticalProperty and to reserve the lookup space

	m_optProp->SetLocation(geopoint, nullptr);

	size_t lineidx = 0;
	ok = ok && m_optProp->CalculateInelasticCrossSections_Reverse(wavenumout, lineidx, &m_wnlookup[lineidx], &m_cdflookup[lineidx]);
	for (lineidx = 1; lineidx < numlines; lineidx++)
	{
		ok = ok && m_optProp->CalculateInelasticCrossSections_Reverse(wavenumout, lineidx, &m_wnlookup[lineidx], &m_cdflookup[lineidx]);
		m_cdflookup[lineidx] += m_cdflookup[lineidx - 1];
	}

	lineidx = std::upper_bound(m_cdflookup.begin(), m_cdflookup.end(), randNum * m_cdflookup.back()) - m_cdflookup.begin();
	ok = ok && lineidx < numlines;
	wavenumin = m_wnlookup[lineidx];

	omp_unset_lock(&m_optPropLock);

	// could add linewidth here if we passed another random number

	incomingwavelength = 1e7 / wavenumin;
	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetIncomingWavelength(const std::vector<double>& outgoingwavelength, const size_t & primary, const HELIODETIC_POINT & point, const double & randNum, std::vector<double>& incomingwavelength, std::vector<double>& xsratios) const
{
	// assumes incomingwavelength and wavelengthfactor are sized already
	bool ok = true;

	ok = ok && primary < outgoingwavelength.size();
	if (ok)
	{
		size_t numlines = m_optProp->NumInelasticLines();
		double wavenumout = 1e7 / outgoingwavelength[primary];
		double wavenumin;
		double primaryxs, secondaryxs;

		GEODETIC_INSTANT geopoint = m_baseTableMC->Coordinates()->PointToGeodetic(point, m_mjd);

		omp_set_lock(&m_optPropLock); // lock to set the location for the skOpticalProperty and to reserve the lookup space

		m_optProp->SetLocation(geopoint, nullptr);

		size_t lineidx = 0;
		ok = ok && m_optProp->CalculateInelasticCrossSections_Reverse(wavenumout, lineidx, &m_wnlookup[lineidx], &m_cdflookup[lineidx]);
		for (lineidx = 1; lineidx < numlines; lineidx++)
		{
			ok = ok && m_optProp->CalculateInelasticCrossSections_Reverse(wavenumout, lineidx, &m_wnlookup[lineidx], &m_cdflookup[lineidx]);
			m_cdflookup[lineidx] += m_cdflookup[lineidx - 1];
		}

		lineidx = std::upper_bound(m_cdflookup.begin(), m_cdflookup.end(), randNum * m_cdflookup.back()) - m_cdflookup.begin();
		ok = ok && lineidx < numlines;
		if (ok) wavenumin = m_wnlookup[lineidx];
		xsratios[primary] = 1.0;
		incomingwavelength[primary] = 1e7 / wavenumin;

		primaryxs = lineidx == 0 ? m_cdflookup[lineidx] : m_cdflookup[lineidx] - m_cdflookup[lineidx - 1];
		for (size_t wlidx = 0; wlidx < outgoingwavelength.size(); wlidx++)
		{
			if (wlidx != primary)
			{
				ok = ok && m_optProp->CalculateInelasticCrossSections_Reverse(1e7 / outgoingwavelength[wlidx], lineidx, &incomingwavelength[wlidx], &secondaryxs);
				xsratios[wlidx] = secondaryxs / primaryxs;
				incomingwavelength[wlidx] = 1e7 / incomingwavelength[wlidx];
			}
		}

		omp_unset_lock(&m_optPropLock);

		// could add linewidth here if we passed another random number

	}
	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetCosScatteringAngle(const double& outgoingwavelength, const HELIODETIC_POINT& point, const double& randNum, double &cosScatAngle, SKTRAN_ScatMat_MIMSNC* pmatrix) const
{
	bool ok = true;

	size_t low, high, mid;

	size_t lowIndex, uppIndex;
	double lowerw, upperw;

	double target = randNum * m_scattercdf.back();

	low = 0;
	high = m_scattercdf.size() - 1;

	while (low < high) {
		mid = (low + high) / 2;
		if (target <= m_scattercdf[mid]) {
			high = mid - 1;
		}
		else {
			low = mid + 1;
		}
	}

	lowIndex = high;
	if ((target < m_scattercdf[lowIndex]) && (lowIndex != 0)) lowIndex--;

	if ((m_scattercdf.size() - 1) == lowIndex) {
		// upper end; no room for uppIndex>lowIndex
		uppIndex = lowIndex;
		lowerw = 1.0;
		upperw = 0.0;
	}
	else {
		// there is room in pdf for (uppIndex = lowIndex+1)
		uppIndex = lowIndex + 1;
		lowerw = (m_scattercdf[uppIndex] - target) / (m_scattercdf[uppIndex] - m_scattercdf[lowIndex]);
		upperw = (target - m_scattercdf[lowIndex]) / (m_scattercdf[uppIndex] - m_scattercdf[lowIndex]);
	}

	cosScatAngle = lowerw * m_baseTableMC->m_scatteranglegrid->At(lowIndex) + upperw * m_baseTableMC->m_scatteranglegrid->At(uppIndex);

	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_3D_UnitSphere::GetWavelengthRange(double & minwavelength, double & maxwavelength) const
{
	bool ok = true;
	const std::vector<double> wavelengths = m_baseTableMC->m_wavelengthgrid->Wavelengths();
	ok = ok && wavelengths.size() > 0;
	if (ok)
	{
		minwavelength = *std::min_element(wavelengths.begin(), wavelengths.end());
		maxwavelength = *std::max_element(wavelengths.begin(), wavelengths.end());
	}
	return ok;
}

bool SKTRAN_TableOpticalProperties_Inelastic_DoNothing::GetWavelengthRange(double & minwavelength, double & maxwavelength) const
{
	minwavelength = 0.0;
	maxwavelength = 0.0;
	return true;
}
