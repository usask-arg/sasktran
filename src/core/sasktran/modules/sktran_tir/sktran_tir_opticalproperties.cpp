/**
 * SASKTRAN TIR Optical Properties Table
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_TableOpticalProperties::SKTRAN_TIR_TableOpticalProperties
 * 2018-09-21
 */
SKTRAN_TIR_TableOpticalProperties::SKTRAN_TIR_TableOpticalProperties()
{
	m_unitsphere = NULL;
	m_alts = NULL;
	m_speedhelper.resize(omp_get_max_threads());
	m_firsttime = true;
	m_forcecacheupdates = false;
	m_calctemperaturederivatives = false;
	m_groundemissivity = 1.0;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::~SKTRAN_TIR_TableOpticalProperties
 * 2018-09-21
 */
SKTRAN_TIR_TableOpticalProperties::~SKTRAN_TIR_TableOpticalProperties()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_TableOpticalProperties::ConfigureGeometry
 * 2018-09-21
 *
 * Initialize the optical properties table from user specifications. Allocates memory for extinction, emission,
 * and inidividual species cross sections. Sets the surface emissivity.
 *
 * @param[in] specs The SKTRAN_TIR_Specs_User object containing user settings
 */
bool SKTRAN_TIR_TableOpticalProperties::ConfigureGeometry(
	const SKTRAN_SpecsUser_Base* specs)
{
	bool ok = true;

	const SKTRAN_TIR_Specs_User* userspecs = dynamic_cast<const SKTRAN_TIR_Specs_User*>(specs);
	ok = ok && nullptr != userspecs;

	size_t numloc = m_unitsphere->NumUnitVectors();
	size_t numalt = m_alts->NumAltitudes();
	size_t numwavel = m_wavelengtharray.size();

	if (0 == numwavel)
	{
		numwavel = 1;
	}

	m_absextinction.resize(numwavel);
	m_emission.resize(numwavel);
	m_groundemission.resize(numwavel);
	for (size_t wavelidx = 0; wavelidx < numwavel; wavelidx++)
	{
		m_absextinction[wavelidx].resize(numloc);
		m_emission[wavelidx].resize(numloc);
		m_groundemission[wavelidx].resize(numloc);
		for (size_t locidx = 0; locidx < numloc; locidx++)
		{
			m_absextinction[wavelidx][locidx].resize(numalt);
			m_emission[wavelidx][locidx].resize(numalt);
			m_groundemission[wavelidx][locidx].resize(1);	// m_groundemission has a third dimension of size=1 so that it can use the same interpolation function as the other tables
		}
	}

	m_airnumberdensity.resize(numloc);
	for (size_t locidx = 0; locidx < numloc; locidx++)
	{
		m_airnumberdensity[locidx].resize(numalt);
	}

	// Store cross sections for WF calculations
	for (size_t speciesidx = 0; speciesidx < userspecs->WeightingFunctionSpecsConst().NumWFSpecies(); speciesidx++)
	{
		m_speciesxs.emplace(userspecs->WeightingFunctionSpecsConst().GetWFSpecies().at(speciesidx), m_absextinction);  // allocates the cross section for this species by copying m_absextinction
		m_speciesnumberdensity_previousrun.emplace(userspecs->WeightingFunctionSpecsConst().GetWFSpecies().at(speciesidx), m_absextinction[0]);  // number density has no wavelength dependence
	}

	m_groundemissivity = userspecs->OpticalPropertiesSpecsConst().GetGroundEmissivity();

	// check if temperature weighting functions are computed
	if (userspecs->WeightingFunctionSpecsConst().GetDoTemperatureWF())
	{
		m_calctemperaturederivatives = true;

		m_dk_dT.resize(numwavel);
		m_dB_dT.resize(numwavel);
		for (size_t wavelidx = 0; wavelidx < numwavel; wavelidx++)
		{
			m_dk_dT[wavelidx].resize(numloc);
			m_dB_dT[wavelidx].resize(numloc);
			for (size_t locidx = 0; locidx < numloc; locidx++)
			{
				m_dk_dT[wavelidx][locidx].resize(numalt);
				m_dB_dT[wavelidx][locidx].resize(numalt);
			}
		}
	}

	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SetWavelengths
 * 2018-09-21
 *
 * Set the wavelengths in nanometers.
 *
 * @param[in] wavelengths should be in order of increasing wavelength
 */
bool SKTRAN_TIR_TableOpticalProperties::SetWavelengths(
	std::vector<double>& wavelengths)
{
	m_wavelengtharray = wavelengths;
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SetAltitudes
 * 2018-09-21
 *
 * Set the altitudes in meters.
 *
 * @param[in] alts
 */
bool SKTRAN_TIR_TableOpticalProperties::SetAltitudes(
	const SKTRAN_GridDefOpticalPropertiesRadii_V21& alts)
{
	bool ok = true;

	alts.AddRef();
	if (m_alts != NULL) m_alts->Release();
	m_alts = &alts;
	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SetUnitSphere
 * 2018-09-21
 *
 * Sets the unitsphere which determines the horizontal points of the optical properties table. The unitsphere will
 * determine the dimensionality of this table. If unitsphere is an SKTRAN_UnitSphere_Dummy object, the table will
 * be one dimensional. If unitsphere is an SKTRAN_UnitSphere_Plane object, the table will be two dimensional.
 *
 * @param[in] unitsphere
 *
 * @post unitsphere reference count is incremented
 */
bool SKTRAN_TIR_TableOpticalProperties::SetUnitSphere(
	const SKTRAN_UnitSphere_V2& unitsphere)
{
	bool ok = true;

	unitsphere.AddRef();
	if (m_unitsphere != NULL) m_unitsphere->Release();
	m_unitsphere = &unitsphere;
	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::ConfigureOptical
 * 2018-09-21
 *
 * Computes extinction and species cross sections at the desired wavelengths, time and locations.
 *
 * @param[in,out] opticalstate
 */
bool SKTRAN_TIR_TableOpticalProperties::ConfigureOptical(
	SKTRAN_TIR_AtmosphericOpticalState& opticalstate)
{
	bool ok = true;

	GEODETIC_INSTANT		geopoint;
	HELIODETIC_VECTOR		heliovector;
	nxVector				geovector;

	if (ok)
	{
		m_mjd = opticalstate.GetTimeAndLocation().mjd;

		geopoint.mjd = m_mjd;
		for (size_t locidx = 0; locidx < m_unitsphere->NumUnitVectors(); locidx++)
		{	// unit sphere is in heliodetic coordinates, so convert to geographic
			heliovector.SetCoords(m_unitsphere->UnitVectorAt(locidx).X(),
									m_unitsphere->UnitVectorAt(locidx).Y(),
									m_unitsphere->UnitVectorAt(locidx).Z());

			geovector = CoordinatesPtr()->HelioVectorToGeographic(heliovector);
			geopoint.latitude = geovector.Latitude();
			geopoint.longitude = geovector.Longitude();
			geopoint.heightm = 0.0;
			ok = ok && opticalstate.SetTimeAndLocation(geopoint, m_forcecacheupdates);
			ok = ok && FillGroundEmissionTableAtIndexMultiWavel(locidx, opticalstate);
			if (m_calctemperaturederivatives)
			{
				ok = ok && opticalstate.ConfigurePerturbedAtmosphere(0.001, m_alts->Altitudes());
			}
			for (size_t altidx = 0; altidx < m_alts->NumAltitudes(); altidx++)
			{
				geopoint.heightm = m_alts->At(altidx);
				ok = ok && opticalstate.SetTimeAndLocation(geopoint, false);
				ok = ok && FillTablesAtIndexMultiWavel(altidx, locidx, opticalstate);
			}
		}
	}
	m_firsttime = false;
	return ok;
}


/**
 * SKTRAN_TIR_TableOpticalProperties::ConfigureOptical
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::ConfigureOpticalFromCache(
	SKTRAN_TIR_AtmosphericOpticalState& opticalstate)
{
	bool ok = true;

	GEODETIC_INSTANT		geopoint;
	HELIODETIC_VECTOR		heliovector;
	nxVector				geovector;

	double					numberdensity_previous;
	double					numberdensity_current;
	CLIMATOLOGY_HANDLE		specieshandle;

	if (ok)
	{
		m_mjd = opticalstate.GetTimeAndLocation().mjd;

		geopoint.mjd = m_mjd;
		if (!m_firsttime)
		{
			for (size_t locidx = 0; locidx < m_unitsphere->NumUnitVectors(); locidx++)
			{	// unit sphere is in heliodetic coordinates, so convert to geographic
				heliovector.SetCoords(m_unitsphere->UnitVectorAt(locidx).X(),
									  m_unitsphere->UnitVectorAt(locidx).Y(),
									  m_unitsphere->UnitVectorAt(locidx).Z());

				geovector = CoordinatesPtr()->HelioVectorToGeographic(heliovector);
				geopoint.latitude = geovector.Latitude();
				geopoint.longitude = geovector.Longitude();
				geopoint.heightm = 0.0;
				ok = ok && opticalstate.SetTimeAndLocation(geopoint, m_forcecacheupdates);
				for (size_t altidx = 0; altidx < m_alts->NumAltitudes(); altidx++)
				{
					geopoint.heightm = m_alts->At(altidx);
					ok = ok && opticalstate.SetTimeAndLocation(geopoint, false);
					ok = ok && opticalstate.UpdateCache();  // needed to update the number density at this time and location
					for (auto& entry : m_speciesnumberdensity_previousrun)
					{
						specieshandle = entry.first;
						numberdensity_previous = entry.second[locidx][altidx];
						ok = ok && opticalstate.GetSpeciesNumberDensity(specieshandle, &numberdensity_current);
						// use change in number density to compute new cross sections
						for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
						{
							m_absextinction[wavelidx][locidx][altidx] += m_speciesxs.at(specieshandle)[wavelidx][locidx][altidx] * (numberdensity_current - numberdensity_previous);
						}
						// set previous numberdensity
						entry.second[locidx][altidx] = numberdensity_current;
					}
				}
			}
		}
		else
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_TableOpticalProperties::ConfigureOpticalFromCache, Can not use cached cross sections on the first radiance calculation because initial cross sections have not been calculated yet.");
		}
	}
	m_firsttime = false;
	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::CalcAltIndices
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::CalcAltIndices(
	const HELIODETIC_POINT& loc,
	double* weights,
	size_t* indices,
	size_t& numindex) const
{
	bool ok = true;

	double alt = loc.Altitude();
	numindex = 2;

	ok = ok && m_alts->FindBoundingIndices(alt, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE,
										   &indices[0], &weights[0], &indices[1], &weights[1]);

	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::CalcSphereIndices
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::CalcSphereIndices(
	const HELIODETIC_POINT& loc,
	double* weights,
	size_t* indices,
	size_t& numindex) const
{
	bool ok = true;
	HELIODETIC_UNITVECTOR locunit = loc.UnitVector();
	nxVector	locvector;

	numindex = m_unitsphere->MaxNumInterpIndex();
	if (numindex == 1)
	{
		weights[0] = 1;
		indices[0] = 0;
		return ok;
	}
	locvector.SetCoords(locunit.X(),
		locunit.Y(),
		locunit.Z());


	ok = ok && m_unitsphere->Triangulate(locvector, indices, weights, 3,
										 m_speedhelper[omp_get_thread_num()]);

	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::InterpTable
 * 2018-09-21
 */
double SKTRAN_TIR_TableOpticalProperties::InterpTable(
	const std::vector<std::vector<double>>& table,
	const HELIODETIC_POINT& loc) const
{
	double result = 0;

	double altweights[2];	// assume only two altitude indices
	double locweights[3];	// assumone three location indicies, for triangular interp
	size_t altindices[2];
	size_t locindices[3];

	size_t numloc;
	size_t numalt;

	CalcAltIndices(loc, altweights, altindices, numalt);
	CalcSphereIndices(loc, locweights, locindices, numloc);

	for (size_t locidx = 0; locidx < numloc; locidx++)
	{
		for (size_t altidx = 0; altidx < numalt; altidx++)
		{
			result += locweights[locidx] * altweights[altidx] * table[locindices[locidx]][altindices[altidx]];
		}
	}
	return result;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::FillTablesAtIndexMultiWavel
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::FillTablesAtIndexMultiWavel(
	size_t altidx,
	size_t locidx,
	SKTRAN_TIR_AtmosphericOpticalState& opticalstate)
{
	bool ok = true;

	std::vector<double> kabs;
	CLIMATOLOGY_HANDLE specieshandle;
	std::map<CLIMATOLOGY_HANDLE, std::vector<double>> speciesxs_temp;
	double speciesnumberdensity_temp;

	std::vector<double> dkabs_temp;	// temperature derivative
	double dB_temp;
	if (m_calctemperaturederivatives)
	{
		dkabs_temp.resize(m_wavelengtharray.size());
	}

	double temperature, pressure, perturbedtemperature;
	double dT = 0;

	for (auto& entry : m_speciesxs)
	{
		specieshandle = entry.first;
		speciesxs_temp.emplace(specieshandle, std::vector<double>(m_wavelengtharray.size()));
	}

	if (m_calctemperaturederivatives)
	{
		ok = ok && opticalstate.CalculateCrossSections(m_wavelengtharray, &kabs, &speciesxs_temp, &dkabs_temp);
	}
	else
	{
		ok = ok && opticalstate.CalculateCrossSections(m_wavelengtharray, &kabs, &speciesxs_temp);
	}

	skClimatology *neutralatmosphere;
	ok = ok && opticalstate.GetAtmosphericStateModel(&neutralatmosphere);
	if (ok)
	{
		ok = ok && neutralatmosphere->GetParameter(SKCLIMATOLOGY_TEMPERATURE_K, opticalstate.GetTimeAndLocation(), &temperature, false);
		ok = ok && neutralatmosphere->GetParameter(SKCLIMATOLOGY_PRESSURE_PA, opticalstate.GetTimeAndLocation(), &pressure, false);
	}
	if (m_calctemperaturederivatives)
	{
		skClimatology *perturbedatmosphere;
		ok = ok && opticalstate.GetPerturbedAtmosphericStateModel(&perturbedatmosphere);
		if (ok)
		{
			ok = ok && perturbedatmosphere->GetParameter(SKCLIMATOLOGY_TEMPERATURE_K, opticalstate.GetTimeAndLocation(), &perturbedtemperature, false);
			if (ok)
			{
				dT = perturbedtemperature - temperature;
			}
		}
	}
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		NXASSERT((kabs[wavelidx] < 10000));
		m_absextinction[wavelidx][locidx][altidx] = kabs[m_wavelengtharray.size() - wavelidx - 1];  // kabs is in order of increasing wavenumber but we want in order of increasing wavelength
		for (auto& entry : m_speciesxs)
		{
			specieshandle = entry.first;
			entry.second[wavelidx][locidx][altidx] = speciesxs_temp.at(specieshandle)[m_wavelengtharray.size() - wavelidx - 1];  // change to increasing wavelength
		}

		// save in order of increasing wavelength
		m_emission[wavelidx][locidx][altidx] = PlanckBlackbody(m_wavelengtharray[wavelidx], temperature);

		if (m_calctemperaturederivatives)
		{
			dB_temp = PlanckBlackbody(m_wavelengtharray[wavelidx], perturbedtemperature) - m_emission[wavelidx][locidx][altidx];
			m_dB_dT[wavelidx][locidx][altidx] = dB_temp / dT;
			m_dk_dT[wavelidx][locidx][altidx] = dkabs_temp[m_wavelengtharray.size() - wavelidx - 1] / dT;
		}
	}

	// Calculate air number density from Ideal Gas Law, n = p / (k T) then convert from m^-3 to cm^-3
	m_airnumberdensity[locidx][altidx] = pressure / (BOLTZMANN * temperature * 1000000.0);

	for (auto& entry : m_speciesnumberdensity_previousrun)
	{
		specieshandle = entry.first;
		ok = ok && opticalstate.GetSpeciesNumberDensity(specieshandle, &speciesnumberdensity_temp);
		entry.second[locidx][altidx] = speciesnumberdensity_temp;
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_TableOpticalProperties::FillTablesAtIndexMultiWavel, Error configuring the Optical State at altidx[%i] locidx[%i]", (int)altidx, (int)locidx);
	}

	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::FillGroundEmissionTableAtIndexMultiWavel
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::FillGroundEmissionTableAtIndexMultiWavel(
	size_t locidx,
	SKTRAN_TIR_AtmosphericOpticalState& opticalstate)
{
	bool ok = true;
	double temperature;

	skClimatology *neutralatmosphere;
	ok = ok && opticalstate.GetAtmosphericStateModel(&neutralatmosphere);
	if (ok)
	{
		ok = ok && neutralatmosphere->GetParameter(SKCLIMATOLOGY_TEMPERATURE_K, opticalstate.GetTimeAndLocation(), &temperature, false);
	}
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		// save in order of increasing wavelength
		m_groundemission[wavelidx][locidx][0] = PlanckBlackbody(m_wavelengtharray[wavelidx], temperature) * m_groundemissivity;
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_TableOpticalProperties::FillGroundEmissionTableAtIndexMultiWavel, Error configuring the Ground Emission at locidx[%i]", (int)locidx);
	}
	
	return ok;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::ReleaseResources
 * 2018-09-21
 */
void SKTRAN_TIR_TableOpticalProperties::ReleaseResources()
{
	if (nullptr != m_unitsphere) m_unitsphere->Release();
	if (nullptr != m_alts) m_alts->Release();
	m_unitsphere = nullptr;
	m_alts = nullptr;

	m_absextinction.clear();
	m_groundemission.clear();
	m_emission.clear();
	m_airnumberdensity.clear();
	m_speciesxs.clear();
}

/**
 * SKTRAN_TIR_TableOpticalProperties::TotalExtinctionPerCM
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::TotalExtinctionPerCM(
	const HELIODETIC_POINT& point,
	std::vector<double>& extinction) const
{
	extinction.resize(m_wavelengtharray.size());
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		extinction[wavelidx] = InterpTable(m_absextinction[wavelidx], point);
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::TotalExtinctionPerCMAtWavel
 * 2018-09-21
 */
double SKTRAN_TIR_TableOpticalProperties::TotalExtinctionPerCMAtWavel(
	const HELIODETIC_POINT& point,
	size_t wavelidx) const
{
	return InterpTable(m_absextinction[wavelidx], point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SpeciesCrossSectionCM2
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::SpeciesCrossSectionCM2(
	const HELIODETIC_POINT& point,
	const CLIMATOLOGY_HANDLE& species,
	std::vector<double>& xs) const
{
	xs.resize(m_wavelengtharray.size());
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		xs[wavelidx] = InterpTable(m_speciesxs.at(species)[wavelidx], point);
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SpeciesCrossSectionCM2AtWavel
 * 2018-09-21
 */
double SKTRAN_TIR_TableOpticalProperties::SpeciesCrossSectionCM2AtWavel(
	const HELIODETIC_POINT& point,
	const CLIMATOLOGY_HANDLE& species,
	size_t wavelidx) const
{
	return InterpTable(m_speciesxs.at(species)[wavelidx], point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::GetEffectiveExtinctionPerCMWithHeight1
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::GetEffectiveExtinctionPerCMWithHeight1(
	const SKTRAN_RayStorage_Base* r,
	size_t startPtIndex,
	double* sigma0,
	double* sigma1,
	size_t wavelidx) const
{
	if (r->ExtinctionAtCellStart(startPtIndex) < 0)
	{
		HELIODETIC_POINT start;
		r->LocationOfPoint(startPtIndex, &start);
		*sigma0 = TotalExtinctionPerCMAtWavel(start, wavelidx);
		r->SetExtinction(startPtIndex, *sigma0);
		NXASSERT(((*sigma0 >= 0) && (*sigma0 < 1.0E5)));
	}
	else
	{
		*sigma0 = r->ExtinctionAtCellStart(startPtIndex);
		NXASSERT(((*sigma0 >= 0) && (*sigma0 < 1.0E5)));
	}
	if (r->ExtinctionAtCellStart(startPtIndex + 1) < 0)
	{
		HELIODETIC_POINT end;
		r->LocationOfPoint(startPtIndex + 1, &end);
		*sigma1 = TotalExtinctionPerCMAtWavel(end, wavelidx);
		r->SetExtinction(startPtIndex + 1, *sigma1);
		NXASSERT(((*sigma1 >= 0) && (*sigma1 < 1.0E5)));
	}
	else
	{
		*sigma1 = r->ExtinctionAtCellStart(startPtIndex + 1);
		NXASSERT(((*sigma1 >= 0) && (*sigma1 < 1.0E5)));
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::GetEffectiveExtinctionPerCMWithHeight1
 * 2018-09-21
 */
bool SKTRAN_TIR_TableOpticalProperties::GetEffectiveExtinctionPerCMWithHeight1(
	const HELIODETIC_POINT& startpoint,
	const HELIODETIC_POINT& endpoint,
	double* sigma0,
	double* sigma1,
	size_t wavelidx) const
{
	*sigma0 = TotalExtinctionPerCMAtWavel(startpoint, wavelidx);
	*sigma1 = TotalExtinctionPerCMAtWavel(endpoint, wavelidx);
	NXASSERT(*sigma0 >= 0.0);
	NXASSERT(*sigma1 >= 0.0);
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::TableSubToInd
 * 2018-09-21
 */
inline SKTRAN_GridIndex SKTRAN_TIR_TableOpticalProperties::TableSubToInd(
	size_t wavidx,
	SKTRAN_GridIndex locidx,
	SKTRAN_GridIndex altidx) const
{
	return (wavidx * m_unitsphere->NumUnitVectors() + locidx) * m_alts->NumAltitudes() + altidx;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SourceTermAtPoint
 * 2019-05-22
 */
bool SKTRAN_TIR_TableOpticalProperties::SourceTermAtPoint(
	const HELIODETIC_POINT& point,
	std::vector<double>& source) const
{
	source.resize(m_wavelengtharray.size());
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		source[wavelidx] = InterpTable(m_emission[wavelidx], point);
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::GroundSourceAtPoint
 * 2019-05-22
 */
bool SKTRAN_TIR_TableOpticalProperties::GroundSourceAtPoint(
	const HELIODETIC_POINT& point,
	std::vector<double>& groundsource) const
{
	groundsource.resize(m_wavelengtharray.size());
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		groundsource[wavelidx] = InterpTable(m_groundemission[wavelidx], point);
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::SourceTermAtPointAndWavel
 * 2019-05-22
 */
double SKTRAN_TIR_TableOpticalProperties::SourceTermAtPointAndWavel(
	const HELIODETIC_POINT& point,
	size_t wavelidx) const
{
	return InterpTable(m_emission[wavelidx], point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::GroundSourceAtPointAndWavel
 * 2019-05-22
 */
double SKTRAN_TIR_TableOpticalProperties::GroundSourceAtPointAndWavel(
	const HELIODETIC_POINT& point,
	size_t wavelidx) const
{
	return InterpTable(m_groundemission[wavelidx], point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::AirNumberDensityAtPoint
 * 2019-07-24
 */
double SKTRAN_TIR_TableOpticalProperties::AirNumberDensityAtPoint(
	const HELIODETIC_POINT& point) const
{
	return InterpTable(m_airnumberdensity, point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::
 * 2020-01-23
 */
bool SKTRAN_TIR_TableOpticalProperties::PlanckFunctionTemperatureDerivativeAtPoint(
	const HELIODETIC_POINT& point,
	std::vector<double>& derivative) const
{
	derivative.resize(m_wavelengtharray.size());
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		derivative[wavelidx] = InterpTable(m_dB_dT[wavelidx], point);
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::
 * 2020-01-23
 */
bool SKTRAN_TIR_TableOpticalProperties::AbsorptionTemperatureDerivativeAtPoint(
	const HELIODETIC_POINT& point,
	std::vector<double>& derivative) const
{
	derivative.resize(m_wavelengtharray.size());
	for (size_t wavelidx = 0; wavelidx < m_wavelengtharray.size(); wavelidx++)
	{
		derivative[wavelidx] = InterpTable(m_dk_dT[wavelidx], point);
	}
	return true;
}

/**
 * SKTRAN_TIR_TableOpticalProperties::
 * 2020-01-23
 */
double SKTRAN_TIR_TableOpticalProperties::PlanckFunctionTemperatureDerivativeAtPointAndWavel(
	const HELIODETIC_POINT& point,
	size_t wavelidx) const
{
	return InterpTable(m_dB_dT[wavelidx], point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::
 * 2020-01-23
 */
double SKTRAN_TIR_TableOpticalProperties::AbsorptionTemperatureDerivativeAtPointAndWavel(
	const HELIODETIC_POINT& point,
	size_t wavelidx) const
{
	return InterpTable(m_dk_dT[wavelidx], point);
}

/**
 * SKTRAN_TIR_TableOpticalProperties::PlanckBlackbody
 * 2019-05-22
 */
double SKTRAN_TIR_TableOpticalProperties::PlanckBlackbody(
	double wavelen_nm,
	double temperature_K)
{
	double       wavelen_m = wavelen_nm / 1.0e9;
	double	     blackbodyradiance;

	// Compute radiance using Planck's Law. Division by 1e13 converts from units of
	// [photons / (m^2 m sr)] to units of [photons / (cm^2 nm sr)].
	blackbodyradiance = (2 * SPEED_OF_LIGHT / pow(wavelen_m, 4)) / (exp(PLANCK * SPEED_OF_LIGHT / wavelen_m / BOLTZMANN / temperature_K) - 1) / 1e13;

	return blackbodyradiance;
}
