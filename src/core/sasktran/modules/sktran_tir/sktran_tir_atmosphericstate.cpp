#include "include/sktran_tir_internals.h"


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::SKTRAN_TIR_AtmosphericOpticalState 2018-06-11 */
/** Default constructor. Sets MSIS90 as the default atmospheric climatology
 *  which is used to get temperatures and pressures for cross section
 *  calculations.
 **/
/*---------------------------------------------------------------------------*/

SKTRAN_TIR_AtmosphericOpticalState::SKTRAN_TIR_AtmosphericOpticalState()
	                               :SKTRAN_AtmosphericOpticalState_V21() 
{
	m_perturbedatmosphericstate = nullptr;
}



/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::~SKTRAN_TIR_AtmosphericOpticalState 2018-06-11 */
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TIR_AtmosphericOpticalState::~SKTRAN_TIR_AtmosphericOpticalState()
{
	if (m_perturbedatmosphericstate != nullptr) m_perturbedatmosphericstate->Release();
}




/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::CalculateCrossSections	2018-06-11 */
/** Calculates the summed absorption (per cm) from all species at the current
 *  location for an array of wavelengths. Additionally, the absorption cross
 *  sections (in cm^2) are returned for the species specified by speciesxs.
 *
 *  Before calling this function, the input argument 'speciesxs' must
 *  contain entries with keys specifying the CLIMATOLOGY_HANDLES of all species
 *  that the user requires. When the function returns 'speciesxs' will contain
 *  these absorption cross sections. If 'speciesxs' contains any species that
 *  have not been added to this object using AddSpecies (or were removed using
 *  RemoveSpecies), that entry in the returned std::map will be all zeros.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::CalculateCrossSections(const std::vector<double>& wavelength, std::vector<double>* kabs, std::map<CLIMATOLOGY_HANDLE, std::vector<double>>* speciesxs, std::vector<double>* dkabs)
{
	iterator									iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*	entry;
	double										n;
	bool										ok = true;
	bool										ok1;
	bool										ok2;
	std::vector<double>							absxs;
	std::vector<double>							extxs;
	std::vector<double>							scattxs;
	std::vector<double>							pertabsxs;	// absorption from perturbed temperature
	CLIMATOLOGY_HANDLE							species;
	std::vector<double>							wavenumber;

	// convert wavelengths to wavenumber
	wavenumber.resize(wavelength.size());
	for (size_t idx = 0; idx < wavelength.size(); idx++)
	{
		// save in order of increasing wavenumber
		wavenumber[wavenumber.size() - idx - 1] = 1e7 / wavelength[idx];
	}

	kabs->assign(wavenumber.size(), 0.0);
	for (auto& entry : *speciesxs)
	{
		entry.second.assign(wavenumber.size(), 0.0);
	}

	ok = CheckClimatologyCacheIsValid(true);	// Check and execute any pending climatology cache updates, warn about bad time settings
	if (ok)
	{
		for (iter = m_species.begin(); !(iter == m_species.end()); ++iter)																	// now update all o fthe species
		{																																	// iterate over all of the species

			entry = &(*iter);																												// get the entry
			ok1 = entry->UpdateNumberDensityPerCM3(m_placeandtime, m_updateclimatologycache);												// Update the number desnity for the climatology
			n = entry->CurrentNumberDensityPerCM3();																						// only update climatologies
			if (ok1)																														// if twe have a nobn zero number density
			{																																// otherwise wavenumber not zero 
				ok1 = CalculateMultiWaveCrossSections( entry,wavenumber, m_atmosphericstate, m_placeandtime, &absxs, &extxs, &scattxs);		// so calculate cross-sections
				if (ok1)
				{
					species = entry->GetSpecies();
					for (size_t iw = 0; iw < wavenumber.size(); iw++)
					{
						if (speciesxs->count(species) > 0)					// if weighting functions are being computed for this species
						{
							speciesxs->at(species)[iw] = absxs.at(iw);		// save the species cross section
						}
						kabs->at(iw) += absxs.at(iw) * n;
					}
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_AtmosphericOpticalState::CalculateCrossSections, Error calculating cross-sections");
				}
				if (dkabs != nullptr)
				{
					ok2 = CalculateMultiWaveCrossSections(entry, wavenumber, m_perturbedatmosphericstate, m_placeandtime, &pertabsxs, &extxs, &scattxs);
					if (ok2)
					{
						for (size_t iw = 0; iw < wavenumber.size(); iw++)
						{
							dkabs->at(iw) += (pertabsxs.at(iw) - absxs.at(iw)) * n;
						}
					}
				}
			}
			ok = ok && ok1;
		}
	}
	m_updateclimatologycache = m_updateclimatologycache && !ok;
	return ok;
}


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalStateEntry::CalculateMultiWaveCrossSections 2018-06-11 */
/** Calculates the extinction, absorption, and scattering cross sections in
 *  cm^2 for an array of wavelengths for one species. Although only absorption
 *  cross sections are used by the TIR engine, scattering and extinction
 *  cross sections are included here because the skOpticalProperties interface
 *  requires them to be there.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::CalculateMultiWaveCrossSections( SKTRAN_AtmosphericOpticalStateEntry_V21* entry,
																		  const std::vector<double>&			   wavenumber,
																		  skClimatology*						   neutralatmosphere,
																		  const GEODETIC_INSTANT&				   placeandtime,
																		  std::vector<double>*					   absxs,
																		  std::vector<double>*					   extxs,
																		  std::vector<double>*					   scattxs
								                                         )
{
	bool ok;
	bool haschanged;

	NXTRACE_ONCEONLY(firsttime, ("SKTRAN_TIR_AtmosphericOpticalStateEntry::CalculateMultiWaveCrossSections, needs re-working so the climatology update is in one loop and the cross-sections is in another\n"));

	absxs->resize(wavenumber.size());
	extxs->resize(wavenumber.size());
	scattxs->resize(wavenumber.size());

	ok = entry->ParticleOpticalProps()->SetAtmosphericState(neutralatmosphere);
	ok = entry->ParticleOpticalProps()->SetLocation(placeandtime, &haschanged);
	ok = ok && entry->UpdateNumberDensityPerCM3(placeandtime, false);
	ok = ok && entry->ParticleOpticalProps()->CalculateCrossSectionsArray(&wavenumber.front(), (int)wavenumber.size(), &absxs->front(), &extxs->front(), &scattxs->front());

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_AtmosphericOpticalStateEntry::CalculateMultiWaveCrossSections, Error updating cross-section for given place and time");
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::SetNumThreads              2019-12-12 */
/** Checks each entry in the atmospheric state and sets the number of threads
 *  to use in the multithreaded cross section calculation, IF that entry uses
 *  a HITRAN optical property.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::SetNumThreads(size_t numthreads)
{
	bool ok = true;
	iterator iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21* entry;
	skOpticalProperties* optpropPtr;
	skOpticalProperties_HitranChemical* hitranPtr;

	for (iter = m_species.begin(); !(iter == m_species.end()); ++iter)
	{
		entry = &(*iter);
		optpropPtr = entry->ParticleOpticalProps();
		hitranPtr = dynamic_cast<skOpticalProperties_HitranChemical*>(optpropPtr);
		if (hitranPtr)  // if hitranPtr is not a null pointer, this is a HITRAN entry
		{
			hitranPtr->HitranChemical()->SetNumThreads(numthreads);
		}
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::GetSpeciesNumberDensity    2018-06-11 */
/** Get the number density of a species. Does not force a cache update if the
 *  time and location has changed.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::GetSpeciesNumberDensity(const CLIMATOLOGY_HANDLE& speciesinlist, double* numberdensity)
{
	iterator iter;
	bool ok = false;
	SKTRAN_AtmosphericOpticalStateEntry_V21* entry;

	for (iter = m_species.begin(); !(iter == m_species.end()); ++iter)
	{
		entry = &(*iter);
		if (entry->GetSpecies() == speciesinlist)
		{
			ok = true;
			*numberdensity = entry->CurrentNumberDensityPerCM3();
			return ok;
		}
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::UpdateCache                2018-06-11 */
/** Forces an update of the explicit particle climatologies and internal
 *  climatologies inside the optical properties.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::UpdateCache()
{
	iterator										iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21*		entry;
	double											n;
	bool											ok = true;
	bool											ok1;

	ok = CheckClimatologyCacheIsValid(true);													// Check and execute any pending climatology cache updates, warn about bad time settings
	if (ok)
	{
		for (iter = m_species.begin(); !(iter == m_species.end()); ++iter)						// now update all o fthe species
		{																						// iterate over all of the species
			entry = &(*iter);																	// get the entry
			ok1 = entry->UpdateNumberDensityPerCM3(m_placeandtime, m_updateclimatologycache);	// Update the number desnity for the climatology
			n = entry->CurrentNumberDensityPerCM3();											// only update climatologies

			ok = ok && ok1;
		}
	}
	m_updateclimatologycache = m_updateclimatologycache && !ok;
	return ok;
}


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::ContainsSpecies            2019-06-26 */
/** Check if this atmospheric state contains a given species.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::ContainsSpecies(const CLIMATOLOGY_HANDLE& species)
{
	bool ok(true);

	iterator								iter;
	SKTRAN_AtmosphericOpticalStateEntry_V21	dummy(species);

	iter = std::find(m_species.begin(), m_species.end(), dummy);
	if (iter == m_species.end())
	{
		ok = false;
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::ConfigurePerturbedAtmosphere 2020-01-23 */
/** Set up an atmospheric state where the temperature has been perturbed by
 *  a specified amount. E.g. for a factor f, the perturbed temperature is
 *  T_perturbed = T_base + T_base * f
 *  Should only be called after the atmospheric state has been finalized
 *
 *  The TIR engine only supports 1D altitude perturbations so the perturbed
 *  atmosphere uses the 1D custom climatology class.
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::ConfigurePerturbedAtmosphere(
	double perturbationfactor,
	const std::vector<double>& altitudes)
{
	bool ok = true;
	std::vector<double> pressure;
	std::vector<double> temperature;
	size_t numbad = 0;
	nx2dArray<double> profiles(altitudes.size(), 3);
	std::vector<CLIMATOLOGY_HANDLE> handles;

	handles.push_back(SKCLIMATOLOGY_PRESSURE_PA);
	handles.push_back(SKCLIMATOLOGY_TEMPERATURE_K);

	if (m_perturbedatmosphericstate != nullptr) m_perturbedatmosphericstate->Release();
	m_perturbedatmosphericstate = new skClimatology_UserDefinedTable;
	m_perturbedatmosphericstate->AddRef();

	pressure.resize(altitudes.size());
	temperature.resize(altitudes.size());

	ok = ok && m_atmosphericstate->GetHeightProfile(SKCLIMATOLOGY_PRESSURE_PA, GetTimeAndLocation(), &altitudes.front(), altitudes.size(), &pressure.front(), true, &numbad);
	ok = ok && m_atmosphericstate->GetHeightProfile(SKCLIMATOLOGY_TEMPERATURE_K, GetTimeAndLocation(), &altitudes.front(), altitudes.size(), &temperature.front(), true, &numbad);

	for (size_t altidx = 0; altidx < altitudes.size(); altidx++)
	{
		profiles.At(altidx, 0) = altitudes[altidx];
		profiles.At(altidx, 1) = pressure[altidx];									// pressure stays the same
		profiles.At(altidx, 2) = temperature[altidx] * (1.0 + perturbationfactor);	// perturb temperature
	}

	m_perturbedatmosphericstate->LoadProfileFrom2dArray(&handles.front(), 2, profiles);

	return ok;
}



/*-----------------------------------------------------------------------------
 * SKTRAN_TIR_AtmosphericOpticalState::GetPerturbedAtmosphericStateModel 2020-01-23 */
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TIR_AtmosphericOpticalState::GetPerturbedAtmosphericStateModel(skClimatology** statemodel)
{
	bool ok = false;
	if (m_perturbedatmosphericstate != nullptr)
	{
		*statemodel = m_perturbedatmosphericstate;
		ok = true;
	}
	return ok;
}

