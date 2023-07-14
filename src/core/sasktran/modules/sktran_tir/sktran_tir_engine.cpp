/**
 * SASKTRAN TIR Engine
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Engine::SKTRAN_TIR_Engine
 * 2018-09-13
 */
SKTRAN_TIR_Engine::SKTRAN_TIR_Engine()
{
	m_coordinates = nullptr;
	m_calcwf = false;
	m_dotemperaturewf = false;
	m_wfspecies.clear();
	m_wavelengths.clear();
}

/**
 * SKTRAN_TIR_Engine::ReferencePoint
 * 2018-09-13
 */
GEODETIC_INSTANT SKTRAN_TIR_Engine::ReferencePoint() const
{
	GEODETIC_INSTANT point;

	point.latitude = m_coordinates->ReferencePtLatitude();
	point.longitude = m_coordinates->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd = m_coordinates->ReferencePointMJD();

	return point;
}

/**
 * SKTRAN_TIR_Engine::~SKTRAN_TIR_Engine
 * 2018-09-13
 */
SKTRAN_TIR_Engine::~SKTRAN_TIR_Engine()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_Engine::ConfigureModel
 * 2018-09-13
 */
bool SKTRAN_TIR_Engine::ConfigureModel(
	SKTRAN_SpecsUser_Base& modelspecifications,
	const std::vector<double>& wavelen,
	const SKTRAN_LineOfSightArray_V21& linesofsight,
	size_t numthreads)
{
	bool ok = true;

	// check for at least one wavelength
	if (wavelen.size() < 1)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Engine::ConfigureModel, no wavelengths specified; radiance will be empty");
	}

	// check that wavelengths are in increasing order
	for (int wavelidx = 0; wavelidx < (int)wavelen.size() - 1; wavelidx++)
	{
		if (wavelen[wavelidx + 1] < wavelen[wavelidx])
		{
			ok = false;
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Engine::ConfigureModel, wavelengths must be in increasing order (calculation results may be incorrect)");
		}
	}
	m_wavelengths = wavelen;

	const SKTRAN_Specifications_Base* modelspecs = dynamic_cast<const SKTRAN_Specifications_Base*>(&modelspecifications);
	const SKTRAN_TIR_Specs_User* userspecs = dynamic_cast<const SKTRAN_TIR_Specs_User*>(modelspecs);
	ok = ok && nullptr != modelspecs && nullptr != userspecs;

	ok = ok && m_internalspecs.CreateCoordinates(&m_coordinates, linesofsight, userspecs->RayTracingSpecsConst().GetTOAHeight(),
												 userspecs->RayTracingSpecsConst().GetGeoidModel(), userspecs->RayTracingSpecsConst().UseManualGeoidModel());
	ok = ok && m_internalspecs.Configure(*modelspecs, linesofsight);
	m_calcwf = m_internalspecs.CalcWf();
	m_dotemperaturewf = m_internalspecs.WeightingFunctionSpecs().DoTemperatureWFCalculation();
	ok = ok && m_linesofsighttable.SetLinesOfSight(linesofsight, *m_coordinates);
	ok = ok && m_internalspecs.OpticalPropertiesSpecs().CreateOpticalTable(m_opticalpropertiestable, *m_coordinates, m_internalspecs.RayTracerSpecs().TOAHeight());
	m_opticalpropertiestable->SetWavelengths(m_wavelengths);
	m_opticalpropertiestable->ConfigureGeometry(&modelspecifications);
	m_opticalpropertiestable->SetCoords(m_coordinates);
	if (m_calcwf || m_dotemperaturewf)
	{
		m_calcwf = true;
		ok = ok && m_internalspecs.WeightingFunctionSpecs().MakePerturbationList();
	}
	ok = ok && m_internalspecs.IntegratorSpecs().CreateIntegrator(*m_opticalpropertiestable, m_integrator);
	ok = ok && m_internalspecs.RayTracerSpecs().CreateLineOfSightRayFactory(m_linesofsightrayfactory, m_coordinates);
	ok = ok && m_linesofsighttable.CreateRays(m_linesofsightrayfactory.get());
	m_docurvedrays = m_internalspecs.RayTracerSpecs().UseCurvedRays();
	m_wfperturbs = m_internalspecs.WeightingFunctionSpecs().PertList();
	m_wfspecies = m_internalspecs.WeightingFunctionSpecs().WFSpecies();
	if (m_dotemperaturewf)
	{
		m_wfspecies.push_back(SKCLIMATOLOGY_TEMPERATURE_K);
	}

	ok = ok && m_threadmanager.Initialize(numthreads, linesofsight, m_coordinates, m_linesofsightrayfactory.get());

	return ok;
}

/**
 * SKTRAN_TIR_Engine::CalculateRadiance
 * 2018-09-13
 */
bool SKTRAN_TIR_Engine::CalculateRadiance(
	std::vector<std::vector<SKTRAN_StokesScalar>>* losradiance,
	SKTRAN_TIR_AtmosphericOpticalState* opticalstate,
	bool usecachedcrosssections)
{
	bool ok = true;
	GEODETIC_INSTANT point;

	// check that the optical state contains the species that the user wanted weighting functions for (as specified when ConfigureModel was called)
	for (auto& species : m_wfspecies)
	{
		if (species == SKCLIMATOLOGY_TEMPERATURE_K) continue;
		if (!opticalstate->ContainsSpecies(species))
		{
			ok = false;
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Engine::CalculateRadiance, atmospheric optical state must contain any species for which weighting functions are to be calculated");
		}
	}

	point.latitude = m_coordinates->ReferencePtLatitude();
	point.longitude = m_coordinates->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd = m_coordinates->ReferencePointMJD();

	opticalstate->SetTimeAndLocation(point, true);
	// Updates the climatology caches of optical property objects
	ok = ok && opticalstate->UpdateCache();

	// set number of threads used in the cross section calculation
	opticalstate->SetNumThreads(m_threadmanager.ThreadStorage().size());

	if (usecachedcrosssections)
	{
		ok = ok && m_opticalpropertiestable->ConfigureOpticalFromCache(*opticalstate);
	}
	else
	{
		ok = ok && m_opticalpropertiestable->ConfigureOptical(*opticalstate);
	}

	losradiance->resize(m_wavelengths.size());
	if (m_calcwf)
	{
		m_wf.resize(m_wavelengths.size());
	}
	for (size_t wavelidx = 0; wavelidx < m_wavelengths.size(); wavelidx++)
	{
		losradiance->at(wavelidx).resize(m_linesofsighttable.NumRays());
		if (m_calcwf)
		{
			m_wf.at(wavelidx).resize(m_linesofsighttable.NumRays());
		}
	}

	ok = ok && CalculateRadianceMultiThreaded(losradiance, opticalstate);

	return ok;
}

/**
 * SKTRAN_TIR_Engine::CalculateRadianceMultiThreaded
 * 2019-06-27
 *
 * Performs the multi-threaded portions of the radiative transfer integration.
 *
 * @param[out] losradiance The result of the scalar radiance calculation, a 2D array with the shape [wavelength][line of sight]
 * @param[in] opticalstate The atmospheric state object
 */
bool SKTRAN_TIR_Engine::CalculateRadianceMultiThreaded(
	std::vector<std::vector<SKTRAN_StokesScalar>>* losradiance,
	SKTRAN_TIR_AtmosphericOpticalState* opticalstate)
{
	bool ok = true;

	GEODETIC_INSTANT point;
	point.latitude = m_coordinates->ReferencePtLatitude();
	point.longitude = m_coordinates->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd = m_coordinates->ReferencePointMJD();

	HELIODETIC_VECTOR observer;
	HELIODETIC_UNITVECTOR look;
	size_t threadindex;

#pragma omp parallel for schedule(dynamic) private (observer, look, threadindex) num_threads ((int)m_threadmanager.ThreadStorage().size())
	for (int wavelidx = 0; wavelidx < (int)m_wavelengths.size(); wavelidx++)
	{
		threadindex = omp_get_thread_num();

		SKTRAN_TIR_Thread_Storage* storage;
		storage = &m_threadmanager.ThreadStorage()[threadindex];

#pragma omp critical
		{
            // The ray factory can only be configured for a single wavelength, so only one thread may access it at a time
			m_linesofsightrayfactory->ConfigureOptical(opticalstate, m_wavelengths[wavelidx], point); // calculates the index of refraction if curved rays are enabled
			storage->Configure(m_linesofsightrayfactory.get());
		}

		for (int losidx = 0; losidx < (int)storage->NumRays(); losidx++)
		{
			ok = ok && storage->TraceRayAt(losidx);
			ok = ok && m_integrator->IntegrateRay(storage->RayAt(losidx), losradiance->at(wavelidx).at(losidx), m_wfspecies, wavelidx);

			if (m_calcwf)
			{
				ok = ok && CalculateWeightingFunctionsForRay(losidx, storage->RayAt(losidx), wavelidx);
			}
		}
	}

	return ok;
}

/**
 * SKTRAN_TIR_Engine::CalculateWeightingFunctionsForRay
 * 2018-09-14
 *
 * Performs the weighting function calculation for a single wavelength and line of sight.
 *
 * @param[in] rayidx The line-of-sight index for this weighting function
 * @param[in] rayptr The optical ray which contains the weighting functions computed at each quadrature point of the ray
 * @param[in] wavelidx The wavelength index; wavelengths are in increasing order
 *
 * @pre m_integrator->IntegrateRay has been called with the ray passed to this function
 * @post The weighting function is calculated and stored in m_wf[wavelidx][rayidx]
 */
bool SKTRAN_TIR_Engine::CalculateWeightingFunctionsForRay(
	const size_t rayidx,
	SKTRAN_RayOptical_Base* rayptr,
	const size_t wavelidx)
{
	bool ok = true;

	HELIODETIC_POINT quadpoint;
	bool isperturb;
	SKTRAN_RayStorage_Base* raystorage = dynamic_cast<SKTRAN_RayStorage_Base*>(rayptr->StorageVar());
	m_wf[wavelidx][rayidx].resize(m_wfspecies.size());
	for (size_t speciesidx = 0; speciesidx < m_wfspecies.size(); speciesidx++)
	{
		m_wf[wavelidx][rayidx][speciesidx].assign(m_wfperturbs.size(), 0.0);
	}
	for (size_t quadidx = 0; quadidx < rayptr->GetNumQuadraturePoints(); quadidx++)
	{
		ok = ok && raystorage->LocationOfPoint(quadidx, &quadpoint);
		for (size_t pertidx = 0; pertidx < m_wfperturbs.size(); pertidx++)
		{
			double pertweight;
			ok = ok && m_wfperturbs[pertidx].PerturbationWeight(quadpoint, &isperturb, &pertweight);
			if (isperturb)
			{
				for (size_t speciesidx = 0; speciesidx < m_wfspecies.size(); speciesidx++)
				{
					m_wf[wavelidx][rayidx][speciesidx][pertidx] += raystorage->WFAtPoint(m_wfspecies[speciesidx], quadidx) * pertweight;
				}
			}
		}
	}

	return ok;
}

/**
 * SKTRAN_TIR_Engine::WFHeights
 * 2018-09-13
 */
std::vector<double> SKTRAN_TIR_Engine::WFHeights() const
{
	std::vector<double> wfheights;
	wfheights.resize(m_wfperturbs.size());
	for (size_t idx = 0; idx < m_wfperturbs.size(); idx++)
	{
		wfheights[idx] = m_wfperturbs[idx].PerturbationLocation(*m_coordinates).Altitude();
	}
	return wfheights;
}

/**
 * SKTRAN_TIR_Engine::ReleaseResources
 * 2018-09-13
 */
void SKTRAN_TIR_Engine::ReleaseResources()
{
	if (nullptr != m_opticalpropertiestable)
	{
		m_opticalpropertiestable->Release();
		m_opticalpropertiestable.release();
	}
	if (nullptr != m_integrator)
	{
		m_integrator->Release();
		m_integrator.release();
	}
}
