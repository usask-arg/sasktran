/**
 * SASKTRAN TIR Thread Manager
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Thread_Manager::SKTRAN_TIR_Thread_Manager
 * 2019-05-22
 */
SKTRAN_TIR_Thread_Manager::SKTRAN_TIR_Thread_Manager()
{
	m_numthreads = 0;
}

/**
 * SKTRAN_TIR_Thread_Manager::~SKTRAN_TIR_Thread_Manager
 * 2019-05-22
 */
SKTRAN_TIR_Thread_Manager::~SKTRAN_TIR_Thread_Manager()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_Thread_Manager::Initialize
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Manager::Initialize(
	size_t numthreads,
	const SKTRAN_LineOfSightArray_V21& linesofsight,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
	const SKTRAN_RayFactory_Base* rayfactory)
{
	bool ok = true;

	if (numthreads > 0)
	{
		m_numthreads = numthreads;
	}
	else
	{
#if defined(NXDEBUG)
		m_numthreads = 1;
#else
		m_numthreads = omp_get_num_procs();
#endif
	}

	m_threadstore.resize(m_numthreads);
	for (size_t idx = 0; idx < m_numthreads; idx++)
	{
		ok = ok && m_threadstore[idx].Initialize(linesofsight, coords);
	}

	omp_set_num_threads((int)m_numthreads);

	return ok;
}

/**
 * SKTRAN_TIR_Thread_Manager::ReleaseResources
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Manager::ReleaseResources()
{
	bool ok = true;

	m_threadstore.clear();

	return ok;
}

/**
 * SKTRAN_TIR_Thread_Manager::CalculateRadianceMultiWavel
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Manager::CalculateRadianceMultiWavel(
	std::vector<std::vector<SKTRAN_StokesScalar>>* losradiance,
	std::vector<double>& wavelen,
	SKTRAN_TIR_AtmosphericOpticalState* opticalstate,
	std::vector<std::vector<skRTStokesVector>>* losvector)
{
	return false;
}

/**
 * SKTRAN_TIR_Thread_Storage::Initialize
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Storage::Initialize(
	const SKTRAN_LineOfSightArray_V21& linesofsight,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords)
{
	bool ok = true;

	ok = ok && m_linesofsighttable.SetLinesOfSight(linesofsight, *coords);

	return ok;
}

/**
 * SKTRAN_TIR_Thread_Storage::Configure
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Storage::Configure(
	const SKTRAN_RayFactory_Base* rayfactory)
{
	bool ok = true;

	ok = ok && m_linesofsighttable.CreateRays(rayfactory);

	return ok;
}

/**
 * SKTRAN_TIR_Thread_Storage::TraceRayAt
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Storage::TraceRayAt(
	size_t losidx)
{
	return m_linesofsighttable.RayAt(losidx)->TraceRay_NewMethod();
}

/**
 * SKTRAN_TIR_Thread_Storage::SKTRAN_TIR_Thread_Storage
 * 2019-05-22
 */
SKTRAN_TIR_Thread_Storage::SKTRAN_TIR_Thread_Storage()
{
}

/**
 * SKTRAN_TIR_Thread_Storage::~SKTRAN_TIR_Thread_Storage
 * 2019-05-22
 */
SKTRAN_TIR_Thread_Storage::~SKTRAN_TIR_Thread_Storage()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_Thread_Storage::ReleaseResources
 * 2019-05-22
 */
bool SKTRAN_TIR_Thread_Storage::ReleaseResources()
{
	return true;
}
