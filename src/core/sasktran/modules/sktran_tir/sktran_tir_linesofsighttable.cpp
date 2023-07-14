/**
 * SASKTRAN TIR Lines of Sight Table
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_LinesOfSightTable::~SKTRAN_TIR_LinesOfSightTable
 * 2019-04-22
 */
SKTRAN_TIR_LinesOfSightTable::~SKTRAN_TIR_LinesOfSightTable()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_LinesOfSightTable::ReleaseResources()
 * 2019-04-22
 */
bool SKTRAN_TIR_LinesOfSightTable::ReleaseResources()
{
	m_opticalrays.clear();
	return true;
}

/**
 * SKTRAN_TIR_LinesOfSightTable::SetLinesOfSight
 * 2019-04-22
 *
 * @param[in] linesofsight Lines of sight to compute radiative transfer along
 * @param[in] coords Coordinate system used in this calculation
 *
 * @pre The coords must be configured prior to calling this function, i.e. by calling the member function
 *      ConfigureCoordinates of the SKTRAN_CoordinateTransform_V2 object, with the reference point used in this
 *      calculation.
 */
bool SKTRAN_TIR_LinesOfSightTable::SetLinesOfSight(
	const SKTRAN_LineOfSightArray_V21& linesofsight,
	const SKTRAN_CoordinateTransform_V2& coords)
{
	bool ok = true;

	ok = ok && m_observerlinesofsight.DeepCopy(linesofsight);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_LinesOfSightTable, There were errors copying the lines of sight");
		return ok;
	}

	// translate observer to osculating sphere
	SKTRAN_LineOfSightEntry_V2* entry;
	nxVector offsetobserver;
	for (size_t idx = 0; idx < m_observerlinesofsight.NumRays(); idx++)
	{
		m_observerlinesofsight.GetRayVar(idx, &entry);
		offsetobserver = coords.TranslateGeoidToOsculatingSphere(entry->Observer());
		entry->Configure(offsetobserver, entry->Look(), entry->Mjd());
	}


	return ok;
}

/**
 * SKTRAN_TIR_LinesOfSightTable::CreateRays
 * 2019-04-22
 *
 * Create the ray objects used in raytracing calculations.
 *
 * @param[in] rayfactory Pointer to the factory object used to create rays.
 *
 * @pre The rayfactory must be properly configured for the desired ray tracing settings, i.e. by using the internal
 *      specifications function SKTRAN_TIR_Specs_Internal_RayTracer::CreateRayFactory
 */
bool SKTRAN_TIR_LinesOfSightTable::CreateRays(
	const SKTRAN_RayFactory_Base* rayfactory)
{
	bool	ok = true;
	bool	ok1;

	const SKTRAN_LineOfSightEntry_V2*	entry;
	size_t								idx;
	size_t								numrays;
	HELIODETIC_VECTOR					observer;
	HELIODETIC_UNITVECTOR				look;

	numrays = m_observerlinesofsight.NumRays();
	m_opticalrays.resize(numrays);
	for (idx = 0; idx < numrays; idx++)
	{
		ok1 = rayfactory->CreateRayObject(&m_opticalrays[idx]);
		ok1 = ok1 && m_observerlinesofsight.GetRay(idx, &entry);
		if (ok1)
		{
			observer = rayfactory->CoordsPtr()->GeographicToHelio(entry->Observer());
			look = rayfactory->CoordsPtr()->GeographicToHelio(entry->Look()).UnitVector();
			ok1 = m_opticalrays[idx]->MoveObserver(observer, look);
		}
		ok = ok && ok1;
	}

	return ok;
}

/**
 * SKTRAN_TIR_LinesOfSightTable::RayAt
 * 2019-04-22
 *
 * @param[in] idx Index of the ray pointer to retrieve
 */
SKTRAN_RayOptical_Base* SKTRAN_TIR_LinesOfSightTable::RayAt(
	size_t idx)
{
	return m_opticalrays.at(idx).get();
}

const SKTRAN_RayOptical_Base* SKTRAN_TIR_LinesOfSightTable::RayAt(
	size_t idx) const
{
	return m_opticalrays.at(idx).get();
}

/**
 * SKTRAN_TIR_LinesOfSightTable::RayEntryAt
 * 2019-04-22
 *
 * @param[in] idx Index of ray to return a unique pointer to
 */
std::unique_ptr<SKTRAN_RayOptical_Base>&  SKTRAN_TIR_LinesOfSightTable::RayEntryAt(
	size_t idx)
{
	return m_opticalrays.at(idx);
}

/**
 * SKTRAN_TIR_LinesOfSightTable::MeanMJD
 * 2019-04-22
 */
double SKTRAN_TIR_LinesOfSightTable::MeanMJD() const
{
	double	mjd = 0;
	size_t	numrays;

	numrays = m_observerlinesofsight.NumRays();
	for (size_t idx = 0; idx < numrays; idx++)
	{
		mjd += m_observerlinesofsight.Entry(idx)->Mjd();
	}
	mjd /= numrays;
	return mjd;
}
