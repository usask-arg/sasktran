/**
 * SASKTRAN TIR Weighting Functions
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Perturbation_Absorption_Linear::Initialize
 * 2018-09-13
 */
bool SKTRAN_TIR_Perturbation_Absorption_Linear::Initialize(
	double center,
	double distancetozerobelow,
	double distancetozeroabove)
{
	m_center = center;
	m_distancetozerobelow = distancetozerobelow;
	m_distancetozeroabove = distancetozeroabove;

	return true;
}

/**
 * SKTRAN_TIR_Perturbation_Absorption_Linear::PerturbationWeight
 * 2018-09-13
 */
bool SKTRAN_TIR_Perturbation_Absorption_Linear::PerturbationWeight(
	const HELIODETIC_POINT& location,
	bool* isperturbation,
	double* value) const
{
	bool ok = true;

	double alt = location.Altitude();
	// TEMPORARY test code
	//alt = floor((alt + 250.0) / 500.0) * 500.0;
	// end test code
	double disttocenter = alt - m_center;

	if (disttocenter < m_distancetozeroabove && disttocenter >= 0.0)
	{
		*isperturbation = true;
		*value = (1 - (disttocenter / m_distancetozeroabove));
	}
	else if (abs(disttocenter) < m_distancetozerobelow && disttocenter < 0.0)
	{
		*isperturbation = true;
		*value = (1 - (abs(disttocenter) / m_distancetozerobelow));
	}
	else
	{
		*isperturbation = false;
		*value = 0.0;
	}

	return ok;
}

/**
 * SKTRAN_TIR_Perturbation_Absorption_Linear::PerturbationWeight
 * 2019-03-05
 */
bool SKTRAN_TIR_Perturbation_Absorption_Linear::PerturbationWeight(
	const double altitude,
	bool* isperturbation,
	double* value) const
{
	bool ok = true;

	// TEMPORARY test code
	//alt = floor((alt + 250.0) / 500.0) * 500.0;
	// end test code
	double disttocenter = altitude - m_center;

	if (disttocenter < m_distancetozeroabove && disttocenter >= 0.0)
	{
		*isperturbation = true;
		*value = (1 - (disttocenter / m_distancetozeroabove));
	}
	else if (abs(disttocenter) < m_distancetozerobelow && disttocenter < 0.0)
	{
		*isperturbation = true;
		*value = (1 - (abs(disttocenter) / m_distancetozerobelow));
	}
	else
	{
		*isperturbation = false;
		*value = 0.0;
	}

	return ok;
}

/**
 * SKTRAN_TIR_Perturbation_Absorption_Linear::PerturbationLocation
 * 2018-09-13
 */
HELIODETIC_POINT SKTRAN_TIR_Perturbation_Absorption_Linear::PerturbationLocation(
	const SKTRAN_CoordinateTransform_V2& coords) const
{
	return coords.ReferencePoint((m_center));
}

/**
 * SKTRAN_TIR_Perturbation_Absorption_Linear::BoundingGeometryObject
 * 2018-09-13
 */
std::unique_ptr<SKTRAN_GeometryObject> SKTRAN_TIR_Perturbation_Absorption_Linear::BoundingGeometryObject(
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
	size_t idx) const
{
	std::unique_ptr<SKTRAN_GeometryObject> ret;

	double groundradius = coords->AltitudeToRadius(0.0);
	if (idx == 0)
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject>(new SKTRAN_GeometryObject_Sphere(m_center - m_distancetozerobelow + groundradius));
	}
	else if (idx == 1)
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject>(new SKTRAN_GeometryObject_Sphere(m_center + groundradius));
	}
	else if (idx == 2)
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject>(new SKTRAN_GeometryObject_Sphere(m_center + m_distancetozeroabove + groundradius));
	}
	else
	{
		ret = nullptr;
	}
	return ret;
}
