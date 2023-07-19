/**
 * SASKTRAN TIR Ray Storage
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::Reserve
 * 2018-09-21
 */
bool SKTRAN_RayStorage_Straight_TIR::Reserve(
	size_t numquadraturepoints)
{
	m_extinction.reserve(numquadraturepoints);
	for (auto& wf : m_wf)
	{
		wf.second.reserve(numquadraturepoints);
	}
	return SKTRAN_RayStorage_Straight::Reserve(numquadraturepoints);
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::Resize
 * 2018-09-21
 */
bool SKTRAN_RayStorage_Straight_TIR::Resize(
	size_t numquadraturepoints)
{
	m_extinction.resize(numquadraturepoints);
	for (auto& wf : m_wf)
	{
		wf.second.resize(numquadraturepoints);
	}
	return SKTRAN_RayStorage_Straight::Resize(numquadraturepoints);
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::TruncateToNumElements
 * 2018-09-21
 */
void SKTRAN_RayStorage_Straight_TIR::TruncateToNumElements(
	size_t numels)
{
	m_extinction.resize(numels);
	for (auto& wf : m_wf)
	{
		wf.second.resize(numels);
	}
	return SKTRAN_RayStorage_Straight::TruncateToNumElements(numels);
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::ClearStorage
 * 2018-09-21
 */
void SKTRAN_RayStorage_Straight_TIR::ClearStorage()
{
	m_extinction.clear();
	for (auto& wf : m_wf)
	{
		wf.second.clear();
	}
	m_wf.clear();
	return SKTRAN_RayStorage_Straight::ClearStorage();
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::PushBack
 * 2018-09-21
 */
bool SKTRAN_RayStorage_Straight_TIR::PushBack(
	SKTRAN_Distance r,
	SKTRAN_Distance distFromTan,
	SKTRAN_Distance s)
{
	m_extinction.push_back(-1);
	for (auto& wf : m_wf)
	{
		wf.second.push_back(-1);
	}
	return SKTRAN_RayStorage_Straight::PushBack(r, distFromTan, s);
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::Insert
 * 2018-09-21
 */
bool SKTRAN_RayStorage_Straight_TIR::Insert(
	SKTRAN_Distance r,
	SKTRAN_Distance distFromTan,
	SKTRAN_Distance s,
	size_t index)
{
	m_extinction.insert(std::begin(m_extinction) + index, -1);
	for (auto& wf : m_wf)
	{
		wf.second.insert(std::begin(wf.second) + index, -1);
	}
	return SKTRAN_RayStorage_Straight::Insert(r, distFromTan, s, index);
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::CellLength
 * 2018-09-21
 */
double SKTRAN_RayStorage_Straight_TIR::CellLength(
	size_t cellindex) const
{
	double	s0;
	double	s1;
	double	length;

	s0 = SKTRAN_RayStorage_Straight::DistanceOfPointFromOrigin(cellindex);
	s1 = SKTRAN_RayStorage_Straight::DistanceOfPointFromOrigin(cellindex + 1);

	length = fabs(s1 - s0);

	return length;
}

/**
 * SKTRAN_TIR_RayStorage_Straight_TIR::AddWFSpecies
 * 2018-09-21
 */
void SKTRAN_RayStorage_Straight_TIR::AddWFSpecies(
	const CLIMATOLOGY_HANDLE& species) const
{
	m_wf.emplace(species, std::vector<double>());
	m_wf[species].clear();
	m_wf[species].resize(SKTRAN_RayStorage_Straight::NumQuadraturePoints());
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::Reserve
 * 2018-10-10
 */
bool SKTRAN_RayStorage_CurvedPiecewise_TIR::Reserve(
	size_t numquadraturepoints)
{
	m_extinction.reserve(numquadraturepoints);
	for (auto& wf : m_wf)
	{
		wf.second.reserve(numquadraturepoints);
	}
	return SKTRAN_RayStorage_CurvedPiecewise::Reserve(numquadraturepoints);
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::Resize
 * 2018-10-10
 */
bool SKTRAN_RayStorage_CurvedPiecewise_TIR::Resize(
	size_t numquadraturepoints)
{
	m_extinction.resize(numquadraturepoints);
	for (auto& wf : m_wf)
	{
		wf.second.resize(numquadraturepoints);
	}
	return SKTRAN_RayStorage_CurvedPiecewise::Resize(numquadraturepoints);
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::TruncateToNumElements
 * 2018-10-10
 */
void SKTRAN_RayStorage_CurvedPiecewise_TIR::TruncateToNumElements(
	size_t numels)
{
	m_extinction.resize(numels);
	for (auto& wf : m_wf)
	{
		wf.second.resize(numels);
	}
	return SKTRAN_RayStorage_CurvedPiecewise::TruncateToNumElements(numels);
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::ClearStorage
 * 2018-10-10
 */
void SKTRAN_RayStorage_CurvedPiecewise_TIR::ClearStorage()
{
	m_extinction.clear();
	for (auto& wf : m_wf)
	{
		wf.second.clear();
	}
	m_wf.clear();
	return SKTRAN_RayStorage_CurvedPiecewise::ClearStorage();
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::PushBack
 * 2018-10-10
 */
bool SKTRAN_RayStorage_CurvedPiecewise_TIR::PushBack(
	HELIODETIC_UNITVECTOR* uv,
	HELIODETIC_POINT* pt,
	SKTRAN_Distance celllength)
{
	m_extinction.push_back(-1);
	for (auto& wf : m_wf)
	{
		wf.second.push_back(-1);
	}
	return SKTRAN_RayStorage_CurvedPiecewise::PushBack(uv, pt, celllength);
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::Insert
 * 2018-10-10
 */
bool SKTRAN_RayStorage_CurvedPiecewise_TIR::Insert(
	HELIODETIC_UNITVECTOR* uv,
	HELIODETIC_POINT* pt,
	SKTRAN_Distance celllength,
	size_t index)
{
	m_extinction.insert(std::begin(m_extinction) + index, -1);
	for (auto& wf : m_wf)
	{
		wf.second.insert(std::begin(wf.second) + index, -1);
	}
	return SKTRAN_RayStorage_CurvedPiecewise::Insert(uv, pt, celllength, index);
}

/**
 * SKTRAN_TIR_RayStorage_CurvedPiecewise_TIR::AddWFSpecies
 * 2018-09-21
 */
void SKTRAN_RayStorage_CurvedPiecewise_TIR::AddWFSpecies(
	const CLIMATOLOGY_HANDLE& species) const
{
	m_wf.emplace(species, std::vector<double>());
	m_wf[species].resize(SKTRAN_RayStorage_CurvedPiecewise::NumQuadraturePoints());
}

