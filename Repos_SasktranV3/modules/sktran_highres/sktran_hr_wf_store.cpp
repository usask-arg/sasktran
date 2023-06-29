#include "include/sktran_hr_internals.h"

template <typename T>
std::pair<size_t, size_t> SKTRAN_HR_WF_TVector<T>::bounding_indicies(double altitude) const
{
	if (!m_issorted)
	{
		return std::pair<size_t, size_t>(0, m_loweraltitudes.size());
	}

	auto lower_it = std::lower_bound(std::begin(m_loweraltitudes), std::end(m_loweraltitudes), altitude);

	int low_idx = std::max((int)0, (int)std::distance(std::begin(m_loweraltitudes), lower_it) - 2);

	auto upper_it = std::lower_bound(std::begin(m_upperaltitudes), std::end(m_upperaltitudes), altitude);

	int upper_idx = std::min((int)m_upperaltitudes.size(), (int)std::distance(std::begin(m_upperaltitudes), upper_it) + 2);

	return std::pair<size_t, size_t>(low_idx, upper_idx);
}


void SKTRAN_HR_WF_Store::ExtinctionPerturbation(const HELIODETIC_POINT & location, std::vector<double>& value) const
{
	value.resize(StoreSize());
	bool temp;
	size_t pertidx = 0;
	const double localt = location.Altitude();
	for (size_t idx = 0; idx < m_boxstore.m_storage.size(); idx++)
	{
		if (abs(localt - m_boxstorealts[idx]) < m_boxstorewidths[idx])
		{
			m_boxstore.m_storage[idx].ExtinctionPerturbation(location, temp, value[pertidx]);
		}
		else
		{
			value[pertidx] = 0.0;
		}
		++pertidx;
	}
	for (size_t idx = 0; idx < m_linearstore.m_storage.size(); idx++)
	{
		m_linearstore.m_storage[idx].ExtinctionPerturbation(location, temp, value[pertidx]);
		++pertidx;
	}
}

void SKTRAN_HR_WF_Store::AddExtinctionPerturbation(const HELIODETIC_POINT & location, std::vector<double>& value, double ds) const
{
	bool temp;
	size_t pertidx = 0;
	const double localt = location.Altitude();
	const size_t numbox = m_boxstore.m_storage.size();
	double kpert;
	for (size_t idx = 0; idx < numbox; idx++)
	{
		if (abs(localt - m_boxstorealts[idx]) < m_boxstorewidths[idx])
		{
			m_boxstore.m_storage[idx].ExtinctionPerturbation(location, temp, kpert);
			value[pertidx] += ds * kpert;
		}
		++pertidx;
	}
	std::pair<size_t, size_t> bounds = m_linearstore.bounding_indicies(location.Altitude());
	for (int idx = bounds.first; idx < bounds.second; idx++)
	{
		m_linearstore.m_storage[idx].ExtinctionPerturbation(location, temp, kpert);
		value[idx + numbox] += ds * kpert;
	}
}

void SKTRAN_HR_WF_Store::PerturbationLocation(const SKTRAN_CoordinateTransform_V2 & coords, std::vector<HELIODETIC_POINT>& data) const
{
	data.resize(StoreSize());

	size_t pertidx = 0;
	for (size_t idx = 0; idx < m_boxstore.m_storage.size(); idx++)
	{
		data[pertidx] = m_boxstore.m_storage[idx].PerturbationLocation(coords);
		++pertidx;
	}
	for (size_t idx = 0; idx < m_linearstore.m_storage.size(); idx++)
	{
		data[pertidx] = m_linearstore.m_storage[idx].PerturbationLocation(coords);
		++pertidx;
	}
}

void SKTRAN_HR_WF_Store::PerturbationAltitudeWidth(std::vector<double>& data) const
{
	data.resize(StoreSize());

	size_t pertidx = 0;
	for (size_t idx = 0; idx < m_boxstore.m_storage.size(); idx++)
	{
		data[pertidx] = m_boxstore.m_storage[idx].PerturbationAltitudeWidth();
		++pertidx;
	}
	for (size_t idx = 0; idx < m_linearstore.m_storage.size(); idx++)
	{
		data[pertidx] = m_linearstore.m_storage[idx].PerturbationAltitudeWidth();
		++pertidx;
	}
}

void SKTRAN_HR_WF_Store::AddGeometryToRayTracer(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_RayTracer_Straight_Generic & raytracer) const
{
	size_t numbound;
	for (size_t idx = 0; idx < m_boxstore.m_storage.size(); idx++)
	{
		numbound = m_boxstore.m_storage[idx].NumBoundingGeometry();
		for (size_t boundidx = 0; boundidx < numbound; boundidx++)
		{
			raytracer.AddGeometryObject(m_boxstore.m_storage[idx].BoundingGeometryObject(coords, boundidx));
		}
	}
	for (size_t idx = 0; idx < m_linearstore.m_storage.size(); idx++)
	{
		numbound = m_linearstore.m_storage[idx].NumBoundingGeometry();
		for (size_t boundidx = 0; boundidx < numbound; boundidx++)
		{
			raytracer.AddGeometryObject(m_linearstore.m_storage[idx].BoundingGeometryObject(coords, boundidx));
		}
	}
}

void SKTRAN_HR_WF_Store::SetPertVal(double value)
{
	size_t pertidx = 0;
	for (size_t idx = 0; idx < m_boxstore.m_storage.size(); idx++)
	{
		m_boxstore.m_storage[idx].SetPertVal(value);
	}
	for (size_t idx = 0; idx < m_linearstore.m_storage.size(); idx++)
	{
		m_linearstore.m_storage[idx].SetPertVal(value);
	}
}

const SKTRAN_HR_Perturbation_Base * SKTRAN_HR_WF_Store::RawAccess(size_t idx) const
{
	if (idx < m_boxstore.m_storage.size())
	{
		return &m_boxstore.m_storage[idx];
	}
	else
	{
		return &m_linearstore.m_storage[idx - m_boxstore.m_storage.size()];
	}
}

SKTRAN_HR_Perturbation_Base * SKTRAN_HR_WF_Store::RawAccess(size_t idx)
{
	if (idx < m_boxstore.m_storage.size())
	{
		return &m_boxstore.m_storage[idx];
	}
	else
	{
		return &m_linearstore.m_storage[idx - m_boxstore.m_storage.size()];
	}
}
