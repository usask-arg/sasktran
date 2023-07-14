/**
 * SASKTRAN TIR Internal Weighting Function Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_Internal_wf::Configure
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_wf::Configure(
	const SKTRAN_TIR_Specs_User& specs)
{
	m_manualwfresolution = specs.WeightingFunctionSpecsConst().m_manualwfresolution;
	m_maxwfheight = specs.WeightingFunctionSpecsConst().m_maxwfheight;
	m_wfinterpwidth = specs.WeightingFunctionSpecsConst().m_wfinterpwidth;
	m_wfspecies = specs.WeightingFunctionSpecsConst().m_wfspecies;
	m_dotemperaturewf = specs.WeightingFunctionSpecsConst().m_dotemperaturewf;

	m_manualwfheights = specs.WeightingFunctionSpecsConst().m_manualwfheights;
	m_manualwfwidths = specs.WeightingFunctionSpecsConst().m_manualwfwidths;

	if (m_manualwfheights.size() == 0)
	{
		size_t numpert = (size_t)(ceil(m_maxwfheight / m_manualwfresolution));
		m_manualwfheights.resize(numpert);
		for (size_t idx = 0; idx < numpert; idx++)
		{
			m_manualwfheights[idx] = (idx + 0.5)*m_manualwfresolution;
		}
	}

	if (m_manualwfwidths.size() == 0)
	{
		m_manualwfwidths.resize(m_manualwfheights.size());
		for (size_t idx = 0; idx < m_manualwfwidths.size(); idx++)
		{
			m_manualwfwidths[idx] = m_wfinterpwidth;
		}
	}

	return true;
}

/**
 * SKTRAN_TIR_Specs_Internal_wf::AddWfGeometryToRayTracer
 * 2019-04-22
 *
 * TODO: make this work with curved raytracer
 */
bool SKTRAN_TIR_Specs_Internal_wf::AddWfGeometryToRayTracer(
	SKTRAN_RayTracer_Straight_Generic& raytracer,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords) const
{
	bool ok = true;
	
	size_t numbound;
	for (size_t idx = 0; idx < m_wfperturbs.size(); idx++)
	{
		numbound = m_wfperturbs[idx].NumBoundingGeometry();
		for (size_t boundidx = 0; boundidx < numbound; boundidx++)
		{
			ok = ok && raytracer.AddGeometryObject(m_wfperturbs[idx].BoundingGeometryObject(coords, boundidx));
		}
	}
	
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_wf::MakeOneDimUniformPerturbations
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_wf::MakeOneDimUniformPerturbations()
{
	bool ok = true;

	bool manualwidthspecified = m_manualwfwidths.size() != 0;
	if (m_manualwfheights.size() == 0)
	{
		size_t numpert = (size_t)(ceil(m_maxwfheight / m_manualwfresolution));
		m_wfperturbs.resize(numpert);
		for (size_t idx = 0; idx < numpert; idx++)
		{
			if (manualwidthspecified)
			{
				m_wfperturbs[idx].Initialize((idx + 0.5) * m_manualwfresolution, m_manualwfwidths[idx], m_manualwfwidths[idx]);
			}
			else
			{
				m_wfperturbs[idx].Initialize((idx + 0.5) * m_manualwfresolution, m_wfinterpwidth, m_wfinterpwidth);
			}
		}
	}
	else
	{
		size_t numpert = m_manualwfheights.size();
		m_wfperturbs.resize(numpert);
		for (size_t idx = 0; idx < numpert; idx++)
		{
			if (manualwidthspecified)
			{
				m_wfperturbs[idx].Initialize(m_manualwfheights[idx], m_manualwfwidths[idx], m_manualwfwidths[idx]);
			}
			else
			{
				m_wfperturbs[idx].Initialize(m_manualwfheights[idx], m_wfinterpwidth, m_wfinterpwidth);
			}
		}
	}

	if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Specs_Internal_wf::Configure, Error occurred");

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_wf::MakePerturbationList
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_wf::MakePerturbationList()
{
	// if different types of perturbations are created (e.g. 2D), the appropriate creation method could be called here
	return MakeOneDimUniformPerturbations();
}
