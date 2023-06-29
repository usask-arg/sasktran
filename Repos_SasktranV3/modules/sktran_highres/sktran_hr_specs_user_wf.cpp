#include "include/sktran_hr_internals.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

SKTRAN_HR_Specs_User_wf::SKTRAN_HR_Specs_User_wf()
{
	m_doextinction = true;
	m_doscatextinction = false;
	m_wfmode = SKTRAN_HR_wf_Mode_None;
	m_wfprecision = SKTRAN_HR_wf_precision_all;

	m_wfinterpwidth = 1000.0;
	m_maxwfheight = 99500;
	m_manualwfresolution = 1000.0;

	m_aerosolsizepercentchange = 5;

	m_aerosolmode = SKTRAN_HR_wf_aerosol_Mode::SKTRAN_HR_wf_aerosol_Mode_numberdensity;
}

void SKTRAN_HR_Specs_User_wf::SetWeightingFunctionSpeciesString(const std::vector<std::string>& species)
{
	m_wfspecies.resize(species.size());
	m_wfspeciesmode.resize(species.size());

	for (int i = 0; i < species.size(); i++)
	{
		auto handle = FindGlobalClimatologyHandle(species[i].c_str(), false);

		if (*handle == SKCLIMATOLOGY_UNDEFINED)
		{
			// Now we have a few options, the wf species could either be wrong, or it could be a special mode of something like
			// aerosol particle size
			if (boost::algorithm::ends_with(species[i].c_str(), "_lognormal_modewidth"))
			{
				m_wfspeciesmode[i] = SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_LogNormal_ModeWidth;
				
				std::string newstr = species[i];
				boost::erase_all(newstr, "_lognormal_modewidth");

				auto handle = FindGlobalClimatologyHandle(newstr.c_str());
				if (*handle == SKCLIMATOLOGY_UNDEFINED)
				{
					nxLog::Record(NXLOG_WARNING, "Error looking up climatology handle: %s", species[i].c_str());
				}
				m_wfspecies[i] = *handle;
			}
			else if (boost::algorithm::ends_with(species[i].c_str(), "_lognormal_medianradius"))
			{
				m_wfspeciesmode[i] = SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_LogNormal_ModeRadius;

				std::string newstr = species[i];
				boost::erase_all(newstr, "_lognormal_medianradius");

				auto handle = FindGlobalClimatologyHandle(newstr.c_str());
				if (*handle == SKCLIMATOLOGY_UNDEFINED)
				{
					nxLog::Record(NXLOG_WARNING, "Error looking up climatology handle: %s", species[i].c_str());
				}
				m_wfspecies[i] = *handle;
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "Error looking up climatology handle: %s", species[i].c_str());
			}
		}
		else
		{
			// Species was in the handle table and so it is a simple numberdensity calculation
			m_wfspecies[i] = *handle;
			m_wfspeciesmode[i] = SKTRAN_HR_wf_Species_Mode::SKTRAN_HR_wf_Species_numberdensity;
		}
	}
}