#include "include/sktran_me_internals.h"
#include "include/sktran_me_atmosphere_interface.h"

namespace sktran_me {
    void WFHandler::set_from_strings(const std::vector<std::string> &wf_handles) {
        m_wf_guids.resize(wf_handles.size());
        m_wf_types.resize(wf_handles.size());

        for (int i = 0; i < wf_handles.size(); i++) {
            auto handle = FindGlobalClimatologyHandle(wf_handles[i].c_str(), false);

            if(*handle == SKCLIMATOLOGY_UNDEFINED) {

            } else {
                // Species was in the handle table and so it is a simple numberdensity calculation
                m_wf_guids[i] = *handle;
                m_wf_types[i] = wftype::NumberDensity;
            }
        }
    }

    void WFHandler::set_scattering_information(AtmosphereInterface &atmo_interface,
                                               const std::vector<GUID>& species_guids
                                               ) {
        m_num_scatterers = 0;

        m_species_index_map.resize(m_wf_types.size());
        m_scattering_map.resize(m_wf_types.size());

        for(int i = 0; i < m_wf_types.size(); ++i) {
            // TODO: Optical property derivatives

            if(m_wf_types[i] == WFHandler::wftype::NumberDensity) {
                auto it = std::find(std::begin(species_guids), std::end(species_guids), m_wf_guids[i]);
                int species_index = std::distance(std::begin(species_guids), it);

                m_species_index_map[i] = species_index;

                if(atmo_interface.species()[species_index].ParticleOpticalProps()->IsScatterer()) {
                    m_scattering_map[i] = m_num_scatterers;
                    ++m_num_scatterers;
                } else {
                    m_scattering_map[i] = -1;
                }
            }
        }
    }

    bool WFHandler::calculating_wf() const {
        return m_wf_guids.size() > 0;
    }

    int WFHandler::num_output_wf() const {
        int num_wf = 0;

        for(int i = 0; i < m_wf_types.size(); ++i) {
            if(m_wf_types[i] == wftype::NumberDensity) {
                num_wf += m_geometry->size();
            }
        }

        return num_wf;

    }

}