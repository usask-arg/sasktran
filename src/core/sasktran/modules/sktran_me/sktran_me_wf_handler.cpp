#include "include/sktran_me_internals.h"

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