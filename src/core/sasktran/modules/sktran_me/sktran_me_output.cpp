#include "include/sktran_me_internals.h"

namespace sktran_me {

    template<int NSTOKES>
    OutputSKIF<NSTOKES>::OutputSKIF(WFHandler &wfhandler,
                                    const std::vector<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>> &species_quantities,
                                    const sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>& total_quantities,
                                    const std::vector<GUID> &species_guids) :
                                    m_wf_handler(wfhandler),
                                    m_species_quantities(species_quantities),
                                    m_total_quantities(total_quantities),
                                    m_species_guids(species_guids)
                                    {

    }

    template<int NSTOKES>
    void OutputSKIF<NSTOKES>::resize(int nlos, int nwavel, int nderiv) {
        sasktran2::Output<NSTOKES>::resize(nlos, nwavel, nderiv);

        m_radiance.resize(nlos*nwavel*NSTOKES);

        m_wf.resize(nlos*nwavel*m_wf_handler.num_output_wf() * NSTOKES);

    }

    template<int NSTOKES>
    void OutputSKIF<NSTOKES>::assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES> &radiance,
                                     int losidx, int wavelidx) {
        int linear_index = NSTOKES * this->m_nlos * wavelidx + NSTOKES * losidx;

        for(int i = 0; i < NSTOKES; ++i) {
            m_radiance(linear_index + i) = radiance.value(i);
        }

        int num_geo = m_wf_handler.geometry()->size();
        int deriv_start = linear_index * m_wf_handler.num_output_wf();

        for(int i = 0; i < m_wf_handler.num_wf_types(); ++i) {
            if(m_wf_handler.wf_type(i) == WFHandler::wftype::NumberDensity) {
                auto it = std::find(std::begin(m_species_guids), std::end(m_species_guids), m_wf_handler.wf_guid(i));
                int species_index = std::distance(std::begin(m_species_guids), it);

                const auto& spec_quantities = m_species_quantities[species_index];

                for(int j = 0; j < num_geo; ++j) {
                    // dk/dn = cross section in cm2
                    double dk = 100*spec_quantities.total_extinction(j, wavelidx);

                    // Scaling factor
                    double f = m_total_quantities.f(j, wavelidx);

                    // Unscaled ssa
                    double ssa = m_total_quantities.ssa(j, wavelidx) / (1 - f + f * m_total_quantities.ssa(j, wavelidx));

                    // Unscaled k
                    double k = m_total_quantities.total_extinction(j, wavelidx) / (1 - ssa*f);

                    // dssa/dn
                    double dssa = dk * (spec_quantities.ssa(j, wavelidx) - ssa) / k;

                    // Model returns back derivatives with respect to k* and w*, have to undo the scalings
                    double dk_scaled_by_dk = (1 - f * ssa);
                    double dssa_scaled_by_dssa = (1 - f) / (1 - f*ssa) + f * (m_total_quantities.ssa(j, wavelidx)) / (1-ssa*f);

                    // Get an extra dssa contribution on dI/dk* from the scaling
                    m_wf(deriv_start) = (dk*dk_scaled_by_dk - f*dssa*k)*radiance.deriv(0, j) + dssa*dssa_scaled_by_dssa*radiance.deriv(0, j + num_geo);

                    ++deriv_start;
                }
            }
        }

    }

    template class OutputSKIF<1>;
    template class OutputSKIF<3>;

}