#pragma once

namespace sktran_me {
    template<int NSTOKES>
    class OutputSKIF : public sasktran2::Output<NSTOKES> {

    private:
        WFHandler& m_wf_handler;
        const std::vector<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>>& m_species_quantities;
        const sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>& m_total_quantities;
        const std::vector<GUID>& m_species_guids;

        Eigen::VectorXd m_radiance;
        Eigen::VectorXd m_wf;

    public:
        OutputSKIF(WFHandler& wfhandler,
                   const std::vector<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>>& species_quantities,
                   const sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>& total_quantities,
                   const std::vector<GUID>& species_guids
                   );

        void resize(int nlos, int nwavel, int nderiv);

        void assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& radiance, int losidx, int wavelidx);


        Eigen::VectorXd& radiance() { return m_radiance; }
        Eigen::VectorXd& wf() { return m_wf; }
    };
}