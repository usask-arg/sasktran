#pragma once


namespace sktran_me {
    class AtmosphereInterface;

    class WFHandler {
    public:
        enum class wftype {
            NumberDensity,

        };
    private:
        std::vector<wftype> m_wf_types;
        std::vector<GUID> m_wf_guids;
        std::vector<int> m_scattering_map; // Map from wf_index to scattering_deriv_index
        std::vector<int> m_species_index_map; // Map from wf_index to species_index

        const sasktran2::Geometry* m_geometry;

        int m_num_scatterers;

    public:
        WFHandler() {};

        void set_from_strings(const std::vector<std::string>& wf_handles);
        void set_geometry(const sasktran2::Geometry& geometry) { m_geometry = &geometry; }
        void set_scattering_information(AtmosphereInterface& atmo_interface, const std::vector<GUID>& species_guids);
        const sasktran2::Geometry* geometry() const { return m_geometry; }

        bool calculating_wf() const;

        int num_output_wf() const;
        int num_wf_types() const { return m_wf_types.size(); }

        wftype wf_type(int i) const { return m_wf_types[i]; }
        GUID wf_guid(int i) const { return m_wf_guids[i]; }
        int scattering_index(int i) const { return m_scattering_map[i]; }
        int species_index(int i) const { return m_species_index_map[i]; }

        int num_scattering_derivatives() const { return m_num_scatterers; }

    };
}