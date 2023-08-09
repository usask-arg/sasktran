#pragma once


namespace sktran_me {
    class WFHandler {
    public:
        enum class wftype {
            NumberDensity,

        };
    private:
        std::vector<wftype> m_wf_types;
        std::vector<GUID> m_wf_guids;

        const sasktran2::Geometry* m_geometry;

    public:
        WFHandler() {};

        void set_from_strings(const std::vector<std::string>& wf_handles);
        void set_geometry(const sasktran2::Geometry& geometry) { m_geometry = &geometry; }
        const sasktran2::Geometry* geometry() const { return m_geometry; }

        bool calculating_wf() const;

        int num_output_wf() const;
        int num_wf_types() const { return m_wf_types.size(); }

        wftype wf_type(int i) const { return m_wf_types[i]; }
        GUID wf_guid(int i) const { return m_wf_guids[i]; }

    };
}