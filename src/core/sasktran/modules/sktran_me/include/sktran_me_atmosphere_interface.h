#pragma once

namespace sktran_me {
    class AtmosphereInterface {
    public:
        // Define a type for each species. PurelyAbsorbing has no scattering componont, ScatteringHeightDependent
        // is a generic scatter, ScatteringHeightIndependent is a scatterer where the legendre coefficients/cross sections
        // are not a function of altitude.
        enum SpeciesType {
            PurelyAbsorbing,
            ScatteringHeightDependent,
            ScatteringHeightIndependent
        };
    private:
        std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21> m_species;
        std::vector<SpeciesType> m_species_types;
        skBRDF* m_brdf;
        skClimatology* m_atmospheric_state;
    public:
        AtmosphereInterface() {
            m_brdf = nullptr;
            m_atmospheric_state = new skClimatology_MSIS90;
            m_atmospheric_state->AddRef();
        }

        ~AtmosphereInterface() {
            if (m_brdf != nullptr)
            {
                m_brdf->Release();
            }
            if (m_atmospheric_state != nullptr)
            {
                m_atmospheric_state->Release();
            }
        }

        void add_species(const CLIMATOLOGY_HANDLE& species, skClimatology* numberdensityclimatology, skOpticalProperties* particleopticalprops) {
            m_species.emplace_back(SKTRAN_AtmosphericOpticalStateEntry_V21(species));
            m_species.back().Configure(species, numberdensityclimatology, particleopticalprops);

            // Set the species type
            if(particleopticalprops->IsScatterer()) {
                if(particleopticalprops->IsHeightDependent()) {
                    m_species_types.push_back(ScatteringHeightDependent);
                } else {
                    m_species_types.push_back(ScatteringHeightIndependent);
                }
            } else {
                m_species_types.push_back(PurelyAbsorbing);
            }

        }

        void set_albedo(double albedo) {
            if(m_brdf != nullptr) {
                m_brdf->Release();
            }

            m_brdf = new SKTRAN_BRDF_Lambertian(albedo);
            m_brdf->AddRef();
        }

        void set_albedo(skBRDF* brdf) {
            if(m_brdf != nullptr) {
                m_brdf->Release();
            }

            m_brdf = brdf;
            m_brdf->AddRef();
        }

        void set_atmospheric_state(skClimatology* clim) {
            m_atmospheric_state->Release();
            m_atmospheric_state = clim;
            m_atmospheric_state->AddRef();
        }

        const skBRDF* brdf() {
            return m_brdf;
        }

        skClimatology* atmospheric_state() {
            return m_atmospheric_state;
        }

        std::vector<SKTRAN_AtmosphericOpticalStateEntry_V21>& species() {
            return m_species;
        }

        const std::vector<SpeciesType>& species_type() const {
            return m_species_types;
        }

    };

    class AtmosphereConstructorBase {
    public:
        virtual ~AtmosphereConstructorBase() {}

        virtual void assign_output_radiance(const double** radiance, int* numwavelens, int* numlinesofsight) = 0;
        virtual void assign_wf_buffer(const double** wf)  = 0;
    };

    template<int NSTOKES>
    class AtmosphereConstructor : public AtmosphereConstructorBase {
    private:
        AtmosphereInterface& m_atmo_interface;
        WFHandler& m_wf_handler;

        std::vector<sasktran2::atmosphere::AtmosphereGridStorageFull<NSTOKES>> m_species_quantities;
        std::vector<Eigen::VectorXd> m_species_nd;
        std::vector<GUID> m_species_guid;

        std::unique_ptr<sktran_me::OutputSKIF<NSTOKES>> m_output;
        std::unique_ptr<sasktran2::atmosphere::Atmosphere<NSTOKES>> m_atmosphere;

        int m_nwavel;
        int m_nlos;
    public:
        AtmosphereConstructor(AtmosphereInterface& atmo_interface, WFHandler& wf_handler) : m_atmo_interface(atmo_interface), m_wf_handler(wf_handler) {};

        void construct_atmospheric_state(const sasktran2::Geometry1D& geometry,
                                         const sasktran2::Config& config,
                                         const sasktran2::viewinggeometry::ViewingGeometryContainer& viewing_geo,
                                         const GEODETIC_INSTANT& refpt,
                                         const std::vector<double>& wavelengths
                                         );

        void assign_output_radiance(const double** radiance, int* numwavelens, int* numlinesofsight) override;
        void assign_wf_buffer(const double** wf) override;

        sasktran2::Output<NSTOKES>& output() { return *m_output; }
        sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere() { return *m_atmosphere; }

    };
}