#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"

namespace sktran_do_detail
{
	// A few forward declarations
	class SphericalRayTracer;
	struct TracedRay;

	// Handles all configuration which is persistent through calls to 
	// calculate radiance (ie. everything that can be done after the 
	// user spec has been given). Also allows objects which need access
	// to different *Properties to be easily constructed (via copy ctor).
    template <int NSTOKES, int CNSTR=-1>
	class PersistentConfiguration: 
		public BasicProperties<NSTOKES>,
		public SolarProperties<NSTOKES>,
		public UserSpecProperties,
		public TestProperties
	{
	public:
		// Configuration members
		typedef SKTRAN_RayFactory<
			SKTRAN_RayOptical_Straight,
			SKTRAN_RayTracer_Straight_Generic,
			SKTRAN_RayStorage_Straight>							SKDO_RayFactory;

		PersistentConfiguration() { m_opticalstate_prefilled = false; }
		~PersistentConfiguration() { }
	public:
		// Configure persistent settings from user spec.
		void configureUserSpec(const SKTRAN_DO_UserSpec& userspec,
							   const SKTRAN_LineOfSightArray_V21& linesofsight,
							   std::mutex* opticalproperties_mutex,
							   std::vector<sktran_do_detail::LOSDiagnostics>* los_diag = nullptr);

		// Configure a test. Overrides standard configuration.
		void configureTest(const SKTRAN_DO_TestSpec<NSTOKES>& userspec);

		// Configure the given objects from the persistent configuration.
		// Essentially sets up the given objects for a radiance calculation.
		void configureRadianceCalculation(std::vector<SKTRAN_StokesScalar>* losradiance,
										  double wavelen,
										  SASKTRANAtmosphereInterface* opticalstate,
										  VectorDim1<LineOfSight>& los,
										  VectorDim2<double>* loswf,
										  std::unique_ptr<sktran_do_detail::BRDF_Base>& brdf,
										  OpticalStateInterface* state = nullptr
			);

        void preConfigureWavelengthTables(const std::vector<double>& wavelengths,
                                          SASKTRANAtmosphereInterface* atmosphereinterface,
                                          OpticalStateInterface* ostate
                                          );

	public: // Configuration getters
		// Returns the current coordinates.
		inline const SKTRAN_CoordinateTransform_V2* coords() const {
			return m_coordinates.get();
		}
		// Returns the number of streams.
		inline const uint nstr() const {
			return this->M_NSTR;
		}
		// Returns the number of layers.
		inline const uint nlyr() const {
			return this->M_NLYR;
		}
		// Returns a pointer to the current user specifications.
		inline const SKTRAN_DO_UserSpec* userSpec() const {
			return m_userspec;
		}
		// Returns true if we are going to use the pseudo-spherical approximation
		inline bool use_pseudo_spherical() const {
			return this->M_USE_PSEUDO_SPHERICAL;
		}

		// Returns true if we are going to use LOS sphericity corrections
		inline bool use_los_spherical() const {
			return this->M_USE_LOS_SPHERICAL;
		}

        inline bool use_greens_function() const {
            return this->M_USE_GREENS_FUNCTION;
        }

		inline bool ss_only() const {
			return this->M_SS_ONLY;
		}

        inline bool opticalstate_prefilled() const {
            return m_opticalstate_prefilled;
        }

		// Returns the method used to construct the layers in the model
		inline SKTRAN_DO_UserSpec::LayerConstructionMethod layer_construction_method() const {
			return this->M_LAYER_CONSTRUCTION;
		}
		// Local solar azimuth angle (relative to north).
		inline const double saz() const {
			return this->M_SAZ;
		}
		// Cosine of the solar zenith angle.
		inline const double csz() const {
			return this->M_CSZ;
		}

		void override_csz(size_t csz_idx) {
			*const_cast<double*>(&this->M_CSZ) = m_spher_cos_sza[csz_idx];
			this->M_LP_CSZ = m_lp_csz_storage[csz_idx].get();
            m_poolindex = csz_idx;
		}

		// Returns whether or not we are currently running a test.
		inline bool thisIsATest() const {
			return m_testing;
		}
		inline const VectorDim3<LegendrePhaseContainer<NSTOKES>>& lpStreamAngles() const {
			return *this->M_LP_MU;
		}
		inline const std::vector<std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>>* perturbation_specs() const {
			return m_perturb_quantities;
		}

		inline std::mutex& opticalPropertiesMutex() {
			return *m_opticalpropertiesmutex;
		}

		inline const std::unique_ptr<SphericalRayTracer>& ray_tracer() const {
			return m_raytracer;
		}

		inline VectorDim1<TracedRay>& traced_rays() {
			return m_traced_rays;
		}

		const std::vector<double>& spherical_cos_sza() {
			return m_spher_cos_sza;
		}

        inline MemoryPool<NSTOKES, CNSTR>& pool() const {
            return m_pool[m_poolindex];
        }


	private:	// Private configuration functions
		void configureModelSpecs(const SKTRAN_DO_UserSpec* userspec);
        void configureLP(const SKTRAN_DO_UserSpec* userspec);
		void configureRayTracing(const SKTRAN_LineOfSightArray_V21& linesofsight);

		void configureDirections(const SKTRAN_LineOfSightArray_V21& linesofsight);
		void configureSolarPosition();

		void correctLineOfSight(const SKTRAN_LineOfSightEntry_V2* los, nxVector& look, nxVector& obs);

		void fillLOSDiagnostics(std::vector<sktran_do_detail::LOSDiagnostics>* los_diag) const;

		std::unique_ptr<SphericalRayTracer>						    m_raytracer;
		VectorDim1<LineOfSight>									    m_unsorted_los;
		VectorDim1<TracedRay>									    m_traced_rays;
		std::mutex*												    m_opticalpropertiesmutex;
		std::mutex												    m_no_mutex_given_mutex;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	    m_coordinates;
        mutable std::vector<MemoryPool<NSTOKES, CNSTR>>             m_pool;
        int                                                         m_poolindex;
		std::vector<std::unique_ptr<LegendrePolynomials<NSTOKES>>>	m_lp_csz_storage;
		const std::vector<std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>>*	m_perturb_quantities;

		double													    m_minsza;
		double													    m_maxsza;

		std::vector<double>										    m_spher_cos_sza;

        bool                                                        m_opticalstate_prefilled;
	};

}
