#pragma once
#include "sktran_disco/sktran_do.h"
#include "sktran_do_specs.h"
#include <sasktran2/config.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/geometry.h>

namespace sasktran_disco
{
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
		PersistentConfiguration() { m_opticalstate_prefilled = false; }
		~PersistentConfiguration() { }
	public:
        void configureLowLevel(SKTRAN_DO_UserSpec& userspecmemory, const sasktran_disco_lowlevel::Config& config, const sasktran_disco_lowlevel::ViewingGeometry& geometry);
        void configure(SKTRAN_DO_UserSpec& userspecmemory,
                       const sasktran2::Config& config,
                       double cos_sza,
                       int nlyr,
                       const std::vector<sasktran2::raytracing::TracedRay>& traced_rays
                       );


	public: // Configuration getters
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
            m_poolindex = (int)csz_idx;
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

		const std::vector<double>& spherical_cos_sza() {
			return m_spher_cos_sza;
		}

        inline MemoryPool<NSTOKES, CNSTR>& pool() const {
            return m_pool;
        }

        inline const VectorDim1<double>* quadrature_weights() const { return this->M_WT; }
        inline const VectorDim1<double>* quadrature_cos_angle() const { return this->M_MU; }


	protected:	// Private configuration functions
		void configureModelSpecs(const SKTRAN_DO_UserSpec* userspec);
        void configureLP(const SKTRAN_DO_UserSpec* userspec);

		VectorDim1<LineOfSight>									    m_unsorted_los;
		std::mutex*												    m_opticalpropertiesmutex;
		std::mutex												    m_no_mutex_given_mutex;
        mutable MemoryPool<NSTOKES, CNSTR>                          m_pool;
        int                                                         m_poolindex;
		std::vector<std::unique_ptr<LegendrePolynomials<NSTOKES>>>	m_lp_csz_storage;
		const std::vector<std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>>*	m_perturb_quantities;

		double													    m_minsza;
		double													    m_maxsza;

		std::vector<double>										    m_spher_cos_sza;

        bool                                                        m_opticalstate_prefilled;
	};

}
