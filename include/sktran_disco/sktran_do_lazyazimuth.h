#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco
{

#pragma region "Lazy-Azimuth Framework"
	// A generic object which has some azimuth dependence. Since not all 
	// order of the azimuth expansion are needed, this object allows for 
	// lazy evaluation of things that depend on the order of the azimuth 
	// expansion.
	class AzimuthDependency
	{
	public:
		// Configures this object for the given azimuth order
		virtual void configureAEOrder(AEOrder m) {};

		virtual void postProcessAEOrder(AEOrder m) {};
	};

	// We need a way of triggering all azimuth dependent objects from a 
	// single call to configureAEOrder (ie. a way of 'cascading' the calls). 
	// This object allows dependencies to be registered, when configureAEOrder
	// is called it will call configureAEOrder on all of its registered 
	// dependencies.
	class AzimuthDependencyCascade: public AzimuthDependency
	{
	public:
		// Configures all registered azimuth dependencies of this object for 
		// the given azimuth expansion order.
		inline virtual void configureAEOrder(AEOrder m) override {
			for(auto ad : m_dependencies)
				ad->configureAEOrder(m);
		}

		inline virtual void postProcessAEOrder(AEOrder m) override {
			for (auto ad : m_dependencies)
				ad->postProcessAEOrder(m);
		}

		// Register an azimuthal dependency with this object. All subsequent
		// calls to ConfigureAEOrder will also call dependency.ConfigureAEOrder 
		inline void registerAzimuthDependency(AzimuthDependency& dependency) {
			m_dependencies.push_back(&dependency);
		}
	private:
		std::list<AzimuthDependency*> m_dependencies;
	};

	// Pure virtual base class which manages the lazy caching of a azimuth dependent
	// function that caches things pre solving the system
	template<class CachedDataType>
	class AzimuthDependentCache: public AzimuthDependency
	{
	public:
		AzimuthDependentCache() = delete;
		AzimuthDependentCache(uint NSTR):
			M_NSTR(NSTR),
			m_cached(M_NSTR, false),
			m_localdata(M_NSTR),
            m_data(m_localdata)
		{
		}
        AzimuthDependentCache(uint NSTR, std::vector<CachedDataType>& data):
                M_NSTR(NSTR),
                m_cached(M_NSTR, false),
                m_data(data)
        {
        }

		// Cache the requested azimuth expansion order if it has not already 
		// been calculated.
		void configureAEOrder(AEOrder m) override {
			if(m_cached[m] == false) {
				calculateAEOrder(m, m_data[m]);
				m_cached[m] = true;
			}
		}

        void reset() {
            for(int i = 0; i < m_cached.size(); ++i) {
                m_cached[i] = false;
            }
        }

		// Cached value accessor.
		inline const CachedDataType& operator[](AEOrder m) const {
			assert(m_cached[m]);
			return m_data[m];
		}
		// Calculation which is performed when an new order of the cache needs
		// to be calculated.
		virtual void calculateAEOrder(AEOrder m, CachedDataType& val) = 0;
	protected:
		const uint M_NSTR;
	private:
        std::vector<CachedDataType> m_localdata;
        std::vector<CachedDataType>& m_data;
		std::vector<bool> m_cached;
	};


	// Pure virtual base class which manages the lazy caching of a azimuth dependent
// function that caches things pre solving the system
	template<class CachedDataType>
	class AzimuthDependentPostCache : public AzimuthDependency
	{
	public:
		AzimuthDependentPostCache() = delete;
		AzimuthDependentPostCache(uint NSTR) :
			M_NSTR(NSTR),
			m_cached(M_NSTR, false),
			m_data(M_NSTR)
		{
		}
		// Cache the requested azimuth expansion order if it has not already 
		// been calculated.
		void postProcessAEOrder(AEOrder m) override {
			if (m_cached[m] == false) {
				calculateAEOrder(m, m_data[m]);
				m_cached[m] = true;
			}
		}
		// Cached value accessor.
		inline const CachedDataType& operator[](AEOrder m) const {
			assert(m_cached[m]);
			return m_data[m];
		}
		// Calculation which is performed when an new order of the cache needs
		// to be calculated.
		virtual void calculateAEOrder(AEOrder m, CachedDataType& val) = 0;
	protected:
		const uint M_NSTR;
	private:
		std::vector<CachedDataType> m_data;
		std::vector<bool> m_cached;
	};

#pragma endregion

#pragma region "Lazy-Azimuth Cache Objects"

	// Caches a Legendre polynomial evaluated at a given value for the required
	// components of the azimuth expansion.
    template <int NSTOKES, int CNSTR=-1>
	class LegendrePolynomials: public AzimuthDependentCache<std::vector<LegendrePhaseContainer<NSTOKES>>>
	{
	public:
		LegendrePolynomials(uint NSTR, double value):
			AzimuthDependentCache<std::vector<LegendrePhaseContainer<NSTOKES>>>(NSTR)
		{
			m_value = value;
		}
		virtual void calculateAEOrder(AEOrder m, std::vector<LegendrePhaseContainer<NSTOKES>>& lepolys) override final;
	private:
		double m_value;
	};

	// Caches multiScatST calculations for sasktran_disco::OpticalLayer.
	// This object build the NSTRxNSTR multiScatST matrix. For the naive 
	// calculation, this calculation is O(NSTR^3). This object improves 
	// performance by ~8x by exploiting the symmetry in the matrix and 
	// in the triple product of the Legendre polynomials.
    template <int NSTOKES, int CNSTR=-1>
	class LegendreSumMatrix: public AzimuthDependentCache<LegendreSumMatrixStorage<NSTOKES>>
	{
	public:
		LegendreSumMatrix(uint NSTR, double ssa, const VectorDim3<LegendrePhaseContainer<NSTOKES>>& lp_mu, const std::vector<LegendreCoefficient<NSTOKES>>& lpe_phasef,
                          std::vector<LegendreSumMatrixStorage<NSTOKES>>& storage):
			AzimuthDependentCache<LegendreSumMatrixStorage<NSTOKES>>(NSTR, storage),
			M_LP_MU(lp_mu),
			M_LPE_PHASEF(&lpe_phasef),
			M_SSA(ssa)
		{
			// empty
		}

        LegendreSumMatrix(uint NSTR, double ssa, const VectorDim3<LegendrePhaseContainer<NSTOKES>>& lp_mu,
                          std::vector<LegendreSumMatrixStorage<NSTOKES>>& storage):
                AzimuthDependentCache<LegendreSumMatrixStorage<NSTOKES>>(NSTR, storage),
                M_LP_MU(lp_mu),
                M_SSA(ssa)
        {
            // empty
        }

		// Override which allows for SSA to be adjusted after ctor (used to 
		// handle SSA dither)
		void adjustSSA(double ssa) {
			const_cast<double&>(M_SSA) = ssa;
		}

        void set_optical(const std::vector<LegendreCoefficient<NSTOKES>>* lpe_phasef, double ssa) {
            const_cast<double&>(M_SSA) = ssa;
            M_LPE_PHASEF = lpe_phasef;
            this->reset();
        }
		
		virtual void calculateAEOrder(AEOrder m, LegendreSumMatrixStorage<NSTOKES>& sum_matrix) override final;
	private:
		const VectorDim3<LegendrePhaseContainer<NSTOKES>>& M_LP_MU;
        const std::vector<LegendreCoefficient<NSTOKES>>* M_LPE_PHASEF;
		const double M_SSA;
        void assign(int linear_index, const TripleProductDerivativeHolder<NSTOKES>& vals, LegendreSumMatrixStorage<NSTOKES>& sum_matrix);
	};

	class AlbedoExpansion: public AzimuthDependentCache<Albedo>
	{
	public:
		AlbedoExpansion(const std::vector<LineOfSight>& los, const std::vector<double>& streams,
						double csz, std::unique_ptr<BRDF_Base> brdf, uint nterms):
			AzimuthDependentCache<Albedo>(static_cast<uint>(streams.size())),
			M_LOS(los),
			M_MU(streams),
			M_CSZ(csz),
			m_brdf(std::move(brdf)),
			m_nterms(nterms)
		{
			// empty
		}

		void injectTestingBRDF(std::unique_ptr<BRDF_Base> brdf) {
			m_brdf = std::move(brdf);
		}

		virtual void calculateAEOrder(AEOrder m, Albedo& sum_matrix) override final;

	private:
		const std::vector<LineOfSight>& M_LOS;
		const std::vector<double>& M_MU;
		const double M_CSZ;
		std::unique_ptr<BRDF_Base> m_brdf;
		uint m_nterms;

	};

}
#pragma endregion


