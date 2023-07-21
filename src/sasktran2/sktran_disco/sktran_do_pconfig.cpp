#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_pconfig.h"
#include "sktran_disco/sktran_do_specs.h"


template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>::configureLowLevel(SKTRAN_DO_UserSpec& userspecmemory, const sasktran_disco_lowlevel::Config& config, const sasktran_disco_lowlevel::ViewingGeometry& geometry) {
    m_opticalpropertiesmutex = &m_no_mutex_given_mutex;

	m_userspec = &userspecmemory;

    const_cast<double&>(this->M_CSZ) = geometry.cos_sza;
    const_cast<double&>(this->M_SAZ) = 0;
    const_cast<double&>(this->M_SOLAR_DIRECT_INTENSITY) = 1.0;

    m_unsorted_los.clear();
    m_unsorted_los.reserve(geometry.nlos);
    for(int i = 0; i < geometry.nlos; ++i) {
        m_unsorted_los.push_back(LineOfSight());
        m_unsorted_los.back().coszenith = geometry.cos_vza[i];
        m_unsorted_los.back().azimuth = geometry.saa[i];
    }

    const_cast<uint&>(this->M_NSTR) = config.nstr;
    const_cast<uint&>(this->M_NLYR) = config.nlyr;

	userspecmemory.configure(config.nstr, config.nlyr);

    this->M_MU = m_userspec->getStreamAbscissae();
    this->M_WT = m_userspec->getStreamWeights();
    configureLP(m_userspec);

	// TODO: Do we have to do anything here for weighting functions?
    // m_perturb_quantities = m_userspec->perturbations();
	const_cast<bool&>(this->M_USE_PSEUDO_SPHERICAL) = config.usepseudospherical;

	const_cast<bool&>(this->M_USE_LOS_SPHERICAL) = false;
	const_cast<bool&>(this->M_SS_ONLY) = false;
	const_cast<size_t&>(this->M_NUM_SZA) = 1;

    m_lp_csz_storage.resize(1);
    m_lp_csz_storage[0] = std::unique_ptr<LegendrePolynomials<NSTOKES>>(new LegendrePolynomials<NSTOKES>(this->M_NSTR, this->M_CSZ));
    this->M_LP_CSZ = m_lp_csz_storage[0].get();
    m_pool.init(this->M_NLYR, this->M_NSTR);
    m_poolindex = 0;

	for (auto& lp_csz : m_lp_csz_storage) {
		for (int m = 0; m < (int)this->M_NSTR; ++m) {
			lp_csz->configureAEOrder(m);
		}
	}

}

template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>::configure(SKTRAN_DO_UserSpec &userspecmemory,
                                                                        const sasktran2::Config &config,
                                                                        double cos_sza,
                                                                        int nlyr,
                                                                        const std::vector<sasktran2::raytracing::TracedRay> &traced_rays) {
    m_opticalpropertiesmutex = &m_no_mutex_given_mutex;

    int nlos = (int)traced_rays.size();
    int nstr = config.num_do_streams();

    m_userspec = &userspecmemory;

    const_cast<double&>(this->M_CSZ) = cos_sza;
    const_cast<double&>(this->M_SAZ) = 0;
    const_cast<double&>(this->M_SOLAR_DIRECT_INTENSITY) = 1.0;

    m_unsorted_los.clear();
    m_unsorted_los.reserve(nlos);
    for(int i = 0; i < nlos; ++i) {
        m_unsorted_los.push_back(LineOfSight());
        m_unsorted_los.back().coszenith = traced_rays[i].observer_and_look.cos_viewing();;
        m_unsorted_los.back().azimuth = 0; // Never used?
    }

    const_cast<uint&>(this->M_NSTR) = config.num_do_streams();
    const_cast<uint&>(this->M_NLYR) = nlyr;

    userspecmemory.configure(nstr, nlyr);

    this->M_MU = m_userspec->getStreamAbscissae();
    this->M_WT = m_userspec->getStreamWeights();
    configureLP(m_userspec);

    // TODO: Do we have to do anything here for weighting functions?
    // m_perturb_quantities = m_userspec->perturbations();
    const_cast<bool&>(this->M_USE_PSEUDO_SPHERICAL) = true; // TODO: Get from config?

    const_cast<bool&>(this->M_USE_LOS_SPHERICAL) = false;
    const_cast<bool&>(this->M_SS_ONLY) = false;
    const_cast<size_t&>(this->M_NUM_SZA) = 1;

    m_lp_csz_storage.resize(1);
    m_lp_csz_storage[0] = std::unique_ptr<LegendrePolynomials<NSTOKES>>(new LegendrePolynomials<NSTOKES>(this->M_NSTR, this->M_CSZ));
    this->M_LP_CSZ = m_lp_csz_storage[0].get();
    m_pool.init(this->M_NLYR, this->M_NSTR);
    m_poolindex = 0;

    for (auto& lp_csz : m_lp_csz_storage) {
        for (int m = 0; m < (int)this->M_NSTR; ++m) {
            lp_csz->configureAEOrder(m);
        }
    }

}

template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>::configureModelSpecs(const SKTRAN_DO_UserSpec* userspec)
{
	m_userspec = userspec;
	const_cast<uint&>(this->M_NSTR) = m_userspec->getNumberOfStreams();
	const_cast<uint&>(this->M_NLYR) = m_userspec->getNumberOfLayers();
	this->M_MU = m_userspec->getStreamAbscissae();
	this->M_WT = m_userspec->getStreamWeights();
    configureLP(userspec);
	m_perturb_quantities = m_userspec->perturbations();
	const_cast<bool&>(this->M_USE_PSEUDO_SPHERICAL) = m_userspec->getUsePseudoSpherical();
	const_cast<SKTRAN_DO_UserSpec::LayerConstructionMethod&>(this->M_LAYER_CONSTRUCTION) = m_userspec->getLayerConstructionMethod();

	const_cast<bool&>(this->M_USE_LOS_SPHERICAL) = m_userspec->getUseLOSSpherical();
	const_cast<bool&>(this->M_SS_ONLY) = m_userspec->getSSOnly();
	const_cast<size_t&>(this->M_NUM_SZA) = m_userspec->getNumSZA();
    const_cast<double&>(this->M_SZA_REL_SEP) = m_userspec->getSZARelSep();
    const_cast<bool&>(this->M_USE_GREENS_FUNCTION) = m_userspec->getUseGreensFunction();


	if(m_userspec->getForcedNumberAzimuthTerms() > this->M_NSTR) {
		throw InvalidConfiguration("Forced number of azimuth terms must be less than or equal to the number of streams!");
	}
}

template <int NSTOKES, int CNSTR>
void sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>::configureLP(const SKTRAN_DO_UserSpec* userspec)
{
    if constexpr(NSTOKES == 1) {
        this->M_LP_MU = userspec->getAbscissaeLegendreP1();
    } else {
        this->M_LP_MU = (const VectorDim3<LegendrePhaseContainer<NSTOKES>>*) userspec->getAbscissaeLegendreP4();
    }
}


SASKTRAN_DISCO_INSTANTIATE_TEMPLATE(sasktran_disco::PersistentConfiguration);

