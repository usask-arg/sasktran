#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_specs.h"
#include "sktran_disco/sktran_do_quadrature.h"


sasktran_disco::SKTRAN_DO_UserSpec::SKTRAN_DO_UserSpec()
{
	configureDefaultDetails();
}

void sasktran_disco::SKTRAN_DO_UserSpec::configureDefaultDetails()
{
	m_ok = true;
	#if 0
	m_disort_fptr = nullptr;
	#endif
	m_nstr = 0;
	m_nlyr = 0;
	setAltitudeGrid(0, 100e3, 201);
	setNumBRDFExpansionTerms();
	setCauchyCriterion();
	setSSAEqual1Dither();
	setWFReturnForm();
	setTOAIntensities();
	setUsePsuedoSpherical();
	setUseUpwellingSpherical();
	setNumSZA();
	setUseLOSSpherical();
	setSSOnly();
	setLayerConstructionMethod();
	setForceNumberAzimuthTerms();
	setLineOfSightAdjustment();
    setSZARelSep();
    setUseGreensFunction();
	m_ptrbs = nullptr;
}

#pragma region "Setters"

void sasktran_disco::SKTRAN_DO_UserSpec::configure(unsigned int num_streams,
                                                   unsigned int num_atmo_layers,
                                                   double solar_direct_intensity)
{
	setNumberOfStreams(num_streams);
	cacheLPOfStreamAngles();
	setNumberOfLayers(num_atmo_layers);
	setTOAIntensities(solar_direct_intensity);
}

void sasktran_disco::SKTRAN_DO_UserSpec::configure(const std::vector<double>& gauss_quad_nodes,
                                                   const std::vector<double>& gauss_quad_weights,
                                                   unsigned int num_atmo_layers,
                                                   double solar_direct_intensity)
{
	setNumberOfStreams(gauss_quad_nodes, gauss_quad_weights);
	cacheLPOfStreamAngles();
	setNumberOfLayers(num_atmo_layers);
	setTOAIntensities(solar_direct_intensity);
}

void sasktran_disco::SKTRAN_DO_UserSpec::configure()
{
	if(!m_nstr) {
		throw sasktran_disco::InvalidConfiguration("Number of streams has not been set!");
	}
	if(!m_nlyr) {
		throw sasktran_disco::InvalidConfiguration("Number of layers has not been set!");
	}
	cacheLPOfStreamAngles();
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setNumberOfStreams(unsigned int num_streams)
{
	using namespace sasktran_disco;
	// check that number of streams is GE 4 and even
	if(num_streams < 2 || num_streams % 2 == 1 || num_streams > 40) {
		throw InvalidConfiguration("Number of streams must be an even number in the range [2, 40]!");
	}
	m_nstr = (uint) num_streams;
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setUsePsuedoSpherical(bool use_ps) {
	m_use_psuedo_spherical = use_ps;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setUseLOSSpherical(bool use_los_spher) {
	m_use_los_spherical = use_los_spher;
	
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setUseUpwellingSpherical(bool use_upwelling_spher) {
	m_use_upwelling_spher = use_upwelling_spher;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setNumSZA(size_t num_sza) {
	m_num_sza = num_sza;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setSSOnly(bool ss_only) {
	m_ss_only = ss_only;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setSurfaceEmissionWavelengths(const std::vector<double>& wavelengths) {
	m_surf_emission_wavelengths = wavelengths;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setSurfaceEmissionValues(const std::vector<double>& emissions) {
	m_surf_emission_values = emissions;

	return this;
}


void sasktran_disco::SKTRAN_DO_UserSpec::cacheLPOfStreamAngles()
{
	using namespace sasktran_disco;
	// get stream angles ang weights from tables
	getStreamsAndWeights(m_nstr, m_abscissae, m_weights);
	// cache legendre polynomials evaluated at stream angles
	m_lp_abscissae4.resize(m_nstr, VectorDim2<sasktran_disco::LegendrePhaseContainer<4>>(m_nstr, VectorDim1<sasktran_disco::LegendrePhaseContainer<4>>(m_nstr)));
    m_lp_abscissae1.resize(m_nstr, VectorDim2<sasktran_disco::LegendrePhaseContainer<1>>(m_nstr, VectorDim1<sasktran_disco::LegendrePhaseContainer<1>>(m_nstr)));


	for(AEOrder m = 0; m < m_nstr; ++m) {
        auto calculatorP = sasktran_disco::WignerDCalculator(m, 0);
        auto calculatorneg = sasktran_disco::WignerDCalculator(m, -2);
        auto calculatorpos = sasktran_disco::WignerDCalculator(m, 2);

        for(LPOrder l = 0; l < m_nstr; ++l) {
		    for(StreamIndex i = 0; i < m_nstr; ++i) {
                double theta = acos(m_abscissae[i]);

                m_lp_abscissae4[m][i][l].P() = calculatorP.d(theta, l);
                m_lp_abscissae1[m][i][l].value = calculatorP.d(theta, l);

                m_lp_abscissae4[m][i][l].R() = -0.5 * (calculatorpos.d(theta, l) + calculatorneg.d(theta, l));
                m_lp_abscissae4[m][i][l].T() = -0.5 * (calculatorpos.d(theta, l) - calculatorneg.d(theta, l));

                if( m % 2 == 1) {
                    //m_lp_abscissae4[m][i][l].R() *= -1;
                    //m_lp_abscissae4[m][i][l].T() *= -1;
                }
			}
		}
	}
}


sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setNumberOfStreams(const std::vector<double>& quadrature_nodes,
                                                                                           const std::vector<double>& quadrature_weights)
{
	using namespace sasktran_disco;
	getStreamsAndWeights(quadrature_nodes, quadrature_weights, m_abscissae, m_weights);
	m_nstr = (uint) m_abscissae.size();
	if(m_nstr == 0) {
		m_ok = false;
	}
    cacheLPOfStreamAngles();

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setNumberOfLayers(unsigned int num_atmo_layers)
{
	using namespace sasktran_disco;
	if(num_atmo_layers < 1) {
		throw InvalidConfiguration("Number of layers must be greater than zero!");
	} else {
		m_nlyr = num_atmo_layers;
	}
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setLayerConstructionMethod(SKTRAN_DO_UserSpec::LayerConstructionMethod method) {
	m_layer_construction_method = method;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setManualLayerAltitudes(const std::vector<double>& altitudes) {
	m_manual_layer_heights = altitudes;

	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setSZARelSep(double rel_sep) {
    m_rel_sza_sep = rel_sep;

    return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setUseGreensFunction(bool use) {
    m_use_greens_function = use;

    return this;
}



sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setTOAIntensities(double direct)
{
	using namespace sasktran_disco;
	if(direct < 0) {
		throw InvalidConfiguration("Radiances at the top of the atmosphere cannot be negative!");
	} else {
		m_itop_direct = direct;
	}
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setAltitudeGrid(const std::vector<double>& altitude_grid)
{
	using namespace sasktran_disco;
	if(altitude_grid.front() >= altitude_grid.back()) {
		throw InvalidConfiguration("The altitude grid must be given in ascending order and the ground altitude must be less than the TOA altitude.");
	}
	m_altitude_grid_storage.resize(altitude_grid.size());
	std::copy(altitude_grid.cbegin(), altitude_grid.cend(), m_altitude_grid_storage.begin());
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setAltitudeGrid(double ground_altitude, double toa_altitude, size_t size_of_grid)
{
	using namespace sasktran_disco;
	if(size_of_grid == 0) {
		throw InvalidConfiguration("Size of optical properties grid must be greater than zero!");
	} else {
		m_altitude_grid_storage.resize(size_of_grid);
		for(size_t i = 0; i < size_of_grid; ++i) {
			m_altitude_grid_storage[i] = i * (toa_altitude - ground_altitude) / (size_of_grid - 1);
		}
	}
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setNumBRDFExpansionTerms(NTERMS_IN_EXPANSION nterms)
{
	m_brdf_quad_terms = static_cast<sasktran_disco::uint>(nterms);
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setForceNumberAzimuthTerms(sasktran_disco::uint num_terms) {
	m_forced_number_azimuth_terms = num_terms;
	return this;
}

#if 0
SKTRAN_DO_UserSpec* SKTRAN_DO_UserSpec::loadDISORT(LPCWSTR filename, LPCSTR procname, double tol)
{
	sasktran_disco::LoadDISORTDLL(filename, procname, &m_disort_fptr);
	m_disort_tol = tol;
	return this;
}

sasktran_disco::DISORT_FPTR SKTRAN_DO_UserSpec::getDISORT_FPTR() const
{
	return m_disort_fptr;
}

bool SKTRAN_DO_UserSpec::compareWithDISORT() const
{
	return m_disort_fptr != nullptr;
}

double SKTRAN_DO_UserSpec::getDISORTTolerance() const
{
	return m_disort_tol;
}

#endif

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setCauchyCriterion(double epsilon)
{
	if(epsilon > 0) {
		m_cc_epsilon = epsilon;
	} else {
		using namespace sasktran_disco;
		throw InvalidConfiguration("Convergence criterion must be greater than zero!");
		m_ok &= false;
	}
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setSSAEqual1Dither(double dither)
{
	if(dither <= 0) {
		using namespace sasktran_disco;
		throw InvalidConfiguration("SSA is 1 dither must greater than zero!");
	}
	m_ssalb1_dither = dither;
	return this;
}

sasktran_disco::SKTRAN_DO_UserSpec* sasktran_disco::SKTRAN_DO_UserSpec::setWFReturnForm(
        sasktran_disco::SKTRAN_DO_UserSpec::WeightingFunctionForm form)
{
	m_return_wf_dI_dLogX = form == WeightingFunctionForm::dI_dLogX;
	return this;
}

const std::vector<std::unique_ptr<sasktran_disco::SKTRAN_DO_UserSpec::WeightingFunctionSpec>>*  sasktran_disco::SKTRAN_DO_UserSpec::perturbations() const
{
	return m_ptrbs;
}

#pragma endregion

#pragma region "Getters"
sasktran_disco::uint sasktran_disco::SKTRAN_DO_UserSpec::getNumberOfStreams() const
{
	return m_nstr;
}

sasktran_disco::uint sasktran_disco::SKTRAN_DO_UserSpec::getNumberOfLayers() const
{
	return m_nlyr;
}

double sasktran_disco::SKTRAN_DO_UserSpec::getTopDirectIntensity() const
{
	return m_itop_direct;
}

double sasktran_disco::SKTRAN_DO_UserSpec::getTopAltitude() const
{
	return m_altitude_grid_storage[m_altitude_grid_storage.size() - 1];
}

double sasktran_disco::SKTRAN_DO_UserSpec::getBottomAltitude() const
{
	return m_altitude_grid_storage[0];
}

const sasktran_disco::VectorDim1<double>* sasktran_disco::SKTRAN_DO_UserSpec::getStreamAbscissae() const
{
	return &m_abscissae;
}
const sasktran_disco::VectorDim1<double>* sasktran_disco::SKTRAN_DO_UserSpec::getStreamWeights() const
{
	return &m_weights;
}
const sasktran_disco::VectorDim3<sasktran_disco::LegendrePhaseContainer<4>>* sasktran_disco::SKTRAN_DO_UserSpec::getAbscissaeLegendreP4() const
{
	return &m_lp_abscissae4;
}

const sasktran_disco::VectorDim3<sasktran_disco::LegendrePhaseContainer<1>>* sasktran_disco::SKTRAN_DO_UserSpec::getAbscissaeLegendreP1() const
{
    return &m_lp_abscissae1;
}


sasktran_disco::uint sasktran_disco::SKTRAN_DO_UserSpec::getNumBRDFQuadratureTerms() const
{
	return m_brdf_quad_terms;
}

double sasktran_disco::SKTRAN_DO_UserSpec::getCCEpsilon() const
{
	return m_cc_epsilon;
}

double sasktran_disco::SKTRAN_DO_UserSpec::getSSAEqual1Dither() const
{
	return m_ssalb1_dither;
}

bool sasktran_disco::SKTRAN_DO_UserSpec::getWFFormIsByLogX() const
{
	return m_return_wf_dI_dLogX;
}

#pragma endregion
