#pragma once
#include "modules/sktran_do_deprecated/include/sktran_do.h"

/**
 * Sasktran Discrete-Ordinates model.
 */
class SKTRAN_DO_EngineInterface {
public:
    virtual void calculateRadiance(std::vector<double>& los_radiances,
                           double wavelength,
                           sktran_do_detail::SASKTRANAtmosphereInterface* opticalstate,
                           std::vector<std::vector<double>>* los_wf = nullptr,
                           sktran_do_detail::RTSDiagnostics* diag = nullptr,
                           std::vector<double>* los_transmission = nullptr,
                           std::vector<double>* layer_od = nullptr,
                           std::vector<double>* layer_ssa = nullptr,
                           int wavelidx = -1
                           ) = 0;
};


template <int NSTOKES, int CNSTR=-1>
class SKTRAN_DO_Engine : public SKTRAN_DO_EngineInterface
{
public:
	/*!
	 * Default ctor.
	 */
	SKTRAN_DO_Engine() { }
	~SKTRAN_DO_Engine() { }

public: // Recommended interface

	/*!
	 * Calculates the scalar radiance at the given wavelength with the given 
	 * optical state. Optionally computes weighting functions as well.
	 * @param[out] los_radiances The calculated scalar results.
	 * @param[in] wavelength The wavelength.
	 * @param[in] opticalstate The current optical state.
	 * @param[out] los_wf The computed weighting functions for each line of 
	 * sight.
	 * 
	 * @pre Engine must be configure. See configureModel.
	 * @post None.
	 */
	void calculateRadiance(std::vector<double>& los_radiances,
						   double wavelength,
						   sktran_do_detail::SASKTRANAtmosphereInterface* opticalstate,
						   std::vector<std::vector<double>>* los_wf = nullptr,
						   sktran_do_detail::RTSDiagnostics* diag = nullptr,
						   std::vector<double>* los_transmission = nullptr,
                           std::vector<double>* layer_od = nullptr,
                           std::vector<double>* layer_ssa = nullptr,
                           int wavelidx = -1
                           );

	/*!
	 * Configures the engine with the user's specifications as well as desired
	 * lines of sight.
	 * @param[in] spec The user's specifications.
	 * @param[in] los The lines of sight the user would like calculations done 
	 * for.
	 * 
	 * @pre None.
	 * @post Calculations can be made. See calculateRadiance.
	 */
	void configureModel(const SKTRAN_DO_UserSpec& spec,
						const SKTRAN_LineOfSightArray_V21& los,
						std::mutex* opticalproperties_mutex = nullptr,
						std::vector<sktran_do_detail::LOSDiagnostics>* los_diag = nullptr);

    void prefillWavelengthTables(const std::vector<double>& wavelengths,
                                 sktran_do_detail::SASKTRANAtmosphereInterface*	atmosphereinterface,
                                 std::mutex* opticalproperties_mutex = nullptr);

	/*!
	 * Recommend a number of threads to be used when performing calculations.
	 */
	static void setNumThreads(int nth) {
		#ifdef SKTRAN_DO_USE_MKL
		mkl_set_num_threads(nth);
		#endif
	}

public: // Testing
	/*! @cond */
	void configureTest(sktran_do_detail::SKTRAN_DO_TestSpec<NSTOKES>& testspec);
	static void runTest(const sktran_do_detail::testing::TestCase<NSTOKES>& testcase, std::vector<double>& abs_diff);
	/*! @endcond */
public: // SASKTRAN required interface
	/*!
	 * See configureModel.
	 */
	virtual bool ConfigureModel(SKTRAN_SpecsUser_Base&				modelspecifications,
								const SKTRAN_LineOfSightArray_V21&  linesofsight,
								size_t								numthreads)
	{
		setNumThreads(static_cast<int>(numthreads));
		configureModel(*dynamic_cast<SKTRAN_DO_UserSpec*>(&modelspecifications), linesofsight);
		return true;
	}
	/*!
	* See calculateRadiance.
	*/
	virtual bool CalculateRadiance(std::vector<SKTRAN_StokesScalar>*	losradiance,
								   double								wavelen,
								   size_t								numordersofscatter,
								   sktran_do_detail::SASKTRANAtmosphereInterface*	opticalstate,
								   std::vector<skRTStokesVector>*		losvector = nullptr,
								   bool									updateclimatology = false,
								   SKTRAN_DiagnosticInterface*			diag = NULL)
	{
		calculateRadiance(*losradiance, wavelen, opticalstate);
		return true;
	}
public: // Dev stuff
	/*! @cond */
	#if 0
	static bool generateTestResults(const sktran_do_detail::testing::TestCase& testcase);
	bool CompareDISORT(std::vector<SKTRAN_StokesScalar>*	losradiance,
					   double								wavelen,
					   size_t								numordersofscatter,
					   SKTRAN_AtmosphericOpticalState_V21*	opticalstate);
	#endif
	/*! @endcond */
private:
	// Configuration of SKDO
	sktran_do_detail::PersistentConfiguration<NSTOKES, CNSTR>		m_config;
    std::unique_ptr<sktran_do_detail::GeometryLayerArray<NSTOKES, CNSTR>> m_geometrylayers;
	std::unique_ptr<sktran_do_detail::SurfaceEmission> m_surfaceemission;
    sktran_do_detail::OpticalState<NSTOKES> m_opticalstate;
    std::vector<std::unique_ptr<sktran_do_detail::LegendrePolynomials<NSTOKES>>>	m_lp_csz_storage;


    void configureRTSDiagnostics(sktran_do_detail::RTSDiagnostics* diag, std::vector<sktran_do_detail::LineOfSight>& linesofsight) const;

	void calculateRadiancePlaneParallel(std::vector<double>& los_radiances,
		double wavelength,
		sktran_do_detail::SASKTRANAtmosphereInterface* opticalstate,
		std::vector<std::vector<double>>* los_wf = nullptr,
		sktran_do_detail::RTSDiagnostics* diag = nullptr,
		std::vector<double>* los_transmission = nullptr,
        std::vector<double>* layer_od = nullptr,
        std::vector<double>* layer_ssa = nullptr,
        int wavelidx = -1);

	void calculateRadianceSpherical(std::vector<double>& los_radiances,
		double wavelength,
		sktran_do_detail::SASKTRANAtmosphereInterface* opticalstate,
		std::vector<std::vector<double>>* los_wf = nullptr,
		sktran_do_detail::RTSDiagnostics* diag = nullptr,
		std::vector<double>* los_transmission = nullptr,
        std::vector<double>* layer_od = nullptr,
        std::vector<double>* layer_ssa = nullptr,
        int wavelidx = -1
        );

	void accumulateStokes(double* ptr, const sktran_do_detail::Radiance<NSTOKES>& radiance);
	void assignStokesDeriv(std::vector<double>* wf, const sktran_do_detail::Radiance<NSTOKES>& radiance, int num_ptrb);
};
