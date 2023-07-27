#pragma once
/*!
 * @file sktran_do_specs.h
 *
 * @author Liam Bindle
 * Contact: liam.bindle@gmail.com
 *
 * @brief This file declares the interface for SASKTRAN-DO's user
 * specifications. Its implementation can be found in
 * '../sktran_do_specs.cpp'.
 *
*/
#include "sktran_disco/sktran_do.h"


namespace sasktran_disco {
    /*!
         * @class SKTRAN_DO_UserSpec
         * @ingroup SASKTRAN-DO
         *
         * @brief An object which stores the users settings for a
         * SKTRAN_DO_Engine instance.
         *
         * @throws An invalid configuration will likely throw std::invalid_argument.
         *
         */
    class SKTRAN_DO_UserSpec
    {
    public: // Constructors and destructor

        /*!
         * Default construct. Configure must be called before this object can be
         * used.
         */
        SKTRAN_DO_UserSpec();

        /*!
        * Default destructor.
        */
        virtual ~SKTRAN_DO_UserSpec() {}

    public: // Setters - typically called by user

        /*!
         * Minimum required user configuration.
         *
         * @param[in] num_streams The number of streams used in the calculation.
         * \p num_streams must be an even number greater than or equal to 4.
         * @param[in] num_atmo_layers The number of homogeneous layers used to
         * approximate the atmosphere. \p num_atmo_layers Must be greater than 0.
         * @param[in] solar_position The solar position vector.
         * @param[in] solar_direct_intensity The intensity of the sun at the TOA.
         * The default value is 1.0.
         *
         * @pre None.
         * @post This specification now meets the absolute minimum required
         * configuration for a SKTRAN_DO_Engine instance. Additional configuration
         * after is fine.
         * @warning Repeated calls can diminish performance. This is because a call
         * to this function generates some cached data.
         */
        void configure(unsigned int num_streams,
                       unsigned int num_atmo_layers,
                       double solar_direct_intensity = 1.0);

        /*!
         * Configure the user spec with the present configuration. This function
         * call requires that the object must be manually configured correctly. If
         * the spec is not configured properly a InvalidConfiguration exception
         * will be thrown. See preconditions for details on manual configuration.
         *
         * @pre The following must be set: number of streams, number of layers,
         * solar position.
         * @post None.
         * @warning Repeated call can diminish performance as some caching occurs.
         */
        void configure();

        /*!
         * Configures the user spec for more than 40 streams. User must provided
         * Gaussian quadrature points and weights.
         *
         * @param[in] gauss_quad_nodes The provided quadrature points.
         * @param[in] gauss_quad_weights The provided quadrature weights.
         *
         * @see configure
         * @see setNumberOfStreams.
         */
        void configure(const std::vector<double>& gauss_quad_nodes,
                       const std::vector<double>& gauss_quad_weights,
                       unsigned int num_atmo_layers,
                       double solar_direct_intensity = 1.0);

        /*!
         * Sets the number of streams used in the calculation. The number of
         * streams corresponds to the number of quadrature angles used to calculate
         * the numeric integral of the source function in the RTE. An increased
         * number of streams can improve precision but is also computationally more
         * expensive.
         *
         * @param[in] num_streams The number of streams used in the calculation.
         * \p num_streams must be an even number greater than or equal to 4 and
         * less than or equal to 40. If you would like to use more than 40 streams
         * you must supply your own. See setNumberOfStreams.
         *
         * @warning Using an insufficient number of streams can cause the
         * calculation to fail. In this event SKTRAN_DO_Engine will safely exit
         * and report the error.
         */
        SKTRAN_DO_UserSpec* setNumberOfStreams(unsigned int num_streams);

        /*!
        * Sets the number of streams from user supplied quadrature tables. This
        * function must be called when more than 40 streams are desired. This is
        * because SASKTRAN-DO has quadrature tables built in to but only supports
        * NSTR between 4 and 40.
        *
        * @param[in] gauss_quad_nodes The user supplied Gaussian quadrature nodes.
        * @param[in] gauss_quad_weights The corresponding Gaussian quadrature
        * weights.
        *
        * @note The number of streams is twice the size of the vectors provided.
        * For example, if you supply a table of 33 nodes and 33 weights then
        * NSTR=66.
        */
        SKTRAN_DO_UserSpec* setNumberOfStreams(const std::vector<double>& gauss_quad_nodes,
                                               const std::vector<double>& gauss_quad_weights);

        /*!
         * Sets the number of homogeneous layers used to approximate the
         * atmosphere. An increased number of layers will improve precision but
         * is also more computationally expensive.
         *
         * @param[in] num_atmo_layers The number of homogeneous layers used to
         * approximate the atmosphere. \p num_atmo_layers Must be greater than 0.
         */
        SKTRAN_DO_UserSpec* setNumberOfLayers(unsigned int num_atmo_layers);

        /*!
         * Sets the incident intensity at the top of the atmosphere.
         *
         * @param[in] solar_direct_intensity The intensity of the sun at the TOA.
         * default value is 0.0.
         */
        SKTRAN_DO_UserSpec* setTOAIntensities(double direct = 1);

        /*!
         * Set the altitude grid where SASKTRAN will poll optical properties
         * before  caching.
         *
         *	@param[in] altitudes The altitudes in meters, and in ascending order.
         *	The first element in this array will become the ground altitude and the
         *	last element will become the TOA altitude.
         */
        SKTRAN_DO_UserSpec* setAltitudeGrid(const std::vector<double>& altitudes);

        /*!
        * Sets the optical properties radii grid size. For the
        * SKTRAN_GridDefOpticalPropertiesRadii_V21 instance.
        *
        *  @param[in] size_of_grid The size of the grid.
        */
        SKTRAN_DO_UserSpec* setAltitudeGrid(double ground_altitude, double altitude, size_t size);

        /*!
         * Valid number of terms used in the numeric integration of the phase
         * function during its expansion.
         */
        enum class NTERMS_IN_EXPANSION
        {
            u64 = 64, u128 = 128, u256 = 256, u512 = 512, u1024 = 1024
        };

        /*!
        * Sets the number of terms used in the numeric integration of the BRDF
        * during its expansion. Default value is 64.
        *
        * @param[in] nterms The number of terms.
        */
        SKTRAN_DO_UserSpec* setNumBRDFExpansionTerms(NTERMS_IN_EXPANSION nterms = NTERMS_IN_EXPANSION::u64);

        /*!
         * Sets the position of the sun.
         *
         * @param[in] sun The solar position vector.
         */

        SKTRAN_DO_UserSpec* setUsePsuedoSpherical(bool use_ps = true);

        SKTRAN_DO_UserSpec* setUseLOSSpherical(bool use_los_spher = false);

        SKTRAN_DO_UserSpec* setUseUpwellingSpherical(bool use_upwelling_spher = false);

        SKTRAN_DO_UserSpec* setNumSZA(size_t num_sza = 2);

        SKTRAN_DO_UserSpec* setSSOnly(bool ss_only = false);

        SKTRAN_DO_UserSpec* setSurfaceEmissionWavelengths(const std::vector<double>& wavelengths);

        SKTRAN_DO_UserSpec* setSurfaceEmissionValues(const std::vector<double>& emissions);


        /*!
         * Sets the calculations convergence criterion. A calculation is said to
         * have converged once it has been observed twice that the relative
         * contribution of an additional order of the azimuth expansion is less
         * than \p epsilon.
         *
         * @param[in] epsilon The relative upper bound of a contribution which is
         * said to have converged.
         */
        SKTRAN_DO_UserSpec* setCauchyCriterion(double epsilon = 0.01);

        /*!
         * Sets the single-scatter dither amount in the event a purely scattering
         * atmospheric layer is encountered. This is a special case which discrete-
         * ordinate algorithms must handle by dithering the SSA.
         *
         * @param[in] dither The amount to dither the SSA by.
         */
        SKTRAN_DO_UserSpec* setSSAEqual1Dither(double dither = 1e-9);

        /*!
         * Forms in which the engine can return weighting functions.
         */
        enum class WeightingFunctionForm
        {
            dI_dLogX,
            dI_dX
        };

        /*!
         * Forms in which the layers can be configured
         */
        enum class LayerConstructionMethod
        {
            uniform_pressure,
            uniform_height,
            match_altitudegrid,
            manual
        };

        enum class LineOfSightAdjustment
        {
            translate_observer,
            match_target_angles
        };

        SKTRAN_DO_UserSpec* setLineOfSightAdjustment(LineOfSightAdjustment method = LineOfSightAdjustment::translate_observer) {
            m_los_adjustment = method;

            return this;
        }

        SKTRAN_DO_UserSpec* setLayerConstructionMethod(LayerConstructionMethod method = LayerConstructionMethod::uniform_pressure);

        SKTRAN_DO_UserSpec* setManualLayerAltitudes(const std::vector<double>& altitudes);

        SKTRAN_DO_UserSpec* setSZARelSep(double rel_sep = 0);

        SKTRAN_DO_UserSpec* setUseGreensFunction(bool use = true);

        /*!
         * Sets the form in which the engine will return the weighting functions.
         *
         * @param[in] form: The form in which to return the weighting function. The
         * default value is dI_dX.
         */
        SKTRAN_DO_UserSpec* setWFReturnForm(WeightingFunctionForm form = WeightingFunctionForm::dI_dX);

        /*!
         * Instantiable base class for weighting function specifications. A vector
         * of this object is what needs to be passed to setWFSpecies. Note that
         * this class is just a concrete base class and doesn't actually do
         * anything.
         */
        struct WeightingFunctionSpec
        {
            // Enum over all types so that we can reinterpret cast instead of dynamic cast the WF
            enum WeightingFunctionType {
                SpeciesWF,
                AlbedoWF,
                TestWF
            };

            virtual ~WeightingFunctionSpec() {};

            virtual WeightingFunctionType type() = 0;
        };

        /*!
         * Specifies a weighting function calculation with respect to ground albedo.
         */
        struct AlbedoWF: WeightingFunctionSpec
        {
            /*!
             * Constructs a ground albedo weighting function.
             *
             * @param[in] p_eps The relative difference to ground albedo used for
             * calculating the numeric derivative.
             */
            AlbedoWF()
            {
                // empty
            }

            WeightingFunctionType type() override { return WeightingFunctionSpec::AlbedoWF; }

        };

        typedef std::vector<std::unique_ptr<WeightingFunctionSpec>> VectorOfUPtrWeightingFunctionSpecs;

        /*!
        * Sets the weighting function specifications. Weighting functions will be
        * calculated for the given vector of weighting function specs.
        *
        * @param[in] wf_spec A vector of weighting function specifications
        */
        SKTRAN_DO_UserSpec* setWFSpecies(const std::vector<std::unique_ptr<WeightingFunctionSpec>>& specs) {
            m_ptrbs = &specs;
            return this;
        }

        /*!
         * Forces the engine to use the specified number of terms in the azimuth
         * expansion. Set to zero to behave normally.
         *
         * @param[in] num_terms Number of terms to use in the azimuth expansion.
         */
        SKTRAN_DO_UserSpec* setForceNumberAzimuthTerms(uint num_terms = 0);


    public: /*! @cond */ // Dev setters - typically called during development
        #if 0
        SKTRAN_DO_UserSpec* loadDISORT(LPCWSTR filename, LPCSTR procname, double tol = 1e-4);
        sasktran_disco::DISORT_FPTR getDISORT_FPTR() const;
        bool compareWithDISORT() const;
        double getDISORTTolerance() const;
            // DISORT Comparison Settings
        double m_disort_tol;
        sasktran_disco::DISORT_FPTR m_disort_fptr;
        #endif

    public: // Getters - typically called by the internals of SASKTRAN-DO
        uint getNumberOfStreams() const;
        const VectorDim1<double>* getStreamAbscissae() const;
        const VectorDim1<double>* getStreamWeights() const;
        const VectorDim3<LegendrePhaseContainer<4>>* getAbscissaeLegendreP4() const;
        const VectorDim3<LegendrePhaseContainer<1>>* getAbscissaeLegendreP1() const;
        uint getNumberOfLayers() const;
        double getTopDirectIntensity() const;
        double getTopAltitude() const;
        double getBottomAltitude() const;
        bool getWFFormIsByLogX() const;
        uint getNumPhaseFQuadratureTerms() const;
        uint getNumBRDFQuadratureTerms() const;
        const std::vector<std::unique_ptr<SKTRAN_DO_UserSpec::WeightingFunctionSpec>>* perturbations() const;
        double getCCEpsilon() const;
        double getSSAEqual1Dither() const;
        void configureDefaultDetails();
        bool ok() const { return m_ok; }
        void cacheLPOfStreamAngles();
        uint getForcedNumberAzimuthTerms() const {
            return m_forced_number_azimuth_terms;
        }
        bool getUsePseudoSpherical() const {
            return m_use_psuedo_spherical;
        }
        bool getUseLOSSpherical() const {
            return m_use_los_spherical;
        }
        bool getUseUpwellingSpherical() const {
            return m_use_upwelling_spher;
        }
        size_t getNumSZA() const {
            return m_num_sza;
        }
        bool getSSOnly() const {
            return m_ss_only;
        }
        LayerConstructionMethod getLayerConstructionMethod() const {
            return m_layer_construction_method;
        }
        const std::vector<double>& getAltitudeGrid() const {
            return m_altitude_grid_storage;
        }
        const std::vector<double>& getManualLayerAltitudes() const {
            return m_manual_layer_heights;
        }
        const std::vector<double>& getSurfaceEmissionWavelengths() const {
            return m_surf_emission_wavelengths;
        }
        const std::vector<double>& getSurfaceEmissionValues() const {
            return m_surf_emission_values;
        }

        const LineOfSightAdjustment getLosAdjustmentMethod() const {
            return m_los_adjustment;
        }

        const double getSZARelSep() const {
            return m_rel_sza_sep;
        }

        const bool getUseGreensFunction() const {
            return m_use_greens_function;
        }

    private: // Members
        bool m_ok;

        // Gaussian Quadrature
        // Just store the full stokes vector components, shouldn't be that much extra calculation
        VectorDim3<LegendrePhaseContainer<4>> m_lp_abscissae4;
        VectorDim3<LegendrePhaseContainer<1>> m_lp_abscissae1;

        VectorDim1<double> m_abscissae;
        VectorDim1<double> m_weights;

        // Incident Intensity Specs
        double m_itop_direct;

        std::vector<double> m_altitude_grid_storage;

        // Discrete-Ordinate Method Configuration
        uint m_nstr;
        uint m_nlyr;

        uint m_forced_number_azimuth_terms;

        LayerConstructionMethod m_layer_construction_method;
        std::vector<double> m_manual_layer_heights;

        // Number of terms used in the expansion of the brdf
        uint m_brdf_quad_terms;

        // Particular solution method
        bool m_use_greens_function;

        // Pseudo Spherical Options
        bool m_use_psuedo_spherical;

        // LOS Sphericity Options
        bool m_use_los_spherical;
        size_t m_num_sza;
        bool m_ss_only;
        bool m_use_upwelling_spher;

        LineOfSightAdjustment m_los_adjustment;

        // Surface emission values
        std::vector<double> m_surf_emission_wavelengths;
        std::vector<double> m_surf_emission_values;

        // Convergence criterion
        double m_cc_epsilon;

        // Dither for handling single-scatter albedo = 1 special case
        double m_ssalb1_dither;

        // Check if SZA Is within some relative fraction of the separation constants
        double m_rel_sza_sep;

        // Flag to indicate whether to return
        bool m_return_wf_dI_dLogX;

        // Weighting-function computations
        const std::vector<std::unique_ptr<WeightingFunctionSpec>>* m_ptrbs;
        /*! @endcond */
    };
}
