#pragma once


namespace sasktran2 {


    /** User Configuration object for the SASKTRAN2 model, all values have default options.
     */
	class Config {

    public:
        /** Enum determining the type of single scatter source used within the model.
         *
         * 'exact' causes solar rays to be traced for each integration point along the line of sight
         *
         * 'solartable' Creates a solar transmission table at a set of grid points.  NOT YET IMPLEMENTED.
         *
         * 'none' Removes the single scatter source from the integration.
         *
         */
        enum class SingleScatterSource {
            exact,
            solartable,
            none
        };

        /** Enum determining the type of multiple scatter source to be included within the model.
         *
         *  'discrete_ordinates' Uses the discrete ordinates method to calculate the multiple scatter source.
         *
         *  'hr' Uses a successive orders of scattering method to calculate the multiple scatter source.
         *
         *  'none' Removes the multiple scatter source from the calculation.
         *
         */
        enum class MultipleScatterSource {
            discrete_ordinates,
            hr,
            none
        };

        /** Enum that determines the accuracy of the weighting function solution within the model.  The exact effect
         *  of each level is determined primarily by the source functions included, however as a rough rule of thumb:
         *
         *  'full` Includes all effects of the weighting function calculation so that the weighting function is exact
         *         to machine precision
         *
         *  'reduced' Allows source terms to neglect hard to calculate weighting function contributions to improve computational speed.
         *
         *
         *  'limited' Includes only core effects of the weighting function, typically just direct line of sight terms and no multiple scatter
         *
         */
        enum class WeightingFunctionPrecision {
            full,
            reduced,
            limited
        };

        /** Enum determining the type of occulation source to include within the model.
         *
         *  'standard' Assumes that each ray looks directly at the sun, including a source of 1 at the end of each LOS that does not intersect the ground.
         *
         *  'none' Removes the occultation source.
         */
        enum class OccultationSource {
            standard,
            none
        };

        Config();

        /**
         *
         * @return The number of threads used for the calculation
         */
        int num_threads() const {return m_nthreads;}

        /** Sets the maximum number of threads used in the calculation.
         *
         * @param nthreads
         */
        void set_num_threads(int nthreads) { m_nthreads = nthreads; }

        /**
         *
         * @return The number of stokes parameters to calculate
         */
        int num_stokes() const { return m_nstokes;}

        /** Sets the number of stokes parameters to calculate, currently only 1 and 3 are supported.
         *
         * @param nstokes
         */
        void set_num_stokes(int nstokes) { m_nstokes = nstokes; }

        /**
         *
         * @return The type of single scatter source to include
         */
        SingleScatterSource single_scatter_source() const { return m_single_scatter_source; }
        /** Sets the type of single scatter source to include
         *
         * @param source
         */
        void set_single_scatter_source(SingleScatterSource source) { m_single_scatter_source = source; }

        /**
         *
         * @return The type of multiple scatter source to include
         */
        MultipleScatterSource multiple_scatter_source() const { return m_multiple_scatter_source; }
        /** Sets the type of multiple scatter source to include
         *
         * @param source
         */
        void set_multiple_scatter_source(MultipleScatterSource source) { m_multiple_scatter_source = source; }

        /**
         *
         * @return The type of occultation source to include
         */
        OccultationSource occultation_source() const { return m_occultation_source; }

        /** Sets the type of occultation source to include
         *
         * @param source
         */
        void set_occultation_source(OccultationSource source) { m_occultation_source = source; }

        /**
         *
         * @return The number of DO streams (full space) to use in the calculation
         */
        int num_do_streams() const { return m_ndostreams; }

        /** Sets the number of DO streams to use in the calculation.  Only used when multiple scatter is activated.
         *  Applies to both the DO and HR multiple scattering sources.
         *
         * @param nstr
         */
        void set_num_do_streams(int nstr)  { m_ndostreams = nstr;}

        /**
         *
         * @return The number of legendre moments to include in the single scatter calculation
         */
        int num_singlescatter_moments() const { return m_nsinglescatter_moments; }

         /**
          *
          * @param moments The number of legendre moments to use in the single scatter calculation
          */
        void set_num_singlescatter_moments(int moments) { m_nsinglescatter_moments = moments; }

        /**
         *
         * @return True if delta-m scaling is enabled
         */
        bool apply_delta_scaling() const { return m_apply_delta_scaling; }

        /**
         *
         * @param scale Set to true to enable delta scaling
         */
        void set_apply_delta_scaling(bool scale) { m_apply_delta_scaling = scale; }

        /**
         *
         * @return
         */
        int num_do_sza() const { return m_ndosza; }

        /** Sets the number of SZA's to use in the DO solution.  For plane parallel this should always be 1.
         *  For spherical mode typically 2 is sufficient for nadir applications, more may be necessary for
         *  limb viewing geometry.
         *
         * @param nsza
         */
        void set_num_do_sza(int nsza) { m_ndosza = nsza; }


        /**
         *
         * @return True if the weighting function calculation is enabled
         */
        bool wf_enabled() const { return m_enable_wfs; }

        /** Enables/disables the calculation of weighting functions
         *
         * @param enable
         */
        void set_wf_enabled(bool enable) { m_enable_wfs = enable; }

        /**
         *
         * @return The number of spherical iterations used in the Discrete Ordinates source (Not currently used)
         */
        int num_do_spherical_iterations() const { return m_ndosphericaliterations; }

        /** Sets the number of spherical iterations used in the discrete ordinates source (not currently used)
         *
         * @param n
         */
        void set_num_do_spherical_iterations(int n) { m_ndosphericaliterations = n; }

        /**
         *
         * @return The number of spherical iterations inside the HR source
         */
        int num_hr_spherical_iterations() const { return m_hr_nspherical_iterations; }

        /** Sets the number of spherical iterations used in the HR source
         *
         * @param n
         */
        void set_num_hr_spherical_iterations(int n) { m_hr_nspherical_iterations = n; }

        /**
         *
         * @return The number of incoming points at each diffuse point in the HR source
         */
        int num_hr_incoming() const { return m_hr_nincoming; }

        /** Sets the number of incoming points at each diffuse point in the HR source
         *
         * @param n
         */
        void set_num_hr_incoming(int n) { m_hr_nincoming = n; }

        /**
         *
         * @return The number of outgoing directions at each diffuse point in the HR source
         */
        int num_hr_outgoing() const { return m_hr_noutgoing; }

        /** Sets the number of outgoing directions at each diffuse point in the HR source
         *
         * @param n
         */
        void set_num_hr_outgoing(int n) { m_hr_noutgoing = n; }

        /**
         *
         * @return  True if the HR source is to be initialized by the DO source
         */
        bool initialize_hr_with_do() const { return m_initialize_hr_with_do_solution; }

        /** Set to True if the HR source is to be intialized by the DO source
         *
         * @param init
         */
        void set_initialize_hr_with_do(bool init) { m_initialize_hr_with_do_solution = init; }

        /**
         *
         * @return The precision of the Weighting Function calculation
         */
        WeightingFunctionPrecision wf_precision() const { return m_wf_precision; }

        /** Sets the precision of the weighting function calculation
         *
         * @param precision
         */
        void set_wf_precision(WeightingFunctionPrecision precision) { m_wf_precision = precision; }

        /**
         *
         * @return Then number of points (per diffuse profile) that we calculate the incoming signal on.  Can be set to -1 to use all the points
         */
        int num_hr_full_incoming_points() const { return m_hr_num_incoming_points; }

        /** Sets the number of points (per diffuse profile) that we calculate the incoming signal at. Can be set to -1
         *  to use all the points.
         *
         * @param points
         */
        void set_num_hr_full_incoming_points (int points) { m_hr_num_incoming_points = points; }



    private:
        // TODO: Refactor these into individual source configs?

        int m_nthreads;
        int m_nstokes;
        int m_ndostreams;
        int m_ndosza;
        int m_ndosphericaliterations;

        int m_nsinglescatter_moments;

        int m_hr_nincoming;
        int m_hr_noutgoing;

        int m_hr_nspherical_iterations;
        int m_hr_num_incoming_points;

        bool m_apply_delta_scaling;

        bool m_enable_wfs;


        SingleScatterSource m_single_scatter_source;
        MultipleScatterSource m_multiple_scatter_source;
        OccultationSource m_occultation_source;

        WeightingFunctionPrecision m_wf_precision;

        bool m_initialize_hr_with_do_solution;
	};
}