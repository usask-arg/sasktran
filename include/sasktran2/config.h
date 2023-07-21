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
         *  'none' Removes the multiple scatter source from the calculation.
         *
         */
        enum class MultipleScatterSource {
            discrete_ordinates,
            hr,
            none
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

        /** Sets the number of DO streams to use in the calculation.  Only used if the multiple scatter source is
         * set to discrete_ordinates
         *
         * @param nstr
         */
        void set_num_do_streams(int nstr)  { m_ndostreams = nstr;}


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


        int num_do_spherical_iterations() const { return m_ndosphericaliterations; }
        void set_num_do_spherical_iterations(int n) { m_ndosphericaliterations = n; }

        int num_hr_spherical_iterations() const { return m_hr_nspherical_iterations; }
        void set_num_hr_spherical_iterations(int n) { m_hr_nspherical_iterations = n; }

        int num_hr_incoming() const { return m_hr_nincoming; }
        void set_num_hr_incoming(int n) { m_hr_nincoming = n; }

        int num_hr_outgoing() const { return m_hr_noutgoing; }
        void set_num_hr_outgoing(int n) { m_hr_noutgoing = n; }

        bool initialize_hr_with_do() const { return m_initialize_hr_with_do_solution; }
        void set_initialize_hr_with_do(bool init) { m_initialize_hr_with_do_solution = init; }

    private:
        int m_nthreads;
        int m_nstokes;
        int m_ndostreams;
        int m_ndosza;
        int m_ndosphericaliterations;

        int m_hr_nincoming;
        int m_hr_noutgoing;

        int m_hr_nspherical_iterations;


        SingleScatterSource m_single_scatter_source;
        MultipleScatterSource m_multiple_scatter_source;
        OccultationSource m_occultation_source;

        bool m_initialize_hr_with_do_solution;
	};
}