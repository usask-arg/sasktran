#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/source_interface.h>
#include <sasktran2/solartransmission.h>
#include <sasktran2/source_integrator.h>
#include <sktran_disco/sktran_do.h>
#include <sasktran2/math/wigner.h>
#include <sasktran2/math/trig.h>
#include <boost/math/quadrature/gauss.hpp>
#include <sktran_disco/sktran_do_specs.h>

namespace sasktran2 {
    /** Utility to store the Legendre functions needed for postprocessing the discrete ordinates solutions.
     *  Main responsibilities are computing the Legendre P R T functions for a given coszen, and then transferring
     *  this information to the container needed by the DO engine.
     *
     *  These values are only needed to post process the DO solution.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template <int NSTOKES, int CNSTR=-1>
    struct LegendrePhaseStorage {
        Eigen::Matrix<double, NSTOKES, -1> storage; /**< Internal storage of P R T functions */
        int nstr;

        /** Constructs the class
         *
         * @param nstr Number of streams
         */
        LegendrePhaseStorage(int nstr);

        /**
         *
         * @param m azimuth expansion order
         * @param l stream index
         * @return Linear index pointing to the storage element for azimuth m, stream l
         */
        int linear_index(int m, int l) const;

        /** Fills the internal storage for a given cosine of zenith angle
         *
         * @param coszen cosine zenith angle
         */
        void fill(double coszen);

        /** Fills a container needed for the DO engine for a given azimuth expansion order
         *
         * @param container Container to fill
         * @param m Azimuth expansion order
         */
        void set_phase_container(std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>>& container, int m) const;
    };


    /** Storage structure to hold all of the necessary objects to perform a single DO calculation for a single
     *  solar zenith angle.
     *
     *  A single calculation requires the internal PersistentConfiguration, a UserSpec object, and the geometry
     *  information defining the layers.
     *
     *  Note that a lot of this information is the same for each SZA calculation, but we do require a separate object
     *  for each.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template<int NSTOKES, int CNSTR>
    struct DOSingleSZACalculator {
        std::unique_ptr<sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>> persistent_config;
        sasktran_disco::SKTRAN_DO_UserSpec userspec;
        std::unique_ptr<sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>> geometry_layers;

        DOSingleSZACalculator() {}
    };


    /** Common storage that every thread requires when computing the DO source, independent of
     *  what is done with it afterwards.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template<int NSTOKES, int CNSTR>
    struct DOSourceThreadStorage {
        std::vector<DOSingleSZACalculator<NSTOKES, CNSTR>> sza_calculators;                   // [sza]
        std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>> legendre_phase_container;
        std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense>> boundary_sources; // [ray, layer_boundary]

        std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense>> layer_sources; // [ray, layer]

        std::vector<sasktran_disco::PostProcessingCache<NSTOKES, CNSTR>> postprocessing_cache; // [layer]

    };

    template<int NSTOKES, int CNSTR=-1>
    class DOSourceDiffuseStorage {
    private:
        struct DOSourceDiffuseThreadStorage {
            sasktran2::Dual<double, sasktran2::dualstorage::denseRowMajor> source_terms_linear; // [nstokes, angle, layer, sza, azi]
            std::vector<sasktran2::LegendrePhaseStorage<NSTOKES, CNSTR>> phase_storage;
            std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>> phase_container;
        };

        std::vector<DOSourceDiffuseThreadStorage> m_storage;

        std::unique_ptr<sasktran2::grids::AltitudeGrid> m_altitude_grid;
        std::unique_ptr<sasktran2::grids::Grid> m_cos_angle_grid;
        const sasktran2::grids::Grid& m_sza_grid;

        Eigen::VectorX<bool> m_need_to_calculate_map; // [source idx] This is a map of source indices that are actually used, some of them may not be required
        Eigen::VectorX<bool> m_converged_map;

        std::vector<Eigen::MatrixXd> m_scattering_matrix_stream_angles;
        std::vector<Eigen::MatrixXd> m_scattering_matrix_interpolation_angles;

        int linear_storage_index(int angleidx, int layeridx, int szaidx, int aziidx) const;
        int ground_start_index() const { return m_ground_start; };
        int ground_storage_index(int angleidx, int szaidx, int aziidx) const;

        const sasktran2::Geometry1D& m_geometry;
        const Config& m_config;

        int m_num_azi;
        int m_ground_start;

        void accumulate_ground_sources(sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
                                       sasktran_disco::AEOrder m,
                                       sasktran2::DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
                                       int szaidx,
                                       int thread_idx);

    public:
        DOSourceDiffuseStorage(const sasktran_disco::GeometryLayerArray<NSTOKES, CNSTR>& layer_geometry,
                               const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>& do_config,
                               const sasktran2::grids::Grid& m_sza_grid,
                               const Config& config,
                               const sasktran2::Geometry1D& geometry
                               );


        void accumulate_sources(sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
                                sasktran_disco::AEOrder m,
                                sasktran2::DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
                                int szaidx,
                                int thread_idx);

        const sasktran2::Dual<double, sasktran2::dualstorage::denseRowMajor>& linear_source(int threadidx) const { return m_storage[threadidx].source_terms_linear;}

        void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmo);


        const sasktran2::grids::AltitudeGrid& altitude_grid() const { return *m_altitude_grid; }
        const sasktran2::grids::Grid cos_angle_grid() const { return *m_cos_angle_grid; }

        std::unique_ptr<sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>> geometry_interpolator(
                const std::vector<sasktran2::raytracing::TracedRay>& rays,
                bool include_azimuth_weights=true
                );

        void create_location_source_interpolator(const std::vector<Eigen::Vector3d>& locations,
                                                 const std::vector<Eigen::Vector3d>& directions,
                                                 const std::vector<bool>& ground_hit_flag,
                                                 Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator
                                                 );

    };

    /** Virtual class that defines common elements to handle the DO source calculation.
     *
     *  This class basically handles the construction of the DO engine, and triggers the radiative transfer solution
     *  The solutions are then passed through to a derived class which handles the storage, postprocessing,
     *  and eventual source term calculations.
     *
     *  The class also determines the set of SZA's to perform calculations on from the input lines of sight
     *  and configuration, and manages the handling of that grid
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template<int NSTOKES, int CNSTR=-1>
    class DOSource : public SourceTermInterface<NSTOKES> {
    private:
        using SInterpolator = std::vector<std::vector<std::vector<std::pair<int, double>>>>;

        std::vector<sasktran_disco::LineOfSight> m_do_los; /**< Lines of sight converted to the LOS objects the DO engine needs */

        int m_nstr; /** Number of streams */

        /** Virtual method that derived classes must implement.  This method is called after the solution
         *  is computed for each SZA for each azimuth order.
         *
         *  Typically derived classes will use this call opportunity to postprocess the solution and to store
         *  the result.
         *
         * @param optical_layer Object containing the solved RTE and layer information
         * @param thread_storage Storage that the derived object may use temporarily
         * @param szaidx Index of the solar zenith angle solution that was calculated
         * @param m Index of the azimuth order that was solved
         * @param threadidx Current thread index
         */
        virtual void accumulate_solved_azimuth(sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
                                               DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
                                               int szaidx,
                                               sasktran_disco::AEOrder m,
                                               int threadidx
                                               ) = 0;

        /** Internal method to generate the SZA grid used by the mode.  Can only be called after m_los_rays
         *  is set. Afterwards m_sza_grid will be initialized.
         *
         */
        void generate_sza_grid();

        /** Internal method which constructs the interpolation weights from the atmosphere grid to the
         *  center of each LOS ray layer
         *
         */
        void construct_los_location_interpolator(const std::vector<sasktran2::raytracing::TracedRay>& los_rays);


    protected:
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere; /**< Reference to the global atmosphere object */
        const sasktran2::Geometry1D& m_geometry; /**< Reference to the global engine geometry */
        const sasktran2::Config* m_config; /**< Reference to the global config */
        const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays; /**< Reference to the LOS rays */
        const sasktran2::raytracing::RayTracerBase& m_raytracer; /**< Reference to the ray tracer */
        std::unique_ptr<sasktran2::grids::Grid> m_sza_grid; /**< SZA Grid values where the DO calculation is performed */
        SInterpolator m_los_source_weights; /** Interpolator from the LOS source query points to the atmosphere table */

        std::vector<DOSourceThreadStorage<NSTOKES, CNSTR>> m_thread_storage; /**< Internal thread storage */
    public:
        /** Initializes the base DOSource object which requires the global geometry and raytracer
         *
         * @param geometry
         * @param raytracer
         */
        DOSource(const sasktran2::Geometry1D& geometry,
                 const sasktran2::raytracing::RayTracerBase& raytracer);

        /** Initializes the configuration
         *
         * @param config
         */
        void initialize_config(const sasktran2::Config& config);

        /** Initializes the geometry
         *
         * @param los_rays
         */
        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays);

        /** Initializes the atmosphere
         *
         * @param atmosphere
         */
        virtual void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /** Triggers calculation of the source term
         *
         * @param wavelidx
         * @param threadidx
         */
        virtual void calculate(int wavelidx, int threadidx);

        void end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;
    };

    template<int NSTOKES, int CNSTR=-1>
    class DOSourceInterpolatedPostProcessing : public DOSource<NSTOKES, CNSTR> {
    private:
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
        bool m_will_integrate_sources;

    protected:
        std::unique_ptr<sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>> m_los_source_interpolator;

        sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>* m_source_interpolator_view;

        std::unique_ptr<DOSourceDiffuseStorage<NSTOKES, CNSTR>> m_diffuse_storage;

        virtual void accumulate_solved_azimuth(sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
                                               DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
                                               int szaidx,
                                               sasktran_disco::AEOrder m,
                                               int threadidx
        );

    public:
        DOSourceInterpolatedPostProcessing(const sasktran2::Geometry1D& geometry,
                                           const sasktran2::raytracing::RayTracerBase& raytracer,
                                           bool will_integrate_sources=true
                                           );

        virtual void calculate(int wavelidx, int threadidx);
        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays) override;
        virtual void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) override;
        virtual void initialize_config(const sasktran2::Config& config) override;


        void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                               const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

        DOSourceDiffuseStorage<NSTOKES, CNSTR>& storage() const { return *m_diffuse_storage; }

    };
}