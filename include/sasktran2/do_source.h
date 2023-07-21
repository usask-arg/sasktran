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
            // TODO: Kill this array storage
            std::array<sasktran2::Dual<double, sasktran2::dualstorage::dense>, NSTOKES> source_terms; // [nstokes(array), angle, layer, sza, azi]

            sasktran2::Dual<double, sasktran2::dualstorage::dense> source_terms_linear; // [nstokes, angle, layer, sza, azi]

            std::array<sasktran2::Dual<double, sasktran2::dualstorage::dense>, NSTOKES> singlescatter_source_terms; // [nstokes(array), angle, layer, sza, azi]
            std::vector<sasktran2::LegendrePhaseStorage<NSTOKES, CNSTR>> phase_storage;
            std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>> phase_container;
        };

        std::vector<DOSourceDiffuseThreadStorage> m_storage;

        std::unique_ptr<sasktran2::grids::AltitudeGrid> m_altitude_grid;
        std::unique_ptr<sasktran2::grids::Grid> m_cos_angle_grid;
        const sasktran2::grids::Grid& m_sza_grid;

        std::vector<Eigen::MatrixXd> m_scattering_matrix_stream_angles;
        std::vector<Eigen::MatrixXd> m_scattering_matrix_interpolation_angles;

        int linear_storage_index(int angleidx, int layeridx, int szaidx, int aziidx) const;
        void generate_scattering_matrices(const sasktran_disco::PersistentConfiguration<NSTOKES, CNSTR>& do_config);

        const sasktran2::Geometry1D& m_geometry;
        const Config& m_config;

        int m_num_azi;
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

        const std::array<sasktran2::Dual<double, sasktran2::dualstorage::dense>, NSTOKES>& source(int threadidx) const { return m_storage[threadidx].source_terms; }
        const std::array<sasktran2::Dual<double, sasktran2::dualstorage::dense>, NSTOKES>& singlescatter_source(int threadidx) const { return m_storage[threadidx].singlescatter_source_terms; }

        std::array<sasktran2::Dual<double, sasktran2::dualstorage::dense>, NSTOKES>& source(int threadidx) { return m_storage[threadidx].source_terms; }

        const sasktran2::Dual<double, sasktran2::dualstorage::dense>& linear_source(int threadidx) const { return m_storage[threadidx].source_terms_linear;}

        int num_azi_terms()const  {
            return (int)(m_storage[0].source_terms[0].value_size() / m_num_azi);
        }


        const sasktran2::grids::AltitudeGrid& altitude_grid() const { return *m_altitude_grid; }
        const sasktran2::grids::Grid cos_angle_grid() const { return *m_cos_angle_grid; }

        std::unique_ptr<sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>> geometry_interpolator(
                const std::vector<sasktran2::raytracing::TracedRay>& rays,
                bool include_azimuth_weights=true
                );

        void create_location_source_interpolator(const std::vector<Eigen::Vector3d>& locations,
                                                 const std::vector<Eigen::Vector3d>& directions,
                                                 Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator
                                                 ) const;

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
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere; /**< Reference to the global atmosphere object */
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


    protected:
        const sasktran2::Geometry1D& m_geometry; /**< Reference to the global engine geometry */
        const sasktran2::Config* m_config; /**< Reference to the global config */
        const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays; /**< Reference to the LOS rays */
        const sasktran2::raytracing::RayTracerBase& m_raytracer; /**< Reference to the ray tracer */
        std::unique_ptr<sasktran2::grids::Grid> m_sza_grid; /**< SZA Grid values where the DO calculation is performed */

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

    /** An implementation of the DO Source term where the Post Processing is done exactly in each layer.
     *  I.e., the DO source from the RTE is taken, and then reintegrated at the exact viewing angles for each layer.
     *
     *  This is a slow process when there are a large number of lines of sight.  Usually this only makes sense to do
     *  when there is a single line of sight.  For example, in a nadir viewing geometry.
     *
     * @tparam NSTOKES
     * @tparam CNSTR
     */
    template<int NSTOKES, int CNSTR=-1>
    class DOSourceExactPostProcessing : public DOSource<NSTOKES, CNSTR> {
    private:

        sasktran_disco::VectorDim2<LegendrePhaseStorage<NSTOKES, CNSTR>> m_lp_segments; /**< Storage for Legendre functions at each layer [ray, layer_boundary] */

        /** Accumulates the source terms on the boundaries of the layers
         *
         * @param layer
         * @param layer_array
         * @param ray_idx
         * @param layer_idx
         * @param thread_storage
         * @param m
         */
        void accumulate_layer_boundary_sources(const sasktran2::raytracing::SphericalLayer& layer,
                                               const sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& layer_array,
                                               int ray_idx,
                                               int layer_idx,
                                               DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage,
                                               sasktran_disco::AEOrder m);

        /** Accumulates the source terms inside the layers
         *
         * @param layer_array
         * @param phase_storage
         * @param lp_segments
         * @param start_location
         * @param end_location
         * @param start_saa
         * @param end_saa
         * @param source_idx
         * @param sza_index
         * @param sza_weight
         * @param source
         * @param m
         * @param thread_storage
         */
        void accumulate_layer_source(const sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& layer_array,
                                     std::vector<sasktran_disco::LegendrePhaseContainer<NSTOKES>>& phase_storage,
                                     const std::vector<LegendrePhaseStorage<NSTOKES, CNSTR>>& lp_segments,
                                     const sasktran2::Location& start_location,
                                     const sasktran2::Location& end_location,
                                     double start_saa,
                                     double end_saa,
                                     int source_idx,
                                     int sza_index,
                                     int sza_weight,
                                     sasktran2::Dual<double, sasktran2::dualstorage::dense>& source,
                                     sasktran_disco::AEOrder m,
                                     DOSourceThreadStorage<NSTOKES, CNSTR> &thread_storage
        );

        /** Triggers the calculation of source terms when an azimuth expansion order has been solved
         *
         * @param optical_layer
         * @param thread_storage
         * @param szaidx
         * @param m
         * @param threadidx
         */
        void accumulate_solved_azimuth(sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
                                       DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
                                       int szaidx,
                                       sasktran_disco::AEOrder m,
                                       int threadidx
        );

    public:
        DOSourceExactPostProcessing(const sasktran2::Geometry1D& geometry,
                                    const sasktran2::raytracing::RayTracerBase& raytracer);

        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays);

        void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                               const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

    };

    template<int NSTOKES, int CNSTR=-1>
    class DOSourceInterpolatedPostProcessingView : public SourceTermInterface<NSTOKES> {
    private:
        const sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>& m_source_interpolator;
        const DOSourceDiffuseStorage<NSTOKES, CNSTR>& m_diffuse_storage;
    public:
        DOSourceInterpolatedPostProcessingView(const sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>& interpolator,
                                               const DOSourceDiffuseStorage<NSTOKES, CNSTR>& source_storage);

        void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                               const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

        void end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const {};
    };

    template<int NSTOKES, int CNSTR=-1>
    class DOSourceInterpolatedPostProcessing : public DOSource<NSTOKES, CNSTR> {
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
                                           const sasktran2::raytracing::RayTracerBase& raytracer);

        virtual void calculate(int wavelidx, int threadidx);
        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays);
        virtual void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                               const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

        const DOSourceDiffuseStorage<NSTOKES, CNSTR>& storage() const { return *m_diffuse_storage; }

    };

    template<int NSTOKES, int CNSTR=-1>
    class DOSourceSphericallyCorrected : public DOSourceInterpolatedPostProcessing<NSTOKES, CNSTR> {
    private:
        int m_nsphericaliterations;

        std::unique_ptr<sasktran_disco::VectorDim2<std::array<Eigen::SparseVector<double>, NSTOKES>>> m_spherical_source_interpolator;
        std::unique_ptr<DOSourceInterpolatedPostProcessingView<NSTOKES, CNSTR>> m_ms_source;
        std::unique_ptr<SourceTermInterface<NSTOKES>> m_singlescatter_source;
        sasktran2::SourceIntegrator<NSTOKES> m_integrator;
        std::vector<SourceTermInterface<NSTOKES>*> m_sources;

        std::vector<std::vector<Eigen::Matrix<double, -1, NSTOKES>>> m_do_solutions;

        std::vector<Eigen::VectorXd> m_altitude_integration_angles;
        std::vector<Eigen::VectorXd> m_altitude_integration_weights;

        const sasktran2::raytracing::RayTracerBase& m_raytracer;

        std::vector<sasktran2::raytracing::TracedRay> m_spherical_correction_rays;

        virtual void accumulate_solved_azimuth(sasktran_disco::OpticalLayerArray<NSTOKES, CNSTR>& optical_layer,
                                               DOSourceThreadStorage<NSTOKES, CNSTR>& thread_storage,
                                               int szaidx,
                                               sasktran_disco::AEOrder m,
                                               int threadidx
        );


    public:
        DOSourceSphericallyCorrected(const sasktran2::Geometry1D& geometry,
                                     const sasktran2::raytracing::RayTracerBase& raytracer);

        virtual void calculate(int wavelidx, int threadidx);
        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays);
        virtual void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);
    };


}