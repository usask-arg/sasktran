#pragma once

#include <sasktran2/source_interface.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/atmosphere/atmosphere.h>
#include <sasktran2/config.h>
#include <sasktran2/dual.h>

namespace sasktran2::solartransmission {
    /** Generates one row of what is known as the geometry matrix.  The optical depth at every necessary integration
     *  point can be written as \f$M k\f$ where \f$k\f$ is the extinction coefficient vector.  In 1D geometry, the
     *  matrix \f$M\f$ is densish, but it may be sparse in higher dimensions.  This function constructs one row
     *  of the matrix \f$M\f$, which corresponds to one traced ray.
     *
     * @param row row index
     * @param traced_ray ray traced to the sun
     * @param geometry global geometry
     * @param result matrix to store the result in
     * @param index_weights buffer to avoid allocs
     */
    inline void assign_dense_matrix_column(int row, const sasktran2::raytracing::TracedRay& traced_ray, const Geometry& geometry, Eigen::MatrixXd& result,
                                           std::vector<std::pair<int, double>>& index_weights) {
        for(int i = 0; i < traced_ray.layers.size(); ++i) {
            const auto& layer = traced_ray.layers[i];

            geometry.assign_interpolation_weights(layer.entrance, index_weights);

            for(const auto& iw : index_weights) {
                result(row, iw.first) += iw.second * layer.od_quad_start;
            }

            geometry.assign_interpolation_weights(layer.exit, index_weights);

            for(const auto& iw : index_weights) {
                result(row, iw.first) += iw.second * layer.od_quad_end;
            }
        }
    }

    class SolarTransmissionBase {
    protected:
        const Geometry1D& m_geometry;
        const sasktran2::raytracing::RayTracerBase& m_raytracer;
    public:
        SolarTransmissionBase(const Geometry1D& geometry, const sasktran2::raytracing::RayTracerBase& raytracer) :
                m_geometry(geometry), m_raytracer(raytracer){}
    };

    class SolarTransmissionExact : public SolarTransmissionBase {
    private:
    public:
        SolarTransmissionExact(const Geometry1D& geometry, const sasktran2::raytracing::RayTracerBase& raytracer) :
                SolarTransmissionBase(geometry, raytracer){}

        virtual void initialize_config(const sasktran2::Config& config) {};

        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& integration_rays) {};

        void generate_geometry_matrix(const std::vector<sasktran2::raytracing::TracedRay>& rays,
                                      Eigen::MatrixXd& od_matrix,
                                      std::vector<bool>& ground_hit_flag
                                      ) const;
    };

    class SolarTransmissionTable : public SolarTransmissionExact {
    private:
        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator> m_location_interpolator;
        const sasktran2::Config* m_config;
        Eigen::MatrixXd m_geometry_matrix;

        std::vector<bool> m_ground_hit_flag;
    public:
        SolarTransmissionTable(const Geometry1D& geometry, const sasktran2::raytracing::RayTracerBase& raytracer) :
                SolarTransmissionExact(geometry, raytracer){}

        void initialize_config(const sasktran2::Config& config) override { m_config = &config; };

        void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& integration_rays) override;

        void generate_interpolation_matrix(const std::vector<sasktran2::raytracing::TracedRay>& rays,
                                      Eigen::SparseMatrix<double, Eigen::RowMajor>& interpolator,
                                      std::vector<bool>& ground_hit_flag
        ) const;

        const Eigen::MatrixXd& geometry_matrix() const { return m_geometry_matrix; }
    };


    template <typename S, int NSTOKES>
    class SingleScatterSource : public SourceTermInterface<NSTOKES> {
    private:
        S m_solar_transmission;
        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;

        Eigen::MatrixXd m_geometry_matrix;
        Eigen::SparseMatrix<double, Eigen::RowMajor> m_geometry_sparse;
        std::vector<bool> m_ground_hit_flag;

        std::vector<Eigen::VectorXd> m_solar_trans;
        std::vector<std::vector<int>> m_index_map;
        std::vector<std::vector<int>> m_phase_index_map;
        std::vector<sasktran2::atmosphere::PhaseInterpolator<NSTOKES, true>> m_phase_interp;

        sasktran2::Dual<double> m_precomputed_sources;

        mutable std::vector<std::vector<std::pair<int, double>>> m_thread_index_cache_one;
        mutable std::vector<std::vector<std::pair<int, double>>> m_thread_index_cache_two;

        mutable std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>> m_start_source_cache;
        mutable std::vector<sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>> m_end_source_cache;

        const Geometry1D& m_geometry;
        const sasktran2::Config* m_config;

        const std::vector<sasktran2::raytracing::TracedRay>* m_los_rays;

        int m_num_cells;

        void integrated_source_quadrature(int wavelidx, int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                                          const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

        void integrated_source_constant(int wavelidx, int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                                              const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

        void integrated_source_linear(int wavelidx, int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                                        const sasktran2::SparseODDualView& shell_od,
                                        sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

    public:
        SingleScatterSource(const Geometry1D& geometry, const sasktran2::raytracing::RayTracerBase& raytracer
                            ) :
        m_solar_transmission(geometry, raytracer), m_geometry(geometry) {};

        void initialize_config(const sasktran2::Config& config) override;

        /** Here the single scatter source term initializes the internal solar transmission object, usually this
         *  involves tracing the required rays and setting up any internal matrices.
         *
         *  @param los_rays The traced line of sight rays
         */
        void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays) override;

        /**
        *
        */
        void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) override;

        /** Triggers an internal calculation of the source term.  This method is called at the beginning of each
         *  'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        void calculate(int wavelidx, int threadidx) override;

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(int wavelidx, int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                               const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const override;

        /** Calculates the source term at the end of the ray.  Common examples of this are ground scattering, ground emission,
         *  or the solar radiance if looking directly at the sun.
         *
         * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
         * @param surface The surface object
         * @param source The returned source term
         */
        void end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const override;
    };

    template <int NSTOKES>
    class OccultationSource : public SourceTermInterface<NSTOKES> {
    private:

    public:
        void initialize_config(const sasktran2::Config& config) override;

        /** Initializes any geometry information that is required for calculating the source term.  This method is called
         *  after the line of sight rays ar traced.
         *
         * @param los_rays The traced line of sight rays
         */
        void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays) override;

        /**
         *
         */
        void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere) override;

        /** Triggers an internal calculation of the source term.  This method is called at the beginning of each
         *  'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        void calculate(int wavelidx, int threadidx) override;

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                               const sasktran2::SparseODDualView& shell_od,
                               sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const override;

        /** Calculates the source term at the end of the ray.  Common examples of this are ground scattering, ground emission,
         *  or the solar radiance if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
         * @param source The returned source term
         */
        void end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const override;
    };

}