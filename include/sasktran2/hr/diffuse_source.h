#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/source_interface.h>
#include <sasktran2/source_integrator.h>
#include <sasktran2/do_source.h>
#include <sasktran2/hr/diffuse_point.h>

namespace sasktran2::hr {
    /** Thread specific storage for the diffuse table
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES>
    struct DiffuseTableThreadStorage {
        sasktran2::Dual<double, sasktran2::dualstorage::dense> m_incoming_radiances; /**< Total incoming stokes vectors, [stokes, direction, diffusegrid] */
        sasktran2::Dual<double, sasktran2::dualstorage::dense> m_firstorder_radiances; /**< Total first order incoming radiances (usaully single scatter + optional thermal/photochemical), [stokes, direction, diffusegrid] */
        sasktran2::Dual<double, sasktran2::dualstorage::dense> m_outgoing_sources; /**< Outgoing source terms, [stokes, direction, diffusegrid] */

        std::vector<Eigen::MatrixXd> point_scattering_matrices; /** For each point, outgoing source = scattering matrix @ incoming */
        Eigen::SparseMatrix<double, Eigen::ColMajor> accumulation_matrix; /** incoming_radiance = accumulation_matrix @ outgoing_sources */
    };

    /** An implementation of the successive orders of scattering technique.  We call this HR mostly for historic
     *  reasons.
     *
     *
     * Note: Currently this source is not perfectly linearized.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES>
    class DiffuseTable : public SourceTermInterface<NSTOKES>  {
        using SInterpolator = std::vector<RaySourceInterpolationWeights>;

    private:
        std::vector<DiffuseTableThreadStorage<NSTOKES>> m_thread_storage; /** Thread (optical) data [nthreads] */

        const sasktran2::Config* m_config; /**< Internal reference to the config */

        sasktran2::SourceIntegrator<NSTOKES> m_integrator; /**< Integrator used to get OD along the rays, as well as integrate the initial sources to get the incoming radiances */
        std::vector<SourceTermInterface<NSTOKES>*> m_initial_sources; /**< Initial sources used to get the first order incoming radiances */

        std::vector<std::unique_ptr<SourceTermInterface<NSTOKES>>> m_initial_owned_sources; /**< Initial sources used to get the first order incoming radiances */

        const sasktran2::raytracing::RayTracerBase& m_raytracer; /**< Raytracer used for the diffuse point incoming rays */
        const sasktran2::Geometry1D& m_geometry; /** Global geometry object */

        std::unique_ptr<sasktran2::grids::SourceLocationInterpolator> m_location_interpolator; /** Interpolates location */

        std::vector<sasktran2::raytracing::TracedRay> m_incoming_traced_rays; /** Traced incoming rays to all diffuse points */

        std::vector<std::unique_ptr<DiffusePoint<NSTOKES>>> m_diffuse_points; /** Stacked vector of all interior diffuse points, including ground, interpolated using m_location_interpolator */
        std::vector<bool> m_diffuse_point_full_calculation; /** True if we are doing the full incoming calculation at this diffuse point, false if it is interpolated from spherical corrections */

        std::vector<std::unique_ptr<sasktran2::hr::IncomingOutgoingSpherePair<NSTOKES>>> m_unit_sphere_pairs; /** Unit sphere ownership for the diffuse points, never accessed after construction by this class */

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere; /** Reference to the atmosphere object */

        std::vector<int> m_diffuse_incoming_index_map; /** element i is the start index of point i in m_incoming_traced_rays */
        std::vector<int> m_diffuse_outgoing_index_map; /** element i is the start index of point i in outgoing directions */

        std::vector<std::vector<std::pair<int, double>>> m_diffuse_point_interpolation_weights; /** Interpolation mapping from the global coordinates to the diffuse locations.  Mostly used to generate the scattering matrices. */

        SInterpolator m_los_source_weights; /** Interpolator mapping from line of sight points to source terms in this table */
        SInterpolator m_diffuse_source_weights; /** Interpolator mapping from incoming rays to source terms in this table */

        int m_total_num_diffuse_weights; /** Total number of diffuse weights, used to help memory allocs */

        Eigen::SparseMatrix<double, Eigen::RowMajor> m_do_to_diffuse_outgoing_interpolator; /** Mapping from the DO source terms to the outgoing sphere sources */
        DOSourceInterpolatedPostProcessing<NSTOKES, -1>* m_do_source; /** Reference to the DO source */

    private:
        sasktran2::grids::Grid generate_cos_sza_grid(double min_cos_sza, double max_cos_sza);
        sasktran2::grids::AltitudeGrid generate_altitude_grid();

        void construct_diffuse_points();
        void trace_incoming_rays();
        void generate_scattering_matrices(int wavelidx, int threadidx);
        void generate_accumulation_matrix(int wavelidx, int threadidx);
        void iterate_to_solution(int wavelidx, int threadidx);
        void interpolate_sources(const Eigen::VectorXd& old_outgoing, sasktran2::Dual<double>& new_outgoing);
        void generate_source_interpolation_weights(const std::vector<sasktran2::raytracing::TracedRay>& rays,
                                                   SInterpolator& interpolator,
                                                   int& total_num_weights
                                                   ) const;

        Eigen::Vector3d rotate_unit_vector(const Eigen::Vector3d& vector,
                                           const Eigen::Vector3d& initial_position,
                                           const Eigen::Vector3d& new_position
                                           ) const;
    public:
        DiffuseTable(const sasktran2::raytracing::RayTracerBase& ray_tracer,
                     const sasktran2::Geometry1D& geometry
                     );

        /** Initializes the config inside the source term
         *
         * @param config
         */
        virtual void initialize_config(const sasktran2::Config& config);

        /** Initializes any geometry information that is required for calculating the source term.  This method is called
         *  after the line of sight rays ar traced.
         *
         * @param los_rays The traced line of sight rays
         */
        virtual void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& los_rays);

        /**
         *
         */
        virtual void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere);

        /** Triggers an internal calculation of the source term.  This method is called at the beginning of each
         *  'wavelength' calculation.
         *
         * @param wavelidx Index of the wavelength being calculated
         */
        virtual void calculate(int wavelidx, int threadidx);

        /** Calculates the integrated source term for a given layer.
         *
         * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
         * @param layeridx Raw index pointing to the layer that was previosuly passed in initialize_geometry
         * @param layer The layer that we are integrating over
         * @param source The returned source term
         */
        virtual void integrated_source(int wavelidx,  int losidx, int layeridx, int threadidx, const sasktran2::raytracing::SphericalLayer& layer,
                                       const sasktran2::SparseODDualView& shell_od,
                                       sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;

        /** Calculates the source term at the end of the ray.  Common examples of this are ground scattering, ground emission,
         *  or the solar radiance if looking directly at the sun.
         *
         * @param wavelidx Raw index for the wavelength we are calculating
         * @param losidx Raw index pointing to the ray that was previously passed in initialize_geometry
         * @param source The returned source term
         */
        virtual void end_of_ray_source(int wavelidx, int losidx, int threadidx, sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& source) const;
    };
}



