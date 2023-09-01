#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/raytracing.h>
#include <sasktran2/source_interface.h>

namespace sasktran2 {

    struct RaySourceInterpolationWeights {
        std::vector<std::pair<std::vector<std::pair<int, double>>, std::vector<std::pair<int, double>>>> interior_weights;
        std::vector<std::pair<int, double>> ground_weights;
        bool ground_is_hit;
    };

    /** Class that integrates source terms along the ray.  Note that in SASKTRAN2, source terms themselves are responsible
     *  for integrating across the layer, this class simply adds the source terms in each layer and attenuates them
     *  by the optical depth.
     *
     *  Integration takes place in three steps.  First initialize_geometry is called with the rays that will be
     *  eventually integrated to set up geometry factors.
     *
     *  Next, initialize_atmosphere is called so that any optical parameters can be pre-calculated, such as the OD
     *  for each layer.
     *
     *  Lastly, integrate is called on each ray individually, summing the overall sources together.
     *
     * @tparam NSTOKES
     */
    template<int NSTOKES>
    class SourceIntegrator {
        using SInterpolator = std::vector<RaySourceInterpolationWeights>;
    private:
        bool m_calculate_derivatives; /**< True if we are calculating derivatives */
        std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> m_traced_ray_od_matrix; /**< Vector of matrices A such that A * atmosphere_extinction = OD for each layer in that ray */

        std::vector<Eigen::MatrixXd> m_shell_od; /**< Vector of matrices that stores the optical depth for each layer */
        std::vector<Eigen::MatrixXd> m_exp_minus_shell_od;

        const std::vector<sasktran2::raytracing::TracedRay>* m_traced_rays; /**< Reference to the rays we are integrating */

        const sasktran2::atmosphere::Atmosphere<NSTOKES>* m_atmosphere;
    public:
        /**
         *
         * @param calculate_derivatives True if the source integrator should calculate derivatives
         */
        SourceIntegrator(bool calculate_derivatives);

        /**
         *
         * @param enable True if the source integrator should calculate derivatives
         */
        void set_calculate_derivatives(bool enable) { m_calculate_derivatives = enable; }

        /** Initializes the geometry of the source integrator
         *
         * @param traced_rays Vector of traced rays
         * @param geometry Global geometry
         */
        void initialize_geometry(const std::vector<sasktran2::raytracing::TracedRay>& traced_rays, const Geometry& geometry);

        /** Initializes the atmosphere
         *
         * @param atmo
         */
        void initialize_atmosphere(const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmo);

        /** Integrates the source terms and stores the result in radiance
         *
         * @param radiance
         * @param source_terms
         * @param wavelidx
         * @param threadidx
         * @param rayidx
         */
        void integrate(
                sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& radiance,
                std::vector<SourceTermInterface<NSTOKES>*> source_terms,
                int wavelidx,
                int rayidx,
                int threadidx
                );

        void integrate_and_emplace_accumulation_triplets(
                sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& radiance,
                std::vector<SourceTermInterface<NSTOKES>*> source_terms,
                int wavelidx,
                int rayidx,
                int threadidx,
                const SInterpolator& source_interpolator,
                std::vector<Eigen::Triplet<double>>& triplets
                );

        /** Calculates the Optical Depth for each ray */
        void integrate_optical_depth(
                sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& radiance,
                int wavelidx,
                int rayidx,
                int threadidx
                );

    };
}