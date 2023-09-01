#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>
#include "../geometry.h"


namespace sasktran2::atmosphere {
    /** Base storage container that defines the scattering phase information at a single point.  Scattering phase
     *  information can either be specified/stored as Legendre coefficients or as the raw phase matrix.
     *  Optionally both can be specified, or conversions can be done to convert between the two representations.
     *  The phase storage information is dependent on the number of stokes parameters being calculated.
     * 
     */
    template <int NSTOKES>
    class PhaseStorage {
    private:
        Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1> m_storage; // (leg coeff (stacked for polarized), geometry location)

        std::vector<Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1>> m_derivatives; // (leg coeff (stacked for polarized), geometry location)
        std::vector<Eigen::Vector<sasktran2::types::leg_coeff, -1>> m_f_derivatives; // (geometry location) derivative of the delta scale factor

        int m_scatderivstart;
    public:
        int max_stored_legendre() const {
            if constexpr(NSTOKES == 1) {
                return (int)m_storage.rows();
            } else if constexpr(NSTOKES == 3) {
                return (int)(m_storage.rows() / 4);
            }
        }
        int num_deriv() const { return (int)m_derivatives.size(); }
        int scattering_deriv_start() const { return m_scatderivstart; }

        const Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1>& storage() const { return m_storage; }
        Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1>& storage() { return m_storage; }

        const Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1>& deriv_storage(int derivindex) const { return m_derivatives[derivindex]; }
        const Eigen::Vector<sasktran2::types::leg_coeff, -1>& f_deriv_storage(int derivindex) const { return m_f_derivatives[derivindex]; }
        Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1>& deriv_storage(int derivindex) { return m_derivatives[derivindex]; }
        Eigen::Vector<sasktran2::types::leg_coeff, -1>& f_deriv_storage(int derivindex) { return m_f_derivatives[derivindex]; }

        void resize(int numgeo, int numlegendre) {
            // TODO: based on NSTOKES
            m_storage.resize(numlegendre, numgeo);

            m_storage.setZero();
        }

        void resize_derivative(int numgeo, int legendre, int numderiv, int derivstart) {
            // TODO: based on NSTOKES
            m_derivatives.resize(numderiv);
            m_f_derivatives.resize(numderiv);

            for(auto& deriv : m_derivatives) {
                deriv.resize(legendre, numgeo);
            }

            for(auto& deriv : m_f_derivatives) {
                deriv.resize(numgeo);
                deriv.setZero();
            }

            m_scatderivstart = derivstart;
        }
    };

    /** Class which performs interpolation over the phase matrix.
     *  This can either mean direct interpolation of phase matrix values, or calculation of the phase
     *  matrix from the Legendre/greek coefficients
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES, bool ssonly=false>
    class PhaseInterpolator {
    using ScatWeightType = typename std::conditional<ssonly, Eigen::Matrix<sasktran2::types::leg_coeff, NSTOKES, -1>, Eigen::Matrix<sasktran2::types::leg_coeff, NSTOKES*NSTOKES, -1>>::type;

    private:
        ScatWeightType m_scattering_weights;
        bool m_geometry_loaded;

    public:
        PhaseInterpolator();

        void load_scattering_angle(
                const PhaseStorage<NSTOKES>& storage,
                const Eigen::Vector3d& incoming_ray, const Eigen::Vector3d& outgoing_ray, bool outgoing_facing_away=true);

        template<sasktran2::dualstorage S>
        void scatter(const PhaseStorage<NSTOKES>& phase_storage,
                     const std::vector<std::pair<int, double>>& index_weights,
                     sasktran2::Dual<double, S, NSTOKES>& source
                     ) const;
    };


    /** Base abstract storage container for specifying the atmospheric constituent parameters on the full
     *  radiative transfer geometry grid, which may be 1dimensional, 2dimensional, or 3dimensional.
     *  The full required information is total extinction, scattering extinction, and phase information for
     *  each geometry grid point.
     */
    class AtmosphereGridStorage {

    };

    class AtmosphereGridStorageAbsorber : public AtmosphereGridStorage {
    private:
        Eigen::MatrixXd m_total_extinction;
    public:
        AtmosphereGridStorageAbsorber(const Eigen::MatrixXd& extinction) { m_total_extinction = extinction; }
    };

    template <int NSTOKES>
    class AtmosphereGridStorageConstantHeightScatterer : public AtmosphereGridStorage {
    private:
        Eigen::MatrixXd m_scat_extinction;
        Eigen::MatrixXd m_total_extinction;

        PhaseStorage<NSTOKES> m_phase;
    public:


    };

    template <int NSTOKES>
    class AtmosphereGridStorageFull : public AtmosphereGridStorage {
    public:
        Eigen::MatrixXd ssa;                // location, wavel
        Eigen::MatrixXd total_extinction;   // location, wavel
        Eigen::Matrix<sasktran2::types::leg_coeff, -1, -1> f; // location, wavel, (Delta scaling factor)

        int applied_f_order;                // Order of the delta_m scaling
        int applied_f_location;             // Index to the phase moment that defines the legendre scaling f

        std::vector<PhaseStorage<NSTOKES>> phase;
    public:
        AtmosphereGridStorageFull(int nwavel, int nlocation, int numlegendre) {
            ssa.resize(nlocation, nwavel);
            total_extinction.resize(nlocation, nwavel);
            f.resize(nlocation, nwavel);
            phase.resize(nwavel);

            for(auto& p : phase) {
                if constexpr(NSTOKES == 1) {
                    p.resize(nlocation, numlegendre);
                } else {
                    p.resize(nlocation, numlegendre*4);
                }
            }

            ssa.setZero();
            total_extinction.setZero();

            applied_f_location = -1;
            applied_f_order = -1;
        }
    };
}