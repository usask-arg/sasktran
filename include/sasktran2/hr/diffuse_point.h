#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/math/unitsphere.h>
#include <sasktran2/atmosphere/grid_storage.h>

namespace sasktran2::hr {
    template <int NSTOKES>
    class IncomingOutgoingSpherePair {
    private:
        std::unique_ptr<const sasktran2::math::UnitSphere> m_incoming_sphere;
        std::unique_ptr<const sasktran2::math::UnitSphere> m_outgoing_sphere;

        // TODO: NSTOKES for something other than 3?
        using ScatMatrices = typename std::conditional<NSTOKES == 1, std::vector<std::array<Eigen::MatrixXd, 1>>, std::vector<std::array<Eigen::MatrixXd, 4>>>::type;

        ScatMatrices m_legendre_scat_mats;

        bool m_is_configured;

        void assign_scat_mat_block(int legendre_idx, int in_idx, int out_idx);
        void configure_geometry();
    public:
        IncomingOutgoingSpherePair(int nlegendre, std::unique_ptr<const sasktran2::math::UnitSphere>&& incoming_sphere, std::unique_ptr<const sasktran2::math::UnitSphere>&& outgoing_sphere);


        void calculate_scattering_matrix(const sasktran2::atmosphere::PhaseStorage<NSTOKES>& phase,
                                         const std::vector<std::pair<int, double>>& index_weights,
                                         double* phase_storage_location) const;

        const sasktran2::math::UnitSphere& incoming_sphere() const { return *m_incoming_sphere; }
        const sasktran2::math::UnitSphere& outgoing_sphere() const { return *m_outgoing_sphere; }
    };


    template <int NSTOKES>
    class DiffusePoint {
    private:
        const IncomingOutgoingSpherePair<NSTOKES>* m_spheres;
        const sasktran2::Location m_location;
        sasktran2::Dual<double, sasktran2::dualstorage::dense> m_phase_matrix_storage;
    public:
        DiffusePoint(const IncomingOutgoingSpherePair<NSTOKES>& spheres,
                     const sasktran2::Location&
                     );

        const sasktran2::Location& location() const { return m_location; }

        int num_incoming() const { return m_spheres->incoming_sphere().num_points(); }
        int num_outgoing() const { return m_spheres->outgoing_sphere().num_points(); }

        const IncomingOutgoingSpherePair<NSTOKES>& sphere_pair() const { return *m_spheres; }

        Eigen::Vector3d incoming_direction(int in_idx) const { return m_spheres->incoming_sphere().get_quad_position(in_idx); }

    };

    template <int NSTOKES>
    class DiffuseGroundPoint : public DiffusePoint<NSTOKES> {

    };
}