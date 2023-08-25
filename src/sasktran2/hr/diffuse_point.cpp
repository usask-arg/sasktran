#include <sasktran2/hr/diffuse_point.h>
#include <sasktran2/math/wigner.h>
#include <sasktran2/math/scattering.h>

namespace sasktran2::hr {
    template<int NSTOKES>
    IncomingOutgoingSpherePair<NSTOKES>::IncomingOutgoingSpherePair(int nlegendre,
                                                                    std::unique_ptr<const sasktran2::math::UnitSphere>&& incoming_sphere,
                                                                    std::unique_ptr<const sasktran2::math::UnitSphere>&& outgoing_sphere) :
                                                                    m_incoming_sphere(std::move(incoming_sphere)),
                                                                    m_outgoing_sphere(std::move(outgoing_sphere)),
                                                                    m_legendre_scat_mats(nlegendre),
                                                                    m_is_configured(false)
                                                                    {
        configure_geometry();
    }

    template<>
    void IncomingOutgoingSpherePair<1>::assign_scat_mat_block(int legendre_idx, int in_idx, int out_idx) {
        auto& scat_mat = m_legendre_scat_mats[legendre_idx][0];

        sasktran2::math::WignerDCalculator wigner(0, 0);

        double cos_theta = (m_incoming_sphere->get_quad_position(in_idx)).dot(m_outgoing_sphere->get_quad_position(out_idx));

        if(cos_theta > 1) {
            cos_theta = 1;
        }
        if(cos_theta < -1) {
            cos_theta = -1;
        }

        scat_mat(out_idx, in_idx) = m_incoming_sphere->quadrature_weight(in_idx) * wigner.d(acos(cos_theta), legendre_idx);
    }

    template<>
    void IncomingOutgoingSpherePair<3>::assign_scat_mat_block(int legendre_idx, int in_idx, int out_idx) {
        auto incoming_ray = m_incoming_sphere->get_quad_position(in_idx);
        auto outgoing_ray = m_outgoing_sphere->get_quad_position(out_idx);

        double C1, C2, S1, S2, theta;
        int negation;

        // Incoming ray is opposite the propagation direction, outgoing ray is in the direction of propagation
        math::stokes_scattering_factors(-1*incoming_ray, outgoing_ray, theta, C1, C2, S1, S2, negation);


        // We have 4 greek coefficients (a1, a2, a3, b1), and we use the expansion (sum over legendre poly assumed)
        // P11 = a1 d_{0, 0}
        // P12 = b1 d_{0, 2}
        // P22 + P33 = (a2 + a3) d_{2, 2}
        // P22 - P33 = (a2 - a3) d_{2,-2}

        // Quick simplification
        // 2 P22 = (a2 + a3) d_{2, 2} + (a2 - a3) d_{2,-2}
        // 2 P33 = (a2 + a3) d_{2, 2} - (a2 - a3) d_{2,-2}

        // 2 P22 = a2 ( d_{2, 2} + d_{2, -2} ) + a3 ( d_{2, 2} - d_{2, -2} )
        // 2 P33 = a2 ( d_{2, 2} - d_{2, -2} ) + a3 ( d_{2, 2} + d_{2, -2} )

        // P22 = 0.5 * (a2 w_add + a3 w_minus)
        // P33 = 0.5 * (a2 w_minus + a3 w_add)


        // Create the wigner functions that we need
        sasktran2::math::WignerDCalculator d00(0, 0);
        sasktran2::math::WignerDCalculator d22(2, 2);
        sasktran2::math::WignerDCalculator d02(0, 2);
        sasktran2::math::WignerDCalculator d2m2(2, -2);

        double w_add = (d22.d(theta, legendre_idx) + d2m2.d(theta, legendre_idx));
        double w_minus = (d22.d(theta, legendre_idx) - d2m2.d(theta, legendre_idx));

        auto& scat_mat = m_legendre_scat_mats[legendre_idx];

        int start_row = 3*out_idx;
        int start_col = 3*in_idx;

        // Scat mat indicies
        // 0 - a1
        // 1 - a2
        // 2 - a3
        // 3 - b1

        // Note that b1 is assumed to be for the "generalized spherical functions" and so it picks up a minus
        // sign when we use the wigner function instead

        // 0, 0 component
        scat_mat[0].block<3, 3>(start_row, start_col)(0, 0) = d00.d(theta, legendre_idx);

        // 1, 0 component
        scat_mat[3].block<3, 3>(start_row, start_col)(1, 0) = -C2 * d02.d(theta, legendre_idx);

        // 2, 0 component
        scat_mat[3].block<3, 3>(start_row, start_col)(2, 0) = -S2 * d02.d(theta, legendre_idx);

        // 0, 1 component
        scat_mat[3].block<3, 3>(start_row, start_col)(0, 1) = -C1 * d02.d(theta, legendre_idx);

        // 1, 1 component
        scat_mat[1].block<3, 3>(start_row, start_col)(1, 1) += 0.5 * C1 * C2 * w_add;
        scat_mat[2].block<3, 3>(start_row, start_col)(1, 1) += 0.5 * C1 * C2 * w_minus;

        scat_mat[1].block<3, 3>(start_row, start_col)(1, 1) -= 0.5 * S1 * S2 * w_minus;
        scat_mat[2].block<3, 3>(start_row, start_col)(1, 1) -= 0.5 * S1 * S2 * w_add;

        // 2, 1 component
        scat_mat[1].block<3, 3>(start_row, start_col)(2, 1) += 0.5 * C1 * S2 * w_add;
        scat_mat[2].block<3, 3>(start_row, start_col)(2, 1) += 0.5 * C1 * S2 * w_minus;

        scat_mat[1].block<3, 3>(start_row, start_col)(2, 1) += 0.5 * S1 * C2 * w_minus;
        scat_mat[2].block<3, 3>(start_row, start_col)(2, 1) += 0.5 * S1 * C2 * w_add;

        // 0, 2 component
        scat_mat[3].block<3, 3>(start_row, start_col)(0, 2) = S1 * d02.d(theta, legendre_idx);

        // 1, 2 component
        scat_mat[1].block<3, 3>(start_row, start_col)(1, 2) -= 0.5 * S1 * C2 * w_add;
        scat_mat[2].block<3, 3>(start_row, start_col)(1, 2) -= 0.5 * S1 * C2 * w_minus;

        scat_mat[1].block<3, 3>(start_row, start_col)(1, 2) -= 0.5 * C1 * S2 * w_minus;
        scat_mat[2].block<3, 3>(start_row, start_col)(1, 2) -= 0.5 * C1 * S2 * w_add;

        // 2, 2 component
        scat_mat[1].block<3, 3>(start_row, start_col)(2, 2) -= 0.5 * S1 * S2 * w_add;
        scat_mat[2].block<3, 3>(start_row, start_col)(2, 2) -= 0.5 * S1 * S2 * w_minus;

        scat_mat[1].block<3, 3>(start_row, start_col)(2, 2) += 0.5 * C1 * C2 * w_minus;
        scat_mat[2].block<3, 3>(start_row, start_col)(2, 2) += 0.5 * C1 * C2 * w_add;

        for(int i = 0; i < 4; ++i) {
            scat_mat[i].block<3, 3>(start_row, start_col) *= m_incoming_sphere->quadrature_weight(in_idx);
        }


    }

    template<int NSTOKES>
    void IncomingOutgoingSpherePair<NSTOKES>::configure_geometry() {
        if(m_is_configured) {
            return;
        }

        int n_incoming = m_incoming_sphere->num_points();
        int n_outgoing = m_outgoing_sphere->num_points();

        // Loop over legendre
        for(int l = 0; l < m_legendre_scat_mats.size(); ++l) {
            for(auto& ele : m_legendre_scat_mats[l]) {
                ele.resize(n_outgoing*NSTOKES, n_incoming*NSTOKES);
                ele.setZero();
            }

            for(int inidx = 0; inidx < n_incoming; ++inidx) {
                for(int outidx = 0; outidx < n_outgoing; ++outidx) {
                    assign_scat_mat_block(l, inidx, outidx);
                }
            }
        }


        m_is_configured = true;
    }

    template<>
    void IncomingOutgoingSpherePair<1>::calculate_scattering_matrix(
            const sasktran2::atmosphere::PhaseStorage<1> &phase,
            const std::vector<std::pair<int, double>>& index_weights,
            double* phase_storage_location) const {
        Eigen::Map<Eigen::MatrixXd> phase_matrix(phase_storage_location, m_legendre_scat_mats[0][0].rows(), m_legendre_scat_mats[0][0].cols());

        phase_matrix.setZero();

        // Interpolate legendre coefficient to the location we are at
        // TODO: Delta scaling
        for(int l = 0; l < m_legendre_scat_mats.size(); ++l) {
            double leg_coeff = 0.0;

            for(const auto& ele : index_weights) {
                leg_coeff += ele.second * (phase.storage()(l, ele.first));
            }
            phase_matrix += leg_coeff * m_legendre_scat_mats[l][0];
        }
    }

    template<>
    void IncomingOutgoingSpherePair<3>::calculate_scattering_matrix(
            const sasktran2::atmosphere::PhaseStorage<3> &phase,
            const std::vector<std::pair<int, double>>& index_weights,
            double* phase_storage_location) const {
        Eigen::Map<Eigen::MatrixXd> phase_matrix(phase_storage_location, m_legendre_scat_mats[0][0].rows(), m_legendre_scat_mats[0][0].cols());

        phase_matrix.setZero();

        // Interpolate legendre coefficient to the location we are at
        // TODO: Delta scaling
        for(int l = 0; l < m_legendre_scat_mats.size(); ++l) {
            std::array<double, 4> leg_coeff({0, 0, 0, 0});

            for(const auto& ele : index_weights) {
                for(int i = 0; i < 4; ++i) {
                    leg_coeff[i] += ele.second * (phase.storage()(l*4 + i, ele.first));

                    #ifdef SASKTRAN_DEBUG_ASSERTS
                    if(leg_coeff[i] != leg_coeff[i]) {
                        BOOST_LOG_TRIVIAL(error) << l << " " << i << " " << ele.second << " " << ele.first << " " << (phase.storage()(l*4 + i, ele.first));
                        BOOST_LOG_TRIVIAL(error) << phase.storage()(Eigen::all, ele.first);

                    }
                    #endif
                }
            }
            for(int i = 0; i < 4; ++i) {
                phase_matrix += leg_coeff[i] * m_legendre_scat_mats[l][i];
            }
        }

    }

    template<int NSTOKES>
    void IncomingOutgoingSpherePair<NSTOKES>::calculate_ground_scattering_matrix(
            const sasktran2::atmosphere::Surface &surface, const std::vector<std::pair<int, double>> &index_weights,
            const sasktran2::Location& loc,
            int wavelidx, double *phase_storage_location) const {
        Eigen::Map<Eigen::MatrixXd> phase_matrix(phase_storage_location, m_legendre_scat_mats[0][0].rows(), m_legendre_scat_mats[0][0].cols());

        // TODO: update for BRDF

        double albedo = surface.albedo()[wavelidx];

        // scattering matrix elements are albedo/pi * cos(theta),
        // but quadrature is sut up for 4pi normalization so we multiply this by 4pi

        phase_matrix.setZero();


        // right now just lambertian, so only scalar
        for(int i = 0; i < phase_matrix.cols(); i += NSTOKES) {
            double cos_theta = loc.cos_zenith_angle(m_incoming_sphere->get_quad_position(i/NSTOKES));

            for(int j = 0; j < phase_matrix.rows(); j += NSTOKES) {
                phase_matrix(j, i) = 4*albedo * cos_theta * m_incoming_sphere->quadrature_weight(i/NSTOKES);
            }
        }
    }

    template<int NSTOKES>
    DiffusePoint<NSTOKES>::DiffusePoint(const IncomingOutgoingSpherePair<NSTOKES> &spheres, const sasktran2::Location& location) :
    m_spheres(&spheres), m_location(location) {

    }

    template class IncomingOutgoingSpherePair<1>;
    template class IncomingOutgoingSpherePair<3>;

    template class DiffusePoint<1>;
    template class DiffusePoint<3>;
}