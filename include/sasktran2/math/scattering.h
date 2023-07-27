#pragma once

#include "sasktran2/internal_common.h"
#include <sasktran2/math/wigner.h>


namespace sasktran2::math {
    #define SCATTERING_EPSILSON 1e-8

    inline void stokes_scattering_factors(const Eigen::Vector3d& incoming,
                                          const Eigen::Vector3d& outgoing,
                                          double& theta,
                                          double& C1,
                                          double& C2,
                                          double& S1,
                                          double& S2,
                                          int& negation
                                          ) {
        // Follows Mishchenko section 4.3, Hovenier 3.2

        // Start by calculating the scattering angle and cos/sin of it
        // Note that if the incoming and outgoing directions are the same, the scattering angle is 0 since
        // These are both propagation directions
        double cos_scatter = incoming.dot(outgoing);

        if(cos_scatter > 1) {
            cos_scatter = 1;
        }
        if(cos_scatter < -1) {
            cos_scatter = -1;
        }

        theta = acos(cos_scatter);

        double sin_scatter = sin(theta);

        negation = 1;

        // Special case if sin of the scattering angle is 0, or 180 degrees
        if(abs(sin_scatter) < SCATTERING_EPSILSON) {
            // then there is no rotation, so cos angles are 1, sin angles ar 0
            C1 = 1;
            C2 = 1;
            S1 = 0;
            S2 = 0;

            return;
        }

        // Mishchenko 2000 equation 4.15 and following

        // inc is angle between incoming and z direction
        double costh_inc = incoming.z();

        // scat is angle between outgoing and z direction
        double costh_scat = outgoing.z();

        // Also can get slight rounding errors here
        if(costh_inc > 1) {
            costh_inc = 1;
        }
        if(costh_inc < -1) {
            costh_inc = -1;
        }

        if(costh_scat > 1) {
            costh_scat = 1;
        }
        if(costh_scat < -1) {
            costh_scat = -1;
        }

        double sinth_inc = sin(acos(costh_inc));
        double sinth_scat = sin(acos(costh_scat));

        double costh1, costh2;
        // Now there are two special cases we have to handle

        // If the incoming direction is the z-direction, then sinth_inc = 0


        double phi_inc, phi_scat, phi_diff;

        if(abs(sinth_inc) < SCATTERING_EPSILSON) {
            // In this case th1 is 180 degrees, so
            costh1 = -1;
            phi_inc = 0;

            C1 = 1;
            C2 = 1;
            S1 = 0;
            S2 = 0;

            return;

        } else {
            costh1 = (costh_scat - costh_inc * cos_scatter) / (sinth_inc * sin_scatter);

            Eigen::Vector3d horiz_inc = incoming;
            horiz_inc.z() = 0;
            horiz_inc = horiz_inc.normalized();

            phi_inc = atan2(horiz_inc.y(), horiz_inc.x());

        }

        // or if the scattering direction is the z-direction then sinth_scat = 0
        if(abs(sinth_scat) < SCATTERING_EPSILSON) {
            // Here theta2 is 180 degrees, so
            costh2 = -1;
            phi_scat = 0;

            C1 = 1;
            C2 = 1;
            S1 = 0;
            S2 = 0;

            return;


        } else {
            costh2 = (costh_inc - costh_scat * cos_scatter) / (sinth_scat * sin_scatter);

            Eigen::Vector3d horiz_scat = outgoing;
            horiz_scat.z() = 0;
            horiz_scat = horiz_scat.normalized();

            phi_scat = atan2(horiz_scat.y(), horiz_scat.x());
        }

        phi_diff = phi_inc - phi_scat;
        if(phi_diff < 0) {
            phi_diff += 2*EIGEN_PI;
        }

        // Can get rounding errors here
        if(costh1 > 1) {
            costh1 = 1;
        }
        if(costh1 < -1) {
            costh1 = -1;
        }

        if(costh2 > 1) {
            costh2 = 1;
        }
        if(costh2 < -1) {
            costh2 = -1;
        }

        C1 = 2*costh1*costh1 - 1;
        C2 = 2*costh2*costh2 - 1;

        S1 = 2 * sqrt(1 - costh1*costh1) * costh1;
        S2 = 2 * sqrt(1 - costh2*costh2) * costh2;

        // We now have to adjust S1,S2 based on the azimuthal differences of the incoming
        if(phi_diff < EIGEN_PI) {
            S1 *= -1;
            S2 *= -1;
        }

        double cos_scat2 = costh_inc*costh_scat + sinth_inc*sinth_scat*cos(phi_inc - phi_scat);

        #ifdef SASKTRAN_DEBUG_ASSERTS
        if(C1 != C1 || C2 != C2 || S1 != S1 || S2 != S2 || abs(C2*C2 + S2*S2 - 1) > 1e-8 || abs(cos_scatter - cos_scat2) > 1e-8) {
            static bool message = true;
            if(message) {
                BOOST_LOG_TRIVIAL(error) << "Error generating scattering matrix elements";

                BOOST_LOG_TRIVIAL(error) << "S1: " << S1;
                BOOST_LOG_TRIVIAL(error) << "S2: " << S2;
                BOOST_LOG_TRIVIAL(error) << "C1: " << C1;
                BOOST_LOG_TRIVIAL(error) << "C2: " << C2;

                BOOST_LOG_TRIVIAL(error) << "incoming: " << incoming;
                BOOST_LOG_TRIVIAL(error) << "outgoing: " << outgoing;

                BOOST_LOG_TRIVIAL(error) << "costh1: " << costh1;
                BOOST_LOG_TRIVIAL(error) << "costh2: " << costh2;

                BOOST_LOG_TRIVIAL(error) << "sin_scatter: " << sin_scatter;
                BOOST_LOG_TRIVIAL(error) << "cos_scatter: " << cos_scatter;

                BOOST_LOG_TRIVIAL(error) << "phi_inc: " << phi_inc;
                BOOST_LOG_TRIVIAL(error) << "phi_scat: " << phi_scat;
            }
            message = false;
        }
        #endif

    }

}