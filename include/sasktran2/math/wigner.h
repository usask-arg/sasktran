#pragma once

#include <sasktran2/internal_common.h>

namespace sasktran2::math {

    /** Calculates the Wigner D functions using recurrence relations found from Mischenko.  The "first" Wigner D
     * function is the standard Legendre polynomial, so this class is used to calculate Legendre polynomials as well.
     * The other Wigner D functions are only used for polarized scattering calculations.
     *
     */
    class WignerDCalculator {
    private:
        int m_m;
        int m_n;
        int m_lmin;

        double m_recurrence_start_factor;

        int m_zeta;

        /** Internal method to determine the start scaling factor of the recurrence
         *
         * @return The start scaling factor of the recurrence relation for the given internal m_m, m_n, and m_lmin
         */
        double recurrence_start_factor() {
            // Start of the recurrence, i.e. d^lmin_mn(theta)
            // Eq. F.4.3, the factor that does not depend on x

            // Start by calculating the factorial chain
            int num = 2*m_lmin;
            int den1 = std::abs(int(m_m) - int(m_n));
            int den2 = std::abs(int(m_m) + int(m_n));

            // We need to calculate num! / (den1! * den2!)
            double factorial = 1;
            // The numerator will always be >= the denominator so we start there
            for(int i = num; i > 1; --i) {
                factorial *= double(i);
                if (i <= den1) {
                    factorial /= double(i);
                }
                if (i <= den2) {
                    factorial /= double(i);
                }
            }
            double otherfactor = std::pow(2.0, -double(m_lmin));

            return m_zeta * otherfactor * sqrt(factorial);
        }

        /** Internal method to calculate the recurrence start value
         *
         * @param theta Angle in radians
         * @return THe recurrence start value
         */
        double recurrence_start(double theta) {
            double x = cos(theta);
            double xfactor = std::pow(1-x, double(std::abs(int(m_m) - int(m_n))) / 2.0) *
                             std::pow(1+x, double(std::abs(int(m_m) + int(m_n))) / 2.0);

            return m_recurrence_start_factor * xfactor;
        }

    public:
        /** Constructs the calculator for \f$d^l_{mn}\f$ for a given m and n
         *
         * @param m
         * @param n
         */
        WignerDCalculator(int m, int n) : m_m(m), m_n(n) {
            // Calculate zeta from F.2.4
            if(m_n >= m_m) {
                m_zeta = 1;
            } else {
                if( (m_m - m_n) % 2 == 0 ) {
                    m_zeta = 1;
                } else {
                    m_zeta = -1;
                }
            }
            m_lmin = std::max(std::abs(int(m_m)), std::abs(int(m_n)));
            m_recurrence_start_factor = recurrence_start_factor();
        };

        /** Calculates \f$d^l_{mn}(\theta)\f$
         *
         * @param theta Angle in radians
         * @param l
         * @return \f$d^l_{mn}(\theta)\f$
         */
        double d(double theta, int l) {
            if( l < m_lmin ) {
                return 0.0;
            } else {
                double x = cos(theta);
                double val_l = recurrence_start(theta); // Value at current l
                double val_lm1 = 0.0; // Value at current l minus 1

                // We use F.4.1 if n != 0, but if n == 0 the recurrence has a 0/0 singularity that we must remove and
                // the recurrence reduces to the standard associated legendre function recurrence relation

                if(m_n == 0) {
                    for(int lidx = m_lmin + 1; lidx <= l; ++lidx) {
                        double multiplier = 1.0 /  (sqrt(lidx*lidx - m_m*m_m) * lidx);
                        double curfactor = (2*lidx - 1) * (lidx * x);
                        double priorfactor = lidx * sqrt((lidx-1)*(lidx-1) - m_m*m_m);

                        double temp = val_l;
                        val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1);
                        val_lm1 = temp;
                    }
                } else {
                    for(int lidx = m_lmin + 1; lidx <= l; ++lidx) {
                        double multiplier = 1.0 /  ((lidx-1) * sqrt(double(lidx*lidx - m_m*m_m)) * sqrt(double(lidx*lidx - m_n*m_n)));
                        double curfactor = (2*lidx - 1) * (lidx * (lidx - 1) * x -  m_n * m_m);
                        double priorfactor = lidx * sqrt(double((lidx-1)*(lidx-1) - m_m*m_m)) * sqrt(double((lidx - 1)*(lidx-1) - m_n*m_n));

                        double temp = val_l;
                        val_l = multiplier * (curfactor * val_l - priorfactor * val_lm1);
                        val_lm1 = temp;
                    }
                }


                return val_l;
            }
        }
    };
}
