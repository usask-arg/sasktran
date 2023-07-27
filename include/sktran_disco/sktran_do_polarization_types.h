#pragma once
#include "sktran_disco/sktran_do.h"

namespace sasktran_disco {
    // Many vector equations use a matrix D = diag(1, 1, -1, -1) that accounts for the symmetry of the Legendre
    // functions.  This implements the diagonal elements of this matrix, templated to not slow down the scalar code
    template <int NSTOKES, int CNSTR=-1>
    double stokes_negation_factor(int linear_index) {
        int stokes_index = linear_index % NSTOKES;

        if (stokes_index >= 2) {
            return -1.0;
        }
        else {
            return 1.0;
        }
    }

    // Class to calculate the Wigner D functions
    // Likely will have to move this into SasktranCore, or copy it, so that we can use it for the optical
    // properties
    // Uses Recurrence relations in Mishchenko (2006) to calculate the functions, all equations refer to this book
    // We use the different notation d^l_mn instead of d^s_mn to match the associate legendre polynomial notations
    class WignerDCalculator {
    private:
        int m_m;
        int m_n;
        int m_lmin;

        double m_recurrence_start_factor;

        int m_zeta;

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

        double recurrence_start(double theta) {
            double x = cos(theta);
            double xfactor = std::pow(1-x, double(std::abs(int(m_m) - int(m_n))) / 2.0) *
                   std::pow(1+x, double(std::abs(int(m_m) + int(m_n))) / 2.0);

            return m_recurrence_start_factor * xfactor;
        }

    public:
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

    template <class T>
    class LegendreCoefficientInterface {

    };

    template <int NSTOKES, int CNSTR=-1>
    class LegendreCoefficientContainer {

        LegendreCoefficientContainer() {
            // Error?
        }
    };

    // Container when operating in Scalar mode, there is only a Legendre expansion of
    // the scalar P11 phase function element and so this container is basically a wrapper
    // around a 1D array
    template<>
    class LegendreCoefficientContainer<1> : LegendreCoefficientInterface<LegendreCoefficientContainer<1>> {
    private:
        Eigen::VectorXd m_storage;
    public:
        LegendreCoefficientContainer(size_t size) {
            resize(size);
        }

        LegendreCoefficientContainer() {};

        void resize(size_t size) {
            m_storage.resize(size);
        }


        LegendreCoefficientContainer<1> operator*(const double& lhs) {
            LegendreCoefficientContainer<1> container(m_storage.size());

            container.data() = lhs * this->data();

            return container;
        }


        Eigen::VectorXd& data() { return m_storage; }
        const Eigen::VectorXd& data() const { return m_storage; }
    };

    // Container when operating in fully vector mode.  The legendre expansion matrix is a 4x4 matrix
    // of coefficients, but only 6 elements are unique, the matrix looks like
    // [ a1 -b1  0   0 ]
    // [-b1  a2  0   0 ]
    // [  0   0 a3 -b2 ]
    // [  0   0 b2  a4 ]
    // where we use the internal storage indexing
    // 0 - a1
    // 1 - a2
    // 2 - a3
    // 3 - a4
    // 4 - b1
    // 5 - b2

    template<>
    class LegendreCoefficientContainer<4> : LegendreCoefficientInterface<LegendreCoefficientContainer<4>> {
    private:
        Eigen::MatrixXd m_storage;
    public:
        LegendreCoefficientContainer(size_t size) {
            resize(size);
        }

        LegendreCoefficientContainer() {};

        void resize(size_t size) {
            m_storage.resize(size, 6);
        }

        Eigen::Ref<Eigen::MatrixXd> a1() {
            return m_storage(Eigen::all, 0);
        }

        Eigen::Ref<Eigen::MatrixXd> a2() {
            return m_storage(Eigen::all, 1);
        }

        Eigen::Ref<Eigen::MatrixXd> a3() {
            return m_storage(Eigen::all, 2);
        }

        Eigen::Ref<Eigen::MatrixXd> a4() {
            return m_storage(Eigen::all, 3);
        }

        Eigen::Ref<Eigen::MatrixXd> b1() {
            return m_storage(Eigen::all, 4);
        }

        Eigen::Ref<Eigen::MatrixXd> b2() {
            return m_storage(Eigen::all, 5);
        }

        Eigen::MatrixXd& data() { return m_storage; }
        const Eigen::MatrixXd& data() const { return m_storage; }

    };

    // For the special case of NSTOKES=3, we have a4=b2=0, and so we only have 4 greek coefficients
    template<>
    class LegendreCoefficientContainer<3> : LegendreCoefficientInterface<LegendreCoefficientContainer<3>> {
    private:
        Eigen::MatrixXd m_storage;
    public:
        LegendreCoefficientContainer(size_t size) {
            resize(size);
        }

        LegendreCoefficientContainer() {};

        void resize(size_t size) {
            m_storage.resize(size, 4);
        }

        Eigen::Ref<Eigen::MatrixXd> a1() {
            return m_storage(Eigen::all, 0);
        }

        Eigen::Ref<Eigen::MatrixXd> a2() {
            return m_storage(Eigen::all, 1);
        }

        Eigen::Ref<Eigen::MatrixXd> a3() {
            return m_storage(Eigen::all, 2);
        }

        Eigen::Ref<Eigen::MatrixXd> b1() {
            return m_storage(Eigen::all, 4);
        }

        Eigen::MatrixXd& data() { return m_storage; }
        const Eigen::MatrixXd& data() const { return m_storage; }

    };


    // Type to store a single Legendre coefficient value
    template <int NSTOKES, int CNSTR=-1>
    struct LegendreCoefficient {

    };

    // NSTOKES = 1 we only have 1 Legendre coefficient
    template <>
    struct LegendreCoefficient<1> {
        double a1;

        LegendreCoefficient() {
            a1 = 0.0;
        }
    };

    // NSTOKES = 4 we have 6 Legendre coefficients
    template <>
    struct LegendreCoefficient<4> {
        double a1;
        double a2;
        double a3;
        double a4;
        double b1;
        double b2;

        LegendreCoefficient() {
            a1 = 0.0;
            a2 = 0.0;
            a3 = 0.0;
            a4 = 0.0;
            b1 = 0.0;
            b2 = 0.0;
        }

    };

    // NSTOKES = 3 we have 4 legendre coefficients
    template <>
    struct LegendreCoefficient<3> {
        double a1;
        double a2;
        double a3;
        double b1;

        LegendreCoefficient() {
            a1 = 0.0;
            a2 = 0.0;
            a3 = 0.0;
            b1 = 0.0;
        }

    };


    // Class to store the Legendre functions evaluated at stream angles, see eq A.6 in Rozanov et. al 2013
    // Storage is 3 elements, [P R T] in the general case, we need these three elements for NSTOKES=3 as well
    // so we leave this as the full general case
    template <int NSTOKES>
    struct LegendrePhaseContainer {
        Eigen::Vector3d values;

        void fill(int m, int l, double coszen) {
            auto calculatorneg = sasktran_disco::WignerDCalculator(m, -2);
            auto calculatorpos = sasktran_disco::WignerDCalculator(m, 2);
            double theta = acos(coszen);


            values(1) = -0.5 * (calculatorpos.d(theta, l) + calculatorneg.d(theta, l));
            values(2) = -0.5 * (calculatorpos.d(theta, l) - calculatorneg.d(theta, l));

            auto calculator = WignerDCalculator(m, 0);

            values(0) = calculator.d(theta, l);
        }

        double& P() {
            return values(0);
        }

        double& R() {
            return values(1);
        }

        double& T() {
            return values(2);
        }

        double P() const {
            return values(0);
        }

        double R() const {
            return values(1);
        }

        double T() const {
            return values(2);
        }
    };

    // Storage is simply a scalar element for the NSTOKES=1 special case
    template <>
    struct LegendrePhaseContainer<1> {
        double value;

        void fill(int m, int l, double coszen) {
            auto calculator = WignerDCalculator(m, 0);
            double theta = acos(coszen);

            value = calculator.d(theta, l);
        }
        double& P() {
            return value;
        }

        double P() const {
            return value;
        }

    };

    template <int NSTOKES, int CNSTR=-1>
    class LayerInputDerivative;

    // Holds the quantities calculated by the LPETripleProduct class and their derivatives
    // The value is a NSTOKES x NSTOKES matrix, where each value in theory has derivatives with respect to every
    // greek parameter.  For NSTOKES=4 this would result in a 4x4 matrix where we need derivatives with respect to
    // 6 quantities for every l.  In reality the derivatives are sparse, alpha1 only affects the 1,1 element of the matrix for
    // example.  This class stores the NSTOKESxNOSTKES matrix values efficiently and the capability to propagate derivatives
    // of the greek parameters efficiently which requires specialized instantiations for every value of NSTOKES
    template <int NSTOKES, int CNSTR=-1>
    class TripleProductDerivativeHolder {
    };

    template <>
    class TripleProductDerivativeHolder<1> {
    public:
        TripleProductDerivativeHolder(){}

        TripleProductDerivativeHolder(int nstr) : nstr(nstr) {
            resize(nstr);
        }

        void resize(int nstr) {
            this->nstr = nstr;
            d_by_legendre_coeff.resize(nstr);
        }

        void calculate(const std::vector<LegendreCoefficient<1>>& coeffs, const std::vector<LegendrePhaseContainer<1>>& lp1s, const std::vector<LegendrePhaseContainer<1>>& lp2s,
                       bool negation, int m) {
            value = 0.0;
            d_by_legendre_coeff.setZero();
            for(int l = m; l < nstr; ++l) {
                const auto &coeff = coeffs[l];
                const auto &lp1 = lp1s[l];
                const auto &lp2 = lp2s[l];

                int negation_factor = 1;
                if (negation) {
                    if ((l - m) % 2 != 0) {
                        negation_factor *= -1;
                    }
                }
                value += coeff.a1 * lp1.P() * lp2.P() * negation_factor;
                d_by_legendre_coeff[l] += lp1.P() * lp2.P() * negation_factor;
            }
        }

        double value;
        double ssa;
        int nstr;
        Eigen::VectorXd d_by_legendre_coeff;
        void reduce(const LayerInputDerivative<1>& layer_deriv, double& deriv) const;

    };

    // The values are a 4x4 matrix, the derivative factors for the 6 greek constants are
    // a1 - 11
    // a2 - 22, 23, 32, 33
    // a3 - 22, 23, 32, 33
    // a4 - 44
    // b1 - 12, 13, 21, 31
    // b2 - 24, 34, 42, 43
    //
    template <>
    class TripleProductDerivativeHolder<4> {
    public:
        Eigen::Matrix<double, 4, 4> value;
        Eigen::VectorXd a1deriv;
        Eigen::MatrixX4d a2deriv;
        Eigen::MatrixX4d a3deriv;
        Eigen::VectorXd a4deriv;
        Eigen::MatrixX4d b1deriv;
        Eigen::MatrixX4d b2deriv;
        int nstr;
        double ssa;

        TripleProductDerivativeHolder(){}

        TripleProductDerivativeHolder(int nstr) : nstr(nstr) {
            resize(nstr);
        }

        void resize(int nstr) {
            this->nstr = nstr;

            a1deriv.resize(nstr);
            a2deriv.resize(nstr, 4);
            a3deriv.resize(nstr, 4);
            a4deriv.resize(nstr);
            b1deriv.resize(nstr, 4);
            b2deriv.resize(nstr, 4);

            // Have to set to 0 because sum only goes from l=m upwards, might be faster to do the full
            // sum and not set to 0?
            a1deriv.setZero();
            a2deriv.setZero();
            a3deriv.setZero();
            a4deriv.setZero();
            b1deriv.setZero();
            b2deriv.setZero();
        }

        void calculate(const std::vector<LegendreCoefficient<4>>& coeffs, const std::vector<LegendrePhaseContainer<4>>& lp1s, const std::vector<LegendrePhaseContainer<4>>& lp2s,
                       bool negation, int m) {
            value.setZero();

            for(int l = m; l < nstr; ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                int negation_factor_upper = 1;
                int negation_factor_lower = 1;
                if( negation ) {
                    negation_factor_lower = -1;
                    if( (l-m) % 2 != 0) {
                        negation_factor_upper *= -1;
                        negation_factor_lower *= -1;
                    }
                }
                // Calculated product by hand
                value(0, 0) += lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                value(0, 1) += -1.0 * lp1.P() * lp2.R() * coeff.b1 * negation_factor_upper;
                value(0, 2) += 1.0 * lp1.P() * lp2.T() * coeff.b1 * negation_factor_upper;
                // 0, 3 is always 0

                value(1, 0) += -1.0 * lp1.R() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(1, 1) += lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper + lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(1, 2) += -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper - lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;
                value(1, 3) += 1.0 * lp1.T() * lp2.P() * coeff.b2 * negation_factor_lower;

                value(2, 0) += lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(2, 1) += -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper - lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(2, 2) += lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper + lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;
                value(2, 3) += -lp1.R() * lp2.P() * coeff.b2 * negation_factor_lower;

                // 3, 0 is always 0
                value(3, 1) += -lp1.P() * lp2.T() * coeff.b2 * negation_factor_lower;
                value(3, 2) += lp1.P() * lp2.R() * coeff.b2 * negation_factor_lower;
                value(3, 3) += lp1.P() * lp2.P() * coeff.a4 * negation_factor_lower;

                // Assign derivatives
                a1deriv(l) = lp1.P() * lp2.P() * negation_factor_upper;

                a2deriv(l, 0) = lp1.R() * lp2.R() * negation_factor_upper;
                a2deriv(l, 1) = -lp1.R() * lp2.T() * negation_factor_upper;
                a2deriv(l, 2) = -lp1.T() * lp2.R() * negation_factor_upper;
                a2deriv(l, 3) = lp1.T() * lp2.T() * negation_factor_upper;

                a3deriv(l, 0) = lp1.T() * lp2.T() * negation_factor_lower;
                a3deriv(l, 1) = -lp1.T() * lp2.R() * negation_factor_lower;
                a3deriv(l, 2) = -lp1.R() * lp2.T() * negation_factor_lower;
                a3deriv(l, 3) = lp1.R() * lp2.R() * negation_factor_lower;

                b1deriv(l, 0) = -1.0 * lp1.P() * lp2.R() * negation_factor_upper;
                b1deriv(l, 1) = lp1.P() * lp2.T() * negation_factor_upper;
                b1deriv(l, 2) = -1.0 * lp1.R() * lp2.P() * negation_factor_upper;
                b1deriv(l, 3) = lp1.T() * lp2.P() * negation_factor_upper;

                b2deriv(l, 0) = lp1.T() * lp2.P() * negation_factor_lower;
                b2deriv(l, 1) = -lp1.R() * lp2.P() * negation_factor_lower;
                b2deriv(l, 2) = -lp1.P() * lp2.T() * negation_factor_lower;
                b2deriv(l, 3) = lp1.P() * lp2.R() * negation_factor_lower;

                a4deriv(l) = lp1.P() * lp2.P() * negation_factor_lower;
            }
        }

        void reduce(const LayerInputDerivative<4>& layer_deriv, Eigen::Matrix<double, 4, 4>& deriv) const;
    };

    // The values are a 3x3 matrix, the derivative factors for the 4 greek constants are
    // a1 - 11
    // a2 - 22, 23, 32, 33
    // a3 - 22, 23, 32, 33
    // b1 - 12, 13, 21, 31
    template <>
    class TripleProductDerivativeHolder<3> {
    public:
        Eigen::Matrix<double, 3, 3> value;
        Eigen::VectorXd a1deriv;
        Eigen::MatrixX4d a2deriv;
        Eigen::MatrixX4d a3deriv;
        Eigen::MatrixX4d b1deriv;
        int nstr;
        double ssa;

        TripleProductDerivativeHolder(){}

        TripleProductDerivativeHolder(int nstr) : nstr(nstr) {
            resize(nstr);
        }

        void resize(int nstr) {
            this->nstr = nstr;

            a1deriv.resize(nstr);
            a2deriv.resize(nstr, 4);
            a3deriv.resize(nstr, 4);
            b1deriv.resize(nstr, 4);

            // Used to set these to 0 and then sum from l=m to NSTR, but now we reuse memory from this class
            // across calculations and so we would have to re-zero everytime calculate is called
            // Instead we just sum from l=0 to nstr everytime and this handles the zeroing
            //a1deriv.setZero();
            //a2deriv.setZero();
            //a3deriv.setZero();
            //b1deriv.setZero();
        }

        void calculate(const std::vector<LegendreCoefficient<3>>& coeffs, const std::vector<LegendrePhaseContainer<3>>& lp1s, const std::vector<LegendrePhaseContainer<3>>& lp2s,
            bool negation, int m) {

            value.setZero();
            for (int l = 0; l < nstr; ++l) {
                const auto& coeff = coeffs[l];
                const auto& lp1 = lp1s[l];
                const auto& lp2 = lp2s[l];

                int negation_factor_upper = 1;
                int negation_factor_lower = 1;
                if (negation) {
                    negation_factor_lower = -1;
                    if ((l - m) % 2 != 0) {
                        negation_factor_upper *= -1;
                        negation_factor_lower *= -1;
                    }
                }

                // Calculated product by hand
                value(0, 0) += lp1.P() * lp2.P() * coeff.a1 * negation_factor_upper;
                value(0, 1) += -1.0 * lp1.P() * lp2.R() * coeff.b1 * negation_factor_upper;
                value(0, 2) += 1.0 * lp1.P() * lp2.T() * coeff.b1 * negation_factor_upper;
                // 0, 3 is always 0

                value(1, 0) += -1.0 * lp1.R() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(1, 1) += lp1.R() * lp2.R() * coeff.a2 * negation_factor_upper + lp1.T() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(1, 2) += -lp1.R() * lp2.T() * coeff.a2 * negation_factor_upper - lp1.T() * lp2.R() * coeff.a3 * negation_factor_lower;

                value(2, 0) += lp1.T() * lp2.P() * coeff.b1 * negation_factor_upper;
                value(2, 1) += -lp1.T() * lp2.R() * coeff.a2 * negation_factor_upper - lp1.R() * lp2.T() * coeff.a3 * negation_factor_lower;
                value(2, 2) += lp1.T() * lp2.T() * coeff.a2 * negation_factor_upper + lp1.R() * lp2.R() * coeff.a3 * negation_factor_lower;


                // Assign derivatives
                a1deriv(l) = lp1.P() * lp2.P() * negation_factor_upper;

                a2deriv(l, 0) = lp1.R() * lp2.R() * negation_factor_upper;
                a2deriv(l, 1) = -lp1.R() * lp2.T() * negation_factor_upper;
                a2deriv(l, 2) = -lp1.T() * lp2.R() * negation_factor_upper;
                a2deriv(l, 3) = lp1.T() * lp2.T() * negation_factor_upper;

                a3deriv(l, 0) = lp1.T() * lp2.T() * negation_factor_lower;
                a3deriv(l, 1) = -lp1.T() * lp2.R() * negation_factor_lower;
                a3deriv(l, 2) = -lp1.R() * lp2.T() * negation_factor_lower;
                a3deriv(l, 3) = lp1.R() * lp2.R() * negation_factor_lower;

                b1deriv(l, 0) = -1.0 * lp1.P() * lp2.R() * negation_factor_upper;
                b1deriv(l, 1) = lp1.P() * lp2.T() * negation_factor_upper;
                b1deriv(l, 2) = -1.0 * lp1.R() * lp2.P() * negation_factor_upper;
                b1deriv(l, 3) = lp1.T() * lp2.P() * negation_factor_upper;
            }
        }

        void reduce(const LayerInputDerivative<3>& layer_deriv, Eigen::Matrix<double, 3, 3>& deriv) const;
    };

    // This class holds the value and derivative for the solar source which is of size NSTOKES
    template <int NSTOKES, int CNSTR=-1>
    class InhomogeneousSourceHolder {
    public:
        Eigen::Vector<double, NSTOKES> value;
        Eigen::VectorXd d_by_a1;
        Eigen::VectorXd d_by_b1_first;
        Eigen::VectorXd d_by_b1_second;
        Eigen::Vector<double, NSTOKES> d_by_ssa;
        int nstr;

        InhomogeneousSourceHolder() {}

        InhomogeneousSourceHolder(int nstr) : nstr(nstr) {
            resize(nstr);
        }

        void resize(int nstr) {
            this->nstr = nstr;
            d_by_a1.resize(nstr);
            d_by_b1_first.resize(nstr);
            d_by_b1_second.resize(nstr);
        }

        void reduce(const LayerInputDerivative<NSTOKES>& layer_deriv, Eigen::Vector<double, NSTOKES>& deriv) const;
    };

    template <>
    class InhomogeneousSourceHolder<1> {
    public:
        double value;
        Eigen::VectorXd d_by_legendre_coeff;
        double d_by_ssa;
        int nstr;

        InhomogeneousSourceHolder() {}

        InhomogeneousSourceHolder(int nstr) : nstr(nstr) {
            resize(nstr);
        }

        void resize(int nstr) {
            this->nstr = nstr;
            d_by_legendre_coeff.resize(nstr);
        }

        void reduce(const LayerInputDerivative<1>& layer_deriv, double& deriv) const;
    };

    // Forward declaration
    template <typename T>
    struct Dual;

    // Class to store a stokes vector (or scalar radiance) and it's derivative
    // We specialize the scalar case to use a double instead of an eigen object
    // (or scalars)
    template <int NSTOKES, int CNSTR=-1>
    struct Radiance {
        Eigen::Vector<double, NSTOKES> value;
        Eigen::Matrix<double, -1, NSTOKES> deriv;

        Radiance() {}

        Radiance(int nderiv, bool zero = true) {
            resize(nderiv, zero);
        }

        void resize(int nderiv, bool zero = true) {
            deriv.resize(nderiv, NSTOKES);

            if (zero) {
                value.setZero();
                deriv.setZero();
            }
        }

        void setzero() {
            value.setZero();
            deriv.setZero();
        }

        void apply_transmission_factor(const Dual<double>& transmission);
        void apply_azimuth_expansion(double angle, int m);
        bool converged(double I, double epsilon);

        double I() const {
            return value(0);
        }
    };

    template <>
    struct Radiance<1> {
        double value;
        Eigen::VectorXd deriv;

        Radiance() {}

        Radiance(int nderiv, bool zero = true) {
            resize(nderiv, zero);
        }

        void resize(int nderiv, bool zero = true) {
            deriv.resize(nderiv);

            if (zero) {
                value = 0.0;
                deriv.setZero();
            }
        }

        void setzero() {
            value = 0.0;
            deriv.setZero();
        }

        void apply_transmission_factor(const Dual<double>& transmission);
        void apply_azimuth_expansion(double angle, int m);
        bool converged(double I, double epsilon);

        double I() const {
            return value;
        }
    };

    template <int NSTOKES, int CNSTR=-1>
    class InputDerivatives;

    // Maps a dual defined with respect to layer quantities to the dual defined with respect to weighting function (
    // engine input) quantities
    template <int NSTOKES, int CNSTR=-1>
    inline Radiance<NSTOKES> convert_dual_to_wf(const Radiance<NSTOKES>& dual, const InputDerivatives<NSTOKES>& in_deriv,
        size_t numwf) {

        Radiance<NSTOKES> result(numwf);

        for (unsigned int l = 0; l < in_deriv.numDerivative(); ++l)
        {
            const auto& qty = in_deriv.layerDerivatives()[l];
            for (unsigned int k = 0; k < qty.group_and_triangle_fraction.size(); ++k) {
                result.deriv(qty.group_and_triangle_fraction[k].first, Eigen::all).noalias() += dual.deriv(l, Eigen::all) * qty.group_and_triangle_fraction[k].second * qty.extinctions[k];
            }
        }

        result.value = dual.value;
        return result;
    }

    template <int NSTOKES, int CNSTR=-1>
    inline void convert_dual_to_wf(const Radiance<NSTOKES>& dual, const InputDerivatives<NSTOKES>& in_deriv, Radiance<NSTOKES>& result) {
        result.deriv.setZero();
        for (unsigned int l = 0; l < in_deriv.numDerivative(); ++l)
        {
            const auto& qty = in_deriv.layerDerivatives()[l];
            for (unsigned int k = 0; k < qty.group_and_triangle_fraction.size(); ++k) {
                result.deriv(qty.group_and_triangle_fraction[k].first, Eigen::all).noalias() += dual.deriv(l, Eigen::all) * qty.group_and_triangle_fraction[k].second * qty.extinctions[k];
            }
        }

        result.value = dual.value;
    }

}