#include "sktran_disco/sktran_do.h"
#include "sktran_disco/sktran_do_linalg.h"


int sasktran_disco::la::dgbsv_pentadiagonal(int N, int NRHS, double* AB, double* B, int LDB) {
    Eigen::Map<Eigen::MatrixXd> eigen_AB(AB, 7, N);

	// Rows 2-7 contain the diagonals of the pentadiagonal matrix, rows 0 and 1 are extra storage
	// we can use

    // c is offset by 1 index
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<7>> c(eigen_AB.data() + 5, N);
    int c_offset = 1;

    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<7>> d(eigen_AB.data() + 4, N);
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<7>> a(eigen_AB.data() + 10, N);

    // e is offset by 2 indicies
    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<7>> e(eigen_AB.data() + 6, N);
    int e_offset = 2;

    Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<7>> b(eigen_AB.data() + 16, N);

    Eigen::VectorXd alpha(N);
	Eigen::VectorXd beta(N);
	Eigen::MatrixXd z(N, NRHS);
	Eigen::VectorXd gamma(N);
	Eigen::VectorXd mu(N);

	Eigen::Map<Eigen::MatrixXd> y(B, N, NRHS);


	mu(0) = d(0);
	alpha(0) = a(0) / mu(0);

	beta(0) = b(0) / mu(0);
	z(0, Eigen::all) = y(0, Eigen::all) / mu(0);

	gamma(1) = c(1 - c_offset);
	mu(1) = d(1) - alpha(0) * gamma(1);
	alpha(1) = (a(1) - beta(0) * gamma(1)) / mu(1);
	beta(1) = b(1) / mu(1);
	z(1, Eigen::all) = (y(1, Eigen::all) - z(0, Eigen::all) * gamma(1)) / mu(1);

	for (int i = 2; i < N - 2; ++i) {
		gamma(i) = c(i - c_offset) - alpha(i - 2) * e(i - e_offset);
		mu(i) = d(i) - beta(i - 2) * e(i - e_offset) - alpha(i - 1) * gamma(i);
		alpha(i) = (a(i) - beta(i - 1) * gamma(i)) / mu(i);
		beta(i) = b(i) / mu(i);
		z(i, Eigen::all) = (y(i, Eigen::all) - z(i - 2, Eigen::all) * e(i - e_offset) - z(i - 1, Eigen::all) * gamma(i)) / mu(i);
	}

    if( N >= 4) {
        gamma(N - 2) = c(N - 2 - c_offset) - alpha(N - 4) * e(N - 2 - e_offset);
        mu(N - 2) = d(N - 2) - beta(N - 4) * e(N - 2 - e_offset) - alpha(N - 3) * gamma(N - 2);
        alpha(N - 2) = (a(N - 2) - beta(N - 3) * gamma(N - 2)) / mu(N - 2);

        gamma(N - 1) = c(N - 1 - c_offset) - alpha(N - 3) * e(N - 1 - e_offset);
        mu(N - 1) = d(N - 1) - beta(N - 3) * e(N - 1 - e_offset) - alpha(N - 2) * gamma(N - 1);
    }

	z(N - 2, Eigen::all) = (y(N - 2, Eigen::all) - z(N - 4, Eigen::all) * e(N - 2 - e_offset) - z(N - 3, Eigen::all) * gamma(N - 2)) / mu(N - 2);
	z(N - 1, Eigen::all) = (y(N - 1, Eigen::all) - z(N - 3, Eigen::all) * e(N - 1 - e_offset) - z(N - 2, Eigen::all) * gamma(N - 1)) / mu(N - 1);

	y(N - 1, Eigen::all) = z(N - 1, Eigen::all);
	y(N - 2, Eigen::all) = z(N - 2, Eigen::all) - alpha(N - 2) * y(N - 1, Eigen::all);

	for (int i = N - 3; i >= 0; --i) {
		y(i, Eigen::all) = z(i, Eigen::all) - alpha(i) * y(i + 1, Eigen::all) - beta(i) * y(i + 2, Eigen::all);
	}

    // There are situations where this linear system is underdetermined, this is usually only the case when we are evaluating
    // an azimuthal moment m, that is greater than the largest phase expansion coefficient.  We could check the
    // determinant of the system in the beginning, but it is more computationally efficient just to solve the system,
    // then if the resulting coefficients are nan just set everything to 0 since there will be 0 contribution from this
    // order anyways

    if(y.hasNaN()) {
        BOOST_LOG_TRIVIAL(warning) << "Pentadiagonal solver failed";
        y.setZero();
    }

	return 0;
}