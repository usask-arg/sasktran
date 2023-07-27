#include <sasktran2/test_helper.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <sasktran2.h>


TEST_CASE("Lebedev Initialization", "[sasktran2][unitsphere]") {
    const std::vector<int> valid_points = {6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890};
    for (const auto& npoints : valid_points) {
        REQUIRE_NOTHROW(sasktran2::math::LebedevSphere(npoints));
    }

    const std::vector<int> invalid_points = {1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20};
    for (const auto& npoints : invalid_points) {
        REQUIRE_THROWS(sasktran2::math::LebedevSphere(npoints));
    }
}


TEST_CASE("Spherical Harmonics Integration", "[sasktran2][unitsphere][spherical_harmonics]") {
    const int max_order = 8;
    const int npoints = 38;
    for (int l = 0; l <= max_order; l++) {
        for (int m = -l; m <= l; m++) {
            sasktran2::math::LebedevSphere sphere(npoints);
            Eigen::Matrix<double, 1, -1> values;
            values.resize(1, npoints);
            for (int i = 0; i < npoints; i++) {
                const auto pos = sphere.get_quad_position(i);
                double theta = atan2(sqrt(pos.x() * pos.x() + pos.y() * pos.y()), pos.z());
                double phi = atan2(pos.y(), pos.x());
                values(0, i) = boost::math::spherical_harmonic(l, m, theta, phi).real();
            }
            Eigen::Vector<double, 1> result;
            sphere.integrate_on_grid<1>(values, result);

            if(l==0 && m==0) {
                REQUIRE(abs(result(0)- (sqrt(4*EIGEN_PI))) < 1e-8);
            } else {
                REQUIRE(abs(result(0)) < 1e-8);
            }
        }
    }
}