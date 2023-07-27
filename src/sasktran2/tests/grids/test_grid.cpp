#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("Grid Interpolation", "[sasktran2][grids]") {
    sasktran2::grids::outofbounds o = GENERATE(sasktran2::grids::outofbounds::extend, sasktran2::grids::outofbounds::setzero);
    sasktran2::grids::gridspacing s = GENERATE(sasktran2::grids::gridspacing::constant, sasktran2::grids::gridspacing::variable);

    // Construct a constant spacing grid

    double dx = GENERATE((take(10, random(374, 1245))));
    int nx = GENERATE((take(10, random(2, 100))));
    double x0 = GENERATE((take(10, random(-5000, 5000))));

    Eigen::VectorXd grid_values(nx);
    for(int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i*dx;
    }

    Eigen::VectorXd copy_grid = grid_values;

    sasktran2::grids::Grid grid = sasktran2::grids::Grid(std::move(copy_grid), s,
                                                         o,
                                                         sasktran2::grids::interpolation::linear);

    std::array<int, 2> index;
    std::array<double, 2> weight;
    int num_contrib;
    for(int i = 0; i < 10; ++i) {

        double x = grid_values(0) - 100 + (grid_values(Eigen::last) - grid_values(0)) / 7 * double(i);

        grid.calculate_interpolation_weights(x, index, weight, num_contrib);
        double interp_to = 0;
        for(int j = 0; j < num_contrib; ++j) {
            interp_to += grid_values(index[j]) * weight[j];
        }

        if(x < grid_values(0)) {
            if(o == sasktran2::grids::outofbounds::extend) {
                x = grid_values(0);
            } else {
                x = 0.0;
            }
        }

        if(x > grid_values(Eigen::last)) {
            if(o == sasktran2::grids::outofbounds::extend) {
                x = grid_values(Eigen::last);
            } else {
                x = 0.0;
            }
        }

        REQUIRE(fabs(interp_to - x) < 1e-8);
    }
}