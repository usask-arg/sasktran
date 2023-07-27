#include <sasktran2/test_helper.h>

#include <sasktran2.h>

TEST_CASE("Spherical Shell Raytracer - Observer Outside limb Viewing", "[sasktran2][raytracing]") {
    // Construct a constant spacing grid

    double dx = 1000;
    int nx = 100;
    double x0 = 0;
    double tangent_altitude = GENERATE(10000.0, 20000.0, 30405.0);

    Eigen::VectorXd grid_values(nx);
    for(int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i*dx;
    }

    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(std::move(grid_values), sasktran2::grids::gridspacing::constant,
                                                         sasktran2::grids::outofbounds::extend,
                                                         sasktran2::grids::interpolation::linear);

    sasktran2::Coordinates coord(1, 0, 6372000);

    sasktran2::Geometry1D geo(std::move(coord), std::move(grid));


    sasktran2::raytracing::SphericalShellRayTracer raytracer(geo);


    sasktran2::viewinggeometry::TangentAltitude ray_policy(tangent_altitude, 0, 600000);

    auto viewing_ray = ray_policy.construct_ray(geo.coordinates());

    sasktran2::raytracing::TracedRay traced_ray;
    raytracer.trace_ray(viewing_ray, traced_ray);

    double min_alt = 1e99;
    for(const auto& layer : traced_ray.layers) {
        double entrance_altitude = layer.entrance.radius() - geo.coordinates().earth_radius();
        double exit_altitude =  layer.exit.radius() - geo.coordinates().earth_radius();
        if(entrance_altitude < min_alt) {
            min_alt = entrance_altitude;
        }
    }

    REQUIRE(fabs(tangent_altitude - min_alt) < 1e-8);

}

TEST_CASE("Spherical Shell Raytracer - Observer Inside Limb Viewing", "[sasktran2][raytracing]") {
    // Construct a constant spacing grid

    double dx = 1000;
    int nx = 100;
    double x0 = 0;
    double tangent_altitude = GENERATE(10000.0, 20000.0, 30405.0);
    double dx_factor = GENERATE(0, 0.01, 0.2, 2, 2.2, 5.7);


    Eigen::VectorXd grid_values(nx);
    for(int i = 0; i < nx; ++i) {
        grid_values(i) = x0 + i*dx;
    }

    double entrance_altitude = grid_values(Eigen::last);
    double exit_altitude = grid_values(Eigen::last-1);

    sasktran2::grids::AltitudeGrid grid = sasktran2::grids::AltitudeGrid(std::move(grid_values), sasktran2::grids::gridspacing::constant,
                                                                         sasktran2::grids::outofbounds::extend,
                                                                         sasktran2::grids::interpolation::linear);

    sasktran2::Coordinates coord(1, 0, 6372000);

    sasktran2::Geometry1D geo(std::move(coord), std::move(grid));


    sasktran2::raytracing::SphericalShellRayTracer raytracer(geo);


    sasktran2::viewinggeometry::TangentAltitude ray_policy(tangent_altitude, 0, tangent_altitude + dx*dx_factor);

    auto viewing_ray = ray_policy.construct_ray(geo.coordinates());

    sasktran2::raytracing::TracedRay traced_ray;
    raytracer.trace_ray(viewing_ray, traced_ray);



    double min_alt = 1e99;
    for(const auto& layer : traced_ray.layers) {
        double diff = layer.exit.radius() - geo.coordinates().earth_radius() - entrance_altitude;
        REQUIRE(abs(diff) < 1e-8);
        entrance_altitude = layer.entrance.radius() - geo.coordinates().earth_radius();
        exit_altitude =  layer.exit.radius() - geo.coordinates().earth_radius();
        if(entrance_altitude < min_alt) {
            min_alt = entrance_altitude;
        }
    }

    REQUIRE(fabs(tangent_altitude - min_alt) < 1e-8);

}

