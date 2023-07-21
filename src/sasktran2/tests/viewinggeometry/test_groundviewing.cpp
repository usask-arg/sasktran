#include <sasktran2/test_helper.h>

#include <sasktran2.h>


TEST_CASE("Ground Viewing Spherical", "[sasktran2][viewinggeometry]") {
    double earth_radius = GENERATE((take(10, random(10000, 63720000))));
    double observer_altitude = GENERATE((take(10, random(0, 100000))));

    double cos_sza_reference = 0.5;
    double cos_sza_viewing = 0.8;
    double cos_viewing = 0.8;

    sasktran2::Coordinates coord(cos_sza_reference, 0, earth_radius);

    sasktran2::viewinggeometry::GroundViewingSolar viewing_geo(cos_sza_viewing, 0, cos_viewing, observer_altitude);

    sasktran2::viewinggeometry::ViewingRay ray = viewing_geo.construct_ray(coord);


    REQUIRE(abs(ray.observer.radius() - earth_radius - observer_altitude) < 1e-6);

}