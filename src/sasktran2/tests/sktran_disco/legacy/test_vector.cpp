#include <sasktran2/test_helper.h>

#include <sktran_disco/sktran_do.h>
#include <sasktran2/testing/do_test_util.h>


/** This file contains all of the legacy vector tests that were included in the original sktran_disco project
 *  They have been transformed to use the new low-level API available
 *  They are primarily included to ensure that the low level interface gives identical answers to the old testing interface
 *  and to ensure that no breaking changes have been done on the core DO model when moving it over to the
 *  sasktran2 project
 *
 */


// Epsilon for the standard tests
#define SKDO_FPC_EPS 1e-8

// Errors when comparing to the Coulsen tables, the tables are only accurate to some number of digits
// and full reproduction requires more than 40 streams
#define SKDO_FPC_COULSEN_EPS 1e-5

// Errors when comparing to the aerosol slab of Siewert.  Slightly larger errors are expected
// since these tests are done neglecting circular polarization
#define SKDO_FPC_SIEWERT_EPS 1e-4

// Used for weighting functions mostly where there is more numerical errors
#define SKDO_FPC_EPS_LOW_PRECISION 1e-6



namespace sasktran_disco
{
    namespace testing
    {
        TestAtmosphere<3> default_atmo_vector = std::vector<TestLayerSpecHG>({
                                                                                     // OD,  SSA, ASYM
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.00 }), //< TOA
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.00 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.10 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.10 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.80, 0.30 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.50 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.65, 0.50 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.40, 0.20 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.20 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.90 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.90 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.10 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.10 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.80, 0.30 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.00 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.65, 0.00 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.40, 0.00 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.20 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.30 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.30 }) //< Ground
                                                                             });

        TestSolarSpec default_sun_vector = {
                0.8, //< CSZ
                0.0, //< SAZ
                {
                        1.0, //< TOA incident
                        0  //< TOA isotropic
                }
        };

        std::vector<TestLOS> default_los_vector = {
                { 1.00, 0 * PI / 6 },
                { 1.00, 1 * PI / 6 },
                { 1.00, 2 * PI / 6 },
                { 1.00, 3 * PI / 6 },
                { 1.00, 4 * PI / 6 },
                { 1.00, 5 * PI / 6 },
                { 1.00, 6 * PI / 6 },
                { 0.80, 0 * PI / 6 },
                { 0.80, 1 * PI / 6 },
                { 0.80, 2 * PI / 6 },
                { 0.80, 3 * PI / 6 },
                { 0.80, 4 * PI / 6 },
                { 0.80, 5 * PI / 6 },
                { 0.80, 6 * PI / 6 },
                { 0.60, 0 * PI / 6 },
                { 0.60, 1 * PI / 6 },
                { 0.60, 2 * PI / 6 },
                { 0.60, 3 * PI / 6 },
                { 0.60, 4 * PI / 6 },
                { 0.60, 5 * PI / 6 },
                { 0.60, 6 * PI / 6 },
                { 0.40, 0 * PI / 6 },
                { 0.40, 1 * PI / 6 },
                { 0.40, 2 * PI / 6 },
                { 0.40, 3 * PI / 6 },
                { 0.40, 4 * PI / 6 },
                { 0.40, 5 * PI / 6 },
                { 0.40, 6 * PI / 6 },
                { 0.20, 0 * PI / 6 },
                { 0.20, 1 * PI / 6 },
                { 0.20, 2 * PI / 6 },
                { 0.20, 3 * PI / 6 },
                { 0.20, 4 * PI / 6 },
                { 0.20, 5 * PI / 6 },
                { 0.20, 6 * PI / 6 }
        };
    }
}

/** Verifies that a polarized calculation with the phase matrix set to only have a (1, 1) component gives an
 *  identical answer to the scalar calculation
 */
TEST_CASE("Polariation Same as Scalar", "[sktran_do_legacy][vector]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.125869120756, 0.125869120756, 0.125869120756, 0.125869120756, 0.125869120756,
            0.125869120756, 0.125869120756, 0.124189551170, 0.122723444410, 0.125059565364,
            0.121080168220, 0.123908864401, 0.122680845542, 0.112038263111, 0.124332351431,
            0.126228656342, 0.124016447685, 0.123546311406, 0.118602297428, 0.115821968757,
            0.121538690817, 0.132601462498, 0.128325866427, 0.123930579732, 0.118279365156,
            0.118239052943, 0.119329185336, 0.114711086377, 0.132223881258, 0.132140081997,
            0.127940792580, 0.125599046239, 0.120860572857, 0.116060123620, 0.117648644198,
    };

    // Convert correct radiances to Stokes vector
    std::vector<double> correct_stokes(correct_radiances.size() * 3);
    for (int i = 0; i < correct_radiances.size(); ++i) {
        correct_stokes[i * 3] = correct_radiances[i];
        correct_stokes[i * 3 + 1] = 0.0;
        correct_stokes[i * 3 + 2] = 0.0;
        //correct_stokes[i * 4 + 3] = 0.0;
    }

    TestAtmosphere<3> atmo_pol = std::vector<TestLayerSpecHG>({
                                                                      // OD,  SSA, ASYM
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.00 }), //< TOA
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.00 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.10 }),
                                                                      TestLayerSpecHG({ 0.04, 0.90, 0.10 }),
                                                                      TestLayerSpecHG({ 0.04, 0.80, 0.30 }),
                                                                      TestLayerSpecHG({ 0.04, 0.90, 0.50 }),
                                                                      TestLayerSpecHG({ 0.04, 0.65, 0.50 }),
                                                                      TestLayerSpecHG({ 0.04, 0.40, 0.20 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.20 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.90 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.90 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.10 }),
                                                                      TestLayerSpecHG({ 0.04, 0.90, 0.10 }),
                                                                      TestLayerSpecHG({ 0.04, 0.80, 0.30 }),
                                                                      TestLayerSpecHG({ 0.04, 0.90, 0.00 }),
                                                                      TestLayerSpecHG({ 0.04, 0.65, 0.00 }),
                                                                      TestLayerSpecHG({ 0.04, 0.40, 0.00 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.20 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.30 }),
                                                                      TestLayerSpecHG({ 0.04, 0.95, 0.30 }) //< Ground
                                                              });

    // Run test
    TestCase<3> testcase(16, default_sun_vector, atmo_pol, 0.7, default_los_vector, correct_stokes);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<3, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}


/** Tests against results given by Siewert
 *
 */
TEST_CASE("Polarization Siewert Aeresol", "[sktran_do_legacy][vector]") {
    using namespace sasktran_disco;
    using namespace testing;


    TestAtmosphere<3> single_layer_atmo = std::vector<TestLayerSpecSiewert>({
                                                                                    TestLayerSpecSiewert({ 1.0, 0.973527 })
                                                                            });

    TestSolarSpec new_sun = {
            0.6, //< CSZ
            0, //< SAZ
            {
                    PI, //< TOA incident
                    0  //< TOA isotropic
            }
    };

    std::vector<TestLOS> new_los = {
            {1.0, 0 * PI / 180 },
            {0.5, 0 * PI / 180 },
            {0.2, 0 * PI / 180 },
            {1.0, 180 * PI / 180 },
            {0.5, 180 * PI / 180 },
            {0.2, 180 * PI / 180 },
            {1.0, 90 * PI / 180 },
            {0.5, 90 * PI / 180 },
            {0.2, 90 * PI / 180 },
    };

    std::vector<double> correct_radiances = {
            0.0506873, -0.00262388, 0,
            0.339136, -0.0282242, 0,
            0.751295, -0.0638561, 0,
            0.0506873, -0.00262388, 0,
            0.0684106, 0.00196215, 0,
            0.0801523, 0.00243740, 0,
            0.0506873, 0.00262388, 0,
            0.124626, 0.00512123, -0.00804140,
            0.169216, 0.00696260, -0.00912219
    };

    // Run test
    TestCase<3> testcase(40, new_sun, single_layer_atmo, 0, new_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<3, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_COULSEN_EPS);
    }
}

/** Validates the Coulsen tables for albedo=0
 *
 */
TEST_CASE("Polarization Coulsen Albedo 0 STOKES=3", "[sktran_do_legacy][vector]") {
    using namespace sasktran_disco;
    using namespace testing;

    // tau=0.5, albedo=0, csz=0.2 table

    TestAtmosphere<3> single_layer_atmo = std::vector<TestLayerSpecRayleigh>({
                                                                                     TestLayerSpecRayleigh({ 0.5, 1.0, 0 })
                                                                             });

    TestSolarSpec new_sun = {
            0.2, //< CSZ
            0, //< SAZ
            {
                    PI, //< TOA incident
                    0  //< TOA isotropic
            }
    };

    std::vector<TestLOS> new_los = {
            {0.02, 0 * PI / 180 },
            {0.4, 0 * PI / 180 },
            {1.00, 0 * PI / 180 },
            {0.02, 60 * PI / 180 },
            {0.4, 60 * PI / 180 },
            {1.00, 60 * PI / 180 },

    };

    std::vector<double> correct_radiances = {
            0.44129802, -0.01753141, 0,
            0.16889020, 0.01119511, 0,
            0.05300496, 0.03755859, 0,
            0.30091208, -0.15965601, 0.07365528,
            0.12752450, -0.06066038, 0.05293867,
            0.05300496, -0.01877930, 0.03252669
    };

    // Run test
    TestCase<3> testcase(40, new_sun, single_layer_atmo, 0, new_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<3, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_COULSEN_EPS);
    }
}

/** Checks the coulsen tables for an albedo value of 0.8
 *
 */
TEST_CASE("Coulsen Albedo 0.8", "[sktran_do_legacy][vector]") {
    using namespace sasktran_disco;
    using namespace testing;

    // tau=0.5, albedo=0.8, csz=0.2 table

    TestAtmosphere<3> single_layer_atmo = std::vector<TestLayerSpecRayleigh>({
                                                                                     TestLayerSpecRayleigh({ 0.5, 1.0, 0 })
                                                                             });

    TestSolarSpec new_sun = {
            0.2, //< CSZ
            0, //< SAZ
            {
                    PI, //< TOA incident
                    0  //< TOA isotropic
            }
    };

    std::vector<TestLOS> new_los = {
            {0.02, 0 * PI / 180 },
            {0.4, 0 * PI / 180 },
            {1.00, 0 * PI / 180 },
            {0.02, 60 * PI / 180 },
            {0.4, 60 * PI / 180 },
            {1.00, 60 * PI / 180 },

    };

    std::vector<double> correct_radiances = {
            0.47382125, -0.01553672, 0,
            0.23059806, 0.01144320, 0,
            0.13280858, 0.03755859, 0,
            0.33343531, -0.15766132, 0.07365528,
            0.18923236, -0.06041229, 0.05293867,
            0.13280858, -0.01877930, 0.03252669,
    };

    // Run test
    TestCase<3> testcase(40, new_sun, single_layer_atmo, 0.8, new_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<3, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        // Errors are slightly larger for this test, but it is just because we are limited to 40 streams (I think)
        REQUIRE(diff < 1.2*SKDO_FPC_COULSEN_EPS);
    }
}


/** Test the boundary conditions by verifying that 3 homogeneous layers give an identical answer to one layer with
 *  3x the optical depth
 */
TEST_CASE("Polarization Boundary Conditions", "[sktran_do_legacy][vector]") {
    using namespace sasktran_disco;
    using namespace testing;


    TestAtmosphere<3> single_layer_atmo = std::vector<TestLayerSpecRayleigh>({
                                                                                     TestLayerSpecRayleigh({ 0.6, 1.0, 0 })
                                                                             });

    TestAtmosphere<3> three_layer_atmo = std::vector<TestLayerSpecRayleigh>({
                                                                                    TestLayerSpecRayleigh({ 0.2, 1.0, 0 }),
                                                                                    TestLayerSpecRayleigh({ 0.2, 1.0, 0 }),
                                                                                    TestLayerSpecRayleigh({ 0.2, 1.0, 0 })
                                                                            });


    std::vector<double> correct_radiances, abs_diff, radiances;
    correct_radiances.resize(default_los_vector.size()*3, 0);

    // Run calculation with single layer atmo, store results in correct_radiances
    TestCase<3> testcase(16, default_sun_vector, single_layer_atmo, 0.8, default_los_vector, correct_radiances);
    SKTRAN_DO_TestSpec<3, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, correct_radiances, abs_diff, nullptr);


    // Run the three layer atmo
    TestCase<3> testcase3(16, default_sun_vector, three_layer_atmo, 0.8, default_los_vector, correct_radiances);
    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase3, spec, radiances, abs_diff, nullptr);

    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}