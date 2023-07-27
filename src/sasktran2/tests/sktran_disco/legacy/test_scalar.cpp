#include <sasktran2/test_helper.h>

#include <sktran_disco/sktran_do.h>
#include <sasktran2/testing/do_test_util.h>


/** This file contains all of the legacy scalar tests that were included in the original sktran_disco project
 *  They have been transformed to use the new low-level API available
 *  They are primarily included to ensure that the low level interface gives identical answers to the old testing interface
 *  and to ensure that no breaking changes have been done on the core DO model when moving it over to the
 *  sasktran2 project
 *
 *  Many of these test cases have been directly verified against the DISORT library
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
        TestAtmosphere<1> default_atmo = std::vector<TestLayerSpecHG>({
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

        TestSolarSpec default_sun = {
                0.8, //< CSZ
                0.0, //< SAZ
                {
                        1.0, //< TOA incident
                        0  //< TOA isotropic
                }
        };

        std::vector<TestLOS> default_los = {
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

/**  Performs a simple calculation with the model using all of the default test settings and ensures the radiances
 *   match
 */
TEST_CASE("Simple test cast", "[sktran_do_legacy][scalar]") {
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

    // Run test
    TestCase<1> testcase(16, default_sun, default_atmo, 0.7, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);

    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}


/**  A special test case whree the SSA is 1 for each layer.  SSA=1 is a special case where singularities occur
 *   in the eigenmatrix solution, and so the SSA has to be slightly dithered lower
 */
TEST_CASE("SSA = 1", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    // Create atmospheric description
    TestAtmosphere<1> ssa1_atmo = std::vector<TestLayerSpecHG>({
                                                                       // OD,  SSA, ASYM
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 }), //< TOA
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.10 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.15 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.20 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.10 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 1, 0.05 })  //< Ground
                                                               });

    std::vector<double> correct_radiances = {
            0.181688961718, 0.181688961718, 0.181688961718, 0.181688961718, 0.181688961718,
            0.181688961718, 0.181688961718, 0.187244455811, 0.186918094867, 0.186054310234,
            0.184935299235, 0.183879843309, 0.183143775974, 0.182881580427, 0.192432819007,
            0.191887519194, 0.190462005248, 0.188650975585, 0.186976727949, 0.185827245307,
            0.185421220374, 0.197391757358, 0.196631781478, 0.194656029074, 0.192167032735,
            0.189885282375, 0.188328794543, 0.187780908105, 0.198608919153, 0.197694477294,
            0.195303303008, 0.192265644351, 0.189460046019, 0.187536735656, 0.186858152629,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, ssa1_atmo, 0.7, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);

    // ** NOTE **
    // ==========
    // SSA = 1 is a special case which is not handled well by the discrete ordinate
    // method. A SSA of 1 causes the eigenmatrix (in the homogeneous solution) to be
    // singular. To handle this SKDO dither's the SSA of any layer with a SSA greater
    // than SSA - userspec.dither by -userspec.dither. While this does prevent the
    // eigenmatrix from being singular, the matrix is poorly conditioned and because
    // of this it is necessary for this test case to have a relaxed tolerance.
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS_LOW_PRECISION);
    }
}

/** Special case where all of the layers are set to have SSA=0, this is essentially a purely absorbing atmosphere,
 *  and we want to make sure no weird things happen
 */
TEST_CASE("SSA = 0", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    // Create atmospheric description
    TestAtmosphere<1> ssa0_atmo = std::vector<TestLayerSpecHG>({
                                                                       // OD,  SSA, ASYM
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 }), //< TOA
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.10 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.15 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.20 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.10 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 }),
                                                                       TestLayerSpecHG({ 0.1, 0, 0.05 })  //< Ground
                                                               });

    std::vector<double> correct_radiances = {
            0.015002350636, 0.015002350636, 0.015002350636, 0.015002350636, 0.015002350636,
            0.015002350636, 0.015002350636, 0.011395367326, 0.011395367326, 0.011395367326,
            0.011395367326, 0.011395367326, 0.011395367326, 0.011395367326, 0.007205708539,
            0.007205708539, 0.007205708539, 0.007205708539, 0.007205708539, 0.007205708539,
            0.007205708539, 0.002881200069, 0.002881200069, 0.002881200069, 0.002881200069,
            0.002881200069, 0.002881200069, 0.002881200069, 0.000184188958, 0.000184188958,
            0.000184188958, 0.000184188958, 0.000184188958, 0.000184188958, 0.000184188958,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, ssa0_atmo, 0.7, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Repeat the same SSA=0, and SSA=1 test with SSA=0.3 for good measure
 *
 */
TEST_CASE("SSA = 0.3", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    // Create atmospheric description
    TestAtmosphere<1> ssa03_atmo = std::vector<TestLayerSpecHG>({
                                                                        // OD,  SSA, ASYM
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 }), //< TOA
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.10 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.15 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.20 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.10 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 }),
                                                                        TestLayerSpecHG({ 0.1, 0.3, 0.05 })  //< Ground
                                                                });

    std::vector<double> correct_radiances = {
            0.031401571879, 0.031401571879, 0.031401571879, 0.031401571879, 0.031401571879,
            0.031401571879, 0.031401571879, 0.029521480102, 0.029427984519, 0.029180833492,
            0.028861335779, 0.028560718933, 0.028351511645, 0.028277079886, 0.027210085574,
            0.027052937148, 0.026642746638, 0.026122983209, 0.025643909180, 0.025315837177,
            0.025200125537, 0.025346909773, 0.025126406136, 0.024553980547, 0.023834643890,
            0.023177073061, 0.022729610211, 0.022572324738, 0.025906164424, 0.025638204916,
            0.024938176013, 0.024050339377, 0.023231894006, 0.022671770295, 0.022474341285,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, ssa03_atmo, 0.7, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Another similar test using an opticall thick atmosphere, where vertical OD > ~1
 */
TEST_CASE("Optically thick atmosphere", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    // Create atmospheric description
    TestAtmosphere<1> thick_atmo = std::vector<TestLayerSpecHG>({
                                                                        // OD,  SSA, ASYM
                                                                        TestLayerSpecHG({ 0.2, 0.94, 0.05 }), //< TOA
                                                                        TestLayerSpecHG({ 0.2, 0.98, 0.05 }),
                                                                        TestLayerSpecHG({ 0.2, 0.93, 0.05 }),
                                                                        TestLayerSpecHG({ 0.2, 0.98, 0.05 }),
                                                                        TestLayerSpecHG({ 0.2, 0.98, 0.05 }),
                                                                        TestLayerSpecHG({ 0.2, 0.97, 0.10 }),
                                                                        TestLayerSpecHG({ 0.2, 0.98, 0.15 }),
                                                                        TestLayerSpecHG({ 0.2, 0.93, 0.20 }),
                                                                        TestLayerSpecHG({ 0.2, 0.98, 0.10 }),
                                                                        TestLayerSpecHG({ 0.2, 0.95, 0.05 }),
                                                                        TestLayerSpecHG({ 0.2, 0.98, 0.05 })  //< Ground
                                                                });

    std::vector<double> correct_radiances = {
            0.151994476090, 0.151994476090, 0.151994476090, 0.151994476090, 0.151994476090,
            0.151994476090, 0.151994476090, 0.158475420729, 0.158187542155, 0.157420557762,
            0.156416394039, 0.155458816672, 0.154785265333, 0.154544234515, 0.164121885345,
            0.163662845115, 0.162449946254, 0.160883477883, 0.159411564142, 0.158388643878,
            0.158025034047, 0.168839259644, 0.168216169585, 0.166574988933, 0.164466875363,
            0.162498310375, 0.161137526364, 0.160655294880, 0.170270208199, 0.169466512864,
            0.167350612239, 0.164636264321, 0.162106876959, 0.160362177787, 0.159744727075,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, thick_atmo, 0.7, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Same test with a very optically thin atmosphere
 */
TEST_CASE("Optically thin atmosphere", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    // Create atmospheric description
    TestAtmosphere<1> thin_atmo = std::vector<TestLayerSpecHG>({
                                                                       // OD,  SSA, ASYM
                                                                       TestLayerSpecHG({ 0.01, 0.94, 0.05 }), //< TOA
                                                                       TestLayerSpecHG({ 0.01, 0.98, 0.05 }),
                                                                       TestLayerSpecHG({ 0.01, 0.93, 0.05 }),
                                                                       TestLayerSpecHG({ 0.01, 0.98, 0.05 }),
                                                                       TestLayerSpecHG({ 0.01, 0.98, 0.05 }),
                                                                       TestLayerSpecHG({ 0.01, 0.97, 0.10 }),
                                                                       TestLayerSpecHG({ 0.01, 0.98, 0.15 }),
                                                                       TestLayerSpecHG({ 0.01, 0.93, 0.01 }),
                                                                       TestLayerSpecHG({ 0.01, 0.98, 0.10 }),
                                                                       TestLayerSpecHG({ 0.01, 0.95, 0.05 }),
                                                                       TestLayerSpecHG({ 0.01, 0.98, 0.05 }) //< Ground
                                                               });

    std::vector<double> correct_radiances = {
            0.176411655415, 0.176411655415, 0.176411655415, 0.176411655415, 0.176411655415,
            0.176411655415, 0.176411655415, 0.176737577342, 0.176662283650, 0.176463239508,
            0.176205990856, 0.175964101861, 0.175795900056, 0.175736089482, 0.176885336248,
            0.176740029763, 0.176361240442, 0.175882590290, 0.175443039848, 0.175143062688,
            0.175037472951, 0.177060426834, 0.176798570869, 0.176122840292, 0.175282428848,
            0.174523037422, 0.174011184550, 0.173832205024, 0.177291898263, 0.176758712378,
            0.175392321280, 0.173710858098, 0.172207283931, 0.171201724101, 0.170851545894,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, thin_atmo, 0.7, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}


/** The default test atmosphere where the surface albedo is adjusted to 1
 */
TEST_CASE("Surface albedo = 1", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.176283725302, 0.176283725302, 0.176283725302, 0.176283725302, 0.176283725302,
            0.176283725302, 0.176283725302, 0.170193142633, 0.168727035872, 0.171063156826,
            0.167083759683, 0.169912455864, 0.168684437005, 0.158041854574, 0.164353057773,
            0.166249362685, 0.164037154027, 0.163567017748, 0.158623003770, 0.155842675099,
            0.161559397159, 0.164518788004, 0.160243191933, 0.155847905238, 0.150196690662,
            0.150156378450, 0.151246510842, 0.146628411883, 0.154749564637, 0.154665765376,
            0.150466475960, 0.148124729619, 0.143386256237, 0.138585807000, 0.140174327578,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, default_atmo, 1.0, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** The default test case where the surface albedo is set to 0
 */
TEST_CASE("Surface albedo = 0", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.034973070891, 0.034973070891, 0.034973070891, 0.034973070891, 0.034973070891,
            0.034973070891, 0.034973070891, 0.041246428216, 0.039780321455, 0.042116442409,
            0.038137045266, 0.040965741446, 0.039737722588, 0.029095140156, 0.052176194459,
            0.054072499370, 0.051860290713, 0.051390154434, 0.046446140456, 0.043665811785,
            0.049382533845, 0.075055462943, 0.070779866872, 0.066384580176, 0.060733365601,
            0.060693053388, 0.061783185780, 0.057165086822, 0.091610736353, 0.091526937092,
            0.087327647676, 0.084985901335, 0.080247427952, 0.075446978716, 0.077035499293,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, default_atmo, 0.0, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Default test case where the surface albedo is adjusted to 0.3
 */
TEST_CASE("Surface albedo = 0.3", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.070127529954, 0.070127529954, 0.070127529954, 0.070127529954, 0.070127529954,
            0.070127529954, 0.070127529954, 0.073325056768, 0.071858950007, 0.074195070961,
            0.070215673818, 0.073044369998, 0.071816351140, 0.061173768708, 0.080082915050,
            0.081979219961, 0.079767011304, 0.079296875025, 0.074352861047, 0.071572532376,
            0.077289254436, 0.097311638966, 0.093036042895, 0.088640756199, 0.082989541624,
            0.082949229411, 0.084039361804, 0.079421262845, 0.107318054131, 0.107234254870,
            0.103034965453, 0.100693219112, 0.095954745730, 0.091154296493, 0.092742817071,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, default_atmo, 0.3, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/**  Default test case modified with the solar zenith angle being set to 0
 *
 */
TEST_CASE("Sun directly overhead", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestSolarSpec overhead_sun = { 1.0, 0,{ 1, 0 } };

    std::vector<double> correct_radiances = {
            0.171860861452, 0.171860861452, 0.171860861452, 0.171860861452, 0.171860861452,
            0.171860861452, 0.171860861452, 0.177175320855, 0.177175320855, 0.177175320855,
            0.177175320855, 0.177175320855, 0.177175320855, 0.177175320855, 0.167067611716,
            0.167067611716, 0.167067611716, 0.167067611716, 0.167067611716, 0.167067611716,
            0.167067611716, 0.156509308960, 0.156509308960, 0.156509308960, 0.156509308960,
            0.156509308960, 0.156509308960, 0.156509308960, 0.151990758394, 0.151990758394,
            0.151990758394, 0.151990758394, 0.151990758394, 0.151990758394, 0.151990758394,
    };

    // Run test
    TestCase<1> testcase(16, overhead_sun, default_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Default test case where the solar zenith angle is set to 65 degrees
 */
TEST_CASE("Solar zenith = 65[deg]", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestSolarSpec overhead_sun = { 0.4226, 0,{ 1, 0 } };

    std::vector<double> correct_radiances = {
            0.067065438230, 0.067065438230, 0.067065438230, 0.067065438230, 0.067065438230,
            0.067065438230, 0.067065438230, 0.075051929402, 0.072449607424, 0.071394554815,
            0.067908269122, 0.067587585853, 0.068647162079, 0.065627974455, 0.079212665688,
            0.077651926458, 0.074734216691, 0.072302753570, 0.069662928277, 0.067651484537,
            0.070850942658, 0.093751067467, 0.091795206438, 0.085322130943, 0.079717990027,
            0.075130246894, 0.071114016385, 0.065632521407, 0.114056393891, 0.105046070881,
            0.095557563707, 0.090298867916, 0.086679023508, 0.083387178498, 0.084498954667,
    };

    // Run test
    TestCase<1> testcase(16, overhead_sun, default_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Test that a single layer atmosphere works correctly
 */
TEST_CASE("Single layer atmosphere", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestAtmosphere<1> single_layer_atmo = std::vector<TestLayerSpecHG>({
                                                                               TestLayerSpecHG({ 0.4, 0.9, 0.03 })
                                                                       });

    std::vector<double> correct_radiances = {
            0.181529186226, 0.181529186226, 0.181529186226, 0.181529186226, 0.181529186226,
            0.181529186226, 0.181529186226, 0.180693155123, 0.180598287723, 0.180342458990,
            0.180000733746, 0.179667657448, 0.179429133958, 0.179342922859, 0.178904078202,
            0.178740210067, 0.178300338702, 0.177717374496, 0.177154186634, 0.176753898424,
            0.176609837113, 0.175806122711, 0.175550632698, 0.174866762944, 0.173964781581,
            0.173098105751, 0.172484919811, 0.172264808127, 0.169811775850, 0.169413087976,
            0.168347711792, 0.166946542512, 0.165604509672, 0.164657548316, 0.164318138934,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, single_layer_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}


/**  A single layer atmosphere where the albedo is set to 0
 */
TEST_CASE("Single layer atmosphere Zero Albedo", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestAtmosphere<1> single_layer_atmo = std::vector<TestLayerSpecHG>({
                                                                               TestLayerSpecHG({ 0.4, 0.9, 0.03 })
                                                                       });

    std::vector<double> correct_radiances = {
            0.027787798461, 0.027787798461, 0.027787798461, 0.027787798461, 0.027787798461,
            0.027787798461, 0.027787798461, 0.034252283230, 0.034157415830, 0.033901587097,
            0.033559861854, 0.033226785555, 0.032988262065, 0.032902050966, 0.043261433948,
            0.043097565813, 0.042657694448, 0.042074730241, 0.041511542380, 0.041111254170,
            0.040967192859, 0.057493940453, 0.057238450441, 0.056554580687, 0.055652599324,
            0.054785923493, 0.054172737554, 0.053952625870, 0.081487682882, 0.081088995008,
            0.080023618824, 0.078622449544, 0.077280416704, 0.076333455348, 0.075994045966,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, single_layer_atmo, 0.0, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** A two layer atmosphere
 */
TEST_CASE("Two layer atmosphere", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestAtmosphere<1> two_layer_atmo = std::vector<TestLayerSpecHG>({
                                                                            TestLayerSpecHG({ 0.4, 0.9, 0.03 }),
                                                                            TestLayerSpecHG({ 0.4, 0.9, 0.03 })
                                                                    });

    std::vector<double> correct_radiances = {
            0.159432472801, 0.159432472801, 0.159432472801, 0.159432472801, 0.159432472801,
            0.159432472801, 0.159432472801, 0.159610506156, 0.159480210417, 0.159128824023,
            0.158659414511, 0.158201837237, 0.157874125377, 0.157755671937, 0.159163804594,
            0.158948104179, 0.158369062160, 0.157601567603, 0.156860010894, 0.156332884521,
            0.156143161924, 0.158327447863, 0.158013926266, 0.157174668418, 0.156067612127,
            0.155003743590, 0.154250950234, 0.153980705431, 0.156620878259, 0.156188433494,
            0.155032792011, 0.153512768027, 0.152056739459, 0.151029242406, 0.150660946768,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, two_layer_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Test where the line of sight angles are modified to be similar to the eigenvalues of the homogeneous solution,
 *  which can lead to a removable singularity in the post processing
 */
TEST_CASE("Intensity singularity", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    // These cos zenith's are such that equal to the reciprocal of eigenvalues found in this solution
    std::vector<TestLOS> conflict_los = {
            { 0.95887396, 0.0 },
            { 0.83467543, 0.0 },
            { 0.65342150, 0.0 },
            { 0.44901703, 0.0 },
            { 0.98274959, 0.0 },
            { 0.91453839, 0.0 },
            { 0.78721704, 0.0 },
            { 0.42366863, 0.0 },
    };

    std::vector<double> correct_radiances = {
            0.140132512805, 0.137803416410, 0.139833537001, 0.140869462076, 0.138027329751,
            0.141855473997, 0.139310941761, 0.142037786346,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, default_atmo, 0.8, conflict_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/**  Test case where the henyey-greenstein asymmetry factor is reduced to 0, i.e., isotropic scattering
 *
 */
TEST_CASE("HG asymmetry factor = 0", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestAtmosphere<1> diffuse_scattering_atmo = std::vector<TestLayerSpecHG>({
                                                                                     // OD,  SSA, ASYM
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }), //< TOA
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.80, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.65, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.40, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.80, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.65, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.40, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0 }) //< Ground
                                                                             });

    std::vector<double> correct_radiances = {
            0.139924516596, 0.139924516596, 0.139924516596, 0.139924516596, 0.139924516596,
            0.139924516596, 0.139924516596, 0.138719007613, 0.138719007613, 0.138719007613,
            0.138719007613, 0.138719007613, 0.138719007613, 0.138719007613, 0.137334516208,
            0.137334516208, 0.137334516208, 0.137334516208, 0.137334516208, 0.137334516208,
            0.137334516208, 0.136139134168, 0.136139134168, 0.136139134168, 0.136139134168,
            0.136139134168, 0.136139134168, 0.136139134168, 0.137466964210, 0.137466964210,
            0.137466964210, 0.137466964210, 0.137466964210, 0.137466964210, 0.137466964210,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, diffuse_scattering_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Test where the Henyey-Greenstein asymmetry factor is modified to 0.9, for more azimuthal dependence
 *
 */
TEST_CASE("HG asymmetry factor = 0.9", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    TestAtmosphere<1> diffuse_scattering_atmo = std::vector<TestLayerSpecHG>({
                                                                                     // OD,  SSA, ASYM
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }), //< TOA
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.80, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.65, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.40, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.80, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.90, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.65, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.40, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }),
                                                                                     TestLayerSpecHG({ 0.04, 0.95, 0.9 }) //< Ground
                                                                             });

    std::vector<double> correct_radiances = {
            0.171566153476, 0.171566153476, 0.171566153476, 0.171566153476, 0.171566153476,
            0.171566153476, 0.171566153476, 0.139241383538, 0.121204090312, 0.156923166781,
            0.126977726494, 0.168112879697, 0.162432569524, 0.044462417067, 0.105783168425,
            0.141843215238, 0.140306010310, 0.157397172862, 0.119223540927, 0.099363637721,
            0.175912150847, 0.180284893222, 0.136838671610, 0.112935795453, 0.081395726254,
            0.116241531918, 0.152254640461, 0.092470416756, 0.027003563384, 0.075594619528,
            0.083424024773, 0.131998265714, 0.098034883021, 0.029254262781, 0.087050645842,
    };

    // Run test
    TestCase<1> testcase(16, default_sun, diffuse_scattering_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Default test case where the number of streams is adjusted to 8
 */
TEST_CASE("8 stream", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.140609793186, 0.140609793186, 0.140609793186, 0.140609793186, 0.140609793186,
            0.140609793186, 0.140609793186, 0.143533282579, 0.142935770865, 0.138883286907,
            0.133500882506, 0.138515998749, 0.136029713060, 0.124764207612, 0.135495672295,
            0.136685001100, 0.140120245794, 0.133867604479, 0.129646206112, 0.136085725994,
            0.129924599973, 0.136899125585, 0.133044967009, 0.131984835599, 0.135266433685,
            0.124129968684, 0.127856226519, 0.131057628314, 0.147924620913, 0.143338640838,
            0.131345750049, 0.132039413651, 0.127932320988, 0.121524775729, 0.122879333394,
    };

    // Run test
    TestCase<1> testcase(8, default_sun, default_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Default test case where the number of streams in adjusted to 32
 */
TEST_CASE("32 stream", "[sktran_do][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.139370233993, 0.139370233993, 0.139370233993, 0.139370233993, 0.139370233993,
            0.139370233993, 0.139370233993, 0.140316425838, 0.139132890709, 0.138949721679,
            0.137555533604, 0.135596186272, 0.136359045313, 0.132079931167, 0.138572530113,
            0.138966170761, 0.136772602886, 0.134512621487, 0.132293507640, 0.131182349726,
            0.131357720470, 0.138251929406, 0.137988231716, 0.134775460723, 0.130749024634,
            0.128644604929, 0.126251726292, 0.126818147767, 0.141756560178, 0.140297260464,
            0.135883223843, 0.130137701591, 0.127667624017, 0.124781537978, 0.124912253141,
    };

    // Run test
    TestCase<1> testcase(32, default_sun, default_atmo, 0.8, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}

/** Test where a functional form of the BRDF is used instead of assuming a lambertian surface
 */
TEST_CASE("BRDF surface", "[sktran_do_legacy][scalar]") {
    using namespace sasktran_disco;
    using namespace testing;

    std::vector<double> correct_radiances = {
            0.103306744162, 0.098348410332, 0.084977805696, 0.068119608493, 0.084977805696,
            0.098348410332, 0.103306744162, 0.117514588892, 0.110815099895, 0.095837312068,
            0.076267074863, 0.092483776660, 0.108041254499, 0.109278094220, 0.110646237926,
            0.107663725734, 0.096468731329, 0.083350034202, 0.091989636084, 0.099879904700,
            0.103431428284, 0.109825639901, 0.108088940054, 0.100999841084, 0.092476846064,
            0.094868985290, 0.096352434630, 0.098391858262, 0.120400918764, 0.118674561085,
            0.113604297251, 0.107247445830, 0.105388697424, 0.103158838599, 0.103556611726,
    };

    auto brdf = [](double out, double in, double caz_diff) {
        // Note: this is a garbage BRDF function. Because BRDF(1,any,any) will
        // be different => correct_radiance(zenith=1,azdiff=any) will be different
        return abs(cos(caz_diff)) * abs(1 - abs(abs(out) - abs(in)));
    };

    // Run test
    TestCase<1> testcase(32, default_sun, default_atmo, brdf, default_los, correct_radiances);
    std::vector<double> abs_diff, radiance;
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    // BRDF NOT CURRENTLY IMPLEMTENTED
    /*
    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, radiance, abs_diff, nullptr);
    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
     */
}

/** Test the boundary conditions by verifying that 3 homogeneous layers give an identical answer to one layer with
 *  3x the optical depth
 */
TEST_CASE("Scalar Boundary Conditions", "[sktran_do_legacy][vector]") {
    using namespace sasktran_disco;
    using namespace testing;


    TestAtmosphere<1> single_layer_atmo = std::vector<TestLayerSpecRayleigh>({
                                                                                     TestLayerSpecRayleigh({ 0.6, 1.0, 0 })
                                                                             });

    TestAtmosphere<1> three_layer_atmo = std::vector<TestLayerSpecRayleigh>({
                                                                                    TestLayerSpecRayleigh({ 0.2, 1.0, 0 }),
                                                                                    TestLayerSpecRayleigh({ 0.2, 1.0, 0 }),
                                                                                    TestLayerSpecRayleigh({ 0.2, 1.0, 0 })
                                                                            });


    std::vector<double> correct_radiances, abs_diff, radiances;
    correct_radiances.resize(default_los.size(), 0);

    // Run calculation with single layer atmo, store results in correct_radiances
    TestCase<1> testcase(16, default_sun, single_layer_atmo, 0.8, default_los, correct_radiances);
    SKTRAN_DO_TestSpec<1, -1>* spec = nullptr;

    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase, spec, correct_radiances, abs_diff, nullptr);


    // Run the three layer atmo
    TestCase<1> testcase3(16, default_sun, three_layer_atmo, 0.8, default_los, correct_radiances);
    sasktran_disco::testing::run_lowlevel_from_old_testspec(testcase3, spec, radiances, abs_diff, nullptr);

    for (double diff : abs_diff) {
        REQUIRE(diff < SKDO_FPC_EPS);
    }
}