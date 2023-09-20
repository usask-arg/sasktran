#ifdef SKTRAN_CATCH2_VERSION3
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_all.hpp>

#else
#include <catch2/catch.hpp>
#endif

#include <sasktran.h>

#include "modules/sasktranv3_impl/enginestubs/skengine_stubs.h"
#include "modules/sasktranv3_impl/climatologystubs/skclimatology_stubs.h"
#include "sources/sasktranif_opticalimpl/skoptprop_stubs.h"
#include "sources/sasktranif_opticalimpl/skbrdf_stubs.h"


TEST_CASE("Simple SKTRAN_ME Stub Test", "[ME]")
{
    ISKEngine_Stub_ME engine;

    ISKClimatology_Stub* msis = new ISKClimatology_Stub_MSIS ( new skClimatology_MSIS90);
    ISKClimatology_Stub* labow = new ISKClimatology_Stub_Base ( new skClimatology_LabowOzoneVMR);

    ISKOpticalProperty_Stub* ray = new ISKOpticalProperty_Stub_Base (new skOpticalProperties_RayleighDryAir );
    ISKOpticalProperty_Stub* dbm = new ISKOpticalProperty_Stub_Base (new skOpticalProperties_O3_DaumontBrionMalicet);

    nxVector obs(367601.31547888496, 1009976.3136400515, -6871601.202127539);
    nxVector look(0.2884568631765662, 0.7925287180643269, 0.5372996083468239);
    int losidx = 0;
    double mjd(54832.5);

    for(int i = 0; i < 10; ++i) {
        engine.AddLineOfSight(mjd, obs, look, &losidx);
    }

    engine.AddSpecies(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis, ray);
    engine.AddSpecies(SKCLIMATOLOGY_O3_CM3, labow, dbm);
    engine.SetAlbedo(0);

    std::vector<double> wavelengths {280, 350, 400};

    engine.SetWavelengths(&wavelengths[0], wavelengths.size());
    engine.SetPolarizationMode(1);

    std::vector<double> sun {0, 0, 1};

    engine.SetPropertyArray("setsun", &sun[0], 3);
    engine.SetPropertyScalar("msmode", 1);

    engine.InitializeModel();

    const double* radiance;
    int numwavel, numlos;
    engine.CalculateRadiance(&radiance, &numwavel, &numlos);



    REQUIRE(true);
}