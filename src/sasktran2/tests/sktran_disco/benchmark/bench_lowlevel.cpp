#include <sasktran2/test_helper.h>

#include <sktran_disco/sktran_do.h>

#ifdef SKTRAN_CATCH2_VERSION3


TEST_CASE("2StreamBenchmark", "[sktran_do][lowlevel][benchmark]") {
    // TODO: 2 stream single layer sometimes fails with derivatives and no full compile? Check this

    int nwavel = 1000;
    int nlyr = 20;
    int nderiv = 4;

    sasktran_disco_lowlevel::CPPApi cppapi(2, nlyr, nwavel, 1, 1, nderiv * nlyr);

    for (int i = 0; i < nwavel; ++i) {
        for (int j = 0; j < nlyr; ++j) {
            cppapi.od(j, i) = 0.2;
            cppapi.ssa(j, i) = 0.8;
            cppapi.a1(0, j, i) = 1;
            //cppapi.a1(2, j, i) = 0.5;

            for (int k = 0; k < nderiv; ++k) {
                cppapi.d_layerindex(k*nlyr + j) = j;
                cppapi.d_ssa(k*nlyr + j, i) = 1;
            }
        }
    }
    // Layer Construction
    double top_alt = 0.001;
    for (int j = 0; j < nlyr; ++j) {
        cppapi.layerboundaryaltitude(j) = top_alt - j*(top_alt / double(nlyr));
    }

    cppapi.earth_radius() = 6372000;

    // Coordinates construction
    cppapi.cos_sza() = 0.8;
    cppapi.cos_vza(0) = 0.7;
    cppapi.saa(0) = 0;


    sasktran_disco_lowlevel::Atmosphere atmosphere;
    sasktran_disco_lowlevel::Config config;
    sasktran_disco_lowlevel::WeightingFunctions weightingfunctions;
    sasktran_disco_lowlevel::ViewingGeometry geometry;
    sasktran_disco_lowlevel::Output output;
    cppapi.initialize_c_api(&atmosphere, &config, &geometry, &weightingfunctions, &output);

    config.nthreads = 1;

    BENCHMARK("Test") {
        return sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);
    };
}


TEST_CASE("8 Stream Benchmark", "[sktran_do][lowlevel][benchmark]") {
    //LAPACKE_set_nancheck(0);

    int nwavel = 50;
    int nlyr = 40;
    int nderiv = 0;

    sasktran_disco_lowlevel::CPPApi cppapi(16, nlyr, nwavel, 1, 1, nderiv * nlyr);

    for (int i = 0; i < nwavel; ++i) {
        for (int j = 0; j < nlyr; ++j) {
            cppapi.od(j, i) = 0.2;
            cppapi.ssa(j, i) = 0.8;
            cppapi.a1(0, j, i) = 1;
            cppapi.a1(2, j, i) = 0.5;

            for (int k = 0; k < nderiv; ++k) {
                cppapi.d_layerindex(k*nlyr + j) = j;
                cppapi.d_od(k*nlyr + j, i) = 1;
            }
        }
    }
    // Layer Construction
    double top_alt = 0.001;
    for (int j = 0; j < nlyr; ++j) {
        cppapi.layerboundaryaltitude(j) = top_alt - j*(top_alt / double(nlyr));
    }

    cppapi.earth_radius() = 6372000;

    // Coordinates construction
    cppapi.cos_sza() = 0.8;
    cppapi.cos_vza(0) = 0.7;
    cppapi.saa(0) = 0;


    sasktran_disco_lowlevel::Atmosphere atmosphere;
    sasktran_disco_lowlevel::Config config;
    sasktran_disco_lowlevel::WeightingFunctions weightingfunctions;
    sasktran_disco_lowlevel::ViewingGeometry geometry;
    sasktran_disco_lowlevel::Output output;

    config.nthreads = 1;

    cppapi.initialize_c_api(&atmosphere, &config, &geometry, &weightingfunctions, &output);
    config.numazimuthexpansion = 1;

    BENCHMARK("Test") {
                          return sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);
                      };
}

#endif