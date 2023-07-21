#include <sasktran2/test_helper.h>

#include <sktran_disco/sktran_do.h>


TEST_CASE("Basic LowLevel", "[sktran_do][lowlevel]") {
	sasktran_disco_lowlevel::CPPApi cppapi(16, 2, 1, 1, 1, 0);

	// Layer Construction
	cppapi.od(0, 0) = 0.2;
	cppapi.od(1, 0) = 0.2;

	cppapi.ssa(0, 0) = 0.8;
	cppapi.ssa(1, 0) = 0.7;

	cppapi.a1(0, 0, 0) = 1;
	cppapi.a1(0, 1, 0) = 1;

	cppapi.a1(2, 0, 0) = 0.5;
	cppapi.a1(2, 1, 0) = 0.5;

	cppapi.layerboundaryaltitude(0) = 100000;
	cppapi.layerboundaryaltitude(1) = 10000;

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


	sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);

}


TEST_CASE("Multiple Layer LowLevel", "[sktran_do][lowlevel]") {
    sasktran_disco_lowlevel::CPPApi cppapi(16, 2, 1, 1, 1, 0);

    // Layer Construction
    cppapi.od(0, 0) = 0.2;
    cppapi.od(1, 0) = 0.2;

    cppapi.ssa(0, 0) = 0.8;
    cppapi.ssa(1, 0) = 0.8;

    cppapi.a1(0, 0, 0) = 1;
    cppapi.a1(0, 1, 0) = 1;

    //cppapi.a1(2, 0, 0) = 0.5;
    //cppapi.a1(2, 1, 0) = 0.5;

    cppapi.layerboundaryaltitude(0) = 0.0001;
    cppapi.layerboundaryaltitude(1) = 0.00001;

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

    sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);

    sasktran_disco_lowlevel::CPPApi cppapi1(16, 1, 1, 1, 1, 0);

    // Layer Construction
    cppapi1.od(0, 0) = 0.4;

    cppapi1.ssa(0, 0) = 0.8;

    cppapi1.a1(0, 0, 0) = 1;

    //cppapi1.a1(2, 0, 0) = 0.5;

    cppapi1.layerboundaryaltitude(0) = 0.0001;

    cppapi1.earth_radius() = 6372000;

    // Coordinates construction
    cppapi1.cos_sza() = 0.8;
    cppapi1.cos_vza(0) = 0.7;
    cppapi1.saa(0) = 0;

    sasktran_disco_lowlevel::Output output2;
    sasktran_disco_lowlevel::Atmosphere atmosphere2;
    sasktran_disco_lowlevel::Config config2;
    cppapi1.initialize_c_api(&atmosphere2, &config2, &geometry, &weightingfunctions, &output2);

    sasktran_disco_lowlevel::calculate(&atmosphere2, &config2, &weightingfunctions, &geometry, &output2);

    int x = 1;

}


TEST_CASE("Profile LowLevel", "[sktran_do][lowlevel]") {
	int nwavel = 5000;
	int nlyr = 20;
	int nderiv = 10;

	sasktran_disco_lowlevel::CPPApi cppapi(2, nlyr, nwavel, 1, 1, nderiv * nlyr);

	for (int i = 0; i < nwavel; ++i) {
		for (int j = 0; j < nlyr; ++j) {
			cppapi.od(j, i) = 0.2;
			cppapi.ssa(j, i) = 0.8;
			cppapi.a1(0, j, i) = 1;
			//cppapi.a1(2, j, i) = 0.5;

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
	cppapi.initialize_c_api(&atmosphere, &config, &geometry, &weightingfunctions, &output);


	sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);
    int x =1 ;

}


TEST_CASE("Check Wavelength Independence Single Thread", "[sktran_do][lowlevel]") {
    int nwavel = 1000;

    sasktran_disco_lowlevel::CPPApi cppapi(16, 2, nwavel, 1, 1, 0);

    // Layer Construction
    for(int i = 0; i < nwavel; ++i) {
        cppapi.od(0, i) = 0.2;
        cppapi.od(1, i) = 0.2;

        cppapi.ssa(0, i) = 0.8;
        cppapi.ssa(1, i) = 0.8;

        cppapi.a1(0, 0, i) = 1;
        cppapi.a1(0, 1, i) = 1;
    }


    //cppapi.a1(2, 0, 0) = 0.5;
    //cppapi.a1(2, 1, 0) = 0.5;

    cppapi.layerboundaryaltitude(0) = 0.0001;
    cppapi.layerboundaryaltitude(1) = 0.00001;

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

    sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);

    Eigen::Map<Eigen::MatrixXd> radiance(output.radiance, 1, nwavel);

    REQUIRE(radiance.isConstant(radiance(0)));
}


TEST_CASE("Check Wavelength Independence Multiple Threads", "[sktran_do][lowlevel]") {
    int nwavel = 1000;

    sasktran_disco_lowlevel::CPPApi cppapi(16, 2, nwavel, 1, 1, 0);

    // Layer Construction
    for(int i = 0; i < nwavel; ++i) {
        cppapi.od(0, i) = 0.2;
        cppapi.od(1, i) = 0.2;

        cppapi.ssa(0, i) = 0.8;
        cppapi.ssa(1, i) = 0.8;

        cppapi.a1(0, 0, i) = 1;
        cppapi.a1(0, 1, i) = 1;
    }


    //cppapi.a1(2, 0, 0) = 0.5;
    //cppapi.a1(2, 1, 0) = 0.5;

    cppapi.layerboundaryaltitude(0) = 0.0001;
    cppapi.layerboundaryaltitude(1) = 0.00001;

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

    config.nthreads = 2;

    sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);

    Eigen::Map<Eigen::MatrixXd> radiance(output.radiance, 1, nwavel);

    REQUIRE(radiance.isConstant(radiance(0)));
}



TEST_CASE("Single Layer Weighting Function", "[sktran_do][lowlevel]") {
    int nwavel = 1;
    int nlyr = 1;
    int nderiv = 1;

    sasktran_disco_lowlevel::CPPApi cppapi(16, nlyr, nwavel, 1, 1, nderiv * nlyr);

    for (int i = 0; i < nwavel; ++i) {
        for (int j = 0; j < nlyr; ++j) {
            cppapi.od(j, i) = 0.5;
            cppapi.ssa(j, i) = 1;
            cppapi.a1(0, j, i) = 1;
            cppapi.a1(2, j, i) = 0.5;

            for (int k = 0; k < nderiv; ++k) {
                cppapi.d_layerindex(k*nlyr + j) = j;
                cppapi.d_od(k*nlyr + j, i) = 1;
            }
        }
    }
    // Layer Construction
    for (int j = 0; j < nlyr; ++j) {
        cppapi.layerboundaryaltitude(j) = 100000 - j*(100000.0 / double(nlyr));

    }

    cppapi.earth_radius() = 6372000;

    // Coordinates construction
    cppapi.cos_sza() = 0.2;
    cppapi.cos_vza(0) = 0.02;
    cppapi.saa(0) = 0;


    sasktran_disco_lowlevel::Atmosphere atmosphere;
    sasktran_disco_lowlevel::Config config;
    sasktran_disco_lowlevel::WeightingFunctions weightingfunctions;
    sasktran_disco_lowlevel::ViewingGeometry geometry;
    sasktran_disco_lowlevel::Output output;
    cppapi.initialize_c_api(&atmosphere, &config, &geometry, &weightingfunctions, &output);


    sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &geometry, &output);

}