#include <sasktran2/test_helper.h>

#include <sasktran2/testing/do_test_util.h>


namespace sasktran_disco::testing {
    /** This is an interface function to run the old SASKTRAN DISCO tests that used the sasktran_disco::testing::TestCase
     *  class along with the standard DO engine.  Here we transform the TestCase input to input used in the lowlevel
     *  interface instead.
     *
     * @tparam NSTOKES
     * @param test_case
     * @param test_spec
     * @param radiance
     * @param wf
     */
    template<int NSTOKES>
    void run_lowlevel_from_old_testspec(sasktran_disco::testing::TestCase<NSTOKES>& test_case,
                                        sasktran_disco::SKTRAN_DO_TestSpec<NSTOKES, -1>* test_spec,
                                        std::vector<double>& radiance,
                                        std::vector<double>& abs_diff,
                                        std::vector<std::vector<double>>* wf
                                        )
    {
        // Set up the low level interface
        int nstr = test_case.nstr;
        int nlyr = test_case.nlyr;
        int nlos = (int)test_case.linesofsight.size();
        int nderiv;

        if(test_spec) {
            int nderiv = (int)test_spec->perturbations()->size();
        }
        else {
            nderiv = 0;
        }

        sasktran_disco_lowlevel::CPPApi cppapi(nstr, nlyr, 1, NSTOKES, nlos, nderiv);

        // Copy the atmosphere
        for(int l = 0; l < nlyr; ++l) {
            cppapi.od(l, 0) = test_case.layers[l].optical_depth;
            cppapi.ssa(l, 0) = test_case.layers[l].ssa;

            for(int k = 0; k < nstr; ++k) {
                cppapi.a1(k, l, 0) = test_case.layers[l].lephasef[k].a1;

                if constexpr (NSTOKES == 3) {
                    cppapi.a2(k, l, 0) = test_case.layers[l].lephasef[k].a2;
                    cppapi.a3(k, l, 0) = test_case.layers[l].lephasef[k].a3;
                    cppapi.b1(k, l, 0) = test_case.layers[l].lephasef[k].b1;
                }
            }
        }

        // Now the albedo
        cppapi.albedo(0) = test_case.lambertian;

        // Copy the solar zenith
        cppapi.cos_sza() = test_case.solar.csz;

        // Set the layer boundaries, not used since we aren't using pseudo spherical anyways
        for(int l = 0; l < nlyr; ++l) {
            cppapi.layerboundaryaltitude(l) = 100000.0 - double((100000.0) * l) / double(nlyr);
        }
        cppapi.earth_radius() = 6372000;

        // And the LOS angles
        for(int i = 0; i < nlos; ++i) {
            cppapi.cos_vza(i) = test_case.linesofsight[i].coszen;
            cppapi.saa(i) = test_case.solar.saz - test_case.linesofsight[i].az;
        }

        sasktran_disco_lowlevel::Config config;
        sasktran_disco_lowlevel::Output output;
        sasktran_disco_lowlevel::Atmosphere atmosphere;
        sasktran_disco_lowlevel::ViewingGeometry viewinggeo;
        sasktran_disco_lowlevel::WeightingFunctions weightingfunctions;

        cppapi.initialize_c_api(&atmosphere, &config, &viewinggeo, &weightingfunctions, &output);

        // Set additional config options
        if(test_spec) {
            config.numazimuthexpansion = test_spec->getForcedNumberAzimuthTerms();
        }
        config.useexactsinglescatter = false;
        config.usepseudospherical = false;

        sasktran_disco_lowlevel::calculate(&atmosphere, &config, &weightingfunctions, &viewinggeo, &output);

        // Copy the results to the old output format
        radiance.resize(nlos*NSTOKES);
        abs_diff.resize(nlos*NSTOKES);
        for(int i = 0; i < nlos*NSTOKES; ++i) {
            radiance[i] = output.radiance[i] * test_case.solar.intensities.direct;
            abs_diff[i] = abs(radiance[i] - (*test_case.correct_radiances)[i]);
        }
    }

    template void run_lowlevel_from_old_testspec(sasktran_disco::testing::TestCase<1>& test_case,
                                                 sasktran_disco::SKTRAN_DO_TestSpec<1, -1>* test_spec,
                                                 std::vector<double>& radiance,
                                                 std::vector<double>& abs_diff,
                                                 std::vector<std::vector<double>>* wf
    );

    template void run_lowlevel_from_old_testspec(sasktran_disco::testing::TestCase<3>& test_case,
                                                 sasktran_disco::SKTRAN_DO_TestSpec<3, -1>* test_spec,
                                                 std::vector<double>& radiance,
                                                 std::vector<double>& abs_diff,
                                                 std::vector<std::vector<double>>* wf
    );

}