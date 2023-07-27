#pragma once

#include <sktran_disco/sktran_do.h>


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
    void run_lowlevel_from_old_testspec(sasktran_disco::testing::TestCase<NSTOKES> &test_case,
                                        sasktran_disco::SKTRAN_DO_TestSpec<NSTOKES, -1> *test_spec,
                                        std::vector<double> &radiance,
                                        std::vector<double> &abs_diff,
                                        std::vector <std::vector<double>>* wf
    );
}