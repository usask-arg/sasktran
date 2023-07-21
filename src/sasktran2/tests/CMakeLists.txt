add_executable(sasktran2_tests
    test_main.cpp
    dual/test_dual.cpp
    viewinggeometry/test_tangentaltitude.cpp
    viewinggeometry/test_groundviewing.cpp
    unitsphere/test_lebedev.cpp
    geometry/test_geometry.cpp
    grids/test_grid.cpp
    raytracing/test_spherical_shell.cpp
    raytracing/test_plane_parallel.cpp
    engine/test_engine_basic.cpp

    engine/hr/test_hr_basic.cpp

    sktran_disco/test_util.cpp
    sktran_disco/tests_lowlevel.cpp

    sktran_disco/legacy/test_scalar.cpp
    sktran_disco/legacy/test_vector.cpp

    sktran_disco/benchmark/bench_lowlevel.cpp
 "engine/benchmark/bench_singlescatter.cpp")

# Some special logic to support both CATCH2 version 2 and 3 for now
find_package(Catch2 2 CONFIG QUIET)

target_link_libraries(sasktran2_tests PRIVATE
                        sasktran2
                        $<$<PLATFORM_ID:Windows>:bcrypt>
                        $<$<PLATFORM_ID:Linux,Darwin>:dl>
                        $<$<PLATFORM_ID:Linux,Darwin>:gfortran>
                        $<$<PLATFORM_ID:Linux>:rt>
                        $<$<PLATFORM_ID:Linux>:z>
                        )
if(Catch2_FOUND)
    # Found Catch2 Version 2
    target_link_libraries(sasktran2_tests PRIVATE Catch2::Catch2)
else()
    # Look for Catch2 Version 3
    message("Could not find Catch2 Version 2, Looking for Version 3")
    find_package(Catch2 3 CONFIG REQUIRED)
    target_link_libraries(sasktran2_tests PRIVATE Catch2::Catch2WithMain)
    target_compile_definitions(sasktran2_tests PUBLIC SKTRAN_CATCH2_VERSION3)
endif()


add_test(NAME TestSASKTRAN2 COMMAND sasktran2_tests)
add_test(NAME TestSasktranDiscoScalar COMMAND sktran_disco_tests [sktran_disco])
#add_test(NAME TestSasktranDiscoVector COMMAND sktran_disco_engine_tests [vector])
#add_test(NAME TestSasktranDiscoLowLevel COMMAND sktran_disco_engine_tests [lowlevel])