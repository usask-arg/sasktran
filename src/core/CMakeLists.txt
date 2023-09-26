# Build properties for the sasktranCore project

# Fortran library settings, windows we add a link directory, linux we build it
if(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/src/core/fortran_libraries/lib/Windows_x64/)
endif()

add_library(sasktranCoreBuildProperties INTERFACE)

target_link_libraries(sasktranCoreBuildProperties INTERFACE sasktranBuildProperties)

target_include_directories(sasktranCoreBuildProperties INTERFACE
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/base>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/sasktran>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/base/nxbase>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/sasktranif/includes>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/base/nxhdf/nxnetcdfio>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/skopticalproperties>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/core/skclimatology>
        $<INSTALL_INTERFACE:include/base>
        $<INSTALL_INTERFACE:include/sasktran>
        $<INSTALL_INTERFACE:include/base/nxbase>
        $<INSTALL_INTERFACE:include/sasktranif/includes>
        $<INSTALL_INTERFACE:include/base/nxhdf/nxnetcdfio>
        $<INSTALL_INTERFACE:include/skopticalproperties>
        $<INSTALL_INTERFACE:include/skclimatology>
        )

# Look for NetCDF, yaml-cpp, and HDF5
find_package(netCDF REQUIRED CONFIG)

find_package(yaml-cpp CONFIG REQUIRED)

message(STATUS "Looking for HDF5")
# hdf5 is annoying because sometimes it doesn't bundle the config
# First try to find config HDF5
find_package(hdf5 CONFIG)

if(NOT hdf5_FOUND)
    # Try to find the module instead and fail if we can't
    find_package(HDF5 REQUIRED)
    set(HDF5_LIBRARY HDF5::HDF5)
else()
    # we found the hdf5 config, but even this does not seem to define proper variables so we have to look for targets
    if(TARGET hdf5::hdf5-static)
        set(HDF5_LIBRARY hdf5::hdf5-static)
    elseif(TARGET hdf5)
        set(HDF5_LIBRARY hdf5)
    else()
        message(ERROR "Could not determine hdf5 target")
    endif()
endif()

message(STATUS "Found HDF5")


target_link_libraries(sasktranCoreBuildProperties INTERFACE yaml-cpp)
target_link_libraries(sasktranCoreBuildProperties INTERFACE ${netCDF_LIBRARIES})
target_link_libraries(sasktranCoreBuildProperties INTERFACE ${HDF5_LIBRARY})


# Platform definitions
if(WIN32)
    target_compile_definitions(sasktranCoreBuildProperties INTERFACE NOMINMAX)
    target_compile_definitions(sasktranCoreBuildProperties INTERFACE _CRT_SECURE_NO_WARNINGS)
    target_compile_definitions(sasktranCoreBuildProperties INTERFACE _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING)
else()
    target_compile_definitions(sasktranCoreBuildProperties INTERFACE BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED)
endif()


# Platform library names
if(WIN32)
    set(CMAKE_RELEASE_POSTFIX _Release)
    set(CMAKE_DEBUG_POSTFIX _Debug)
endif()


# Project naming
set(NXBASE_NAME nxbase)
set(NXNETCDFIO_NAME nxnetcdfio)

set(SASKTRANIF_NAME sasktranif)
set(SKCLIMATOLOGY_NAME skclimatology21)
set(SKOPTICALPROPERTIES_NAME skopticalproperties21)
set(SASKTRAN_BASE_NAME sktran_base)

set(SASKTRAN_CORE_NAME _sasktran_core_internals)
set(SASKTRAN_STUB_NAME sktran_stubs)

set(SKTRAN_TIR_NAME sktran_tir_engine)

write_basic_package_version_file( "../../sasktran-coreconfigVersion.cmake"
        VERSION ${sasktran-core_VERSION}
        COMPATIBILITY SameMajorVersion
        )

install (TARGETS sasktranCoreBuildProperties
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR})


# Project installation
install(EXPORT sasktran-coreTargets
        FILE sasktran-coreTargets.cmake
        NAMESPACE sasktran-core::
        DESTINATION lib/cmake/sasktran-core
        )

install( FILES "../../sasktran-coreConfig.cmake" DESTINATION lib/cmake/sasktran-core)

if(WIN32)
    # Bundle the fortran libraries along with the installed libraries
    install(DIRECTORY "fortran_libraries/lib/Windows_x64/"
            DESTINATION lib
            FILES_MATCHING
            PATTERN "*.dll"
            )

    install(DIRECTORY "fortran_libraries/lib/Windows_x64/"
            DESTINATION lib
            FILES_MATCHING
            PATTERN "*.lib"
            )
endif()


# Add subprojects
add_subdirectory (base)
add_subdirectory (sasktranif)
add_subdirectory (skclimatology)
add_subdirectory (skopticalproperties)
add_subdirectory (sasktran)
add_subdirectory (fortran_libraries)