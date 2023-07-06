cmake_minimum_required      (VERSION 3.16)                                                          # Version 3.16 provides support for pre-compiled headers
include                     (CMakePackageConfigHelpers )                                            # Provids function/macro write_basic_package_version_file

set(CMAKE_VERBOSE_MAKEFILE ON)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

add_definitions("-Wno-inconsistent-missing-override -Wno-suggest-override")

execute_process(COMMAND git describe
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DESCRIBE_TEMP
)

if(EXISTS ${GIT_DESCRIBE_TEMP})
    string(REGEX MATCH "v([0-9\\.]+)" SASKTRAN_CORE_VERSION ${GIT_DESCRIBE_TEMP})
    string(REPLACE "v" "" SASKTRAN_CORE_VERSION ${GIT_DESCRIBE_TEMP})
else()
    set(SASKTRAN_CORE_VERSION 1.5.0)
endif()

project (sasktran-core VERSION ${SASKTRAN_CORE_VERSION})

# Fortran library settings, windows we add a link directory, linux we build it
if(WIN32)
    link_directories(${CMAKE_SOURCE_DIR}/fortran_libraries/lib/Windows_x64/)
else()
    enable_language(Fortran)
endif()

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(FORCE_BOOST_STATIC "Force the use of boost static libraries" ON)
# Try to find boost 1.6 or greater
set(Boost_USE_STATIC_LIBS ${FORCE_BOOST_STATIC})
find_package(Boost REQUIRED COMPONENTS filesystem system regex timer random thread)

# Try to find OpenMP
find_package(OpenMP REQUIRED)

# Look for NetCDF, yaml-cpp, and HDF5
find_package(netCDF REQUIRED CONFIG)

find_package(yaml-cpp CONFIG REQUIRED)

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

# add the additional compiler flags
set (CMAKE_CXX_FLAGS " ${CompilerFlags} ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Directory to install the static library files to
set (STATIC_LIB_INSTALL_DIR lib)
set (INCLUDE_INSTALL_DIR include)
set (EXPORT_NAME sasktran-coreTargets)

# Platform independent including
include_directories(
    ${CMAKE_SOURCE_DIR}/Repos_SasktranV3
    ${CMAKE_SOURCE_DIR}/Repos_skclimatology
    ${CMAKE_SOURCE_DIR}/Repos_BaseCode/nxbase
    ${CMAKE_SOURCE_DIR}/Repos_SasktranIF/includes
    ${CMAKE_SOURCE_DIR}/Repos_BaseCode/nxhdf/nxnetcdfio
    ${CMAKE_SOURCE_DIR}/Repos_skopticalproperties

    ${Boost_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    ${netCDF_INCLUDE_DIR}
)

# Platform definitions
if(WIN32)
    add_definitions(-DNOMINMAX)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else()
    add_definitions(-DBOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED)
endif()

# Platform library names 
if(WIN32)
    set(CMAKE_RELEASE_POSTFIX _Release)
    set(CMAKE_DEBUG_POSTFIX _Debug)
endif()
set(NXBASE_NAME nxbase)
set(NXNETCDFIO_NAME nxnetcdfio)

set(SASKTRANIF_NAME sasktranif)
set(SKCLIMATOLOGY_NAME skclimatology21)
set(SKOPTICALPROPERTIES_NAME skopticalproperties21)
set(SASKTRAN_BASE_NAME sktran_base)

set(SASKTRAN_CORE_NAME _sasktran_core_internals)
set(SASKTRAN_STUB_NAME sktran_stubs)

write_basic_package_version_file( "sasktran-coreconfigVersion.cmake"
        VERSION ${sasktran-core_VERSION}
        COMPATIBILITY SameMajorVersion
        )


# Add subprojects
add_subdirectory (Repos_BaseCode)
add_subdirectory (Repos_SasktranIF)
add_subdirectory (Repos_skclimatology)
add_subdirectory (Repos_skopticalproperties)
add_subdirectory (Repos_SasktranV3)
add_subdirectory (fortran_libraries)

install(EXPORT sasktran-coreTargets
        FILE sasktran-coreTargets.cmake
        NAMESPACE sasktran-core::
        DESTINATION lib/cmake/sasktran-core
        )

install( FILES "sasktran-coreConfig.cmake" DESTINATION lib/cmake/sasktran-core)

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
