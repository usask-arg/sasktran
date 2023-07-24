cmake_minimum_required(VERSION 3.20.0)
include                     (CMakePackageConfigHelpers )                                            # Provids function/macro write_basic_package_version_file

set(CMAKE_VERBOSE_MAKEFILE ON)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(UNIX)
    add_definitions("-Wno-inconsistent-missing-override -Wno-suggest-override -Wno-comment -Wno-deprecated-declarations")
endif()

execute_process(COMMAND git describe --tags
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

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    # This is nominally a windows definition but for some reason BOOST < 1.8.3 has problems on clang without this?
    add_definitions(-D_HAS_AUTO_PTR_ETC 0)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
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

find_package(Eigen3 REQUIRED)

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
    ${CMAKE_SOURCE_DIR}/src/core/base
    ${CMAKE_SOURCE_DIR}/src/core/sasktran
    ${CMAKE_SOURCE_DIR}/src/core/base/nxbase
    ${CMAKE_SOURCE_DIR}/src/core/sasktranif/includes
    ${CMAKE_SOURCE_DIR}/src/core/base/nxhdf/nxnetcdfio
    ${CMAKE_SOURCE_DIR}/src/core/skopticalproperties
    ${CMAKE_SOURCE_DIR}/src/core/skclimatology
    ${CMAKE_SOURCE_DIR}/include


    ${Boost_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    ${netCDF_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
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

set ( SKTRAN_TIR_NAME   "sktran_tir_engine")

write_basic_package_version_file( "sasktran-coreconfigVersion.cmake"
        VERSION ${sasktran-core_VERSION}
        COMPATIBILITY SameMajorVersion
        )


# Add subprojects
add_subdirectory (src/core/base)
add_subdirectory (src/core/sasktranif)
add_subdirectory (src/core/skclimatology)
add_subdirectory (src/core/skopticalproperties)
add_subdirectory (src/core/sasktran)
add_subdirectory (fortran_libraries)

add_subdirectory (src/sasktran2)

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
