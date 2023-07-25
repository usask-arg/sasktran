cmake_minimum_required(VERSION 3.20.0)
include                     (CMakePackageConfigHelpers )                                            # Provids function/macro write_basic_package_version_file

set(CMAKE_VERBOSE_MAKEFILE ON)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

if(UNIX)
    add_definitions("-Wno-inconsistent-missing-override -Wno-suggest-override -Wno-comment -Wno-deprecated-declarations -Wno-psabi")
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

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    # This is nominally a windows definition but for some reason BOOST < 1.8.3 has problems on clang without this?
    add_compile_definitions(_HAS_AUTO_PTR_ETC=0)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(FORCE_BOOST_STATIC "Force the use of boost static libraries" ON)
# Try to find boost 1.6 or greater
set(Boost_USE_STATIC_LIBS ${FORCE_BOOST_STATIC})
find_package(Boost REQUIRED COMPONENTS filesystem system regex timer random thread log)

# Try to find OpenMP
find_package(OpenMP REQUIRED)


find_package(Eigen3 REQUIRED)

# add the additional compiler flags
set (CMAKE_CXX_FLAGS " ${CompilerFlags} ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Directory to install the static library files to
set (STATIC_LIB_INSTALL_DIR lib)
set (INCLUDE_INSTALL_DIR include)
set (EXPORT_NAME sasktran-coreTargets)


# Add an interface target for storing build properties common to all sasktran modules
add_library(sasktranBuildProperties INTERFACE)

target_link_libraries(sasktranBuildProperties INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(sasktranBuildProperties INTERFACE Eigen3::Eigen)
target_link_libraries(sasktranBuildProperties INTERFACE
        Boost::system
        Boost::regex
        Boost::timer
        Boost::log
        Boost::filesystem
        Boost::random
        Boost::thread
        )

# -------- LAPACK FINDING -----------
# We support multiple blas vendors, but these can require some funny options
set(SKTRAN_BLAS_VENDOR "OpenBLAS" CACHE STRING
        "Sets the BLAS/LAPACK vendor. See https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors.")
set_property(CACHE SKTRAN_BLAS_VENDOR PROPERTY STRINGS OpenBLAS Intel10_64lp Apple)

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
    set(BLA_STATIC FALSE)
    set(BLA_SIZEOF_INTEGER 4)
elseif(SKTRAN_BLAS_VENDOR MATCHES "Apple")
    set(BLA_STATIC TRUE)
    set(BLA_SIZEOF_INTEGER 4)
    target_include_directories(sasktranBuildProperties INTERFACE /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers)
    target_compile_definitions(sasktranBuildProperties INTERFACE SKTRAN_USE_ACCELERATE)
else()
    set(BLA_STATIC FALSE)
    set(BLA_SIZEOF_INTEGER 4)
endif()
set(BLA_VENDOR ${SKTRAN_BLAS_VENDOR})
set(CMAKE_FIND_LIBRARY_PREFIXES "" lib)  # openblas include a "lib" prefix in it's names
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(sasktranBuildProperties INTERFACE BLAS::BLAS LAPACK::LAPACK)
target_compile_definitions(sasktranBuildProperties INTERFACE
        $<$<AND:$<PLATFORM_ID:Windows>,$<STREQUAL:${BLA_VENDOR},OpenBLAS>>:__WIN64__>   # for compatibility with Eigen/src/misc/blas.h
        )

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
    find_path(LAPACKE_H_INCLUDE_DIR mkl_lapacke.h REQUIRED)  # include directory for <lapacke.h>
    target_include_directories(sasktranBuildProperties INTERFACE ${LAPACKE_H_INCLUDE_DIR})
    target_compile_definitions(sasktranBuildProperties INTERFACE SKTRAN_USE_MKL)
else()
    find_path(LAPACKE_H_INCLUDE_DIR lapacke.h REQUIRED)  # include directory for <lapacke.h>
    find_path(CBLAS_H_INCLUDE_DIR cblas.h REQUIRED)    # include directory for <cblas.h>
    target_include_directories(sasktranBuildProperties INTERFACE ${LAPACKE_H_INCLUDE_DIR} ${CBLAS_H_INCLUDE_DIR})

    if(WIN32)
        # OpenBlas can give warnings about std::complex linkage on windows/visual studio, but we don't use any complex lapack functions anyway so just disable the warning
        target_compile_options(sasktranBuildProperties INTERFACE /wd4190)
    endif()

endif()

# -------- END LAPACK FINDING -------

install (TARGETS sasktranBuildProperties
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR})

add_subdirectory (src/sasktran2)
add_subdirectory (src/core)
