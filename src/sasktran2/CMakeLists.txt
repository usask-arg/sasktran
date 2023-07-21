
# Add an interface target for storing build properties to be inherited by the sasktran2 object
add_library(sasktran2BuildProperties INTERFACE)

# Set the include directories
target_include_directories(sasktran2BuildProperties INTERFACE
		${CMAKE_SOURCE_DIR}/include)

# Find dependencies and "link" them to sasktran2
set(Boost_USE_STATIC_LIBS ${FORCE_BOOST_STATIC})
find_package(Boost 1.60.0 REQUIRED COMPONENTS system regex timer log)
target_link_libraries(sasktran2BuildProperties INTERFACE Boost::system Boost::regex Boost::timer Boost::log)

find_package(OpenMP REQUIRED)
target_link_libraries(sasktran2BuildProperties INTERFACE OpenMP::OpenMP_CXX)

find_package(Eigen3 REQUIRED)
target_link_libraries(sasktran2BuildProperties INTERFACE Eigen3::Eigen)


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
	target_include_directories(sasktran2BuildProperties INTERFACE /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX12.3.sdk/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers)
	target_compile_definitions(sasktran2BuildProperties INTERFACE SKTRAN_USE_ACCELERATE)
else()
	set(BLA_STATIC FALSE)
	set(BLA_SIZEOF_INTEGER 4)
endif()
set(BLA_VENDOR ${SKTRAN_BLAS_VENDOR})
set(CMAKE_FIND_LIBRARY_PREFIXES "" lib)  # openblas include a "lib" prefix in it's names
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_link_libraries(sasktran2BuildProperties INTERFACE BLAS::BLAS LAPACK::LAPACK)
target_compile_definitions(sasktran2BuildProperties INTERFACE
		$<$<AND:$<PLATFORM_ID:Windows>,$<STREQUAL:${BLA_VENDOR},OpenBLAS>>:__WIN64__>   # for compatibility with Eigen/src/misc/blas.h
		)

if(SKTRAN_BLAS_VENDOR MATCHES "Intel")
	find_path(LAPACKE_H_INCLUDE_DIR mkl_lapacke.h REQUIRED)  # include directory for <lapacke.h>
	target_include_directories(sasktran2BuildProperties INTERFACE ${LAPACKE_H_INCLUDE_DIR})
	target_compile_definitions(sasktran2BuildProperties INTERFACE SKTRAN_USE_MKL)
else()
	find_path(LAPACKE_H_INCLUDE_DIR lapacke.h REQUIRED)  # include directory for <lapacke.h>
	find_path(CBLAS_H_INCLUDE_DIR cblas.h REQUIRED)    # include directory for <cblas.h>
	target_include_directories(sasktran2BuildProperties INTERFACE ${LAPACKE_H_INCLUDE_DIR} ${CBLAS_H_INCLUDE_DIR})

	if(WIN32)
		# OpenBlas can give warnings about std::complex linkage on windows/visual studio, but we don't use any complex lapack functions anyway so just disable the warning
		target_compile_options(sasktran2BuildProperties INTERFACE /wd4190)
	endif()

endif()


# Add sub libraries
add_library(sasktran2 
			STATIC

			coordinates/geodetic.cpp
		    coordinates/coordinates.cpp

		    config/config.cpp

		    unitsphere/lebedev_autogen.cpp
		    unitsphere/lebedev.cpp

			grids/grid.cpp
			grids/source_location_interpolator.cpp
		    grids/altitude_sza_source_location_interpolator.cpp
			dual/dual.cpp
			atmosphere/atmosphere.cpp
			atmosphere/surface.cpp
			atmosphere/userdefined1d_constituent.cpp
			phasefunction/interpolator.cpp
		    sourceintegrator/sourceintegrator.cpp

			geometry/geometry.cpp
			geometry/geometry1d.cpp
			viewinggeometry/tangentaltitude.cpp
		    viewinggeometry/groundviewing.cpp
			raytracing/spherical_shell.cpp
			raytracing/plane_parallel.cpp
		    solar/solartransmissionexact.cpp
		    solar/solartransmissiontable.cpp
		    solar/singlescattersource.cpp
		    solar/occultation.cpp

			output/output.cpp
		    output/outputidealdense.cpp

		    engine/engine.cpp

		    hr/diffuse_point.cpp
		    hr/diffuse_table.cpp

			sktran_disco/sktran_do.cpp
			sktran_disco/sktran_do_geometrylayerarray.cpp
			sktran_disco/sktran_do_layerarray.cpp
			sktran_disco/sktran_do_lazyazimuth.cpp
			sktran_disco/sktran_do_linearization_types.cpp
			sktran_disco/sktran_do_lowlevelinterface.cpp
			sktran_disco/sktran_do_misc.cpp
			sktran_disco/sktran_do_opticallayer.cpp
			sktran_disco/sktran_do_pconfig.cpp
			sktran_disco/sktran_do_pentadiagonal.cpp
			sktran_disco/sktran_do_polarization_types.cpp
			sktran_disco/sktran_do_quadrature.cpp
			sktran_disco/sktran_do_rte.cpp
			sktran_disco/sktran_do_specs.cpp
			sktran_disco/sktran_do_surface.cpp
			sktran_disco/sktran_do_testing.cpp
			sktran_disco/sktran_do_types.cpp
		sktran_disco/source_term/do_source_legendrestorage.cpp
		sktran_disco/source_term/do_source_diffuse_storage.cpp
        sktran_disco/source_term/do_source.cpp
		sktran_disco/source_term/do_source_exact_pp.cpp
		sktran_disco/source_term/do_source_interpolated_pp.cpp
		sktran_disco/source_term/do_source_spherical_corrected.cpp
		sktran_disco/sktran_do_dll_expose.cpp)

target_link_libraries(sasktran2 PUBLIC sasktran2BuildProperties)

include(CTest)
add_subdirectory(tests)


# On install, copy sasktran2/include directory to install_directory/include
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sasktran2)
