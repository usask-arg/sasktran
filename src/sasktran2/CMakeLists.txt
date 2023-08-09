
# Add an interface target for storing build properties to be inherited by the sasktran2 object
add_library(sasktran2BuildProperties INTERFACE)

# Set the include directories
target_include_directories(sasktran2BuildProperties INTERFACE
		$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:include>)

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
		    viewinggeometry/tangentaltitudesolar.cpp
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

target_link_libraries(sasktran2 PUBLIC sasktran2BuildProperties sasktranBuildProperties)

include(CTest)
add_subdirectory(tests)

install (TARGETS sasktran2BuildProperties
		EXPORT ${EXPORT_NAME}
		LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
		INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR})

install (TARGETS sasktran2
		EXPORT ${EXPORT_NAME}
		LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
		INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR})