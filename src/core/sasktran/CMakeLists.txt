# ssasktranv3 source code

set(sasktran_common_src
    modules/sktran_common/curvedrays/sktran_raystorage_curved.cpp
	modules/sktran_common/curvedrays/sktran_raytracer_curved_shells.cpp
    modules/sktran_common/curvedrays/sktran_ray_curved.cpp
	modules/sktran_common/engine/sktran_enginebase.cxx
	modules/sktran_common/diagnostics/sktran_diagnosticinterface.cxx
	modules/sktran_common/emissiontables/sktran_emissiontable_1d.cpp
	modules/sktran_common/emissiontables/sktran_emissiontable_base.cpp
	modules/sktran_common/emissiontables/sktran_emissiontable_notable.cpp
	modules/sktran_common/grids/sktran_griddefscatterangle.cpp
	modules/sktran_common/grids/sktran_raytracingshells.cpp
	modules/sktran_common/miscellaneous/sktran_codetimer.cxx
	modules/sktran_common/miscellaneous/sktran_sourceterm.cpp
	modules/sktran_common/opticalpropertytables/sktran_opticalproperties_1d_heightv3.cpp
	modules/sktran_common/opticalpropertytables/sktran_opticalproperties_3d_unitsphere.cpp
	modules/sktran_common/opticalpropertytables/sktran_polarizationprops.cpp
	modules/sktran_common/opticalpropertytables/sktran_scattermatrix.cpp
	modules/sktran_common/quadrature/opticaldepthcalculator_linearwithheight.cxx
	modules/sktran_common/quadrature/sktran_integrator_adaptive.cpp
	modules/sktran_common/quadrature/sktran_integrator_base.cpp
	modules/sktran_common/quadrature/sktran_integrator_straight.cpp
	modules/sktran_common/quadrature/sktran_srcintegrator.cpp
	modules/sktran_common/solartransmissiontables/sktran_solartransmission_2d.cpp
	modules/sktran_common/solartransmissiontables/sktran_solartransmission_3d.cpp
	modules/sktran_common/solartransmissiontables/sktran_solartransmission_base.cpp
	modules/sktran_common/solartransmissiontables/sktran_solartransmission_notable.cpp
	modules/sktran_common/solartransmissiontables/sktran_sun.cpp
	modules/sktran_common/specifications/raytracingregionmanager.cxx
	modules/sktran_common/specifications/sktran_lineofsightarray.cxx
	modules/sktran_common/straightrays/sktran_indexofrefraction.cpp
	modules/sktran_common/straightrays/sktran_rayminimumcontainer.cpp
	modules/sktran_common/straightrays/sktran_raytracer_shells.cpp
	modules/sktran_common/straightrays/sktran_raytracer_straight_generic.cpp
	modules/sktran_common/straightrays/sktran_ray_base.cpp
	modules/sktran_common/straightrays/sktran_ray_straight.cpp
	modules/sktran_common/straightrays/sktran_storage_curved.cpp
	modules/sktran_common/unitspheres/sktran_unitsphere.cxx
	modules/sktran_common/unitspheres/sktran_unitsphere_bucky.cxx
	modules/sktran_common/unitspheres/sktran_unitsphere_delaunay.cpp
	modules/sktran_common/unitspheres/sktran_unitsphere_dummy.cpp
	modules/sktran_common/unitspheres/sktran_unitsphere_latlongrid.cxx
	modules/sktran_common/unitspheres/sktran_unitsphere_me.cpp
	modules/sktran_common/unitspheres/sktran_unitsphere_plane.cpp
)

set(sasktran_mc_src
    modules/monte_carlo/sktran_averagingkernel.cpp
	modules/monte_carlo/sktran_configurationmanager_mc.cpp
	modules/monte_carlo/sktran_engine_mc_v21.cpp
	modules/monte_carlo/sktran_erptmanager.cpp
	modules/monte_carlo/sktran_hpflipop.cpp
	modules/monte_carlo/sktran_integrator_mc.cpp
    modules/monte_carlo/sktran_mc_optimal_scatter_sequence.cpp
	modules/monte_carlo/sktran_mc.cpp
	modules/monte_carlo/sktran_mcamfcalculator.cpp
	modules/monte_carlo/sktran_mcobservermanager.cpp
	modules/monte_carlo/sktran_mcphoton.cpp
	modules/monte_carlo/sktran_mcscatteroperator.cpp
    modules/monte_carlo/sktran_mcsolartransmission.cpp
	modules/monte_carlo/sktran_random_number_generator.cpp
	modules/monte_carlo/sktran_specifications_mc.cpp
	modules/monte_carlo/sktran_sun_randomdisc.cpp
    modules/monte_carlo/sktran_tableopticalproperties_inelastic_mc.cpp
	modules/monte_carlo/sktran_tableopticalproperties_mc.cpp
)
set(sasktran_hr_src
    modules/sktran_highres/sktran_diffuse_table_cpu_eatcalculator_dorotation.cpp
	modules/sktran_highres/sktran_diffuse_table_cpu_etacalculator_norotation.cpp
	modules/sktran_highres/sktran_hr.cpp
	modules/sktran_highres/sktran_hr_debug.cpp
	modules/sktran_highres/sktran_hr_diffuse_aval_base.cpp
	modules/sktran_highres/sktran_hr_diffuse_aval_matrixstore.cpp
	modules/sktran_highres/sktran_hr_diffuse_aval_matrixtable.cpp
	modules/sktran_highres/sktran_hr_diffuse_aval_scalarstore.cpp
	modules/sktran_highres/sktran_hr_diffuse_point.cpp
	modules/sktran_highres/sktran_hr_diffuse_radstore_base.cpp
	modules/sktran_highres/sktran_hr_diffuse_radstore_polarized.cpp
	modules/sktran_highres/sktran_hr_diffuse_radstore_pseudovec1.cpp
	modules/sktran_highres/sktran_hr_diffuse_radstore_scalar.cpp
	modules/sktran_highres/sktran_hr_diffuse_second_order_source.cpp
	modules/sktran_highres/sktran_hr_diffuse_source.cpp
	modules/sktran_highres/sktran_hr_diffuse_table_cpu.cpp
	modules/sktran_highres/sktran_hr_diffuse_table_sza.cpp
	modules/sktran_highres/sktran_hr_engine.cpp
	modules/sktran_highres/sktran_hr_linesofsighttable.cpp
	modules/sktran_highres/sktran_hr_perturbation_absorption.cpp
	modules/sktran_highres/sktran_hr_rayregionmanager.cpp
	modules/sktran_highres/sktran_hr_specs_internal_core.cpp
	modules/sktran_highres/sktran_hr_specs_internal_diffuse.cpp
	modules/sktran_highres/sktran_hr_specs_internal_integrator.cpp
	modules/sktran_highres/sktran_hr_specs_internal_opticalpropertiestable.cpp
	modules/sktran_highres/sktran_hr_specs_internal_raytracer.cpp
	modules/sktran_highres/sktran_hr_specs_internal_wf.cpp
	modules/sktran_highres/sktran_hr_specs_user_diffuse.cpp
	modules/sktran_highres/sktran_hr_specs_user_integrator.cpp
	modules/sktran_highres/sktran_hr_specs_user_opticalpropertiestable.cpp
	modules/sktran_highres/sktran_hr_specs_user_raytracer.cpp
	modules/sktran_highres/sktran_hr_specs_user_wf.cpp
	modules/sktran_highres/sktran_hr_thread_manager.cpp
	modules/sktran_highres/sktran_hr_wf_speciesinformation.cpp
	modules/sktran_highres/sktran_hr_wf_extinction_table.cpp
	modules/sktran_highres/sktran_hr_wf_integrator.cpp
	modules/sktran_highres/sktran_hr_wf_ray.cpp
	modules/sktran_highres/sktran_hr_wf_store.cpp
	modules/sktran_highres/sktran_outgoingsphereobject.cpp
)
set(sasktran_so_src
    modules/sktran_so/jindex/sktran_jindex.cpp
	modules/sktran_so/sktran_legacy/sktran_opticalproperties_v2.cpp
	modules/sktran_so/sktran_legacy/sktran_raytracespecs_general.cpp
)
set(sasktran_occ_src
	modules/sktran_occ/skoccult_engine.cpp 
    modules/sktran_occ/skoccult_opticalproperties1d_heightwavelength.cpp 
    modules/sktran_occ/sktran_ray_curved_piecewise.cpp
	modules/sktran_occ/refrac.cpp  
)

set(sasktran_tir_src
	modules/sktran_tir/sktran_indexofrefraction_tir.cpp
	modules/sktran_tir/sktran_tir.cpp
	modules/sktran_tir/sktran_tir_atmosphericstate.cpp
	modules/sktran_tir/sktran_tir_engine.cpp
	modules/sktran_tir/sktran_tir_integrator.cpp
	modules/sktran_tir/sktran_tir_linesofsighttable.cpp
	modules/sktran_tir/sktran_tir_opticalproperties.cpp
	modules/sktran_tir/sktran_tir_rayregionmanager.cpp
	modules/sktran_tir/sktran_tir_raystorage.cpp
	modules/sktran_tir/sktran_tir_specs_internal_core.cpp
	modules/sktran_tir/sktran_tir_specs_internal_integrator.cpp
	modules/sktran_tir/sktran_tir_specs_internal_opticalpropertiestable.cpp
	modules/sktran_tir/sktran_tir_specs_internal_raytracer.cpp
	modules/sktran_tir/sktran_tir_specs_internal_wf.cpp
	modules/sktran_tir/sktran_tir_specs_user_integrator.cpp
	modules/sktran_tir/sktran_tir_specs_user_opticalpropertiestable.cpp
	modules/sktran_tir/sktran_tir_specs_user_raytracer.cpp
	modules/sktran_tir/sktran_tir_specs_user_wf.cpp
	modules/sktran_tir/sktran_tir_thread_manager.cpp
	modules/sktran_tir/sktran_tir_wf.cpp
		)

set(sasktran_do_src
		modules/sktran_do_deprecated/sktran_do.cpp
		modules/sktran_do_deprecated/sktran_do_engine.cpp
		modules/sktran_do_deprecated/sktran_do_geometrylayerarray.cpp
		modules/sktran_do_deprecated/sktran_do_layerarray.cpp
		modules/sktran_do_deprecated/sktran_do_lazyazimuth.cpp
		modules/sktran_do_deprecated/sktran_do_linearization_types.cpp
		modules/sktran_do_deprecated/sktran_do_misc.cpp
		modules/sktran_do_deprecated/sktran_do_opticallayer.cpp
		modules/sktran_do_deprecated/sktran_do_opticalstate.cpp
		modules/sktran_do_deprecated/sktran_do_pconfig.cpp
		modules/sktran_do_deprecated/sktran_do_pentadiagonal.cpp
		modules/sktran_do_deprecated/sktran_do_polarization_types.cpp
		modules/sktran_do_deprecated/sktran_do_postprocessing.cpp
		modules/sktran_do_deprecated/sktran_do_quadrature.cpp
		modules/sktran_do_deprecated/sktran_do_rte.cpp
		modules/sktran_do_deprecated/sktran_do_specs.cpp
		modules/sktran_do_deprecated/sktran_do_spherical.cpp
		modules/sktran_do_deprecated/sktran_do_surface.cpp
		modules/sktran_do_deprecated/sktran_do_testing.cpp
		modules/sktran_do_deprecated/sktran_do_types.cpp
		)

set(sasktran_me_src
	modules/sktran_me/sktran_me_atmosphere_interface.cpp
	modules/sktran_me/sktran_me_los_converter.cpp
	modules/sktran_me/sktran_me_refpoint_estimator.cpp
	modules/sktran_me/sktran_me_wf_handler.cpp
	modules/sktran_me/sktran_me_output.cpp
	)

set(sasktran_core_src
    modules/sasktranv3_impl/climatologystubs/iskclimatology_stub_base.cpp
	modules/sasktranv3_impl/climatologystubs/iskclimatology_stub_base_specifics.cpp
	modules/sasktranv3_impl/climatologystubs/iskclimatology_stub_linearcombo.cpp
	modules/sasktranv3_impl/climatologystubs/iskclimatology_stub_msis.cpp
	modules/sasktranv3_impl/enginestubs/iskengine_stub_hr.cpp
	modules/sasktranv3_impl/enginestubs/iskengine_stub_mc.cpp
	modules/sasktranv3_impl/enginestubs/iskengine_stub_occ.cpp
	modules/sasktranv3_impl/enginestubs/iskengine_stub_tir.cpp
	modules/sasktranv3_impl/enginestubs/iskengine_stub_do.cpp
	modules/sasktranv3_impl/enginestubs/iskengine_stub_me.cpp
	modules/sasktranv3_impl/enginestubs/iskgeodetic_stub_std.cpp
		 "modules/sasktranv3_impl/miestubs/iskmie_stub_wiscombe.cpp")

set(sasktran_dll_src
		modules/sasktranv3_impl/dllimplementation/configuresasktranregistry.cxx
		modules/sasktranv3_impl/dllimplementation/dllmain.cpp
		modules/sasktranv3_impl/dllimplementation/stdafx.cpp
		modules/sasktranv3_impl/dllimplementation/SasktranImpl.cpp
		)

set(sasktran_test_src
	Application/TestInstallation/runmc.cpp
	Application/TestInstallation/sasktran_unit_tests.cpp
	Application/TestInstallation/short_test_suite.cpp
	Application/TestInstallation/simple_ray_test.cpp
	Application/TestInstallation/test_curvedrays.cpp
	Application/TestInstallation/test_hitran.cpp
	Application/TestInstallation/test_hr.cpp 
	Application/TestInstallation/test_linearcombo.cpp 
	Application/TestInstallation/test_mc.cpp 
	Application/TestInstallation/test_mcamf.cpp 
	Application/TestInstallation/test_msisatomico.cpp 
	Application/TestInstallation/test_o2_o2.cpp 
	Application/TestInstallation/test_so2.cpp
	Application/TestInstallation/testdelaunaysphere.cpp 
	Application/TestInstallation/testscattermatrices.cpp
	Application/TestInstallation/test_sktran_me_stub.cpp
)

# Build the library
add_library(${SASKTRAN_BASE_NAME} ${sasktran_common_src} ${sasktran_hr_src} ${sasktran_mc_src} ${sasktran_so_src} ${sasktran_occ_src} ${sasktran_tir_src} ${sasktran_do_src} ${sasktran_me_src})

target_link_libraries(${SASKTRAN_BASE_NAME} PUBLIC sasktranCoreBuildProperties)
target_link_libraries(${SASKTRAN_BASE_NAME} PUBLIC sasktran2)


# Set the header files
install (TARGETS ${SASKTRAN_BASE_NAME}
		EXPORT ${EXPORT_NAME}
		LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
		INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${SASKTRAN_BASE_NAME})

install(DIRECTORY "."
		DESTINATION ${INCLUDE_INSTALL_DIR}/${SASKTRAN_BASE_NAME}
		FILES_MATCHING
		PATTERN "*.h"
		)


add_library(${SASKTRAN_STUB_NAME} ${sasktran_core_src})
add_library(${SASKTRAN_CORE_NAME} SHARED ${sasktran_dll_src})

target_link_libraries(${SASKTRAN_STUB_NAME} PUBLIC sasktranCoreBuildProperties)
target_link_libraries(${SASKTRAN_CORE_NAME} PUBLIC sasktranCoreBuildProperties)

target_link_libraries(${SASKTRAN_CORE_NAME} PUBLIC sasktran2)
target_link_libraries(${SASKTRAN_STUB_NAME} PUBLIC sasktran2)



# Set the header files
install (TARGETS ${SASKTRAN_STUB_NAME}
		EXPORT ${EXPORT_NAME}
		LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
		INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${SASKTRAN_BASE_NAME})


add_executable(sasktran_tests ${sasktran_test_src})

set_target_properties(${SASKTRAN_CORE_NAME} PROPERTIES RELEASE_POSTFIX "")
set_target_properties(${SASKTRAN_CORE_NAME} PROPERTIES DEBUG_POSTFIX "")

# Windows specific links
set(sasktran_win_links
	hitran_tips
	wiscombemie
	msis90e
	tmatrixrandomep
	bcrypt
)
# linux specific links
set(sasktran_unix_links
	dl
	rt
	z
)
# mac specific links
set(sasktran_darwin_links
		z
		gfortran
		)

set(sasktran_links
		${SASKTRAN_STUB_NAME}
		${SASKTRAN_BASE_NAME}
		${SKCLIMATOLOGY_NAME}
		${SKOPTICALPROPERTIES_NAME}
		${SASKTRANIF_NAME}
		${NXBASE_NAME}
		${NXNETCDFIO_NAME}
)

target_link_libraries(${SASKTRAN_CORE_NAME}
		PUBLIC
	${sasktran_links}
	)
target_link_libraries(sasktran_tests
		PUBLIC
		${sasktran_links}
		)

if(WIN32)
	target_link_libraries(${SASKTRAN_CORE_NAME} PUBLIC ${sasktran_win_links})
	target_link_libraries(sasktran_tests PUBLIC ${sasktran_win_links})
elseif(LINUX)
	target_link_libraries(${SASKTRAN_CORE_NAME} PUBLIC ${sasktran_unix_links})
	target_link_libraries(sasktran_tests PUBLIC ${sasktran_unix_links})
elseif(APPLE)
	target_link_libraries(${SASKTRAN_CORE_NAME} PUBLIC ${sasktran_darwin_links})
	target_link_libraries(sasktran_tests PUBLIC ${sasktran_darwin_links})
endif()

# Some special logic to support both CATCH2 version 2 and 3 for now
find_package(Catch2 2 CONFIG QUIET)

if(Catch2_FOUND)
	# Found Catch2 Version 2
	target_link_libraries(sasktran_tests PUBLIC Catch2::Catch2)
else()
	# Look for Catch2 Version 3
	message("Could not find Catch2 Version 2, Looking for Version 3")
	find_package(Catch2 3 CONFIG REQUIRED)
	target_link_libraries(sasktran_tests PUBLIC Catch2::Catch2WithMain)
	target_compile_definitions(sasktran_tests PUBLIC SKTRAN_CATCH2_VERSION3)
endif()


# Build wheels for each python environment
if(TRUE)
	if(DEFINED PythonEnvironments)
		message(STATUS "Building wheels for ${PythonEnvironments}")

		add_custom_command(TARGET ${SASKTRAN_CORE_NAME}
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${SASKTRAN_CORE_NAME}> ${CMAKE_SOURCE_DIR}/src/core/sasktran/python_sasktrancore/sasktran_core/
		)

		get_target_property(YAML_LIBRARY_DIR yaml-cpp LOCATION_RELEASE)
		get_target_property(YAML_INCLUDER_DIRS yaml-cpp INTERFACE_INCLUDE_DIRECTORIES)


		# This is the proper way to do it but for some reason this segfaults in the linux build container
		# So instead we allow passing of pure files to the setup.py and remove it there
		# cmake_path(REMOVE_FILENAME YAML_LIBRARY_DIR)

		set(SASKTRAN_IF_INCLUDE_PATH "${Boost_INCLUDE_DIRS}\;${YAML_CPP_INCLUDE_DIR}\;${YAML_INCLUDER_DIRS} ")
		set(SASKTRAN_IF_LIBRARY_PATH "$<TARGET_FILE_DIR:${SASKTRANIF_NAME}>\;$<TARGET_FILE_DIR:${NXBASE_NAME}>\;${Boost_LIBRARY_DIRS}\;${YAML_LIBRARY_DIR};${YAML_CPP_LIBRARY_DIR} ")

		foreach(env ${PythonEnvironments})
			if(WIN32)
				add_custom_command(TARGET ${SASKTRAN_CORE_NAME}
					POST_BUILD
					COMMAND ./buildscripts/build_whl_sasktran_win.bat ${env} "${SASKTRAN_IF_INCLUDE_PATH}" "${SASKTRAN_IF_LIBRARY_PATH}"
					WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
				)
			else()
				add_custom_command(TARGET ${SASKTRAN_CORE_NAME}
					POST_BUILD
					COMMAND bash buildscripts/build_whl_sasktran_linux ${env} "${SASKTRAN_IF_INCLUDE_PATH}" "${SASKTRAN_IF_LIBRARY_PATH}" ${CMAKE_C_COMPILER} ${CMAKE_CXX_COMPILER}
					WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
				)
			endif()
		endforeach(env)
	endif()
endif()

