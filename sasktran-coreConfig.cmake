include( "${CMAKE_CURRENT_LIST_DIR}/sasktran-coreTargets.cmake")

if(WIN32)
    get_target_property(SASKTRAN_LIB_DIRECTORY sasktran-core::sktran_base LOCATION_RELEASE)
    cmake_path(REMOVE_FILENAME SASKTRAN_LIB_DIRECTORY)
    set(SASKTRAN_CORE_LIBRARIES
            sasktran-core::sktran_stubs
            sasktran-core::sktran_base
            sasktran-core::skopticalproperties21
            sasktran-core::skclimatology21
            sasktran-core::sasktranif
            sasktran-core::nxbase
            sasktran-core::nxnetcdfio
            msis90e
            hitran_tips
            )
else()
    set(SASKTRAN_CORE_LIBRARIES
            sasktran-core::sktran_stubs
            sasktran-core::sktran_base
            sasktran-core::skopticalproperties21
            sasktran-core::skclimatology21
            sasktran-core::sasktranif
            sasktran-core::nxbase
            sasktran-core::nxnetcdfio
            )
endif()