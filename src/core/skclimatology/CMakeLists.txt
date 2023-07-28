# skclimatology source code

set(skclimatology_src
    sources/nxmsis90.cxx
    sources/skclimatology.cxx 
    sources/skclimatology_msis90.cxx
    sources/skclimatology_labow.cxx
    sources/skclimatology_one.cxx 
    sources/skclimatology_pratmo.cpp
    sources/skclimatology_textfiledbase.cxx 
    sources/skclimatology_userdefined_latlon_table.cpp
    sources/skclimatology_userdefinedplane.cpp
    sources/skclimatology_userdefinedtable.cxx 
    sources/skclimatology_usertablespline.cxx 
    sources/skclimatology_zero.cxx 
    sources/skclimatology_linearcombo.cxx
)

if(WIN32)
    set(skclimatology_platform_src)
else()
    set(skclimatology_platform_src
            ../fortran_libraries/msis90/sources/msise90.for
    )
endif()

set(CMAKE_Fortran_FLAGS "-std=legacy -fPIC -fdefault-real-8 -fdefault-double-8 -fno-automatic -ffixed-line-length-132")

# Build the library
add_library(${SKCLIMATOLOGY_NAME} ${skclimatology_src} ${skclimatology_platform_src})

target_link_libraries(${SKCLIMATOLOGY_NAME} PUBLIC sasktranCoreBuildProperties)

install (TARGETS ${SKCLIMATOLOGY_NAME}
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${SKCLIMATOLOGY_NAME})

install(DIRECTORY "."
        DESTINATION ${INCLUDE_INSTALL_DIR}/${SKCLIMATOLOGY_NAME}
        FILES_MATCHING
        PATTERN "*.h"
        )
