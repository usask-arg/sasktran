# sasktranif source code

set(sasktranif_src
    sources/climatology_constants.cxx 
    sources/iskclimatology.cpp
    sources/iskemission.cpp
    sources/iskgeodetic.cpp 
    sources/iskoptprop.cpp
    sources/iskbrdf.cpp
    sources/sasktranif.cpp
    sources/stokesvector.cpp 
    sources/sasktranif_classfactory.cpp
        ../sasktran/modules/sktran_common/miscellaneous/iskmodule_base.cpp

         "sources/iskmie.cpp")

# Build the library
add_library(${SASKTRANIF_NAME} ${sasktranif_src})

set(sasktranif_headers
        includes/climatology_handles.h
        includes/sasktran_interfaces.h
        includes/sasktran_polarization.h
        includes/sasktranif.h
        includes/sasktranif_internals.h
        )

install (TARGETS ${SASKTRANIF_NAME}
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${SASKTRANIF_NAME})

install(DIRECTORY "includes"
        DESTINATION ${INCLUDE_INSTALL_DIR}/${SASKTRANIF_NAME}
        FILES_MATCHING
        PATTERN "*.h"
        )
