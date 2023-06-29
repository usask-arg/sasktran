# nxhdf source code

set(nxhdfeos5_src
#    nxhdfeos5/hdfeos5reader.cxx
#    nxhdfeos5/nxhdfeos5sw_base.cxx
)
set(nxhdfeos_src
#    nxhdfeos/nxhdfeospt.cxx 
#    nxhdfeos/nxhdfeospt_base.cxx 
#    nxhdfeos/nxhdfeossw_base.cxx
)
set(nxnetcdfio_src
    nxnetcdfio/netcdfio.cpp
)

# Build the library
add_library(${NXNETCDFIO_NAME} ${nxnetcdfio_src})

install(DIRECTORY "nxnetcdfio/"
        DESTINATION ${INCLUDE_INSTALL_DIR}/${NXNETCDFIO_NAME}
        FILES_MATCHING
        PATTERN "*.h"
        )

install (TARGETS ${NXNETCDFIO_NAME}
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${NXNETCDFIO_NAME})
