# nxbase source code

set(nxbase_core_src
    module/system/strings/nxfuncs.cxx
    module/system/strings/nxstring.cxx
    module/system/strings/nxstringw.cxx
    module/system/win32/nxregistryconfig.cxx
    module/system/win32/yaml_registrykey.cxx
    module/system/multithread/nxworkerthread.cxx
    module/system/loggers/nxlog.cxx
    module/system/loggers/nxlogbase.cxx
    module/system/loggers/nxlogcon.cxx
#    module/system/getopt/nxgetopt.cxx
    module/system/fileio/nxbinfil.cxx
    module/system/fileio/nxdir.cxx
    module/system/fileio/nxfile.cxx
    module/system/fileio/nxfilesessionname.cxx
    module/system/fileio/nxfilesp.cxx
    module/system/fileio/nxfilestrucfragments.cxx
    module/system/fileio/wildcard.cxx
    module/system/debug/nxtrace.cxx
    module/system/cominterface/nxunknown.cxx
    module/system/cominterface/nxlibCOMIF.cxx
    module/science/physicalconstants/nxcgsconstants.cxx
    module/science/physicalconstants/nxsiconstants.cxx
    module/science/physicalconstants/refractiveindexdryairatstp.cpp
    module/sktran_core/sktran_grid_definition.cpp
    module/science/geodesy/timestmp.cxx
    module/science/geodesy/timestmp_consts.cxx
    module/CommonUnits/commonphysicalunits.cxx
)

set(nxbase_geodesy_src
    module/science/geodesy/keplerorbit.cxx
    module/science/geodesy/moon.cxx
    module/science/geodesy/nxgeodetic.cxx
    module/science/geodesy/nxgeodetic_shellheight.cxx
    module/science/geodesy/planets.cxx
    module/science/geodesy/sun.cxx
    module/science/geodesy/sgp8.cxx
    module/science/geodesy/tangentpoint.cxx
    module/science/meteorology/nxmeteorology.cpp
    module/science/novas/novas.c 
    module/science/novas/novascon.c 
    module/science/novas/readeph0.c 
    module/science/novas/solsys2.c 
)

if (WIN32)
    set(nxbase_platform_src
        module/system/win32/nxregistryconfig.cxx
        module/system/win32/nxregistrykey.cxx 
        module/system/win32/nxregistrykeyunix.cxx 
        module/system/win32/nxwin32_helperfuncs.cxx 
        module/system/win32/yaml_registrykey.cxx
    )
else()
    set(nxbase_platform_src
    )
endif(WIN32)

set(nxbase_math_src
    module/math/crc16.cxx
    module/math/nxconicellipse.cxx
    module/math/nxcrc_ccitt.cxx
    module/math/nxgaussquadrature.cxx
    module/math/nxlinearinterpolate.cxx
    module/math/nxlinearinterpolate2.cxx
    module/math/nxmath.cxx
    module/math/nxpiecewiselinear.cxx
    module/math/nxquaternion.cxx
    module/math/nxspline.cxx
    module/math/nxspline2.cxx
    module/math/nxvector.cxx
    module/math/nxvector2d.cxx
    module/math/nxbase36.cxx
    module/math/arrays/nxdblarr.cxx
    module/math/arrays/nxarray_imp.cxx
    module/math/arrays/rankspecification.cxx
    module/math/data_analysis/lsqanalytic.cxx
    module/math/data_analysis/nxfpitophatfit.cxx
    module/math/data_analysis/nxrunningstatistic.cxx
    module/math/data_analysis/runningstraightlinefit.cxx
)


set(nxbase_other_src
    module/system/unix/unix_com_clsid.cxx
    module/system/unix/unix_com_dllentry.cxx
    module/system/unix/unix_com_emmulate.cxx
    module/system/unix/unix_com_registrykey.cxx
    module/system/automation/registerdll.cxx
)


# Build the library
add_library(${NXBASE_NAME} ${nxbase_core_src}
        ${nxbase_platform_src}
        ${nxbase_other_src}
        ${nxbase_geodesy_src}
        ${nxbase_math_src}
        )

target_link_libraries(${NXBASE_NAME} PUBLIC sasktranCoreBuildProperties)

install(DIRECTORY "."
        DESTINATION ${INCLUDE_INSTALL_DIR}/${NXBASE_NAME}
        FILES_MATCHING
        PATTERN "*.h"
        )
install(DIRECTORY "."
        DESTINATION ${INCLUDE_INSTALL_DIR}/${NXBASE_NAME}
        FILES_MATCHING
        PATTERN "*.hpp"
        )


install (TARGETS ${NXBASE_NAME}
        EXPORT ${EXPORT_NAME}
        LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
        INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${NXBASE_NAME})
