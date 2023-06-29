# skopticalproperties source code

set(skopticalproperties_src
    sources/skabsorptiontable.cxx
	sources/skabsorptiontablepressure.cpp
    sources/skopticalpropertylistentry.cxx
    sources/skmiesphericalparticle.cxx
	sources/sknonsphericalparticle.cxx
	sources/sknonsphericalparticlediscretedipole.cxx
	sources/sknonsphericalparticletmatrixoriented.cxx
	sources/sknonsphericalparticletmatrixrandom.cxx
	sources/skparticlesize.cxx
	sources/skrefractiveindex.cxx
	sources/skrefractiveindex_ice.cxx
	sources/skrtextinction_tabulatedphasematrix.cxx
	sources/skrtextinctionaerosolprofile.cxx
	sources/skrthenyeygreenstein.cxx
    sources/skrtscattericecrystal.cxx
	sources/skrtscattericecrystalcached.cxx
	sources/skrtscatteringcrosssectionbase.cxx
	sources/skrtscattermieaerosol.cxx
	sources/skrtscattermieaerosolcached.cxx
	sources/skconvolvedopticalpropentry.cxx
	sources/skconvolvedopticalproperties.cxx
	sources/skrtextinction_convolved.cxx
	sources/skrtscatterrayleighdryair.cxx
    sources/skrtscatterrayleighdryairinelastic.cpp
	sources/skrtsimplerayleigh.cpp
	sources/skrtextinction_multipleoverlappingspectra.cxx
	sources/skphasematrix.cxx
	sources/skrtalbedo.cxx
	sources/skstokesmatrix.cxx
	sources/sktran_atmopshericstate.cxx
	sources/o3/sko3xsection.cxx
	sources/o3/sko3xsection_basspaur.cxx
	sources/o3/sko3xsection_basspaurtable.cxx
	sources/o3/sko3xsection_daumont_brion_malicet.cxx
	sources/o3/sko3xsection_ftsvoigt.cxx
	sources/o3/sko3xsection_osiris.cxx
	sources/o3/sko3xsection_sciabogumil.cxx
	sources/o3/sko3xsection_sciabogumilv4.cxx
	sources/o3/sko3xsection_serdyuchenko.cxx
	sources/no2/Vandaele_1998/vandaele_no21998.cpp
	sources/no2/skno2xsection.cxx
	sources/no2/skno2xsection_osiris.cxx
	sources/so2/vandaele_so2_2009.cpp
	sources/o2-o2/Hermans/o4_herman.cxx
	sources/o2-o2/Thalman_2013/o4_thalman2013.cxx
	sources/o2-o2/Hitran2016/o4_hitran2016.cxx
	sources/baum_bulk_icecloudmodel/baumicecrystals_database2014.cxx
	sources/baum_bulk_icecloudmodel/skrtscatter_tabulatedcirrusproperties.cxx
    sources/baum_bulk_icecloudmodel/skopticalproperties_baumicecrystals2014.cxx
	sources/skspectralline/hitran/hitranlinestruct.cxx
	sources/skspectralline/hitran/hitran_spectrallineio.cpp
	sources/skspectralline/hitran/hitran_xs_cache.cpp
	sources/skspectralline/hitran/hitran_partition_cache.cpp 
	sources/skspectralline/hitran/skhitranmoleculemanager.cxx
	sources/skspectralline/skspectrallinecollection_molecule.cxx
	sources/skspectralline/voigthumlicek_kuntz.cxx
	sources/skspectralline/voigthumlicek_wells.cxx
	sources/skspectralline/voigttabulated.cxx
	sources/solarspectrum/FontelaUVIS_3micron/fontela_uvis3micron.cpp
	sources/solarspectrum/FontelaUVIS_100micron/fontela_uvis100micron.cpp
	sources/solarspectrum/sksolarspectrum.cpp
	sources/solarspectrum/soa2010/solarspectrum_soa2010.cpp
	sources/emissions/hitran_upperstates.cpp
	sources/emissions/skemission.cpp
	sources/emissions/skemission_hitranchemical.cpp
	sources/emissions/spectralemission_collection.cpp
	sources/emissions/skemissiontabulatedheightwavelength.cpp
	sources/emissions/skemissionthermal.cpp
    sources/emissions/sktranatmosphericemission.cxx
	sources/sasktranif_opticalimpl/iskbrdf_stub_base.cpp
	sources/sasktranif_opticalimpl/iskemission_stub_base.cpp
	sources/sasktranif_opticalimpl/iskopticalproperty_stub_base.cpp
	sources/sasktranif_opticalimpl/isksolarspectrum_stub_base.cpp
	sources/brdf/sktran_brdf.cpp
	sources/brdf/sktran_brdf_combination.cpp
	sources/brdf/sktran_brdf_cox_munk.cpp
	sources/brdf/sktran_brdf_hapke.cpp
	sources/brdf/sktran_brdf_lambertian.cpp
	sources/brdf/sktran_brdf_li.cpp
	sources/brdf/sktran_brdf_modis.cpp
	sources/brdf/sktran_brdf_rahman.cpp
	sources/brdf/sktran_brdf_ross.cpp
	sources/brdf/sktran_brdf_roujean.cpp
	sources/brdf/sktran_brdf_snow_kokhanovsky2012.cpp
	sources/brdf/sktran_brdf_userdefinedlatlon.cpp
	sources/brdf/sktran_brdf_spectralvarying.cpp
 	sources/skuserscatterconstantheight.cpp
)

set(CMAKE_Fortran_FLAGS "-std=legacy -fPIC -fdefault-real-8 -fdefault-double-8 -fno-automatic -ffixed-line-length-132")


if(WIN32)
	set(skopticalproperties_platform_src)
else()
	set(skopticalproperties_platform_src 
		../fortran_libraries/wiscombemie/sources/MIEV0.f
		../fortran_libraries/wiscombemie/sources/ErrPack.f
		../fortran_libraries/wiscombemie/sources/RDI1MACHmodule.F90
		../fortran_libraries/hitran/bd_tips_2017/sources/BD_TIPS_2017_v1p0.f
		../fortran_libraries/nonspherescatter/TMatrixRandom/tmd.lp.f
	)
endif()

# Build the library
add_library(${SKOPTICALPROPERTIES_NAME} ${skopticalproperties_src} ${skopticalproperties_platform_src})

install(DIRECTORY "./"
		DESTINATION ${INCLUDE_INSTALL_DIR}/${SKOPTICALPROPERTIES_NAME}
		FILES_MATCHING
		PATTERN "*.h"
		)


install (TARGETS ${SKOPTICALPROPERTIES_NAME}
		EXPORT ${EXPORT_NAME}
		LIBRARY DESTINATION ${STATIC_LIB_INSTALL_DIR}
		INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}/${SKOPTICALPROPERTIES_NAME})
