#undef NOMINMAX
#define NOMINMAX 1 // Prevent a macro definition from overriding the STL implementation of max, min

#pragma once
#if !defined(SKTRAN_COMMON_H_INCLUDED)
#define SKTRAN_COMMON_H_INCLUDED


#if !defined(SKTRANV2_MINIMIZE_MEMORY)			// if the option is not set on the compiler line then
	#define SKTRANV2_MINIMIZE_MEMORY  0			// True if we wish to minimize memory, eg. use floats for storage instead of doubles
#endif

#if defined(_MSC_VER) && defined(NXDEBUG)
#include <crtdbg.h>
#endif

#include "nxbase_core.h"
#include "nxbase_math.h"
#include "nxbase_geodesy.h"
#include <vector>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <skclimatology21.h>
#include <skopticalproperties21.h>

#if !defined(SKCLIMATOLOGY_VERSION)       || (SKCLIMATOLOGY_VERSION < 210)
#error The skclimatology version does not meet the requirements of sasktranv21. Please upgrade your version of skclimatology from SVN.
#endif

#if !defined(SKOPTICALPROPERTIES_VERSION) || (SKOPTICALPROPERTIES_VERSION < 210)
#error The skopticalproperties version does not meet the requirements of sasktranv21. Please upgrade your version of skopticalproperties from SVN.
#endif


#include "sktran_typedefs.h"
#include "opticalpropertytables/sktran_scattermatrix.h"
#include "grids/sktran_grid_definition.h"
#include "miscellaneous/codetimer.h"
#include "specifications/sktran_specs_base.h"
#include "specifications/sktran_lineofsightarray.h"
#include "specifications/raytracingregionmanager.h"
#include "specifications/sktran_userspecificationsclasses.h"
#include "unitspheres/sktran_unitsphere.h"
#include "straightrays/sktran_rayminimumcontainer.h"
#include "opticalpropertytables/sktran_polarizationprops_base.h"
#include "opticalpropertytables/sktran_opticalpropertiestable_pointcache.h"
#include "opticalpropertytables/sktran_opticalpropertiestable_base.h"
#include "quadrature/sktran_quadrature_tls_base.h"
#include "straightrays/sktran_ray_basev3.h"
#include "straightrays/sktran_raytracer_base.h"
#include "solartransmissiontables/sktran_sun.h"
#include "solartransmissiontables/sktran_source_term.h"
//#include "thermalemission/sktran_thermalemission.h"
//#include "miscellaneous/sktran_brdf.h"
#include "quadrature/sktran_optpropintegrator_base.h"
#include "quadrature/opticaldepthcalculator_linearwithheight.h"
#include "solartransmissiontables/sktran_solartransmission_base.h"
#include "straightrays/sktran_ray_straight.h"
#include "curvedrays/sktran_ray_curved.h"
//#include "straightrays/sktran_ray_shells.h"
#include "quadrature/sktran_optpropintegrator_straight.h"
#include "quadrature/sktran_srcintegrator_base.h"
#include "quadrature/sktran_srcintegrator.h"
//#include "quadrature/sktran_integrator_shells.h"
#include "quadrature/sktran_optpropintegrator_adaptive.h"
#include "straightrays/sktran_indexofrefraction.h"
#include "straightrays/sktran_raytracer_shells.h"
#include "straightrays/sktran_raytracer_straight_generic.h"
#include "straightrays/sktran_rayfactory_base.h"
#include "curvedrays/sktran_raytracer_curved_shells.h"
#include "diagnostics/sktran_diagnostics.h"
#include "engine/sktran_enginebase.h"
#include "opticalpropertytables/sktran_polarizationprops.h"
#include "opticalpropertytables/sktran_opticalproperties_1d_heightv3.h"
#include "opticalpropertytables/sktran_opticalproperties_3d_unitsphere.h"
#include "solartransmissiontables/sktran_solartransmission_2d.h"
#include "solartransmissiontables/sktran_solartransmission_3d.h"
#include "solartransmissiontables/sktran_solartransmission_notable.h"
#include "emissiontables/sktran_emissiontables.h"
//#include "modules/sktran_engine/include/sktran_engine.h"


#endif
