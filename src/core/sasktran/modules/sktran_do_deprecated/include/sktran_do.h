#pragma once

// STL dependencies 
#include <algorithm>
#include <numeric>
#include <mutex>
#include <queue>
#include <math.h>
#include <iostream>

// Boost dependencies
#include <boost/math/special_functions/legendre.hpp>

// Eigen settings

// Almost always better to thread over wavelength (or even layer) instead of using linear algebra parallelization
#define EIGEN_DONT_PARALLELIZE 1

// Whether or not to use the Eigensolver from Eigen instead of dgeev call directly
// I think the Eigen EigenSolver is faster for small matrices but slower for large ones, and sometimes the Eigen
// EigenSolver can have precision problems so right now we just disable it
#define SKTRAN_DO_USE_EIGEN_EIGENSOLVER false

// SKTRAN_DO is templated over two main parameters, the first is NSTOKES which we always explicitly instantiate
// every class over.  DO is also templated over the number of streams, the default is -1 which is "dynamic" and allows
// for any number of streams.  We also can specialize specific values of the number of streams but this takes a long
// time to compile so we typically only do it on release.  Useful values for speed are 2, 4, and 16

#define SKTRAN_DO_FULL_COMPILE

#ifdef SKTRAN_DO_FULL_COMPILE
    #define INSTANTIATE_TEMPLATE(classname) \
     template class classname<1>;\
     template class classname<3>;\
     template class classname<4>;\
     template class classname<1, 2>;\
     template class classname<3, 2>;\
     template class classname<4, 2>;\
     template class classname<1, 4>;\
     template class classname<3, 4>;\
     template class classname<4, 4>;\
     template class classname<1, 16>;\
     template class classname<3, 16>;\
     template class classname<4, 16>;
#else
    #define INSTANTIATE_TEMPLATE(classname) \
     template class classname<1>;\
     template class classname<3>;\
     template class classname<4>;
#endif


// Setup dependening on what linear algebra package is being linked
#ifdef SKTRAN_DO_USE_MKL
    // Using MKL for linear algebdra
    #define EIGEN_USE_MKL_ALL 1

    #include <mkl_lapacke.h>
#else
    // Unsure if this is faster or not
    #define EIGEN_USE_BLAS 1
    #ifdef SKTRAN_DO_USE_ACCELERATE
        // Using apple Accelerate for linear algebra, which doesn't have a LAPACKE interface

        #define lapack_int int
        #include <clapack.h>
        #include <cblas.h>
    #else
        // Using a standard LAPACKE compatible package

        #define lapack_complex_float std::complex<float>
        #define lapack_complex_double std::complex<double>
        #include <lapacke.h>
        #include <cblas.h>
    #endif
#endif

#include <Eigen/Dense>

// SASKTRAN dependencies
#include <modules/sktran_common/sktran_common.h>

// SASKTRAN-DO internal dependencies
#include "modules/sktran_do_deprecated/include/sktran_do_linearization_types.h"
#include "modules/sktran_do_deprecated/include/sktran_do_polarization_types.h"
#include "modules/sktran_do_deprecated/include/sktran_do_types.h"
#include "modules/sktran_do_deprecated/include/sktran_do_linalg.h"
#include "modules/sktran_do_deprecated/include/sktran_do_memory.h"
#include "modules/sktran_do_deprecated/include/sktran_do_quadrature.h"
#include "modules/sktran_do_deprecated/include/sktran_do_specs.h"
#include "modules/sktran_do_deprecated/include/sktran_do_opticalstate.h"
#include "modules/sktran_do_deprecated/include/sktran_do_spherical.h"
#include "modules/sktran_do_deprecated/include/sktran_do_surface.h"
#include "modules/sktran_do_deprecated/include/sktran_do_testing.h"
#include "modules/sktran_do_deprecated/include/sktran_do_lazyazimuth.h"
#include "modules/sktran_do_deprecated/include/sktran_do_properties.h"
#include "modules/sktran_do_deprecated/include/sktran_do_pconfig.h"
#include "modules/sktran_do_deprecated/include/sktran_do_opticallayer.h"
#include "modules/sktran_do_deprecated/include/sktran_do_geometrylayerarray.h"
#include "modules/sktran_do_deprecated/include/sktran_do_layerarray.h"
#include "modules/sktran_do_deprecated/include/sktran_do_rte.h"
#include "modules/sktran_do_deprecated/include/sktran_do_postprocessing.h"
#include "modules/sktran_do_deprecated/include/sktran_do_engine.h"

