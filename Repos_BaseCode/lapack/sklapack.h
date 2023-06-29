#ifndef SKLIBRARY_LAPACK
#define SKLIBRARY_LAPACK

#include "skblas.h"
using namespace blas;

namespace lapack
{
typedef f77logical (*L_fp)(...);
#if defined (NX_UNIX_VER)
	#define LAPACK(X) lapack::X					// Macro to access Lapack namespace. Redefine if acessing Intel MKL
	#include "lapackinterface_lowcase.h"
#else
	#define LAPACK(X) lapack::X					// Macro to access Lapack namespace. Redefine if acessing Intel MKL
	#include "lapackinterface_upcase.h"
#endif
};


#endif
