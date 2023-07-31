#ifndef FILE_NXBASE_LINEARALGEBRA_H
#define FILE_NXBASE_LINEARALGEBRA_H
//#pragma once

#include "nxbase_core.h"
#include "nxbase_math.h"


#if defined(_MSC_VER)
															// Include the intel math kernel library
 #if !defined(NXBASE_LINEARALGEBRA_LIB)
	#pragma message("Including the nxbase_linearalgebra libraries into this project")

	#if defined (_MANAGED)									// Is Visual C++ generating MSIL with /CLR
		#if defined (NXDEBUG)
			#pragma comment( lib, "nxbase_linearalgebradclr")
		#else
			#pragma comment( lib, "nxbase_linearalgebrarclr")
		#endif
	#else
		#if defined (NXDEBUG)
			#pragma comment( lib, "nxbase_linearalgebra_Debug")
		#else
			#pragma comment( lib, "nxbase_linearalgebra_Release")
		#endif
		//#pragma comment( lib, "blaslapack")
	#endif
	#pragma message("Including the Fortran LAPACK libraries (nxbase_linearalgebra.lib and Fortran dll library blaslapack) into this project")


  #else
	#pragma message("Link projects with nxbase_linearalgebra.lib and Fortran library LapackIntelFortranLib into this project")
  #endif
#endif

#endif //FILE_NXBASE_MATH_H


