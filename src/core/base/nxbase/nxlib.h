
/* ------------------- start of System documentation ---------------------*/

/** \defgroup System
    Documentation System
*/

/** \defgroup system_multithread	Multi-threading
	\ingroup System
    Multithreading functions
*/

/** \defgroup System_win32	Win32 classes
	\ingroup System
    A few Win32 helper functions
*/

/** \defgroup nxwin32registry Registry Class
	\ingroup System
	Registry class and functions. Most people will use #nxRegistryConfiguration to save and read
	configuration information for a specific application
**/

/** \defgroup system_os9	OS9 helper classes
	\ingroup System
    OS9 Helper classes
*/

/** \defgroup system_loggers	Logging Classes
	\ingroup System
    Logging classes
*/

/** \defgroup system_fileio		File I/O
	\ingroup System
    File I/O classes and functions
*/

/** \defgroup system_strings		Strings
	\ingroup System
    String Classes and functions
*/

/** \defgroup system_getopt		getopt
	\ingroup System
   	 See #nxGetOpt for a description of usage.
*/

/** \defgroup system_debug		Debug functions
	\ingroup System
   	Debugging classes and functions
*/

/** \defgroup system_com		COM and Automation
	\ingroup System
   	COM and Automation classes
*/

/** \defgroup system_onyx		Onyx helper functions
	\ingroup System
   	Onyx helper classes and functions
*/
/* ------------------- start of Dougs Common Functions ---------------------*/

/**   \defgroup Math_common Dougs Common Functions
   	  Documentation for Dougs common functions
*/

/* ------------------- start of Math documentation ---------------------*/

/**   \defgroup Math
   	  Documentation for Math
*/


/**   \defgroup Math_netlib Netlib
	  \ingroup Math
	  Functions from the extensive netlib.org list of fortran routines
*/

/**   \defgroup Math_fitting Curve Fitting
	  \ingroup Math
	  Classes for fitting curves
	  - #nxConicEllipse: Fits an ellipse to a curve
	  - #nxCubic : finds the real roots of a cubic equation
	  - #nxSpline : fits a spline to data points, interpolates and integrates
*/

/**   \defgroup MATH_functions General Functions
	  \ingroup Math
	  The nxmath namespace holds a range of general purpose functions
*/

/**   \defgroup Math_Matrix  Matrices 
	  \ingroup Math
	  Classes for matrix mathematics
*/
/**   \defgroup Math_Matrix_deprecated Deprecated
	  \ingroup Math_Matrix
	  Deprecated matrix mathematics
*/

/**   \defgroup nxLinearArray Arrays
	  \ingroup Math
	  
*/
/**   \defgroup nxLinearArray_Algorithms Algorithms
	  \ingroup nxLinearArray
	  Array algorithms
*/


/**   \defgroup nxLinearArray_deprecated Deprecated
	  \ingroup nxLinearArray
	  Deprecated arrays
*/

/**   \defgroup nxMath_Quaternion Quaternions
	  \ingroup Math
	  Quaternion Class
*/
/**   \defgroup nxMath_Vector nxVector
	  \ingroup Math
	  Vector Class
*/

/**   \defgroup  Math_quadrature Quadrature 
	  \ingroup Math
	  Quadrature classes
   	  
*/
/**   \defgroup Math_quadrature_internals Internals
	  \ingroup Math_quadrature
	  Quadrature internals
   	  
*/

/* ------------------- start of Science documentation ---------------------*/

/**  \defgroup Science
   Documentation for Science
*/

/**   \defgroup Geodesy   Space, Time and Geoids
	  \ingroup Science
	  Classes for managing space time and geodesy
*/

#include "nxbase.h"
#include "nxbase_math.h"
#include "nxbase_geodesy.h"
#if defined(NX_WINDOWS)
#include "dadcommon.h"
#endif


