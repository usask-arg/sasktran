
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

/*---------------------------------------------------------------------------
 *		The NX library include file
 *
 *	Include this header file to get access to all of the NX library.
 *
 * @const MACRO | NXDEBUG							| Defined if compiling DEBUG versions of code. Possibly used to changes behaviour of trace and assert macros
 * @const MACRO | NXASSERT							| Used in Debug compilations to make various tests and assertions in the code
 * @const MACRO | NX_WINDOWS						| Defined if this is a Microsoft WIN32 system
 * @const MACRO | NX_UNIX_VER						| Defined if this is a Unix/ANSI compiler
 * @const MACRO | STL(x)							| Defines the namespace for STL library for various compilers.
 * @const MACRO | SUPPORTS_ANSI_STL					| Compiler supports ANSI STL, ie version using ansi STL filenames
 * @const MACRO | SUPPORT_ANSI_FUNCTIONS_ONLY		| Only support genuine ANSI C runtime library functions.
 *  												  This mostly effects nxfile.cxx since fstat, stat and
 *													  fileno are NOT ANSI functions.
 * @const MACRO | NX_INITGUID						| If defined then define the GUID's for various interfaces. eg IID_InxLog.
 *													  If required then it will only need defiing when comipling one source file
 *													  in a project.
 * @const MACRO  | DIRECTORY_CHAR					| The character used to delimit directories, '\' on Windows, '/' on Unix
 * @const MACRO	 | DIRECTORY_SEARCH_PATH_DELIMITER	| Charcters used to delimit a list of directories ';' on Window, ':' on Unix
 * @const MACRO  | N_ELEMENTS(x)					| Returns number of elements in a fixed array sizeof(x)/sizeof(x[0])
 * @const MACRO  | nxTRUE							| Boolean TRUE. Defines to "true" on systems that support it otherwise maps to "1"
 * @const MACRO  | nxFALSE							| Boolean FALSE. Defines to "false" on systems that support it otherwise maps to "1"
 * @type nxBOOL										| Implements Boolean variables. Maps to bool on systems that support it otherwise maps to int
 * @type nxBYTE										| 8 bit unsigned integer
 * @type nxCHAR										| 8 bit signed integer
 * @type nxWORD										| 16 bit unsigned integer
 * @type nxSHORT									| 16 bit signed integer
 * @type nxDWORD									| 32 bit unsigned integer
 * @type nxLONG										| 32 bit signed integer
 * @type nxHYPER									| 64 bit signed integer, if supported
 *
 *--------------------------------------------------------------------------- */


#if !defined(NXLIB_NXCORE_H)									/* Is this the first time through here */
#define NXLIB_NXCORE_H 1										/* Yip so define so we dont come through again */

#define     N_ELEMENTS(x)		(sizeof(x)/sizeof(x[0]))		/*!<  Calculate the number of elements in an array structure */

#if defined(__cplusplus)
	#define NX_START_STRUCT										/*!<  Structure definition, provides for differences between C and C++ compilers. On C it maps to typedef */
	#define NX_END_STRUCT(x)									/*!<  Structure definition, provides for differences between C and C++ compilers. */
#else															/*!defined _cplusplus */
	#define NX_START_STRUCT		typedef							/* Nope so define typedef struct <name> { ... }<name>; */
	#define NX_END_STRUCT(x)		x
#endif


/*---------------------------------------------------------------------------
 *					Microsoft Visual C++ compiler					2003-6-26
 *-------------------------------------------------------------------------*/


#if defined(_MSC_VER)										/* is this a Microsoft Visual C++ compiler	*/
	#define NX_WINDOWS	1									/*!<  1 if this is a WINDOWS system */

	/* ---- Check for the _DEBUG flag set by the compiler ---- */

	#if defined(_DEBUG)
		#if !defined(NXDEBUG)
			#define NXDEBUG									/*!< Defined if this is a DEBUG build */
		#endif
	#endif

	/* ---- check for the exception, bool and STL compatibility, */

	#if (_MSC_VER < 1200)						/* If this is earlier than Visual C++ 2.0  */
    	#error "This version of Visual C++ is too old"
	#elif (_MSC_VER < 1200)						/* if this is VC++ 5.0 */
		#define VC50COMPILER					/* then define we are using VC+ 5.0 */
		#define STL(x) std::x					/* Put the STL library in the std namespace */
		#define STD(x) std::x					/* define the std namespace macro */
		#define SUPPORTS_ANSI_STL				/* It uses ANSO STL filenames */
		#define intptr_t long
		#pragma warning( disable : 4786)		/* This is neededto suppress STL warning (about const const declarations) and also variable line length exceeds 255 chars. */
	#elif (_MSC_VER < 1300)						/* if this is VC++ 6.0 */
		#define VC60COMPILER
		#pragma warning( disable : 4786)		/* This is neededto suppress STL warning (about const const declarations) and also variable line length exceeds 255 chars. */
		#define STL(x) std::x					/* Put the STL library in the std namespace */
		#define STD(x) std::x					/* define the std namespace macro */
		#define intptr_t long

		#define SUPPORTS_ANSI_STL
		#include <crtdbg.h>						/* load in the definitions for the _ASSERTE macro */
	#else										/* if this is VC++ 6.0 */
		#define VC70COMPILER
		#define STL(x) std::x					/* Put the STL library in the std namespace */
		#define STD(x) std::x					/* define the std namespace macro */
		#define SUPPORTS_ANSI_STL
		#if !defined(_CRT_SECURE_NO_WARNINGS)
		#define _CRT_SECURE_NO_WARNINGS
		#endif
		#include <crtdbg.h>						/* load in the definitions for the _ASSERTE macro */



	#endif
	#if defined (__AFX_H__)
		#define NX_USING_AFX
	#endif

	#if ( _MSC_VER < 1900 )
		#define  NXFINITE(X)			(_finite((X)) != 0)
	#else
		#define  NXFINITE(X)			(std::isfinite((X)))
	#endif

	
	typedef		unsigned char		nxBYTE;		/* 8 bit unsigned variable */
	typedef     signed char			nxCHAR;		/* 8 bit signed  integer */
	typedef		wchar_t				WCHAR;		/* 16 bit Unicode character */
	typedef		unsigned short		nxWORD;		/* 16 bit unsigned variable */
	typedef     signed short		nxSHORT;	/* 16 bit signed integer */
	typedef		unsigned int		nxDWORD;	/* 32 bit unsigned variable */
	typedef		signed int			nxLONG;		/* 32 bit signed integer */
	#define		nxTRUE				true		/*!<  define a value for logical true. Its 1 on some systems */
	#define 	nxFALSE				false		/*!<  define a value for logical false. Its 0 on some systems */
	typedef		bool				nxBOOL;		/*!<  Define a value for logical variables. Its bool on some systems and int on others*/
	typedef		__int64				nxHYPER;	/*!<  64 bit signed integer */

/*---------------------------------------------------------------------------
 *					Borland C++ compiler					2003-6-26
 *-------------------------------------------------------------------------*/

#elif defined(__BORLANDC__)					/* is this a Borland C++ compiler */
	#define STL(x) std::x					/* Put the STL library in the appropriate namespace */
	#define STD(x) x						/* define the std namespace macro */
	#define NX_WINDOWS 1
	#if (__BORLANDC__ < 0x460)				/* is this before version 5.0 */
		#define BC40COMPILER				/* flag this as BC 4.0 */
	#else
		#define BC50COMPILER				/* otherwise we are using BC 5.0 orhigher */
		#define SUPPORTS_ANSI_STL
	#endif
	typedef		unsigned short		nxWORD;		/* 16 bit unsigned variable */
	typedef		unsigned int		nxDWORD;	/* 32 bit unsigned variable */
	typedef     signed char			nxCHAR;		/* 8 bit signed  integer */
	typedef     signed short		nxSHORT;	/* 16 bit signed integer */
	typedef		signed int			nxLONG;		/* 32 bit signed integer */
	#define		nxTRUE				true		/* as well as true and false */
	#define 	nxFALSE				false
	typedef		bool				nxBOOL;		/* Use bool as the definition of boolean either implicit on compiler or done by STL*/
	typedef		unsigned char		nxBYTE;		/* 8 bit unsigned variable */
	typedef		__int64				nxHYPER;

/*---------------------------------------------------------------------------
 *					GNU compiler (Assume its UNIX) C++ compiler		2003-6-26
 *-------------------------------------------------------------------------*/

#elif defined( __GNUG__ )							/* is this the GNU compiler */

	#if !(__GNUC__ > 2) && !(__GNUC__ ==2 && __GNUC_MINOR__ > 8)
	   #error This version of GNU g+= is too old (__VERSION__). You need 2.95 or higher
	#elif (__GNUC__ == 2)
		#define NX_UNIX_VER							/* assume this is unix (needs adjustment for DOS based GNU C++ DJGPP) */
		#define STL(x) x							/* put the Standard template library in default namespace */
		#define STD(x) x							/* put some of the standard functtions in the default namespace. */
		#define SUPPORTS_ANSI_STL					/* this compiler does use the ANSI filename notations */
		typedef long intptr_t;						/* Define the intptr_t datatype */
	#elif (__GNUC__ > 2)							/* g++ 3.x */
		#define NX_UNIX_VER							/* assume this is unix (needs adjustment for DOS based GNU C++ DJGPP) */
		#define STL(x) std::x						/* put the Standard template library in default namespace */
		#define STD(x) std::x						/* put some of the standard functtions in the default namespace. */
		#define SUPPORTS_ANSI_STL
		#if (__GNUC__ == 3 && __GNUC_MINOR__ < 3)
			typedef long intptr_t;
		#endif
	#endif
	#define		NXFINITE(X)			(std::isfinite((X)))

	typedef		unsigned char		nxBYTE;		/* 8 bit unsigned variable */
	typedef     signed char			nxCHAR;		/* 8 bit signed  integer */
	typedef     signed short		nxSHORT;	/* 16 bit signed integer */
	typedef		unsigned short		nxWORD;		/* 16 bit unsigned variable */
	typedef		wchar_t				WCHAR;		/* 16 bit Unicode character */
	typedef		unsigned int		nxDWORD;	/* 32 bit unsigned variable */
	typedef		signed int			nxLONG;		/* 32 bit signed integer */
	#define		nxTRUE				true		/* as well as true and false */
	#define 	nxFALSE				false
	typedef		bool				nxBOOL;		/* Use bool as the definition of boolean either implicit on compiler or done by STL*/
//	typedef		__int64				nxHYPER;	/* Does GNU support 64 bit ? */
	#define		_finite				finite
	#define		_isnan				!finite


/*---------------------------------------------------------------------------
 *					DEC Unix C++ compiler		2003-6-26
 *-------------------------------------------------------------------------*/

#elif defined(__DECCXX)							/* Is this the DEC C++ cxx compiler */
	#define NX_UNIX_VER							/* this is unix. */
	#define STL(x) x							/* put the Standard template library in default namespace */
	#define STD(x) x							/* put some of the standard functtions in the default namespace. */
	#define SUPPORTS_ANSI_STL					/* this compiler does use the ANSI filename notations */

	typedef		unsigned char		nxBYTE;		/* 8 bit unsigned variable */
	typedef     signed char			nxCHAR;		/* 8 bit signed  integer */
	typedef		wchar_t				WCHAR;		/* 16 bit Unicode character */
	typedef		unsigned short		nxWORD;		/* 16 bit unsigned variable */
	typedef     signed short		nxSHORT;	/* 16 bit signed integer */
	typedef		unsigned int		nxDWORD;	/* 32 bit unsigned variable */
	typedef		signed int			nxLONG;		/* 32 bit signed integer */
	typedef     signed char			nxCHAR;		/* 8 bit signed  integer */
//	typedef		__int64				nxHYPER;	/* Does GNU support 64 bit ? */
	typedef		int					nxBOOL;		/* Use bool as the definition of boolean either implicit on compiler or done by STL*/

/*---------------------------------------------------------------------------
 *					SGI Unix C++ compiler		2003-6-26
 *-------------------------------------------------------------------------*/
#elif defined(__sgi__) || defined (__sgi)		/* Is this the sgi C++ cxx compiler */

	#define		NX_UNIX_VER						/* this is unix. */
	#define		STL(x) std::x					/* put the Standard template library in default namespace */
	#define		STD(x) std::x					/* put some of the standard functtions in the default namespace. */
	#define		SUPPORTS_ANSI_STL				/* this compiler does use the ANSI filename notations */

	typedef		unsigned char		nxBYTE;		/* 8 bit unsigned variable */
	typedef     signed char			nxCHAR;		/* 8 bit signed  integer */
	typedef		wchar_t				WCHAR;		/* 16 bit Unicode character */
	typedef		unsigned short		nxWORD;		/* 16 bit unsigned variable */
	typedef     signed short		nxSHORT;	/* 16 bit signed integer */
	typedef		unsigned int		nxDWORD;	/* 32 bit unsigned variable */
	typedef		signed int			nxLONG;		/* 32 bit signed integer */
	typedef		bool				nxBOOL;		/* Use bool as the definition of boolean either implicit on compiler or done by STL*/
#else
	#error This C++ compiler has not yet been included in the nxlib system
#endif



/*------------------------------------------------------------------------ */
/*	Standard Includes from the C/C++ system								  */
/*------------------------------------------------------------------------ */

//#include <stdio.h>
//#include <stdarg.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>

#if defined(SUPPORTS_ANSI_THROW)														/*    if the compiler supports the C++ throw */
#define NXTHROW(x) throw x																/*!<  A macro for implementing the C++ throw statement. On some systems it just prints a message and stops the program.*/
#else																					/*    otherwise */
#define NXTHROW(x) {printf x; printf("\n***ABORTING PROGRAM ***\n"); abort();}			/*    simply print the message, program will probably bomb */
#endif


   //------------------------------------------------------------------------
   //  If your compiler chokes on these lines its not a C++ compiler.
   // You must undefine the predefined preprocessor macro __cplusplus
   // the C++ code from your compilation.
   //
   // Standard NX class libray includes
   //------------------------------------------------------------------------

#if defined(SUPPORTS_ANSI_STL)					//Does compiler support ANSI standard STL filenames
	#include <list>
	#include <limits>
	#include <vector>
	#include <set>
	#include <map>
	#include <iterator>
	#include <algorithm>
	#include <istream>
	#include <ostream>
	#include <fstream>
	#include <sstream>
//	#include <strstream>				// strstream is deprecated (May 2009). It generates warnings in the GNU C++ compiler

#else										// While older implementations use these files for STL
	#include <list.h>
	#define stlvector_h <vector.h>
    #include <multiset.h>
	#include <set.h>
	#include <map.h>
	#include <iterator.h>
	#include <algo.h>
#endif

/*-----------------------------------------------------------------------------------
 *			DIRECTORY_CHAR stuff
 *------------------------------------------------------------------------------------*/

#if defined(NX_WINDOWS)
    #include <io.h>
    #include <direct.h>
	#define DIRECTORY_CHAR		'\\'				/*!<  The directory seperator specifier '\' on windows, '/' on Unix */
	#define DIRECTORY_SEARCH_PATH_DELIMITER ';'		/*!<  The delimiter when speciffying multiple search paths in an environment variable ';' on windows, ':' on Unix */
	#include <ole2.h>								// Include OLE definitions etc. etc.
#endif

#if defined(NX_UNIX_VER) 
//	#define NX_USEBOOST_THREADLIBS						// 2012-01-22 All versions now use Boost. Unix versions must use Boost Multi-threading Libs
    #include /* */ <dirent.h>						// unix specifics
	#define DIRECTORY_CHAR		'/'					// Unix Directory specifier
	#define DIRECTORY_SEARCH_PATH_DELIMITER ':'		// delimiter when speciffying multiple search paths in an envirnment variable
	#define MAX_PATH 256							// max path name is 256 chars
	#include "module/system/unix/unixole2.h"		// Include Nick's Unix OLE definitions.
#endif

#include <sys/types.h>
#include <sys/stat.h>

#if !defined(nxmax)
template <class T> T nxmax( const T& x1, const T& x2)	{ return (x1 > x2)? x1 :x2;}	/*!< Get the maximum of two objects */
template <class T> T nxmin( const T& x1, const T& x2)	{ return (x1 < x2)? x1 :x2;}	/*!< Get the minimum of two objects */
#endif

/*-----------------------------------------------------------------------------------
 *			NX_INITGUID
 *------------------------------------------------------------------------------------*/

#if defined(NX_INITGUID)												/*!< \def NX_INITGUID Used when importing COM objects such as ONYX */
	#ifndef INITGUID
		#define INITGUID
	#endif
	#include "module/system/cominterface/nxinitguid.h"			// Is this a WIndows System.
#endif

#include <stdint.h>
#include <string.h>
#include "module/system/debug/nxtrace.h"
#include "module/system/nxbitmanip.h"
#include "module/science/physicalconstants/nxsiconstants.h"
#include "module/science/physicalconstants/nxcgsconstants.h"
#include "module/science/physicalconstants/refractiveindexdryairatstp.h"
#include "module/system/strings/nxstring.h"
#include "module/science/geodesy/timestmp.h"
#include "module/system/loggers/nxlog.h"
#include "module/system/fileio/nxfile.h"
#include "module/system/fileio/wildcard.h"
#include "module/system/fileio/dirscan.h"
#include "module/system/fileio/nxdir.h"
#include "module/system/fileio/nxfilestrucfragments.h"
#include "module/system/fileio/nxfilesessionname.h"
#include "module/system/automation/nxvarianttype.h"
#include "module/system/win32/nxregistrykey.h"


#ifdef NX_WINDOWS
  //#include "module/system/multithread/nxsyncobjects.h"		// 2013-01-17, Now include nxbase_threads.h after nxbase_core.h in code that needs it
  #include "module/system/automation/registerdll.h"
  //#include "module/system/automation/nxidispatchenum.h"
  //#include "module/system/os9/os9module.h"
  #include "module/system/win32/nxwin32.h"
  //#include "module/system/automation/variantstuff.h"
  //#include "module/math/arrays/nxsafearray.hpp"
#endif

//#include "module/system/multithread/nxworkerthread.h"			// 2013-01-17, Now include nxbase_threads.h nxworkerthread.h, added May 28, 2009, support for boost::threads, (developed with boost 1.39.0 )
//#include "module/system/getopt/nxgetopt.h"

#endif  //NXLIB_NXCORE_H
