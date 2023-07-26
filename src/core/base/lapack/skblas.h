
#ifndef SKLIBRARY_BLAS
#define SKLIBRARY_BLAS

namespace blas
{
#if true
	#define F77CALL 						// Method for calling Intel Fortran functions. Redefine with __pascal or __stdcall or whatever.
	#define F77EXTERNC extern "C"			// Disable Name Mangling for Intel Fortran
	#define BLAS(X) blas::X					// Macro to access Lapack namespace. Redefine if acessing Intel MKL
	typedef char							f77char;
	typedef int 							f77integer;
//	typedef unsigned long int 				f77uinteger;
	typedef short int 						f77shortint;
	typedef float 							f77real;
	typedef double 							f77doublereal;
	typedef struct { f77real r, i; } 		f77complex;
	typedef struct { f77doublereal r, i; } 	f77doublecomplex;
	typedef int 							f77logical;
	typedef short int 						f77shortlogical;
//	typedef char 							f77logical1;
//	typedef char 							f77integer1;
	typedef long long 						f77longint;
	typedef unsigned long long 				f77ulongint;
	#include "blasinterface_lowcase.h"		// Include the BLAS definitions that use lowercase Fortran library calls (with appended underscroe)

#else

	//#define F77CALL __stdcall				// Method for calling Lahey Fortran functions. Redefine if accessing INtel MKL
	#define F77CALL 						/* Calling Method for Fortran functions. Redefine if accessing Lahey Fortran */
	#define F77EXTERNC extern "C"			// Disable Name Mangling for calls to Fortran compiled functions
	#define BLAS(X) blas::X					// Macro to access Blas namespace.
	typedef char							f77char;
	typedef int 							f77integer;
//	typedef unsigned long int 				f77uinteger;
	typedef short int 						f77shortint;
	typedef float 							f77real;
	typedef double 							f77doublereal;
	typedef struct { f77real r, i; } 		f77complex;
	typedef struct { f77doublereal r, i; } 	f77doublecomplex;
	typedef int 							f77logical;
	typedef short int 						f77shortlogical;
//	typedef char 							f77logical1;
//	typedef char 							f77integer1;
	typedef long long 						f77longint;
	typedef unsigned long long 				f77ulongint;
	#include "blasinterface_upcase.h"		// Include the BLAS definitions that use Uppercase Fortran library calls (with appended underscroe)
#endif
};



/**
 *  \ingroup Lapack
    \page Page1 Outline
    \par Introduction
    This is an installation of the Lapack and BLAS linear algebra routines available at www.netlib.org.
    There are over 1400 various functions within this FORTRAN library and represent an industry standard. Other
	optimized implementations are available from Intel and AMD, as the Math Kernel Library for example.
	This implemetation is written in fortran and has been compiled with the Lahey lf95 fortran compiler
	and linked into the dll, lapack_win32.dll.  Function stubs for the Microsoft Visual C++ have been written and
	also included in the DLL. The lapack_win32 dll only depends upon the microsoft operating system and does not
	require any other Lahey support DLL's, a major bonus for redistribution.

	\par Data Types
	The fortran was compiled with the following assumptions:
	- Fortran INTEGERS are C++ int which are 32 bit signed integers.
	- Fortran REAL  are C++ float which are 32 bit floating point
	- Fortran DOUBLE PRECISION are C++ double, which are 64 bit floating point

	\par Lapack Include file
	You can get access to the Lapack and BLAS routines by include "lapack.h" into your code. Please
	note that the LAPACK and BLAS functions are declared in C++ namespace called lapack. This was done to avoid namespace conflicts especially
	with some of the type definitions.


	\par LAPACK macro
	A macro is defined , LAPACK(X), which will wrap the functiion call with the
	lapack namespace definition. The macro can be redefined for other Lapack and BLAS
	implementations and permits porting of code to these implementations with only modification to the header file.

	\par Example
	Here is an example that performs a least squares fit of a 2nd order polynomial to x,y data.
	The 2nd order polynomial is predefined with coeffiecients a,b and c and the x array array defined
	as the integeres 0 to 99. The y array is then computed.  The least squares array 100x3 matrix is defined
	relating the coeffiecients a,b,c to the signal y.  This requires computation of matrix A with 3 columns
	of x**2, x and 1.  This array is fed into the Lapack routine DGELS which can be used to solve least
	squares problems.  The fitted coefficients a,b,c are returned as the first 3 elements of y.
	computed.

	\code

	double a = 0.004;
	double b = 0.3;
	double c = 5;
	nx1dArray<double>	x;
	nx1dArray<double>	y;
	LapackMatrix		A;
	LapackMatrix		column;
	nx1dArray<double>	coeffs;
	bool				ok;

	x.Indgen(100);
	y = (a*x + b)*x + c;
	A.SetSize(100,3);
	A.Slice( -1, -1, 1,1, &column ) = x*x;
	A.Slice( -1, -1, 2,2, &column ) = x;
	A.Slice( -1, -1, 3,3, &column ) = 1;

	int M      = A.NumRows();
	int N      = A.NumColumns();
	int LDA    = M;
	int INFO   = -1;
	int LWORK  = 100*M;
	char TRANS = 'N';
	int NRHS   = 1;
	int LDB    = M;
	nx1dArray<double>	work;

	work.SetSize( LWORK );
	LAPACK(DGELS) ( &TRANS, &M, &N, &NRHS, A.UnsafeArrayBasePtr(), &LDA, y.UnsafeArrayBasePtr(), &LDB, work.UnsafeArrayBasePtr(), &LWORK, &INFO );
	for (int i = 1; i <= 3; i++)
	{
		printf("%d %f\n", (int)i, (double)y.At(i));
	}
	\endcode
*/


//namespace lapack
//{

	/* The following are typedefs defined by f2c. I have placed them into the
	   lapack namespace to avoid corrup[ting the main namespace. However
	   users should be aware that they use symbols that other programmers may
	   have #define'd elsewhere in other includes.

	   Nick Lloyd July 19, 2007
	 */

/*
typedef int integer;
typedef unsigned int uinteger;
typedef char *address;
typedef short  shortint;
typedef float real;
typedef double doublereal;
struct complex{ real r, i; } ;
struct doublecomplex { doublereal r, i; } ;
typedef int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;


typedef int // Unknown procedure type  (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
//typedef // Complex  VOID (*C_fp)(...);
//typedef // Double Complex  VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
//typedef // Character  VOID (*H_fp)(...);
//typedef // Subroutine  int (*S_fp)(...);

// E_fp is for real functions when -R is not specified
//typedef VOID C_f;	// complex function
//typedef VOID H_f;	// character function
//typedef VOID Z_f;	// double complex function
typedef doublereal E_f;	// real function with -R not specified
#include "msvclapackinterface.h"
#include "msvcspecials.h"
};
*/


#endif
