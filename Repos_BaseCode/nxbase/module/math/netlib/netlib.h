#ifndef FILE_NXBASE_NETLIB_H
#define FILE_NXBASE_NETLIB_H

/*---------------------------------------------------------------------------
 *					Namespace nxnetlib								2003-10-10*/
/** \ingroup Math
 *	\namespace nxnetlib
 *	This is where functions from the vast array of math functions stored at
 *	netlib.org are stored in the nxbase library. The functions are compiled
 *	with a Fortran compiler (typically LAHEY Fortran) into a netlib.dll
 *	and corresponding library netlib.lib
 **/
 /*-------------------------------------------------------------------------*/

namespace nxnetlib
{
/** \ingroup Math_netlib **/
	extern double(* dlgama)( double *x );		//!< The Log(Gamma) Function

};

#endif


