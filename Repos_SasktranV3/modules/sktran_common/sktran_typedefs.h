

/*-----------------------------------------------------------------------------
 *					typedef SKTRAN_GridIndexShort								*/
/** A typedef used to index grids. This type is synonomous with \c unsigned \c short.
 *	It was introduced to reduce memory storage requirements when storing
 *	indexes into the grids. If we start to use table indexes that require more
 *	than 65536 elements then change this type to unsigned int or size_t
**/
/*---------------------------------------------------------------------------*/

#if (SKTRANV2_MINIMIZE_MEMORY)
	typedef float			SKTRAN_GridIndexFloat;
	typedef unsigned short	SKTRAN_GridIndexShort;
	typedef unsigned char	SKTRAN_GridIndexByte;
	#define					SKTRAN_MAX_GridIndexByte	256
	#define					SKTRAN_MAX_GridIndexShort	65536
	#define					SKTRAN_DBL_TO_GridIndexFloat(x)	(SKTRAN_GridIndexFloat)((x)) 

#else
	typedef double			SKTRAN_GridIndexFloat;
	typedef unsigned short	SKTRAN_GridIndexShort;
	typedef unsigned short	SKTRAN_GridIndexByte;
	#define					SKTRAN_MAX_GridIndexByte	65536
	#define					SKTRAN_MAX_GridIndexShort	65536
	#define					SKTRAN_DBL_TO_GridIndexFloat(x)	(x) 
#endif

/*-----------------------------------------------------------------------------
 *					typedef SKTRAN_StokesScalar								*/
/** A typedef used to values which are radiances. This will allow future
 *	versions of SASKTRAN to identify where polarization adjustments need to be
 *	made.
**/
/*---------------------------------------------------------------------------*/

#if (SKTRANV2_MINIMIZE_MEMORY)
	typedef float	SKTRAN_StokesScalar;
	#define			SKTRAN_DBL_TO_STOKES_SCALAR(x)	(SKTRAN_StokesScalar)((x)) 
#else
	typedef double	SKTRAN_StokesScalar;
	#define			SKTRAN_DBL_TO_STOKES_SCALAR(x)	 (x) 
#endif

/*-----------------------------------------------------------------------------
 *					typedef SKTRAN_PhaseMatrixScalar								*/
/** A typedef used to values which are Phase Matrices used for scattering
 *	calculations. This will allow future versions of SASKTRAN to identify
 *	where polarization adjustments need to be made.
**/
/*---------------------------------------------------------------------------*/

typedef SKTRAN_StokesScalar SKTRAN_PhaseMatrixScalar;

/*-----------------------------------------------------------------------------
 *					SKTRAN_Distance						2008-1-29*/
/** A class used to represent the intercept of a ray with a ray tracing shell
 *	I only store the distance of the intercept from the ray origin. There are
 *	a relatively large number of ray intercepts, typically 1-5 million, so we
 *	don't want to store too much information.
 **/
/*---------------------------------------------------------------------------*/

#if (SKTRANV2_MINIMIZE_MEMORY)
	typedef float	SKTRAN_Distance;														// This MACRO is used in V3 but other than that they are identical
	#define			SKTRAN_DBL_TO_DISTANCE(x)	(SKTRAN_Distance)((x)) 
#else
	typedef double	SKTRAN_Distance;
	#define			SKTRAN_DBL_TO_DISTANCE(x)	 (x) 
#endif



/*-----------------------------------------------------------------------------
 *					enum ENUM_SKTRAN_JSOURCE		2008-2-5*/
/** Specifies the different type of source functions and radiances
 *	encountered in the SASKTRAN model. Used by the code to clarify what
 *	function is doing what.
 **/
/*---------------------------------------------------------------------------*/

enum	ENUM_SKTRAN_JSOURCE { SKTRAN_JUNDEFINED,				//!< Reserved to place enumeration in "known" undefined state.
							  SKTRAN_JDIFFUSE,					//!< Used by rays to calculate MS source function along rays
							  SKTRAN_JSINGLESCATTER,			//!< Used by rays to calculate single scatter terms along ray
							  SKTRAN_JGROUNDMS,					//!< Used by rays and albedo to get the multiple scatter from the ground
							  SKTRAN_JGROUNDSS,					//!< Used by rays and albedo to get the single scatter from the ground
							  SKTRAN_JGROUND_UPFLUX,			//!< Used by ground points to calculate upward flux for albedo scattering
							  SKTRAN_JEMISSION };				//!< Reserved for future use, used by rays to get thermal emission along ray



