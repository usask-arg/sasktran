#if !defined(NXVARIANTTYPE_H)

#define NXVARIANTTYPE_H

/*--------------------------------------------------------------------------
 *					nxVariantType											*/
/**	\ingroup system_com
 *	A function to determine the Variant type of different data types. Used by
 *	template overloaded functions to get the correct VARTYPE of different
 *	data types.
 *
 *	\param argument
 *			The object for which the variant type is required. It can be
 *			double*, float*, int*, long*, short*, char*, nxWORD*, nxDWORD*, nxBYTE*
 *			
 *	\return 
 *		the variant data type.
**/
/*------------------------------------------------------------------------*/

inline VARTYPE nxVariantType( double*  /*argument*/  ) { return VT_R8;}
inline VARTYPE nxVariantType( float*   /*argument*/ ) { return VT_R4;}
inline VARTYPE nxVariantType( int*     /*argument*/ ) { return VT_I4;}
inline VARTYPE nxVariantType( long*    /*argument*/ ) { return VT_I4;}
inline VARTYPE nxVariantType( short*   /*argument*/ ) { return VT_I2;}
inline VARTYPE nxVariantType( char*    /*argument*/ ) { return VT_I1;}
inline VARTYPE nxVariantType( nxWORD*  /*argument*/ ) { return VT_UI2;}
inline VARTYPE nxVariantType( nxDWORD* /*argument*/) { return VT_UI4;}
inline VARTYPE nxVariantType( nxBYTE*  /*argument*/) { return VT_UI1;}

#endif

