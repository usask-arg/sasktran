

/*--------------------------------------------------------------------------
 *					nxcrc_ccitt									2002-9-30*/
/**	\ingroup Math_Functions
 *	Calculates the 32 bit CRC for a chunk of memory
 *
 *	\param a_ptr
 *		Points to the start of the memory.
 *
 *	\param bytelength
 *		The length of the chunk of memory in bytes			
 *
 *	\return
 *		the 32 bit CRC
**/
 /*------------------------------------------------------------------------*/

nxDWORD  nxcrc_ccitt( void *a_ptr, nxDWORD bytelength);

/*--------------------------------------------------------------------------
 *					nxadditional_crc_ccitt						2002-9-30*/
/** \ingroup Math_Functions
 *	Calculates the 32 bit CRC for a partial chunk of memory.  The CRC of
 *	several chunks of memory can be evaluated by calling this function
 *	several times.
 *
 *	\param oldcrc32
 *		The current grand sum CRC of all chunks of memory checked.  Set to
 *		0 if this is the start of a new  calculation
 *			
 *	\param a_ptr
 *		Points to the start of the memory.
 *			
 *	\param bytelength
 *		The length of the chunk of memory in bytes			
 *
 *	\return 
 *		the new 32 bit CRC
**/
/*------------------------------------------------------------------------*/

nxDWORD  nxadditional_crc_ccitt( nxDWORD oldcrc32, void *a_ptr, nxDWORD bytelength);
