
#ifdef NX_NXFUNCS_H
   #error Multiple includes of nxfuncs.h, Only use nxlib.h!
#else

/*--------------------------------------------------------------------------
 *					nxStrtok												*/
/**	\ingroup system_strings
 * Tokenizes a string and returns the tokens in a nxStringArray object
 *
 *	\param utstring
 *		The string to be tokenized
 *			
 *	\param List
 *		The list of tokens to be returned
 *			
 *	\param separator
 *		The token separator string
 *			
 *
 *	\return
 *		the number of tokens generated
**/
/*------------------------------------------------------------------------*/
int nxStrtok( const char *utstring, nxStringArray *List, const char *separator = NULL );

#endif

