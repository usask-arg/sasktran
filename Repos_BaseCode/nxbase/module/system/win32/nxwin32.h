

/*----------------------------------------------------------------------------
 * 					nxGetClsid												*/
/** \ingroup system_com
 *	Gets the CLSID associated with a string based program name (eg "excel.Spreadsheet")
 *
 *	\param cls
 *		Returns the cclsid if successful
 *
 *	\param pgmname
 *		The name of the program to find, eg. "excel.Spreadsheet"
 *
 *	\return
 *		Returns nxTRUE if successful
**/
/*---------------------------------------------------------------------------*/

extern  nxBOOL nxGetClsid( CLSID* cls, const char *pgmname );

/*--------------------------------------------------------------------------
 *					YieldToSystem											*/
/**	\ingroup System_win32
 *	Yields control to the system.  It will empty the message queue at least
 *	N times before returning.
 *
 *	\param N
 *		Number of times to completely empty the message queue before returning
**/
/*------------------------------------------------------------------------*/

extern  void  YieldToSystem( int N );


/*--------------------------------------------------------------------------
 *					Wait													*/
/** \ingroup System_win32
 * Yields control to the system for at least N milliseconds.
 * It may wait longer as it will be processing the message queue while waiting.
 * Designed to provide message processing in the middle of long tasks and/or simple
 * throttling of CPU usage by a single thread. Not designed to provide accurate timing.
 *
 *	\param n
 *		Number of milliseconds to wait
**/
/*------------------------------------------------------------------------*/

extern  void  Wait(int n);



