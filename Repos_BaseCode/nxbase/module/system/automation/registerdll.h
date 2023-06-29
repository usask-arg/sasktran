
/*--------------------------------------------------------------------------
 *					EncodeClassID											*/
/**	\ingroup system_com
 *	Encode a CLSID as a hexadecimal string representation as used by the
 *	Registry
 *
 *	\param clsid
 *		The CLSID to be converted to string representation
 *
 *	\return
 *		The string representation
**/
/*------------------------------------------------------------------------*/

extern nxString EncodeClassID( REFCLSID clsid );

/*--------------------------------------------------------------------------
 *					ConfigureCLSIDEntries									*/
/**	\ingroup system_com
 *	Configure the registry for the automation object of the given clsid 
 *
 *	\param clsid
 *		The CLSID of the automation object
 *			
 *	\param modulename
 *		The name of the module (normally a DLL name)
 *			
 *	\param classname
 *		The user friendly "text" representation of this clsid object (eg "excel.Spreadsheet")
 *			
 *	\param threadmodel
 *		Value of the ThreadingModel keyword, NULL, Apartment or Free
 *
 *	\return
 *		nxTRUE is success.
 *
 *	\par HISTORY
 *	2002-9-30
**/
/*------------------------------------------------------------------------*/

extern nxBOOL ConfigureCLSIDEntries( REFCLSID clsid, const char* modulename, const char *classname, const char* threadmodel );
