#include "nxbase_core.h"

//---------------------------------------------------------------------------
//						DecodeHex
//	quicky helper function for DecodeClassIDstring
//---------------------------------------------------------------------------

static unsigned long DecodeHex( const char * str, int l1, int l2 )
{
	char		buffer[20];
	char*		endptr;
	int			n;
	const char*	bufptr = str + l1;

	n = l2-l1 + 1;
	for (int i=0; i < n; i++) buffer[i] = *bufptr++;
	buffer[n] = '\0';
	return strtoul( buffer, &endptr, 16 );
}
	

//---------------------------------------------------------------------------
//						DecodeClassIDString
//---------------------------------------------------------------------------

CLSID DecodeClassIDString( const char * str )
{
	CLSID	clsid;

	if (strlen(str) == 38)
	{
		clsid.Data1    = (unsigned long) DecodeHex( str, 1, 8);
		clsid.Data2    = (unsigned short)DecodeHex( str,10,13);
		clsid.Data3    = (unsigned short)DecodeHex( str,15,18);
		clsid.Data4[0] = (unsigned char) DecodeHex( str,20,21);
		clsid.Data4[1] = (unsigned char) DecodeHex(	str,22,23);
		clsid.Data4[2] = (unsigned char) DecodeHex( str,25,26);
		clsid.Data4[3] = (unsigned char) DecodeHex( str,27,28);
		clsid.Data4[4] = (unsigned char) DecodeHex( str,29,30);
		clsid.Data4[5] = (unsigned char) DecodeHex( str,31,32);
		clsid.Data4[6] = (unsigned char) DecodeHex( str,33,34);
		clsid.Data4[7] = (unsigned char) DecodeHex( str,35,36);
	}
	else memset( &clsid, 0, sizeof(clsid) );
	return clsid;
}

//---------------------------------------------------------------------------
//						EncodeClassID
//	Encodes the COM object class ID as a HEX string in brace {} brackets.
/*
          1111111111222222222233333333
01234567890123456789012345678901234567
{00000000-0000-0000-0000-000000000000}
*/
//---------------------------------------------------------------------------
nxString EncodeClassID( REFCLSID clsid )
{
	nxString	buffer;

	buffer.sprintf( "{%08X-%04X-%04X-%02X%02X-%02X%02X%02X%02X%02X%02X}",  (int)clsid.Data1,
																		   (int)clsid.Data2,
																		   (int)clsid.Data3,
																		   (int)clsid.Data4[0],
																		   (int)clsid.Data4[1],
																		   (int)clsid.Data4[2],
																		   (int)clsid.Data4[3],
																		   (int)clsid.Data4[4],
																		   (int)clsid.Data4[5],
																		   (int)clsid.Data4[6],
																		   (int)clsid.Data4[7]);
	return buffer;
}

//----------------------------------------------------------------------------
//						ConfigureCLSIDEntries
//----------------------------------------------------------------------------
#if defined(NX_WINDOWS)
// for windows systems write entries to the system regsitry.
nxBOOL ConfigureCLSIDEntries( REFCLSID clsid, const char* modulename, const char *classname, const char* threadmodel )
{
	long		status;
	HKEY		CLSIDKey;
	HKEY		key;
	HKEY		prockey;
	DWORD		disposition;
	nxString	classid = EncodeClassID( clsid ); 
	TCHAR		buffer[512];
	nxBOOL		ok;
	HMODULE		hmod;
	DWORD		threadmodelsize;


	hmod = GetModuleHandle( modulename );
	ok  = (hmod != NULL);
	ok  = ok && (::GetModuleFileName( hmod, buffer, 512 ) > 0);										// Get the module filename
	if (ok)																					// if we got one
	{																						// then
		status = RegOpenKeyEx( HKEY_CLASSES_ROOT, "CLSID", 0, KEY_ALL_ACCESS, &CLSIDKey );	// open the CLSID key
  		ok = (status == ERROR_SUCCESS);														// and see if it worked
		if (ok)																				// if it did 
		{																					// then 
		status = RegCreateKeyEx( CLSIDKey, (const char *)classid, 0, "Onyx",			// create/open the key for this COM object
									 REG_OPTION_NON_VOLATILE, KEY_ALL_ACCESS,NULL,	
									 &key, &disposition); 
			ok = (status == ERROR_SUCCESS);													// see if it as ok
			if (ok)																				// if ok
			{																					// then set the values for this key

				ok  = RegSetValueEx( key, NULL, 0, REG_SZ, (const BYTE *)classname, (DWORD)strlen(classname)+1 ) == ERROR_SUCCESS;

				status = RegCreateKeyEx( key, "InProcServer32", 0, "Onyx",					// create/open the key for InProcServer32
									 REG_OPTION_NON_VOLATILE, KEY_ALL_ACCESS,NULL,	
									 &prockey, &disposition); 
				ok = ok && (status == ERROR_SUCCESS);
				if (ok)
				{
					ok =       RegSetValueEx( prockey, NULL, 0, REG_SZ, (const BYTE *)buffer, (DWORD)strlen(buffer)+1    ) == ERROR_SUCCESS;
					if (threadmodel != NULL)
					{
						threadmodelsize = (DWORD)strlen(threadmodel) + 1;
						ok = ok && (RegSetValueEx( prockey, "ThreadingModel", 0, REG_SZ, (CONST BYTE *)threadmodel, threadmodelsize ) == ERROR_SUCCESS);
					}
				}
			}
			RegCloseKey( key );
		}
		RegCloseKey(CLSIDKey);
	}
	else
	{
		strcpy( buffer, "????");
	}
	if (!ok)
	{
		nxString str;
		str.sprintf("Error registering <%s>", (const char *)buffer);
		::MessageBox( NULL, (const char *)str, "COM registration", MB_OK | MB_ICONSTOP);
	}
	return ok;
}
#else
// ---- for NON windows systems (aka UNIX) do nothing.

nxBOOL ConfigureCLSIDEntries( REFCLSID , const char* , const char *)
{
	return nxTRUE;
}
#endif
