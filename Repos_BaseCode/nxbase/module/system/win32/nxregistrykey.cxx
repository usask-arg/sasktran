#include "nxbase_core.h"


//---------------------------------------------------------------------------
//						nxRegistryKeyWin32::CreateKey
//	Create a key that takes a keyname in the form xxx\yyy\zzz\xxx and makes
//	a key to 
//---------------------------------------------------------------------------

nxRegistryKeyWin32* nxRegistryKeyWin32::CreateKey( const char * keyname, const char* filekey, nxRegistryConfiguration::INIRootLocation rootlocation, nxRegistryConfiguration::INIAccessRights accessrights)
{
	nxStringArray		tokens;
	nxRegistryKeyWin32*	lastkey = NULL;
	nxRegistryKeyWin32*	key = NULL;
	HKEY				basekey;
	REGSAM				access;
	nxString			fullkey(keyname);

	fullkey += filekey;
	basekey = RootLocation(rootlocation);
	access  = AccessRights( accessrights );

	tokens.Strtok( fullkey, "/" );
	for (int i= 0; i < tokens.GetSize(); i++)
	{
		if (tokens[i].GetLength() > 0 )
		{
			key = new nxRegistryKeyWin32( lastkey, tokens[i], basekey, access );
			lastkey = key;
		}
	}
	return key;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::Root		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

HKEY nxRegistryKeyWin32::RootLocation( nxRegistryConfiguration::INIRootLocation rootlocation )
{
	HKEY	handle = NULL;

	switch (rootlocation)
	{
	case nxRegistryConfiguration::GLOBAL_APPL_INI	: handle =  HKEY_LOCAL_MACHINE; break;
	case nxRegistryConfiguration::GLOBAL_INI		: handle =  HKEY_LOCAL_MACHINE; break; 
	case nxRegistryConfiguration::USER_APPL_INI		: handle =  HKEY_CURRENT_USER;  break;
	case nxRegistryConfiguration::USER_INI		    : handle =  HKEY_CURRENT_USER;  break; 
	default                                         : nxLog::Record(NXLOG_WARNING, "nxRegistryKeyWin32::RootLocation, Uknown rootlocatioon. Thats a problem.");
	};
	return handle;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::AccessRights		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

REGSAM nxRegistryKeyWin32::AccessRights( nxRegistryConfiguration::INIAccessRights accessrights)
{
	REGSAM rights;

	switch (accessrights)
	{
		case nxRegistryConfiguration::READ_INI   : rights = KEY_READ; break;
		case nxRegistryConfiguration::WRITE_INI  : rights = KEY_WRITE; break;
		case nxRegistryConfiguration::FULLIO_INI : rights = KEY_ALL_ACCESS; break;
		default								     : rights = KEY_QUERY_VALUE; break;
	}
	return rights;
}


//---------------------------------------------------------------------------
//						nxRegistryKeyWin32::Constructior
//---------------------------------------------------------------------------

nxRegistryKeyWin32::nxRegistryKeyWin32( nxRegistryKeyWin32* parent, const char * name, HKEY basekey, REGSAM access )
{
//	DWORD	disposition;
	LONG	status;
	bool	ok;

	m_parent = parent;
	m_keyname = name;
//	NXTRACE(("Creating registry key that is connected to <%s>\n", (const char*)m_keyname ));

	if (m_parent != NULL) basekey = m_parent->Key();
	status = RegOpenKeyEx( basekey, name, 0, access, &m_key);

/*
	status = RegCreateKeyEx( basekey, name, 0, "nx",							// create/open the key for this COM object
									 REG_OPTION_NON_VOLATILE,
									 access,
									 NULL,	
									 &m_key,
									 &disposition); 
*/
 	 ok = (status == ERROR_SUCCESS);														// and see if it worked
	 if (!ok)
	 {
		 nxLog::Record( NXLOG_WARNING, "nxRegistryKeyWin32::Constructor, Error creating key <%s>", (const char *)name );
		 m_key = NULL;
	 }
}

//---------------------------------------------------------------------------
//						nxRegistryKeyWin32::Destructor
//	Follow down the chain until we hit the key with no parent.
//---------------------------------------------------------------------------

nxRegistryKeyWin32::~nxRegistryKeyWin32()
{
	RegCloseKey( m_key );
//	NXTRACE(("Destroying registry key that is connected to <%s>\n", (const char*)m_keyname ));
	m_keyname.Empty(nxTRUE);
	if (m_parent != NULL) delete m_parent;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::SetValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::SetValue( const char* valuename,  DWORD dwType, const BYTE* lpData, DWORD cbData)
{
	bool	ok;
	LONG	status;

	ok = (m_key != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxRegistryKeyWin32::SetValue, the Windows key is NULL. Cannot set the value");
	}
	else
	{
		status = ::RegSetValueEx( m_key, valuename, 0, dwType, lpData, cbData);
		ok = (status == ERROR_SUCCESS);
		if (!ok)
		{
			char	buffer[1000];
			::FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM, NULL, status, 0, &buffer[0], N_ELEMENTS(buffer), NULL );
			nxLog::Record(NXLOG_WARNING, "nxRegistryKeyWin32::SetValue, Error setting value of key <%s>, Message =  %s", (const char *)valuename, (const char *)buffer );
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::GetValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::GetValue( const char* valuename,  DWORD dwType,  BYTE* lpData, DWORD* cbData, DWORD mandatorysize)
{
	bool	ok;
	LONG	status;
	DWORD	type;

	type = dwType;
	ok = (m_key != NULL);
	if (ok)
	{
		status = ::RegQueryValueEx( m_key, valuename, NULL, &type, lpData, cbData );
		ok = (status == ERROR_SUCCESS);
		if (ok)
		{
			ok = ok && (type == dwType);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "nxRegistryKeyWin32::GetValue, Key <%s> is not of the correct registry type", (const char *)valuename );
			}
			else
			{
				if (mandatorysize != 0)
				{
					ok = ok && (*cbData == mandatorysize);
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING, "nxRegistryKeyWin32::GetValue, Key <%s> is not of the correct size type, expected (%d) but got (%d)", (const char *)valuename, (int)mandatorysize, (int)(*cbData) );
					}
				}
			}
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::SetStringValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::SetString( const char * valuename, const char *str )
{
	return SetValue( valuename, REG_SZ, (const BYTE *)str, (DWORD)strlen(str)+1 );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::GetStringValue		2005-3-30*/
/**	Fetch the value of a string associated with this key in the registry
 *	Return NULL if the string is not found
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::GetString( const char * valuename, nxString* str )
{
	char		buffer[1024];
	DWORD		bufsize = sizeof(buffer)-1;
	bool		ok;

	ok = GetValue( valuename, REG_SZ, (LPBYTE) buffer, &bufsize, 0 );
	if (ok)
	{
		buffer[bufsize+1] = '\0';			//Make sure we have a NULL terminated string to avoid buffer overrun (Common registry problem apparently)
		*str = (const char *)buffer;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::SetDoubleValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::SetDouble( const char * valuename, double value)
{
	return SetValue( valuename, REG_BINARY, (const BYTE *)&value, (DWORD)sizeof(value) );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::GetDoubleValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::GetDouble( const char * valuename, double* value)
{
	double		buffer;
	DWORD		bufsize = sizeof(buffer);
	bool		ok;

	ok = GetValue( valuename, REG_BINARY, (LPBYTE)&buffer, &bufsize, sizeof(buffer) );
	if (ok)
	{
		*value = buffer;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::SetIntValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::SetInteger( const char * valuename, int value)
{
	DWORD	buffer = (DWORD)value;
	return SetValue( valuename, REG_DWORD, (const BYTE *)&buffer, (DWORD)sizeof(buffer) );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::GetIntValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::GetInteger( const char * valuename, int* value)
{
	DWORD		buffer;
	DWORD		bufsize = sizeof(buffer);
	bool		ok;

	ok = GetValue( valuename, REG_DWORD, (LPBYTE)&buffer, &bufsize, sizeof(buffer) );
	if (ok)
	{
		*value = (int)buffer;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::SetBoolValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::SetBool( const char * valuename, bool value)
{
	DWORD	buffer = value ? 1 : 0;
	return SetValue( valuename, REG_DWORD, (const BYTE *)&buffer, (DWORD)sizeof(buffer) );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyWin32::GetBoolValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyWin32::GetBool( const char * valuename, bool* value)
{
	DWORD		buffer;
	DWORD		bufsize = sizeof(buffer);
	bool		ok;

	ok = GetValue( valuename, REG_DWORD, (LPBYTE)&buffer, &bufsize, sizeof(buffer) );
	if (ok)
	{
		*value =  (buffer != 0);
	}
	return ok;
}




