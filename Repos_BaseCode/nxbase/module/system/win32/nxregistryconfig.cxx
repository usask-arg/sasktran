#include "nxbase_core.h"
#include "nxbase_threads.h"


static std::recursive_mutex	g_threadlock;

/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::nxRegistryConfiguration		2005-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryConfiguration::nxRegistryConfiguration(const char*  companyname, const char* filekey, nxRegistryConfiguration::INIRootLocation rootlocation, bool usenativeregistry )
{
	m_usenativeregistry = usenativeregistry && nxRegistryKey::RegistryLocation().UseNativeRegistry();
	SetRootLocation  ( rootlocation );
	SetCompanyKeyName( companyname);
	SetFileKeyName   ( filekey);
}

/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::ShallowCopy		2008-7-29*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::ShallowCopy( const nxRegistryConfiguration& other)
{
	m_usenativeregistry = other.m_usenativeregistry;
	m_rootlocation = other.m_rootlocation;
	m_filekey      = other.m_filekey;
	m_companyname  = other.m_companyname;
	return true;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::~nxRegistryConfiguration		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryConfiguration::~nxRegistryConfiguration()
{
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetRootLocation		2005-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetRootLocation( nxRegistryConfiguration::INIRootLocation sharemode)
{
	m_rootlocation = sharemode;
	return nxTRUE;
}




/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetFileKeyName		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetFileKeyName( const char* filekey)
{
	m_filekey = filekey;
	return true;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::RemoveLeadingAndTrailingDirChars		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::RemoveLeadingAndTrailingDirChars( nxString& basekey )
{
	bool	ok;

	ok = !(basekey.IsEmpty());
	if (ok)
	{
		if (basekey[basekey.GetLength()-1] == '/')
		{
			basekey = basekey.Left(basekey.GetLength()-1);
		}

		if (basekey[0] == '/')
		{
			basekey = basekey.Right(basekey.GetLength()-1);
		}
		ok = !basekey.IsEmpty();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::GenerateKeyName		2005-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

nxString nxRegistryConfiguration::BaseKeyName()
{
	nxString	fullkeyname("/Software/");


	fullkeyname += m_companyname;
	fullkeyname += "/";
	if (( m_rootlocation == USER_APPL_INI) ||( m_rootlocation == GLOBAL_APPL_INI)) fullkeyname	+=  "ApplicationSettings/";
	return fullkeyname;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::CloseKey		2005-3-30*/
/** Closes the current key.  Not normally required by the user as this is
 *	normally done in the destructor.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::CloseKey(nxRegistryKey*	key)
{
	bool	ok;

	ok = (key == NULL);
	if (!ok) ok = key->DestroyKeyHierarchy();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::CheckDirtyAndReopen		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::OpenKey( nxRegistryKey**	key, INIAccessRights accessmode)
{
	nxString		basename;
	bool			ok;

	basename  = BaseKeyName();
#if defined(NX_WINDOWS)
	if ( m_usenativeregistry ) *key  = nxRegistryKeyWin32::CreateKey(basename, m_filekey, m_rootlocation, accessmode);
	else					   *key  = nxRegistryKeyYaml::CreateKey(basename, m_filekey, m_rootlocation, accessmode);
#else
	*key  = nxRegistryKeyYaml::CreateKey(basename, m_filekey, m_rootlocation, accessmode);
#endif

	ok        = (*key != NULL);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::FlushRegistry		 2016- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::FlushRegistry( INIRootLocation location)
{
	bool	ok = true;

	ok = nxRegistryKeyYaml::FlushRegistry(location);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetCompanyKeyName		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetCompanyKeyName( const char* companyname)
{
	m_companyname = companyname;
	RemoveLeadingAndTrailingDirChars( m_companyname );
	return true;
}



/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetString		2005-3-31*/
/** Sets the string of a named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	The value of the null terminated string stored in the named value.
 *
 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetString( const char * name, const char *  str )
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey( &key, FULLIO_INI);
	ok = ok && key->SetString( name, str );
	ok = ok && CloseKey( key );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::GetString		2005-3-31*/
/** retrieves the string of a named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	 Returns the string stored in the named value.
 *
 *	\param defaultstr
 *	This string is returned if there is an error retrieving a string from the named value.

 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::GetString( const char * name, nxString* str)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key, READ_INI);
	ok = ok && key->GetString( name, str );
	CloseKey( key);
	if (!ok) str->Empty(false);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetString		2005-3-31*/
/** Sets the string of a named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	The value of the null terminated string stored in the named value.
 *
 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetPath( const char * name, const char *  str )
{
	bool			ok;
	nxRegistryKey*	key;
	nxString		vstr(str);

	vstr.MakeDirectorySeparatorsOSConsistent('/');
	{
		std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

		ok = OpenKey( &key, FULLIO_INI);
		ok = ok && key->SetString( name, vstr );
		ok = ok && CloseKey( key );
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::GetPath		2005-3-31*/
/** retrieves the string of a named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	 Returns the string stored in the named value.
 *
 *	\param defaultstr
 *	This string is returned if there is an error retrieving a string from the named value.

 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::GetPath( const char * name, nxString* str)
{
	bool			ok;
	nxRegistryKey*	key;

	{
		std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

		ok = OpenKey(&key, READ_INI);
		ok = ok && key->GetString( name, str );
		CloseKey( key);
	}
	if (!ok) str->Empty(false);
	str->MakeDirectorySeparatorsOSConsistent();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetDouble		2005-3-31*/
/** Sets a named value belonging to this key to the binary representation of a double. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param value
 *	The double value to stored in the named value.
 *
 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetDouble( const char * name, double value)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key, FULLIO_INI);
	ok = ok && key->SetDouble( name, value );
	ok = ok && CloseKey( key);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::GetDouble		2005-3-31*/
/** retrieves a double from the named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	 Returns the double stored in the named value.
 *
 *	\param defaultstr
 *	This value is returned if there is an error retrieving the value from the named value.

 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::GetDouble( const char * name, double* value, double defaultvalue)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key, READ_INI);
	ok = ok && key->GetDouble( name, value );
	if (!ok) *value = defaultvalue;
	CloseKey( key);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetInteger		2005-3-31*/
/** Sets a named value belonging to this key to the 32 bit representation of an integer.
 *
 *	\param name
 *	The name of the key value
 *
 *	\param value
 *	The integer value to stored in the named value.
 *
 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetInteger( const char * name, int value)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key,FULLIO_INI);
	ok = ok && key->SetInteger( name, value );
	ok = ok && CloseKey( key);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::GetInteger		2005-3-31*/
/** retrieves a 32 bit integer from the named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	 Returns the integer stored in the named value.
 *
 *	\param defaultstr
 *	This value is returned if there is an error retrieving the value from the named value.

 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::GetInteger( const char * name, int* value, int defaultvalue)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key, READ_INI);
	ok = ok && key->GetInteger( name, value );
	CloseKey( key);
	if (!ok) *value = defaultvalue;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::SetBool		2005-3-31*/
/** Sets a named value belonging to this key to the "1 or 0" representation of a bool. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param value
 *	The boolean value to stored in the named value.
 *
 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::SetBool( const char * name, bool value)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key, FULLIO_INI);
	ok = ok && key->SetBool( name, value );
	ok = ok && CloseKey( key);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::GetBool		2005-3-31*/
/** retrieves a boolean from the named value belonging to this key. 
 *
 *	\param name
 *	The name of the key value
 *
 *	\param str
 *	 Returns the boolean stored in the named value.
 *
 *	\param defaultstr
 *	This value is returned if there is an error retrieving the value from the named value.

 *	\return
 *	nxTRUE if successful.
 **/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::GetBool( const char * name, bool* value, bool defaultvalue)
{
	bool			ok;
	nxRegistryKey*	key;
	std::lock_guard<std::recursive_mutex>	lock( g_threadlock );

	ok = OpenKey(&key, READ_INI);
	ok = ok && key->GetBool( name, value );
	CloseKey( key);
	if (!ok) *value = defaultvalue;
	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::LocateDirectoryFromKey		2005-3-31*/
/** Locates a directory from the specified registry key and returns the
 *	directory as a string. There is no directory delimiter on the end of the directory name. 
 *	The code provides the option to create the directory if it does not already and
 *	will create an entire directory tree as necessary. If the named directory cannot
 *	be created or the key does not exist then
 *	there is an option to browse for a new directory and update/create the key entry.
 *
 *	\param	keyname
 *	The name of the registry subkey.
 *
 *	\param dirname
 *	Returns the name of the directory.  It will be empty if the procedure fails.
 *	Cannot be NULL.
 *
 *	\param createifnotexist
 *	If TRUE then create the directory specified by the key if it does not exist.
 *	
 *	\param browseifnokey
 *	If there is no value for this key or we cannot create the directory then if
 *	this \b browseifnokey is TRUE then pop up a dialog for the user to browse
 *	for a new directory.
 *
 *	\param browsepromptstring
 *	If the code decides to browse for a new directory then it will display
 *	this string (if its not null) inside the dialog box
 *
 *	\return
 *	nxTRUE if we have success, nxFALSE otehrwise.
**/
/*---------------------------------------------------------------------------*/

bool nxRegistryConfiguration::LocateDirectoryFromKey( const char *	keyname,
														nxString*		dirname,
														bool			createifnotexist,
														bool			browseifnokey,
														const char*		browsepromptstring)
{
	nxString	name;
	bool		ok;
	bool		exists;
	static		std::recursive_mutex	threadlock;

	ok = GetString(keyname, &name);					// Get the directory name from the key.
	if (ok)
	{
		exists = nxFile::Exists(name);
		if (!exists && createifnotexist) 
		{
			std::lock_guard<std::recursive_mutex> lock(threadlock);
			ok = nxFile::Exists(name);								// Once we get the thread back see if the directory exists as another thread may have created it
			if (!ok) ok = nxDirectory::CreateADirectory(name);		// If it does not then we can create it
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"nxRegistryConfiguration::LocateDirectoryFromKey: There were errors trying to create directory <%s>", (const char*)name);
			}
		}
	}
	if (!ok && browseifnokey)
	{
		ok = BrowseForDirectory( &name, browsepromptstring );
		if (ok)
		{
			ok = SetString( keyname, name );
		}
	}
	if ( ok ) *dirname = name;
	else       dirname->Empty(nxFALSE);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxRegistryConfiguration::LocateDirectoryFromKey, There were errors looking up directory setting from registry key [%s]", (const char*)keyname);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::BrowseForDirectory		2005-3-31*/
/** Browse for a new directory. On linux use the console to get the information
 **/
/*---------------------------------------------------------------------------*/

//#if defined(NX_UNIX_VER)
bool nxRegistryConfiguration::BrowseForDirectory( nxString* /*dirname*/, const char * promptstr )
{
	nxLog::Record( NXLOG_WARNING, "nxRegistryConfiguration::BrowseForDirectory, Your registry settings for [%s] are either improperly installed or there is an access issue. You need to refer to installation documentation or see someone who knows how to fix the problem.", (const char*)promptstr);
	return false;
}
//#else
/*-----------------------------------------------------------------------------
 *					nxRegistryConfiguration::BrowseForDirectory		2005-3-31*/
/** Windows version: uses the Shell functions to browse for a new directory.**/
/*---------------------------------------------------------------------------*/

#if 0
#include "shlobj.h"
bool nxRegistryConfiguration::BrowseForDirectory( nxString* dirname, const char * promptstr )
{
	BROWSEINFO	 info;
	char		 displayname[MAX_PATH+2];
	LPITEMIDLIST pidlSelected = NULL;
	IMalloc*	 pmalloc;
	bool		ok;


	CoInitialize(NULL);
	ok = (SHGetMalloc( &pmalloc ) == NOERROR);
	if (ok)
	{
		info.hwndOwner      = NULL;
		info.pidlRoot       = NULL;
		info.pszDisplayName = &displayname[0];
		info.lpszTitle      = (promptstr != NULL) ? promptstr : "Browse for a directory";
		info.ulFlags        = BIF_NEWDIALOGSTYLE;
		info.lpfn           = NULL;
		info.lParam         = 0;
		info.iImage         = 0;

		pidlSelected = SHBrowseForFolder( &info );
		ok = (pidlSelected != NULL);
		if (ok)
		{
			ok = (SHGetPathFromIDList( pidlSelected, &displayname[0] ) == TRUE);
			if (ok)
			{
				*dirname = displayname;
			}
			pmalloc->Free( pidlSelected);
		}
		pmalloc->Release();
	}
	if (!ok) dirname->Empty( nxFALSE);
	CoUninitialize();
	return ok;
}
#endif
