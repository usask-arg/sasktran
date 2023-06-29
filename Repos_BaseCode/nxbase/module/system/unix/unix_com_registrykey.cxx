/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_core.h"
#include <boost/thread.hpp>

nxRegistryKey_RegistryLocation nxRegistryKey::g_registrylocation;

//#define NXMAKESTRG(str) # str								// The stringizing macro, places quites around the parameter (without macro expansion)
//#define NXMAKESTRGA(str) NXMAKESTRG(str)					// Expand the str macro into its base for

/*-----------------------------------------------------------------------------
 *					nxRegistryKey_RegistryLocation::nxRegistryKey_RegistryLocation		 2016- 11- 1*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKey_RegistryLocation::nxRegistryKey_RegistryLocation()
{
#if defined(NX_WINDOWS)
	m_usenativeregistry = true;
#else
	m_usenativeregistry = false;
#endif
	Default_BaseDirectory();
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::GetSystem_BaseDirectory		2009-12-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKey_RegistryLocation::Default_BaseDirectory(  )
{
	const char* envptr;
	nxString	envvalue;

	m_basedirectory.Empty(true);
	envptr      = getenv("NXGLOBALREGISTRYDIR");		// On windows  or Linux you can use an environment variable
	if ( envptr != NULL)
	{
		m_basedirectory = envptr;
	}
	else
	{
		envptr = getenv("HOME");
		if (envptr != NULL)
		{
			m_basedirectory = envptr;
		}
		else
		{
			m_basedirectory = getenv("HOMEDRIVE");
			m_basedirectory += getenv("HOMEPATH");
		}
		if (!m_basedirectory.IsEmpty())
		{
			m_basedirectory.EnsureLastCharIsDirectoryChar();
			m_basedirectory += ".argregistry/";
			m_basedirectory.MakeDirectorySeparatorsOSConsistent();
		}
	}
//	printf("The default registry setting in this software module is <%s>\n", (const char*)m_basedirectory);
	return true;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKey_RegistryLocation::Set_BaseDirectory		 2016- 11- 1*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKey_RegistryLocation::Set_BaseDirectory( const char* dirname )
{
	m_basedirectory = dirname;
	if (dirname[0] != '\0')
	{
		m_usenativeregistry = false;
	}
	else
	{
		m_usenativeregistry = true;
	}
//	nxLog::Record(NXLOG_INFO,"Registry Location: The USASK-ARG registry for this software module is located at <%s>", (const char*)m_basedirectory);
//	printf("Registry Location: The USASK-ARG registry for this software module is located at <%s>", (const char*)m_basedirectory);
	return true;
}



/****************************************************************************
*
*	Section Names are case insensitive.
*   Values  Names are case insensitive
*   Values are case sensitive
*   Exclamation marks (even inside strings!!!) mark the start of a comment.
*
*    [1stSectionName]
*	           = avalue				! Optional Default Value
*      Key1    = value2
*	   Key2    = value1
*	   [subkey1]
*        Key1  = value2
*	     Key2  = value1
*	   []
*	 []
*
*	 [2ndSectionName]
*      Key1  = value2
*	   Key2  = value1
*	 []
*
****************************************************************************/



/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::nxRegistryKeyUnix		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyUnix::nxRegistryKeyUnix( nxRegistryKeyUnix* parent, nxRegistryConfiguration::INIAccessRights accessmode )
{
	m_parent      = parent;
	m_isdirty     = false;
	m_accessmode  = accessmode;
}

//---------------------------------------------------------------------------
//						nxRegistryKeyUnix::destructor
//--------------------------------------------------------------------------

nxRegistryKeyUnix::~nxRegistryKeyUnix()
{
	erase();
}

//---------------------------------------------------------------------------
//						nxRegistryKeyUnix::SetName
//--------------------------------------------------------------------------

void nxRegistryKeyUnix::SetSectionName( const nxString& name )
{
	m_sectionname = name;
}

//---------------------------------------------------------------------------
//						nxRegistryKeyUnix::erase
//--------------------------------------------------------------------------

void nxRegistryKeyUnix::erase()
{
	subkeyIterator ptr;
	nxRegistryKeyUnix* keyptr;

	m_values.erase( m_values.begin(), m_values.end() );
	for (ptr = m_subkeys.begin(); !(ptr == m_subkeys.end()); ptr++ )
	{
		keyptr = (*ptr);
		if (keyptr != NULL)
		{
			NXTRACE(( "nxRegistryKeyUnix::erase for key [%s], deleting subkey [%s]\n", (const char *)SectionName(), (const char *)keyptr->SectionName() ));
			delete keyptr;
		}
	}
	m_subkeys.erase( m_subkeys.begin(), m_subkeys.end() );
}

//---------------------------------------------------------------------------
//						nxRegistryKeyUnix::ParseLine
//	Parse a line of text from an ini file and return the state of the line.
//---------------------------------------------------------------------------

nxRegistryKeyUnix::KEYSTATE nxRegistryKeyUnix::ParseLine( nxString &line, nxString* name, nxString* value )
{
	int idx;
	KEYSTATE state;
	nxString nextline;

	idx = line.Find('!');									// Find any comment fields
	if (idx >= 0) line = line.Left(idx+1 );					// and remove those guys
	line.RemoveWhiteSpace();								// remove white space

	if (line .GetLength() < 1 )								// does the string have any length
	{														// nope
		state = NOP;										// so return straight away
	}														// otherwise do the parsing
	else													// so
	{
		idx = line.Find( '[' );								// Is this a new subkey of the end of the current key
		if (idx >= 0)										// Yup
		{													// so decide if its a new subkey or the end of the current key
			int endidx = line.Find( ']' );					// find the terminating bracket
			int nchar = endidx - idx - 1;					// get the number of characters between the brackets
			if (nchar >= 0 )								// if the number is zero or positive
			{												// then
				*name = line.Mid(idx+1, nchar );			// get the subkey name
				name->RemoveWhiteSpace();					// remove any white space
				if (name->GetLength() < 1 )					// if its empty
				{											// then
					state = ENDOFKEY;						// this is the end of the current key
				}											// otherwise
				else										// this
				{											// is the start of
					state    = NEWSUBKEY;						// a new subkey.
					nextline = line.Right(line.GetLength()-endidx-1);
				}											// and we have
			}												// finished the processing of valid '[' key
			else											// otherwise we failed to find proper placement of ']'
			{												// so
				state = NOP;								// invalidate this line
			}
		}
		else												// This is not a new subkey or end of this key so see if its a value.
		{													// So this is a key value
			idx = line.Find( '=' );							// all values must have an "=" character
			if (idx >= 0 )									// is there an equals sign
			{												// yup
				*name  = line.Left(idx);					// Get the value name
				*value = line.Right( line.GetLength() - idx -1 );	// Get the value;
				name->RemoveWhiteSpace();					// Remove white space from the name
				value->RemoveWhiteSpace();					// remove white space from the value
				state  = NEWVALUE;							// and set the state to a new value
			}
			else											// otherwise we found no "=" sign
			{												// so this is an invalid line
				state = NOP;								// so set the line to invalid
			}
		}
	}
	line = nextline;
	return state;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::AddNewValue		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::AddNewValue( const char* name, const char* str )
{
	nxRegistryValueUnix* val;
	nxRegistryValueUnix	 dummyvalue;
	nxString			 astr(str);

	m_values.push_back( dummyvalue );
	val = &m_values.back();
	val->SetName(name);
	val->SetValue(astr);
	SetDirty();
	return true;
}

//---------------------------------------------------------------------------
//						nxRegistryKeyUnix
//--------------------------------------------------------------------------

void  nxRegistryKeyUnix::ReadKey( nxFile& infile )
{
	nxString		line;
	nxString		name;
	nxString		value;
	KEYSTATE		state;
	nxBOOL				 more = nxTRUE;
	nxRegistryKeyUnix*   key;
//	nxRegistryValueUnix	 dummyvalue;
//	nxRegistryValueUnix* val;

	key = this;
   	while (!infile.eof() && more )										//
	{

		if (line.IsEmpty())line = infile.ReadALine();			// Read a line from a file
		state = ParseLine( line, &name, &value );				// Parse the line

		switch (state)
		{
		case NEWSUBKEY :	key = new nxRegistryKeyUnix(key, m_accessmode);
							if (key == NULL) NXTHROW(("Memory Allocation error in nxRegistryKeyUnix::ReadKey\n"));
							key->SetSectionName( name );
							key->ReadKey( infile );
							m_subkeys.push_back( key );
							break;
		case NEWVALUE  :	AddNewValue( name, value );
							break;
		case ENDOFKEY  : 	more = nxFALSE;
							break;
		case NOP       :	break;
		};
	}
}

static boost::recursive_mutex unixkey_iolock;

//--------------------------------------------------------------------------
//                            nxRegistry::ReadFile
//--------------------------------------------------------------------------

void nxRegistryKeyUnix::ReadFile( const char * FileName )
{
	nxFile		infile;
	boost::lock_guard<boost::recursive_mutex> lock(unixkey_iolock);

	m_fullpathname = FileName;
	erase();									// Erase this key
	infile.Open( FileName,  "rt");				// Open the file
	ReadKey( infile );							// Read this key
	infile.Close();      						// then close the file
	m_isdirty  = false;							// And after reading teh file, we are are clean
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::CheckAndSaveFile		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::CheckAndSaveFile()
{
	bool	ok;

	ok = (this == RootParent());
	if (m_accessmode != nxRegistryConfiguration::READ_INI)
	{
		ok = !m_isdirty;
		if (!ok)
		{
			ok = WriteFile();
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::DestroyKeyHierarchy		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::DestroyKeyHierarchy()
{
	bool	ok;

	nxRegistryKeyUnix*	root;

	root = RootParent();
	ok = root->CheckAndSaveFile();
	delete root;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::WriteEntriesToFile		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::WriteEntriesToFile( nxFile& f )
{
	valueIterator			ptr;
	nxRegistryValueUnix*	value;
	subkeyIterator			subkeyptr;
	bool					ok;
	bool					ok2;
	bool					ok1;
	nxString				str;
	bool					writesection;

	ok = true;

	writesection = ( !m_sectionname.IsEmpty() );
	if (writesection)
	{
		str.sprintf("[%s]\n",(const char*)m_sectionname);
		ok1 = f.WriteString( str );
		ok = ok && ok1;
	}

	for (ptr = m_values.begin(); !(ptr == m_values.end()); ++ptr)
	{
		value = &(*ptr);
		str.sprintf( "%s = %s\n", (const char *)value->Name(), (const char*) value->Value() );
		ok1 = f.WriteString( str );
		ok = ok && ok1;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," nxRegistryKeyUnix::WriteEntriesToFile, Error writing value entries to file <%s>", (const char*)f.Filename() );
	}

	ok2 = true;
	for (subkeyptr = m_subkeys.begin(); !(subkeyptr == m_subkeys.end()); ++subkeyptr)
	{
		ok1 = (*subkeyptr)->WriteEntriesToFile( f );
		ok2 = ok2 && ok1;
	}
	if (!ok2)
	{
		nxLog::Record(NXLOG_WARNING," nxRegistryKeyUnix::WriteEntriesToFile, Error writing subkey entries to file <%s>", (const char*)f.Filename() );
	}

	if (writesection)
	{
		str.sprintf("[ ]\n");
		ok1 = f.WriteString( str );
		ok = ok && ok1;
	}
	ok = ok && ok2;

	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::WriteFile		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::WriteFile( )
{
	nxFile	f;
	bool	ok;
	boost::lock_guard<boost::recursive_mutex> lock(unixkey_iolock);

	f.Open(m_fullpathname, "wt");
	ok = f.IsOpen();
	ok = ok && WriteEntriesToFile(f);
	f.Close();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxRegistryKeyUnix::WriteFile, there were errors writing registry settings to file <%s>", (const char*)m_fullpathname);
	}
	m_isdirty = !ok;
	return ok;
}

//--------------------------------------------------------------------------
//                            nxRegistry::GetKey
//--------------------------------------------------------------------------


bool nxRegistryKeyUnix::FindKey( nxString& name, nxRegistryKeyUnix** userkey )
{
	bool						ok;
	nxRegistryKeyUnix*			key;
	EqualKeys					compare(name );
	subkeyIterator				ptr;

	ptr = std::find_if(	m_subkeys.begin(), m_subkeys.end(),  compare );
	ok =  !(ptr == m_subkeys.end());
	if (ok)
	{
		key = (*ptr);
	}
	else
	{
		if (m_accessmode == nxRegistryConfiguration::READ_INI)
		{
			nxLog::Record(NXLOG_WARNING,"nxRegistryKeyUnix::FindKey, Cannot find the key with name <%s> in the section", (const char*) name );
			key = NULL;
		}
		else
		{
			key = new nxRegistryKeyUnix(this, m_accessmode);
			ok = (key != NULL);
			if (!ok)
			{
				nxLog::Record(NXLOG_ERROR, "Memory Allocation error in nxRegistryKeyUnix::FindKey");
			}
			else
			{
				key->SetSectionName( name );
				m_subkeys.push_back( key );
			}
		}
	}
	*userkey = key;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::FindValue		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::FindValue( nxString& name, nxRegistryValueUnix** userkey )
{
	nxRegistryValueUnix*		key;
	Equalvalues					compare(name );
	valueIterator				ptr;

	ptr = std::find_if(	m_values.begin(), m_values.end(),  compare );
	if (!(ptr == m_values.end())) key = &(*ptr);
	else                          key = NULL;
	*userkey = key;
	return (key != NULL);
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::RootParent		2009-6-2*/
/** Return the key which is the root parent.**/
/*---------------------------------------------------------------------------*/

nxRegistryKeyUnix* nxRegistryKeyUnix::RootParent( )
{
	nxRegistryKeyUnix*	thisobject = this;
	nxRegistryKeyUnix*	nextparent = m_parent;
	while( nextparent != NULL )
	{
		thisobject = nextparent;
		nextparent = thisobject->m_parent;
	}
	return thisobject;
}


//--------------------------------------------------------------------------
//                            nxUnixRegistry::GetValue
//--------------------------------------------------------------------------

bool nxRegistryKeyUnix::GetValue( const char* name, nxString* uservalue ) const
{
	const nxRegistryValueUnix*	key;
	Equalvalues 			compare(name );
	const_valueIterator		ptr;
	bool					ok;

	ptr = STL(find_if)(	m_values.begin(), m_values.end(),  compare );
	if (!(ptr == m_values.end())) key = &(*ptr);
	else                           key = NULL;
	ok = (key != NULL);
	if (ok) *uservalue = key->Value();
	else    uservalue->Empty(false);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::SetDoubleValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::SetDouble( const char * valuename, double value)
{
	nxString	str;
	str.sprintf( "%g", (int)value);
	return SetString(valuename, str );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::GetDoubleValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::GetDouble( const char * valuename, double* value)
{
	nxString	str;
	bool		ok;

	ok = GetString( valuename, &str );
//	str.MakeUpper();
	*value = atof( str );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::SetIntValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::SetInteger( const char * valuename, int value)
{
	nxString	str;
	str.sprintf( "%d", (int)value);
	return SetString(valuename, str );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::GetIntValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::GetInteger( const char * valuename, int* value)
{
	nxString	str;
	bool		ok;

	ok = GetString( valuename, &str );
//	str.MakeUpper();
	*value = atoi( str );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::SetBoolValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::SetBool( const char * valuename, bool value)
{
	const char*	strtrue = "True";
	const char* strfalse = "False";
	const char* str;

	str = value ? strtrue : strfalse;
	return SetString( valuename, str );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::GetBoolValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::GetBool( const char * valuename, bool* value)
{
	nxString	str;
	bool		ok;

	ok = GetString( valuename, &str );
	str.MakeUpper();
	*value = ok &&  ( (str == "1") || (str == "TRUE") || (str == "T") );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::SetString		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::SetString( const char * valuename, const char* str  )
{
	nxRegistryValueUnix*	value;
	bool					ok = true;
	nxString				vstr;
	nxString				bstr;

	ok = (m_accessmode != nxRegistryConfiguration::READ_INI);				// MAke sure we are not in readonly mode
	if (!ok)																// if we are in readonly mode
	{																		// then log a message
		nxLog::Record(NXLOG_WARNING, "nxRegistryKeyUnix::SetString, You cannot set values in keys that opened in read only mode");
	}
	else																	// otherwsie we are in write mode
	{																		// so
		RootParent()->SetDirty();											// notify the root we are dirty
		vstr = valuename;
		vstr.MakeLower();
		ok = FindValue( vstr, &value );										// find the value we wish to change
		if (ok)
		{
			bstr = str;
			value->SetValue( bstr );
		}
		else
		{
			ok = AddNewValue( vstr, str );
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::SetString		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyUnix::GetString( const char * valuename, nxString* str  )
{
	nxRegistryValueUnix*	value;
	bool					ok = true;
	nxString				vstr;

	vstr = valuename;
	vstr.MakeLower();
	ok = FindValue( vstr, &value );								// find the value we wish to change
	if (ok)
	{
		*str = value->Value();
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"nxRegistryKeyUnix::GetString, Error returning value for <%s> from key <%s>", (const char*)valuename, (const char*)m_sectionname);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix::CreateKey		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyUnix* nxRegistryKeyUnix::CreateKey(  const char *									basekey,
												  const char*									filekeyname,
												  nxRegistryConfiguration::INIRootLocation		rootlocation,
												  nxRegistryConfiguration::INIAccessRights		access)
{
	nxRegistryKeyUnix*	key    = NULL;
	nxStringArray		filetokens;
	bool				ok = true;
	nxString			filename;
	nxString			thename;
	nxString			dirname;
	nxString			name;
	nxString			userbase;
	nxString			regdir;
	nxString			basekeyname(basekey);
	nxString			filekey( filekeyname);
	int	i;

	// ---- Encode the access mode into a base directory location

	switch (rootlocation)
	{
	case  nxRegistryConfiguration::USER_APPL_INI       : userbase = nxRegistryKey::RegistryLocation().BaseDirectory();
														 regdir  = ".nxuserregistry";
														 thename = "userapplkey.ini";
														 break;

		case  nxRegistryConfiguration::USER_INI        : userbase = nxRegistryKey::RegistryLocation().BaseDirectory();
														 regdir  = ".nxuserregistry";
														 thename = "userkey.ini";
														 break;

		case  nxRegistryConfiguration::GLOBAL_INI      : userbase = nxRegistryKey::RegistryLocation().BaseDirectory();
														 regdir  = ".nxglobalregistry";
														 thename = "globalkey.ini";
														 break;

		case  nxRegistryConfiguration::GLOBAL_APPL_INI : userbase = nxRegistryKey::RegistryLocation().BaseDirectory();
														 regdir  = ".nxglobalregistry";
														 thename = "globalapplkey.ini";
														 break;

		default										   : ok = false;
														 break;
	}

	if (ok)
	{
		basekeyname.MakeLower();
		filekey.MakeLower();
		dirname.sprintf("%s/%s/%s/%s/", (const char*)userbase, (const char*)regdir, (const char*)basekeyname, (const char*) filekey);
		filetokens.Strtok( dirname, "\\/" );							// Parse the filekey. The first part is the filename
		if (dirname.GetAt(0) == DIRECTORY_CHAR) filename = DIRECTORY_CHAR;	// if the first char is "/" then make sure that goes across
		for (i=0; i < filetokens.GetSize(); i++)						// and reconstruct removing any "doubling of folder chars"
		{
			filename += filetokens.GetAt(i);
			filename += DIRECTORY_CHAR;
		}

		// ---- For write access make sure we can create the output directory

		if ( (access == nxRegistryConfiguration::FULLIO_INI) || ( access == nxRegistryConfiguration::WRITE_INI) )
		{
			ok = nxDirectory::CreateADirectory( filename );
			if (!ok) nxLog::Record(NXLOG_WARNING,"nxRegistryKeyUnix::CreateKey, the key cannot be opened for modification as directory <%s> cannot be created", (const char*)dirname);
		}

		// ---- now generate the full path to the filename

		if (ok)
		{
			filename += thename;
			key = new nxRegistryKeyUnix(NULL, access );						// Create the parent key
			key->ReadFile(filename);										// Read in the file keys (this reads in everything)
		}
	}
	return key;
}
