
#include "nxbase_core.h"
#include <boost/thread.hpp>

static nxRegistryKeyYamlRoot	g_globalapplkey ("globalapplkey.yaml");
static nxRegistryKeyYamlRoot	g_globalkey		("globalkey.yaml");
static nxRegistryKeyYamlRoot	g_userkey		("userkey.yaml");
static nxRegistryKeyYamlRoot	g_userapplkey	("userapplkey.yaml");
boost::mutex					g_yamlmutex;
boost::mutex					g_yamloutputmutex;

nxRegistryKeyYamlRoot*	nxRegistryKeyYamlRoot::GlobalApplKey()  { return &g_globalapplkey;}
nxRegistryKeyYamlRoot*	nxRegistryKeyYamlRoot::GlobalKey()		{ return &g_globalkey;}
nxRegistryKeyYamlRoot*	nxRegistryKeyYamlRoot::UserKey()		{ return &g_userkey;}
nxRegistryKeyYamlRoot*	nxRegistryKeyYamlRoot::UserApplKey()	{ return &g_userapplkey;}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYamlRoot::nxRegistryKeyYamlRoot		 2016- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyYamlRoot::nxRegistryKeyYamlRoot( const char* filename )
	                  :m_filename(filename)
{
	m_isdirty = false;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYamlRoot::~nxRegistryKeyYamlRoot		 2016- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyYamlRoot::~nxRegistryKeyYamlRoot()
{
	NXASSERT( (!m_isdirty) );					
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYamlRoot::SetYAMLOutputFolder		 2016- 11- 9*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYamlRoot::SetYAMLOutputFolder	( const char* dirname )
{
	nxString	name(dirname);

	if (name.GetLength() > 0)
	{
		name.MakeDirectorySeparatorsOSConsistent();
		name.EnsureLastCharIsDirectoryChar();
		name += m_filename.c_str();
		name.MakeDirectorySeparatorsOSConsistent();
		m_fullfilename = (const char*)name;
	}
	else
	{
		m_fullfilename.clear();
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYamlRoot::CheckDirtyAndSaveIfRequired		 2016- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYamlRoot::CheckDirtyAndSaveIfRequired	()
{
	bool ok;

	ok = !m_isdirty;
	if (!ok)
	{
		boost::lock_guard<boost::mutex> lock(g_yamloutputmutex);			// Lock the threads

		if (m_fullfilename.size() == 0)
		{
			nxLog::Record(NXLOG_WARNING,"nxRegistryKeyYamlRoot::CheckDirtyAndSaveIfRequired, the output yaml filename is empty. Thats should not happen. It needs debugging");
		}
		else
		{
	//		printf("Writing YAML registry to file <%s>\n",(const char*)m_fullfilename.c_str());

			nxFileSpec	spec(m_fullfilename.c_str());
			ok = nxDirectory::CreateADirectory(spec.FullDirSpec() );
			if (ok)
			{
				YAML::Emitter	out;
				out << m_root;
				FILE* f = fopen(m_fullfilename.c_str(), "wt");
				ok =       (f != nullptr) ;
				ok = ok && (fprintf(f, "%s\n", (const char*) out.c_str()) > 0);
				if (f != nullptr) fclose(f);
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"nxRegistryKeyYamlRoot::CheckDirtyAndSaveIfRequired, there were errors writing the yaml text to file <%s>.  Thats not good. It needs debugging", (const char*)m_fullfilename.c_str() );
			}
		}
		m_isdirty = false;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::CheckYamlLoaded		 2016- 11- 7*/
/** LOads an existing yaml file into the root registry. **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYamlRoot::CheckYamlLoaded()
{
	bool ok;
	bool goterrors;

	ok = (!m_root.IsNull()) || m_isdirty;
	if (!ok)
	{
		boost::lock_guard<boost::mutex> lock(g_yamlmutex);			// Lock the threads
		ok = (!m_root.IsNull());									// and check that another thread did not beat us to the load.
		if (!ok)													// If it did not 
		{															// we will do the load
		
			nxString	fullname = nxRegistryKey::RegistryLocation().BaseDirectory();
			fullname.MakeDirectorySeparatorsOSConsistent();
			fullname.EnsureLastCharIsDirectoryChar();
			fullname += m_filename.c_str();
			fullname.MakeDirectorySeparatorsOSConsistent();
			m_fullfilename = (const char*)fullname;
	//		printf("nxRegistryKeyYaml::CheckYamlLoaded, reading file %s\n", (const char*)fullname);


			ok = (!nxDirectory::FileExists(fullname));								// See if the file exists, We are good if it does not
			if (!ok)																// if file does exist
			{																		// then we should load the file
				goterrors = false;													// prepare for the load
				try																	// try an load the YAML file
				{
					m_root = YAML::LoadFile(m_fullfilename);
				}
				catch (...)															// and cactch any errors that may occur
				{	
					goterrors = true;
				}				
				ok = (!goterrors);													// and see if we loaded the YAML file ok (note it may be empty)
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"nxRegistryKeyYaml::CheckYamlLoaded, Errors occurred reading yaml registry file <%s>. Thats is probably going to cause problems later on. You need to run installation scripts to reset the registry values", (const char*)fullname);
				}
			}
		}
		if (!ok)
		{
			m_fullfilename.clear();
			m_root.reset();
		}
		m_isdirty = false;
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::FindRootRegistry		 2016- 11- 10*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyYamlRoot* nxRegistryKeyYaml::FindRootRegistry(nxRegistryConfiguration::INIRootLocation		rootlocation)
{
	nxRegistryKeyYamlRoot*	parentnode = nullptr;

	switch (rootlocation)
	{
		case  nxRegistryConfiguration::USER_APPL_INI   : parentnode = &g_userapplkey;
														 break;

		case  nxRegistryConfiguration::USER_INI        : parentnode = &g_userkey;
														 break;

		case  nxRegistryConfiguration::GLOBAL_INI      : parentnode = &g_globalkey;
														 break;

		case  nxRegistryConfiguration::GLOBAL_APPL_INI : parentnode = &g_globalapplkey;
														 break;

		default										   : parentnode = nullptr;
														 break;
	};
	return parentnode;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::FlushRegistry		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::FlushRegistry( nxRegistryConfiguration::INIRootLocation		rootlocation )
{
	bool ok = true;
	nxRegistryKeyYamlRoot*	parentnode;

	parentnode = FindRootRegistry( rootlocation );
	ok = (parentnode != nullptr);
	if (ok)
	{
		ok = parentnode->CheckDirtyAndSaveIfRequired();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"nxRegistryKeyYaml::FlushRegistry, There were errors flushing the YAML registry to file. The file may be corrupt and that is never good.");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::CreateKey		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyYaml* nxRegistryKeyYaml::CreateKey(  const char *									basekey,
												  const char*									filekeyname,
												  nxRegistryConfiguration::INIRootLocation		rootlocation,
												  nxRegistryConfiguration::INIAccessRights		access)
{
	nxRegistryKeyYaml*		key    = NULL;
	nxStringArray			filetokens;
	bool					ok = true;
	bool					ok1;
	std::string				dirname;
	nxString				basekeyname(basekey);
	nxString				filekey( filekeyname);
	int						i;
	YAML::Node*				thisnode;
	std::list<YAML::Node>	nodelist;
	nxRegistryKeyYamlRoot*	parentnode;



	basekeyname += filekey;

	// ---- Encode the access mode into a base directory location

	parentnode = FindRootRegistry( rootlocation );
	ok = (parentnode != nullptr);
	ok = ok && parentnode->CheckYamlLoaded();
	if (ok)
	{
		basekeyname.MakeLower();
		filetokens.Strtok( basekeyname, "/" );										// Parse the filekey. The first part is the filename
		thisnode = parentnode->Root();
		for (i=0; i < filetokens.GetSize(); i++)									// and reconstruct removing any "doubling of folder chars"
		{
			dirname = (const char*)filetokens.GetAt(i);								// Get the Registry folder name. This becomes a mapped  Node in YAML
			if (dirname.size() > 0)													// if we have any length
			{																		// then see if the  NODE already exists
				ok1 = true;
				YAML::Node	testnode( (*thisnode)[dirname]);						// if it does then save the node on our list of nodes

//				YAML::Emitter	out;
//				out << *(parentnode->Root());
//				printf( "%s\n", (const char*) out.c_str());
																					// and get a pointer to it
				if (testnode)														// if the test node is valid then we are good to go
				{																	// so
					nodelist.push_back( testnode );									// push the node on our list so we have permanent values we can take addresses of)
					thisnode = &nodelist.back();									// and get the addres so fthe new node
				}
				else																// otherwise the node does not exist
				{
					ok1 = (    (access == nxRegistryConfiguration::WRITE_INI)		// we have write access to the values
					        || (access == nxRegistryConfiguration::FULLIO_INI) );	// as we dont create keys in read only mode
					if (ok1)														// if we are good
					{																// the
						(*thisnode)[dirname]["unused"]="x";							// Make a dummy entry to make sure that Node dirname is a map object
						nodelist.push_back( (*thisnode)[dirname] );
						thisnode = &nodelist.back();								// and get the addres so fthe new node
						thisnode->remove("unused");									// and remove the dummy entry

	//					YAML::Emitter	outy;
	//					outy << *(parentnode->Root());
	//					printf( "%s\n", (const char*) outy.c_str());

					}
				}
				ok = ok && ok1;
			}
		}

		if (ok)
		{
			key = new nxRegistryKeyYaml( *thisnode, parentnode, access );						// Create the parent key
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"nxRegistryKeyYaml::CreateKey, There were errors creating the YAML registry key entry <%s>. This usually means you are trying to create a new key in READ only mode. This usually indicates an installation error/problem", (const char*)basekeyname);
		}
	}
//	(*thisnode)["nicksvalue"]="Rain In Spain";
//
//	YAML::Emitter	out;
//	out << *(parentnode->Root());
//	printf("%s\n", (const char*) out.c_str());

	return key;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::nxRegistryKeyYaml		2009-6-2*/
/** **/
/*---------------------------------------------------------------------------*/

nxRegistryKeyYaml::nxRegistryKeyYaml( YAML::Node& node,  nxRegistryKeyYamlRoot* basenode, nxRegistryConfiguration::INIAccessRights accessmode )
                  : m_node( node)

{
	m_accessmode = accessmode;
	m_basenode   = basenode;
}

//---------------------------------------------------------------------------
//						nxRegistryKeyYaml::destructor
//--------------------------------------------------------------------------

nxRegistryKeyYaml::~nxRegistryKeyYaml()
{
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::DestroyKeyHierarchy		 2016- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::DestroyKeyHierarchy	()
{
//	if (    (m_accessmode == nxRegistryConfiguration::INIAccessRights::WRITE_INI) 
//		 || (m_accessmode == nxRegistryConfiguration::INIAccessRights::FULLIO_INI) )
//	{
//		m_basenode->CheckDirtyAndSaveIfRequired();
//	}
	delete this;
	return true;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::SetDoubleValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::SetDouble( const char * valuename, double value)
{
	nxString	str;
	str.sprintf( "%g", (int)value);
	return SetStringInternal(valuename, str );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::GetDoubleValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::GetDouble( const char * valuename, double* value)
{
	nxString	str;
	bool		ok;

	ok = GetString( valuename, &str );
//	str.MakeUpper();
	*value = atof( str );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::SetIntValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::SetInteger( const char * valuename, int value)
{
	nxString	str;
	str.sprintf( "%d", (int)value);
	return SetStringInternal(valuename, str );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::GetIntValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::GetInteger( const char * valuename, int* value)
{
	nxString	str;
	bool		ok;

	ok = GetString( valuename, &str );
//	str.MakeUpper();
	*value = atoi( str );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::SetBoolValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::SetBool( const char * valuename, bool value)
{
	const char*	strtrue = "True";
	const char* strfalse = "False";
	const char* str;

	str = value ? strtrue : strfalse;
	return SetStringInternal( valuename, str );
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::GetBoolValue		2005-3-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::GetBool( const char * valuename, bool* value)
{
	nxString	str;
	bool		ok;

	ok = GetString( valuename, &str );
	str.MakeUpper();
	*value = ok &&  ( (str == "1") || (str == "TRUE") || (str == "T") );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::SetString		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::SetStringInternal( const char * valuename, const char* str  )
{
	bool					ok ;
	nxString				vstr;
	nxString				val(valuename);

	ok = (m_accessmode != nxRegistryConfiguration::READ_INI);				// MAke sure we are not in readonly mode
	if (!ok)																// if we are in readonly mode
	{																		// then log a message
		nxLog::Record(NXLOG_WARNING, "nxRegistryKeyYaml::SetString, You cannot set values in keys that opened in read only mode");
	}
	else																	// otherwsie we are in write mode
	{																		// so
		RootParent()->SetDirty();											// notify the root we are dirty
		val.MakeLower();
		m_node[(const char*)val] = str;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::SetString		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::SetString( const char * valuename, const char* str  )
{
	nxString				vstr;

//	vstr = "'";
	vstr = str;
//	vstr += "'";
	return SetStringInternal(valuename, vstr);
}


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml::SetString		2009-6-3*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxRegistryKeyYaml::GetString( const char * valuename, nxString* str  )
{
	bool					ok = true;
	nxString				vstr;
	std::string				bstr;


	vstr = valuename;
	vstr.MakeLower();
	if (m_node[(const char*)vstr])
	{
		ok   = true;
		bstr = (m_node[(const char*)vstr]).as<std::string>();
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING,"nxRegistryKeyYaml::GetString, Error returning value for <%s> ", (const char*)valuename);
	}
	*str = bstr.c_str();
	return ok;
}


