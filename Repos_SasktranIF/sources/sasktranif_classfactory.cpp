#include "sasktranif_internals.h"
#include <boost/thread/mutex.hpp>
#include <boost/filesystem.hpp>

#if defined(NX_WINDOWS)
#define LOADLIBRARY(name)				LoadLibrary(name)
#define GETPROCADDRESS( handle, name)	GetProcAddress( handle, name)
#define FREELIBRARY(handle)				FreeLibrary(handle)

#else
#include <dlfcn.h>
#define HMODULE	 void *
#define FARPROC  void *
#define LOADLIBRARY(name)				dlopen(name, RTLD_LAZY)
#define GETPROCADDRESS( handle, name)   dlsym( handle, name)
#define FREELIBRARY(handle)				dlclose(handle)
#endif


/*-----------------------------------------------------------------------------
 *					class SasktranIF_DllEntry		 2015- 8- 26*/
/** A small internal class to hold the lockcount and handle of each DLL loaded
 *	into memory.
 **/
/*---------------------------------------------------------------------------*/

class SasktranIF_DllEntry
{
	public:
		size_t			m_lockcount;
		HMODULE			m_dllhandle;

	public:
		bool			LoadFunctionFromDLL	( const char* functionname, void** funcptr );
};

static  std::map< std::string, SasktranIF_DllEntry>					g_dllmap;			// A static map to hold the names, handles and locakcounts of each DLL
static boost::mutex													g_mutex;			// A mutex to ensure thread safe loading and unloading of DLL's and the associated map.
typedef std::map< std::string, SasktranIF_DllEntry>::iterator		mapiterator;
typedef std::map< std::string, SasktranIF_DllEntry>::value_type		mapvalue_type;


/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::InitializeDLLLogger		2014-3-6*/
/** INitializ ethe logger inside the DLL so it uses the current default logger.
 *	This is really handy to ensure the loaded dll/shareable objects modules
 *	report back to the centralized logging system.
 **/
/*---------------------------------------------------------------------------*/

static bool InitializeDLLLogger( void* funcptr )
{
	typedef				bool  (* logfuncptrtype)( InxLog* );
	logfuncptrtype		logfunc;
	bool				ok;

	ok = (funcptr == nullptr);
	if (!ok)
	{
		logfunc = (logfuncptrtype)funcptr;
		ok = (*logfunc)( nxLogBase::DefaultLog() );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_IFSetRegistryDirectory		 2016- 11- 14*/
/** This is a function exported  by SasktranIF to external programs so that Python (and others)
 *	can initialize the registry directory of this DLL when they load the SasktranIF module from within python.
 *  The SasktranIF module, in turn, initializes any SasktranIF component DLL's that it loads.
 **/
/*---------------------------------------------------------------------------*/

extern "C" bool SKTRAN_IFSetRegistryDirectory(const char* registrydirname)
{
	bool	ok = true;
	
	if (registrydirname == nullptr) 
	{
		ok = nxRegistryKey::RegistryLocation().UseNativeRegistry();
		if (!ok)
		{
			printf("SasktranIF Internal Registry Initialization::SKTRAN_IFSetRegistryDirectory, the caller has requested using the native registry but that is not available on this build");
		}
	}
	else
	{

		nxString   regdir(registrydirname);
		ok =       regdir.MakeDirectorySeparatorsOSConsistent();
		ok = ok && nxRegistryKey::RegistryLocationVar()->Set_BaseDirectory(registrydirname);
		if (!ok)
		{
			printf("SasktranIF Internal Registry Initialization::SKTRAN_IFSetRegistryDirectory, there were errors setting the registry to use directory <%s>", (const char*)registrydirname);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					CreateRegistryEntries		 2016- 11- 14*/
/** **/
/*---------------------------------------------------------------------------*/

static bool CreateRegistryEntries( void* funcptr, const char* paramstr )
{
	typedef				bool  (* regfuncptrtype)( const char*);
	regfuncptrtype		logfunc;
	bool				ok;


	ok = (funcptr == nullptr);
	if (!ok)
	{
		logfunc = (regfuncptrtype)funcptr;
		ok = (*logfunc)(paramstr);

	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::InitializeDLLLogger		2014-3-6*/
/** INitializ ethe logger inside the DLL so it uses the current default logger.
 *	This is really handy to ensure the loaded dll/shareable objects modules
 *	report back to the centralized logging system.
 **/
/*---------------------------------------------------------------------------*/

static bool InitializeDLLRegistryDir( void* funcptr )
{
	typedef				bool  (* logfuncptrtype)( const char* );
	logfuncptrtype		logfunc;
	bool				ok;
	const char*			dirname;


	ok = (funcptr == nullptr);
	if (!ok)
	{
		if (nxRegistryKey::RegistryLocation().UseNativeRegistry())
		{
			// We are using the native registry, so we don't have to pass anything into 
			// the sub object
			ok = true;
		}
		else
		{
			// We are using a local sandboxed registry, pass its location into the sub object
			dirname =  (const char*)nxRegistryKey::RegistryLocation().BaseDirectory();
			logfunc = (logfuncptrtype)funcptr;
			ok = (*logfunc)( dirname);
		}

	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					InitializeGlobalHandleTable		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

static bool InitializeGlobalHandleTable( void* funcptr )
{
	typedef								bool  (* handlefuncptrtype)( GlobalClimatologyHandleTable** entry, int* numpoints );
	handlefuncptrtype					logfunc;
	GlobalClimatologyHandleTable*		entry;
	int									npts;
	bool								ok;

	ok = (funcptr == nullptr);
	if (!ok)
	{
		logfunc = (handlefuncptrtype)funcptr;
		ok = (*logfunc)( &entry, &npts);
		if (ok && (entry != NULL))
		{
			for ( int i = 0; i < npts; i++)
			{
				AddGlobalClimatologyHandle ( entry[i].name, entry[i].handle);
			}
		}
	}
	return ok;
}

static bool InitializeParentHandleTable(void* funcptr)
{
	typedef								bool(*handlefuncptrtype)(std::map<nxString, CLIMATOLOGY_HANDLE>* parenttable);
	handlefuncptrtype					logfunc = (handlefuncptrtype)funcptr;

	return logfunc(InternalGlobalClimatologyHandleTable());
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_IFCreateRegistryEntriesForDLL		 2015- 8- 26*/
/** This is a function typically called during installation of a Sasktran 
 *	component. The process is
 *	- python component finds out is dll file name using standard sys and glob python modules. Note the component DLL is typically in the modules site-package folder
 *  - python component then imports sasktranif and calls sasktranif.sasktranif.SKTRAN_IFCreateRegistryEntriesForDLL( modulename ) to create its registry entries
 **/
/*---------------------------------------------------------------------------*/

extern "C" bool SKTRAN_IFCreateRegistryEntriesForDLL( const char* dllname, const char* userparamstr )
{
	bool				ok;
	HMODULE				dllhandle;
	nxString			paramstr(userparamstr);
	nxFileSpec			dllpath(dllname);

#if defined(NX_WINDOWS)
	SetDllDirectory( (const char*)(const char*)dllpath.FullDirSpec() );
#else	
	boost::filesystem::path cwd = boost::filesystem::current_path();
	boost::filesystem::path newdirectory( (const char*)dllpath.FullDirSpec() );
	//printf("Changing directory from <%s> to <%s>\n",(const char*)newdirectory.string().c_str(), (const char*)newdirectory.string().c_str());
	boost::filesystem ::current_path( newdirectory );
#endif
	dllhandle = LOADLIBRARY(dllname);			// Now lets load the library
	ok = (dllhandle != NULL);					// and make sure it worked
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF::CreateRegistryEntriesForDLL, Cannot find DLL/Shareable object <%s>. This probably indicates an incorrect installation or try adjusting the PATH (or LD_LIBRARY_PATH) to include this directory", (const char *)dllname);
#if defined(NX_WINDOWS)
		nxLog::Record(NXLOG_WARNING,"SasktranIF::CreateRegistryEntriesForDLL, Cannot find DLL/Shareable object <%s>. This probably indicates an incorrect installation or try adjusting the PATH (or LD_LIBRARY_PATH) to include this directory", (const char *)dllname);
#else
		const char* errmes = dlerror();
		nxLog::Record(NXLOG_WARNING,"SasktranIF::CreateRegistryEntriesForDLL, Cannot find DLL/Shareable object <%s>. dlopen reports error <%s>.", (const char *)dllname, (const char*)errmes);
#endif
	}
	else																					// We have the DLL load
	{																						// So
		FARPROC logproc;

		logproc = GETPROCADDRESS( dllhandle, "SKTRAN_IFSetRegistryDirectoryInChildDLL");		// Do we have a hook to pass in the registry information
		ok = (logproc != NULL) && InitializeDLLRegistryDir( logproc);								// If we do then initialize the registry directory
		if (!ok)
		{
			if (logproc  == NULL)
			{
				nxLog::Record(NXLOG_WARNING,"SasktranIF Creating Registry Entries: Error loooking up function SKTRAN_IFSetRegistryDirectoryInChildDLL in the DLL or shareable object. That probably means its not implemented in the C++ code, which is not good");
			}
			else
			{
				nxLog::Record(NXLOG_WARNING,"SasktranIF Creating Registry Entries: Error setting the registry directory to defaults. Entry point SKTRAN_IFSetRegistryDirectoryInChildDLL exists but failed during execution. Thats not good");
			}
		}
		else		
		{
			
			logproc = GETPROCADDRESS( dllhandle, "SKTRAN_IFCreateRegistryEntriesForChildDLL");			// Do we have a hook to pass in the registry creatinon folder
			ok = (logproc != NULL) && CreateRegistryEntries( logproc, paramstr);						// and then pass in the registry creation parameters

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SasktranIF Creating Registry Entries: there were errors creating the registry setting for DLL <%s> with parameters <%s>", (const char*)dllname, (const char*)paramstr);
			}
		}
	}
#if defined(NX_WINDOWS)
	FreeLibrary(dllhandle);
#else
	boost::filesystem ::current_path( cwd );
#endif
	return ok;
}




/*-----------------------------------------------------------------------------
 *					LoadFunctionFromDLL		 2015- 8- 26*/
/** **/
/*---------------------------------------------------------------------------*/

static bool LoadFunctionFromDLL( const std::string& dllname, const char* functionname, void** funcptr )
{
	bool				ok;
	FARPROC				proc = NULL;

	g_mutex.lock();
	auto iter = g_dllmap.find(dllname);
	ok = ! (iter == g_dllmap.end());						// HAve we already loaded this DLL into our code
	if (!ok)												// Nope
	{														// so
		SasktranIF_DllEntry	entry;							// Lets create an entry for this DLL
		std::pair<mapiterator, bool>	answer;

		nxFileSpec			dllpath(dllname.c_str());
#if defined(NX_WINDOWS)
		SetDllDirectory( (const char*)(const char*)dllpath.FullDirSpec() );
#else	
		boost::filesystem::path cwd = boost::filesystem::current_path();
		boost::filesystem::path newdirectory( (const char*)dllpath.FullDirSpec() );
		//printf("Changing directory from <%s> to <%s>\n",(const char*)newdirectory.string().c_str(), (const char*)newdirectory.string().c_str());
		boost::filesystem ::current_path( newdirectory );
#endif

		entry.m_lockcount = 0;								// The lock count on the DLL will be 0
		entry.m_dllhandle = LOADLIBRARY(dllname.c_str());			// Now lets load the library
		ok = (entry.m_dllhandle != NULL);					// and make sure it worked
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SasktranIF::LoadFunctionFromDLL, Cannot find DLL/Shareable object <%s>. This probably indicates an incorrect installation or try adjusting the PATH (or LD_LIBRARY_PATH) to include this directory", (const char *)dllname.c_str());
		}
		else																					// We have the DLL load
		{																						// So
			if (g_dllmap.find(dllname) == g_dllmap.end())
			{
				FARPROC logproc;
			
				logproc = GETPROCADDRESS( entry.m_dllhandle, "SKTRAN_IFSetRegistryDirectoryInChildDLL");	// Do we have a hook to pass in the registry information
				if (logproc != NULL) InitializeDLLRegistryDir( logproc);									// If we do then initialize the registry directory
				logproc = GETPROCADDRESS( entry.m_dllhandle, "SKTRAN_IFInitializeLogger");					// Do we have a hook to pass in the logger information
				if (logproc != NULL) InitializeDLLLogger( logproc);											// If we do then initialize our logger
				logproc = GETPROCADDRESS( entry.m_dllhandle, "SKTRAN_IFGlobalHandleTable");
				if (logproc != NULL) InitializeGlobalHandleTable( logproc);
				logproc = GETPROCADDRESS(entry.m_dllhandle, "SKTRAN_IFSetParentHandleTable");
				if (logproc != NULL)
				{
					InitializeParentHandleTable(logproc);
				}
				answer = g_dllmap.insert( mapvalue_type( dllname, entry) );								// Insert this DLL and its handle into the map of stored DLL's
				ok   = answer.second;																	// See if it worked
				iter = answer.first;																	// and get the pointer to the map entry.
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"SasktranIF::LoadFunctionFromDLL, Error inserting DLL/Shareable object <%s> into our internal map", (const char*)dllname.c_str());
				}
			}
		}
#if !defined(NX_WINDOWS)
		boost::filesystem ::current_path( cwd );
#endif


	}
	if (ok)
	{
		proc = GETPROCADDRESS( (*iter).second.m_dllhandle, functionname);
		ok = (proc != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SasktranIF::LoadFunctionFromDLL, Error locating function <%s> inside DLL/shareable objects. This probably indicates an installation error or deficiency in the DLL implementation", (const char*)functionname, (const char*)dllname.c_str() );
		}
		else
		{
			(*iter).second.m_lockcount++;							// Increment the lock count 
		}
	}
	*funcptr = ok ? (void*)proc : NULL;
	g_mutex.unlock();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ReleaseDLL_LockCount		 2015- 8- 26*/
/** **/
/*---------------------------------------------------------------------------*/

static bool ReleaseDLL_LockCount( const nxString& userdllname )
{
	bool				ok;
	std::string			dllname( (const char*)userdllname);

	g_mutex.lock();
	auto iter = g_dllmap.find(dllname);

	ok = ! (iter == g_dllmap.end());						// HAve we already loaded this DLL into our code
	if (ok)												// Nope
	{	
		(*iter).second.m_lockcount--;
		if ((*iter).second.m_lockcount == 0)
		{
			FREELIBRARY((*iter).second.m_dllhandle);
			g_dllmap.erase(iter);
		}
	}
	g_mutex.unlock();
	return ok;
}


/*---------------------------------------------------------------------------
 *                  ISKModuleBase::ISKModuleBase                  2019-06-13 */
/** **/
/*---------------------------------------------------------------------------*/

ISKModuleBase::ISKModuleBase()
{
	m_dllname = nullptr;
}

/*-----------------------------------------------------------------------------
 *					ISKModuleBase::~ISKModuleBase		 2015- 8- 26*/
/** Whenever we destroy a SasktranIF object we must decrement the reference
 *	count on the DLL/shareable object. A call is made to unload the DLL when the
 *	reference count goes to zero.
 **/
/*---------------------------------------------------------------------------*/

ISKModuleBase::~ISKModuleBase()
{
	if (m_dllname != nullptr) delete [] m_dllname;
	//ReleaseDLL_LockCount(m_dllname); THe DLL unload functionality has been commented out until it is better implemented. CReates a lot of loading/unloading issues
}



/*-----------------------------------------------------------------------------
 *					ISKModuleBase::SetProperty		 2015- 11- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase::SetProperty( const char* propertyname, void* valueorobject, int numpoints_or_type  )
{
	bool ok;

	if      (numpoints_or_type >=  0) ok = SetPropertyArray ( propertyname,   (const double*)valueorobject, numpoints_or_type);
	else if (numpoints_or_type == -1) ok = SetPropertyScalar( propertyname, *((const double*)valueorobject) );
	else if (numpoints_or_type == -2) ok = SetPropertyObject( propertyname, (ISKModuleBase*)valueorobject );
	else if (numpoints_or_type == -3) 
	{
		nxString	str( (const char*)valueorobject );
		ok = SetPropertyString( propertyname, str);
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"ISKModuleBase::SetProperty, invalid value (%d) for numpoints_or_type", (int)numpoints_or_type);
		ok = false;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKModuleBase::GetProperty		 2015- 11- 4*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase::GetProperty( const char* propertyname, const double** value, int* numpoints )
{
	nxLog::Record(NXLOG_WARNING,"This SasktranIF object does not support property %s", (const char*)propertyname);
	*value = nullptr;
	*numpoints = 0;
	return false;
}


/*-----------------------------------------------------------------------------
 *					SasktranIF_FindRegistrySetting		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::FindRegistrySetting( const char*classname, const char* entityname, std::string* userdllname )
{
	nxString	dllname;
	nxString	str;
	bool		ok;

	str.sprintf( "/SasktranIF/%s/%s/", (const char*)classname, (const char*)entityname );								// Get the full registry name
	nxRegistryConfiguration		register_object		( "USask-ARG", str, nxRegistryConfiguration::GLOBAL_INI, true);		// Connect to the entry

	ok = register_object.GetString("DLLName", &dllname);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF C++ Interface, Cannot find the DLLName registry entry for %s", (const char*)str);
		dllname.Empty(false);
	}
	userdllname->assign( (const char*)dllname );
	return ok;
}


/*---------------------------------------------------------------------------
 *         SasktranIF_ClassFactoryLocator::AssignDLLname          2019-06-13 */
/** Assigns the DLLName to the  interface object we are creating
**/
/*---------------------------------------------------------------------------*/

void SasktranIF_ClassFactoryLocator::AssignDLLname( char**userdllname, const std::string& dllname)
{
	if (*userdllname != nullptr) delete [] * userdllname;
	size_t n = dllname.size() + 1;
	*userdllname = new char[n + 1];
	strncpy(*userdllname, dllname.c_str(), n);
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKEngine		2014-2-8*/
/** Finds the DLL name from the registry /SasktranIF/Engines/enginename/
 *	Finds a pointer to the function SKTRANIF_CreateEngine in the DLL.
 *	Calls the function as SKTRANIF_CreateEngine( const char* enginename, ISKEngine_Stub** engine)
 *	If successful returns the engine with a reference count of 1.
**/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKEngine( const char* enginename, ISKEngine_Stub** engine, char** userdllname )
{ 
	typedef			bool (* funcptrtype)( const char*, ISKEngine_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok =       FindRegistrySetting( "Engines", enginename, &dllname );
	ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateEngine2", &voidfuncptr );
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(enginename, engine);
	}
	if (!ok)
	{
		*engine = NULL;
		nxLog::Record(NXLOG_WARNING,"SasktranIF_ClassFactoryLocator::CreateISKEngine, Error creating SasktranIF engine [%s]. This usually indicates a configuration issue", (const char*)enginename);
	}
	else (*engine)->AddRef();
	AssignDLLname( userdllname, dllname);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKEngine		2014-2-8*/
/** Finds the DLL name from the registry /SasktranIF/Climatology/climatename/
 *	Finds a pointer to the function SKTRANIF_CreateClimatology in the DLL.
 *	Calls the function as SKTRANIF_CreateClimatology( const char* climatename, ISKClimatology_Stub** engine)
 *	If successful returns the engine with a reference count of 1.
**/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKClimatology( const char* climatename, ISKClimatology_Stub** climatology, char** userdllname)
{
	typedef			bool (* funcptrtype)( const char*, ISKClimatology_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok =       FindRegistrySetting( "Climatology", climatename, &dllname );
	ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateClimatology2", &voidfuncptr );
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(climatename, climatology);
	}
	if (!ok)
	{
		*climatology = NULL;
		nxLog::Record(NXLOG_WARNING,"SasktranIF_ClassFactoryLocator::CreateISKClimatology, Error creating climatology [%s]. This usually indicates a configuration issue", (const char*)climatename);
	}
	else (*climatology)->AddRef();
	AssignDLLname(userdllname, dllname);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKEngine		2014-2-8*/
/** Finds the DLL name from the registry /SasktranIF/OpticalProperty/optpropname/
 *	Finds a pointer to the function SKTRANIF_CreateOptProp in the DLL.
 *	Calls the function as SKTRANIF_CreateOpticalProperty( const char* enginename, ISKOpticalProperty_Stub** engine)
 *	If successful returns the engine with a reference count of 1.
**/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKOpticalProperty( const char* optpropname, ISKOpticalProperty_Stub** optprop, char** userdllname)
{
	typedef			bool  (* funcptrtype)( const char*, ISKOpticalProperty_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	if ( strlen(optpropname) == 0 )			
	{
		ok = true;
		*optprop = nullptr;
		AssignDLLname(userdllname, std::string("") );
	}
	else
	{
		ok =       FindRegistrySetting( "OpticalProperty", optpropname, &dllname );
		ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateOpticalProperty2", &voidfuncptr );
		if (ok)
		{
			func = (funcptrtype)voidfuncptr;
			ok = (*func)(optpropname, optprop);
		}
		if (!ok) 
		{
			nxLog::Record(NXLOG_WARNING,"SasktranIF_ClassFactoryLocator::CreateISKOpticalProperty, Error creating optical property [%s]. This usually indicates a configuration issue", (const char*)optpropname);
			*optprop = NULL;
		}
		else (*optprop)->AddRef();
		AssignDLLname(userdllname, dllname);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKEmission		 2015- 3- 11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKEmission( const char* optpropname, ISKEmission_Stub** optprop, char** userdllname)
{
	typedef			bool  (* funcptrtype)( const char*, ISKEmission_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok =       FindRegistrySetting( "Emission", optpropname, &dllname );
	ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateEmission2", &voidfuncptr );
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(optpropname, optprop);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF_ClassFactoryLocator::CreateISKEmission, Error creating emission object [%s]. This usually indicates a configuration issue", (const char*)optpropname);
		*optprop = NULL;
	}
	else (*optprop)->AddRef();
	AssignDLLname(userdllname, dllname);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKBrdf		 2015- 3- 11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKBrdf( const char* brdfname, ISKBrdf_Stub** brdf, char** userdllname )
{
	typedef			bool  (* funcptrtype)( const char*, ISKBrdf_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok =       FindRegistrySetting( "BRDF", brdfname, &dllname );
	ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateBRDF2", &voidfuncptr );
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(brdfname, brdf);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF_ClassFactoryLocator::CreateISKBrdf, Error creating BRDF object [%s]. This usually indicates a configuration issue", (const char*)brdfname);
		*brdf = NULL;
	}
	else (*brdf)->AddRef();
	AssignDLLname(userdllname, dllname);
	return ok;
}

bool SasktranIF_ClassFactoryLocator::CreateISKMie(const char* miename, ISKMie_Stub** mie, char** userdllname) {
	typedef			bool  (*funcptrtype)(const char*, ISKMie_Stub**);
	void* voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok = FindRegistrySetting("MIE", miename, &dllname);
	ok = ok && LoadFunctionFromDLL(dllname, "SKTRANIF_CreateMie2", &voidfuncptr);
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(miename, mie);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SasktranIF_ClassFactoryLocator::CreateISKMie, Error creating Mie object [%s]. This usually indicates a configuration issue", (const char*)miename);
		*mie = NULL;
	}
	else (*mie)->AddRef();
	AssignDLLname(userdllname, dllname);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKSolarSpectrum		 2015- 1- 7*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKSolarSpectrum( const char* solarname, ISKSolarSpectrum_Stub** solar, char** userdllname)
{
	typedef			bool  (* funcptrtype)( const char*, ISKSolarSpectrum_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok =       FindRegistrySetting( "SolarSpectrum", solarname, &dllname );
	ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateSolarSpectrum2", &voidfuncptr );
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(solarname, solar);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF_ClassFactoryLocator::CreateISKSolarSpectrum, Error creating solar spectrum [%s]. This usually indicates a configuration issue", (const char*)solarname);
		*solar = NULL;
	}
	else (*solar)->AddRef();
	AssignDLLname(userdllname, dllname);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SasktranIF_ClassFactoryLocator::CreateISKGeodetic		 2014- 5- 6*/
/** Finds the DLL name from the registry /SasktranIF/Geodetic/optpropname/
 *	Finds a pointer to the function SKTRANIF_CreateOptProp in the DLL.
 *	Calls the function as SKTRANIF_CreateOpticalProperty( const char* enginename, ISKOpticalProperty_Stub** engine)
 *	If successful returns the engine with a reference count of 1.
 **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_ClassFactoryLocator::CreateISKGeodetic( const char* propname, ISKGeodetic_Stub** geoid, char** userdllname)
{
	typedef			bool  (* funcptrtype)( const char*, ISKGeodetic_Stub** );
	void*			voidfuncptr;
	funcptrtype		func;
	bool			ok;
	std::string		dllname;

	ok =       FindRegistrySetting( "Geodetic", propname, &dllname );
	ok = ok && LoadFunctionFromDLL( dllname, "SKTRANIF_CreateGeodetic2", &voidfuncptr );
	if (ok)
	{
		func = (funcptrtype)voidfuncptr;
		ok = (*func)(propname, geoid);
	}
	if (!ok) *geoid = nullptr;
	else (*geoid)->AddRef();
	AssignDLLname(userdllname, dllname);
	return ok;
}
