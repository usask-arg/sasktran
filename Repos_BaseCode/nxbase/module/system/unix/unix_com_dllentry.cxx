/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_core.h"

#if defined(NX_UNIX_VER)
	#include <dlfcn.h>
	#include <sys/stat.h>
	#include <unistd.h>
#else
	extern "C" void * dlopen( const char* m_dllname, int mode );
	extern "C" void * dlsym ( void * handle, const char* functionptr);
	extern "C" const char * dlerror();
	#define RTLD_NOW 1
#endif


#include "unix_com_emmulate.h"


nxComDllEntry::nxComDllEntry()
{
	m_handle         = NULL;
	Close();
}

nxComDllEntry::~nxComDllEntry()
{
//	NXTRACE(( "nxComDllEntry::destructor invoked\n"));
	Close();
//	NXTRACE(( "nxComDllEntry::destructor succeeded\n"));
}

void nxComDllEntry::Close()
{
//	NXTRACE(( "nxComDllEntry::Close, dlclose skipped because we get Segmentation Problems if we unload\n"));
/*
	if (m_handle != NULL)
	{
		FPRINTF(( stderr(("nxComDllEntry::Close, closing DLL <%s> %d\n", (const char *)m_dllname, (int)m_handle));
		dlclose( m_handle );
		FPRINTF(( stderr(("nxComDllEntry::Close, closed down ok\n"));
	}
	m_unload         = NULL;
	m_getclassobject = NULL;
	m_handle         = NULL;
	m_dllname.Empty( nxFALSE);
	memset(&m_class, 0, sizeof(m_class));
*/
}


nxBOOL nxComDllEntry::Load(REFCLSID cls, const char* dllname)
{
	nxBOOL ok;
	const char *	noerror  = "No error";
 	const char *	errorsym = NULL;
	DLLSETREGISTRY	reginitfunc;

	m_class          = cls;
	m_dllname        = dllname;
	m_unload         = NULL;
	m_getclassobject = NULL;

//	printf( "nxComDllEntry::Load, Opening DLL <%s>\n", (const char *)dllname );
	m_handle  = dlopen( m_dllname, RTLD_NOW );
	ok = (m_handle != NULL);
	if (!ok)
	{
		errorsym = dlerror();
		if (errorsym == NULL) errorsym = noerror;
		nxLog::Record(NXLOG_ERROR, "nxComDllEntry::Load, failed to load DLL <%s>, error is <%s>\n", (const char *)m_dllname, (const char *)errorsym);
	}
	else
	{
		reginitfunc = (DLLSETREGISTRY)(dlsym( m_handle, "DllSetRegistryDirectory" ));
		ok = (reginitfunc != NULL);
		if (ok)
		{
//			printf("nxComDllEntry::Load,, base registry = %s\n", (const char*)nxRegistryKey::RegistryLocation().BaseDirectory() );
			reginitfunc(nxRegistryKey::RegistryLocation().BaseDirectory());
		}
		else
		{
 			errorsym = dlerror();
			if (errorsym == NULL) errorsym = noerror;
			nxLog::Record(NXLOG_ERROR,"nxComDllEntry::Load, failed to find symbol DllSetRegistryDirectory in DLL <%s>, error is <%s>\n", (const char *)m_dllname, (const char *)errorsym);
		}


		//NXTRACE(( "nxComDllEntry::Load, returned handle is %p. Now get entry point DllGetClassObject\n", (void *)m_handle));
		m_getclassobject = (DLLGETCLASSOBJECT)(dlsym( m_handle, "DllGetClassObject" ));
		ok = (m_getclassobject != NULL);
		if (!ok)
		{
 			errorsym = dlerror();
			if (errorsym == NULL) errorsym = noerror;
			nxLog::Record(NXLOG_ERROR,"nxComDllEntry::Load, failed to find symbol DllGetClassObject in DLL <%s>, error is <%s>\n", (const char *)m_dllname, (const char *)errorsym);
		}

		else
		{
			//NXTRACE(( "nxComDllEntry::Load, Now look for entry point DllCanUnloadNow\n");)
			m_unload = (DLLCANUNLOADNOW)  (dlsym( m_handle, "DllCanUnloadNow"));
			ok = (m_unload != NULL);
			if (!ok)
			{
 				errorsym = dlerror();
				if (errorsym == NULL) errorsym = noerror;
				nxLog::Record(NXLOG_ERROR,"nxComDllEntry::Load, failed to find symbol DllCanUnloadNow in DLL <%s>, error is <%s>\n", (const char *)m_dllname, (const char *)errorsym);
			}
		}
	}
	if (!ok)
	{
		fprintf(stderr, "nxComDllEntry::Load, error opening <%s>\n", (const char *) m_dllname);
		Close();
	}
	else
	{
		//NXTRACE(("nxComDllEntry::Load, success opening <%s>, the pointers are 0x%p, and 0x%08p\n", (const char *) m_dllname, (void*)m_getclassobject, (void *)m_unload ));
	}
//	printf("nxComDllEntry::Load, finished\n");
	return ok;
}


HRESULT nxComDllEntry::GetClassObject( REFCLSID rclsid, REFIID riid, void **ppv )
{
	HRESULT result;

	if (m_getclassobject != NULL)
	{
		//NXTRACE(("nxComDllEntry::GetClassObject, calling DllGetClassObject (0x%08p)\n", (void *)m_getclassobject));
		result = (*m_getclassobject)(rclsid, riid, ppv );
		//NXTRACE(("nxComDllEntry::GetClassObject, back from calling DllGetClassObject (result = 0x%08x, pointer = 0x%08p)\n", (int)result, (void *)(*ppv) ));
	}
	else
	{
		*ppv = NULL;
		result = E_NOINTERFACE;
	}
	return result;
}



HRESULT nxComDllEntry::UnloadIfUnused()
{
	HRESULT result;

	NXTRACE(( "nxComDllEntry::UnloadIfUnused,  THIS FUNCTION MAY CAUSE PROBLEMS\n"));
	if (m_unload != NULL)
	{
		result = (*m_unload)();
	}
	else
	{
		result = S_FALSE;
	}
	return result;
}
