/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_core.h"
#include "unix_com_emmulate.h"


//---------------------------------------------------------------------------
//						CoInitialize
//	Unix emulation of the Win32 CoInitialize function.  This simply
//	allocates the CLSID registry object and increases the lockcount.
//---------------------------------------------------------------------------

static size_t 	g_lockcount			= 0;

#if defined(NX_UNIX_VER)
const GUID 		IID_IUnknown		= { 0x00000000L, 0x0000, 0x0000, {0xc0,0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46}};
const GUID 		IID_IClassFactory	= { 0x00000001L, 0x0000, 0x0000, {0xc0,0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46}};
#define nxCoInitialize     CoInitialize
#define nxCoUninitialize   CoUninitialize
#define nxCoGetClassObject CoGetClassObject
#define nxCoCreateInstance CoCreateInstance
#define nxVariantInit      VariantInit
#define nxVariantClear     VariantClear
#endif

/*-----------------------------------------------------------------------------
 *					CoInitialize		2009-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" HRESULT nxCoInitialize( LPVOID )
{
	g_lockcount++;
	return S_OK;
}

//---------------------------------------------------------------------------
//						CoUninitialize
//	Unix emulation of the Win32 CoUninitialize function.  This simply
//	decrements the lockcount and deletes the registry object when it hits zero
//---------------------------------------------------------------------------

extern "C" void nxCoUninitialize()
{
	if (g_lockcount == 0)
	{
		nxLog::Record(NXLOG_WARNING,"nxCoUninitialize(), You do dont have balanced calls to CoInitialize and CoUninitialize. Thats a problem.");
	}
	else
	{
		--g_lockcount;
	}
}

/*-----------------------------------------------------------------------------
 *					CoGetClassObject		2009-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" HRESULT nxCoGetClassObject( REFCLSID rclsid, DWORD,	LPVOID,	REFIID riid,	LPVOID * ppv)
{
	HRESULT			status = E_FAIL;
	nxComDllEntry	dll;
	nxString		dllname;
	nxUnixCLSID		clsid;
	bool			ok;

	ok =      clsid.GetDllName( rclsid, &dllname );
	ok = ok && dll.Load( rclsid, dllname );
	if (ok)
	{
		status = dll.GetClassObject( rclsid, riid, ppv );
		ok     = (status == S_OK);
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "CoGetClassObject, USask-ARG Unix implementation of CoGetClassObject failed. Thats a problem.");
		*ppv = NULL;
	}
	return status;
}


/*-----------------------------------------------------------------------------
 *					nxCoCreateInstance		2009-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" HRESULT nxCoCreateInstance( REFCLSID rclsid, LPVOID aggregate, DWORD clsctx, REFIID riid,	LPVOID * ppv)
{
	IClassFactory*	ptr;
	HRESULT status;

	*ppv = NULL;
	status = nxCoGetClassObject( rclsid, clsctx, NULL, IID_IClassFactory, (LPVOID *)&ptr);
	if (SUCCEEDED(status))
	{
		NXTRACE(( "CoCreateInstance, got the factory, now create an instance\n"));
		status = ptr->CreateInstance( (IUnknown*)aggregate, riid, ppv );
		NXTRACE(( "CoCreateInstance, got the object now release the factory\n"));
		ptr->Release();
		NXTRACE(( "CoCreateInstance, done that now back to the caller\n"));
	}
	return status;
}



/*-----------------------------------------------------------------------------
 *					nxVariantInit		2009-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" void nxVariantInit( VARIANT*  pvarg  )
{
	memset( pvarg, 0, sizeof(VARIANT) );
	pvarg->vt = VT_EMPTY;
}


/*-----------------------------------------------------------------------------
 *					nxVariantClear		2009-6-8*/
/** **/
/*---------------------------------------------------------------------------*/

extern "C" HRESULT nxVariantClear( VARIANT *  pvarg )
{
	VariantInit( pvarg );
	return S_OK;
}

//#endif


